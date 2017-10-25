#ifndef FAKE_SIMPLEX_TREE_H
#define FAKE_SIMPLEX_TREE_H

#include <gudhi/Simplex_tree.h>
#include <gudhi/Filtered_toplex_map.h>

#include <boost/graph/adjacency_list.hpp>


namespace Gudhi {

class Fake_simplex_tree : public Filtered_toplex_map {

public:

    typedef Vertex Vertex_handle;

    typedef Simplex Simplex_handle;

    typedef void Insertion_result_type;

    /** \brief Inserts a given range `Gudhi::rips_complex::Rips_complex::OneSkeletonGraph` in the simplicial
     * complex. */
    template<class OneSkeletonGraph>
    void insert_graph(const OneSkeletonGraph& skel_graph);

    /** \brief Expands the simplicial complex containing only its one skeleton until a given maximal dimension as
     * explained in \ref ripsdefinition. */
    void expansion(int max_dim);

    /** \brief Returns the number of vertices in the simplicial complex. */
    std::size_t num_vertices() const;

    Simplex_ptr_set candidates() const;

    std::size_t dimension() const;

    std::size_t num_simplices() const;

    void set_dimension(int d);

    Simplex simplex_vertex_range(const Simplex& s) const;

    std::vector<Simplex> max_simplices() const;

    std::vector<Simplex> filtration_simplex_range() const;

    std::vector<Simplex> skeleton_simplex_range(int d=std::numeric_limits<int>::max()) const;

    std::size_t dimension(Simplex_ptr& sptr) const;

protected:

    /** \internal Does all the facets of the given simplex belong to the complex ?
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    bool all_facets_inside(const Input_vertex_range &vertex_range) const;

};

void Fake_simplex_tree::set_dimension(int d){

}

template<class OneSkeletonGraph>
void Fake_simplex_tree::insert_graph(const OneSkeletonGraph& skel_graph){
    if (boost::num_vertices(skel_graph) == 0) return;
    typename boost::graph_traits<OneSkeletonGraph>::vertex_iterator v_it, v_it_end;
    for (std::tie(v_it, v_it_end) = boost::vertices(skel_graph); v_it != v_it_end; ++v_it){
        Simplex s;
        s.insert(*v_it);
        insert_simplex_and_subfaces(s, boost::get(vertex_filtration_t(), skel_graph, *v_it));
    }

    typename boost::graph_traits<OneSkeletonGraph>::edge_iterator e_it, e_it_end;
    for (std::tie(e_it, e_it_end) = boost::edges(skel_graph); e_it != e_it_end; ++e_it) {
        Vertex u = source(*e_it, skel_graph);
        Vertex v = target(*e_it, skel_graph);
        if (u < v) {
            Simplex s;
            s.insert(u);
            s.insert(v);
            insert_simplex_and_subfaces(s, boost::get(edge_filtration_t(), skel_graph, *e_it));
        }
    }
}

void Fake_simplex_tree::expansion(int max_dim){
    for(int d=2; d <= max_dim; d++){
        Simplex_ptr_set cs = candidates(); //dimension ?
        if(cs.empty()) std::cout << d << std::endl;
        if(cs.empty()) return;
        for(const Simplex_ptr& sptr: cs)
            insert_simplex_and_subfaces(*sptr); //filtration ?
    }
}

template <typename Input_vertex_range>
bool Fake_simplex_tree::all_facets_inside(const Input_vertex_range &vertex_range) const{
    Simplex sigma(vertex_range);
    for(const Simplex& s : facets(sigma))
        if(!filtrations.count(get_key(s))) return false;
    return true;
}

Simplex_ptr_set Fake_simplex_tree::candidates() const{
    Simplex_ptr_set c;
    std::unordered_map<Simplex_ptr, std::vector<Vertex>, Sptr_hash, Sptr_equal> facets_to_max;
    for(const auto& kv : filtrations){
        Simplex sigma (*(kv.first));
        if(sigma.size()>1)
            for(Vertex v : *(kv.first)){
                sigma.erase(v);
                auto sptr = get_key(sigma);
                if(!facets_to_max.count(sptr))
                    facets_to_max.emplace(sptr, std::vector<Vertex>());
                facets_to_max.at(sptr).emplace_back(v);
                sigma.insert(v);
            }
    }
    for(const auto& kv : facets_to_max){
        std::unordered_set<Vertex> facets(kv.second.begin(), kv.second.end());
        for(Vertex v : kv.second){
            facets.erase(v);
            for(Vertex w : facets){
                Simplex sigma(*(kv.first));
                sigma.insert(v);
                sigma.insert(w);
                if(all_facets_inside(sigma))
                    c.emplace(get_key(sigma));
            }
            facets.emplace(v);
        }
    }
    return c;
}

std::size_t Fake_simplex_tree::dimension() const {
    std::size_t max = 0;
    for(auto kv : filtrations)
        max = std::max(max, kv.first->size());
    return max-1;
}

std::size_t Fake_simplex_tree::num_simplices() const {
    //return filtration_simplex_range().size();
    return max_simplices().size();
}

std::size_t Fake_simplex_tree::num_vertices() const {
    std::unordered_set<Vertex> vertices;
    for(auto kv : filtrations)
        for (Vertex v : *(kv.first))
            vertices.emplace(v);
    return vertices.size();
}

Simplex Fake_simplex_tree::simplex_vertex_range(const Simplex& s) const {
    return s;
}

std::vector<Simplex> Fake_simplex_tree::filtration_simplex_range() const{
    std::vector<Simplex> m = max_simplices();
    std::vector<Simplex> seen1;
    Simplex_ptr_set seen2;
    while(m.begin()!=m.end()){
        Simplex s(m.back());
        m.pop_back();
        if(seen2.find(get_key(s))==seen2.end()){
            seen1.emplace_back(s);
            seen2.emplace(get_key(s));
            if(s.size()>0)
                for(Simplex& sigma : facets(s))
                    m.emplace_back(sigma);
        }
    }
    return seen1;
}

std::vector<Simplex> Fake_simplex_tree::skeleton_simplex_range(int d) const{
    std::vector<Simplex> simplices;
    for(auto s: filtration_simplex_range())
        if(s.size()<=d)
            simplices.emplace_back(s);
    return simplices;
}

std::vector<Simplex> Fake_simplex_tree::max_simplices() const{
    std::vector<Simplex> s;
    for(auto kv : filtrations)
        s.emplace_back(*(kv.first));
    return s;
}

std::size_t Fake_simplex_tree::dimension(Simplex_ptr& sptr) const{
    return sptr->size();
}

} //namespace Gudhi

#endif /* FAKE_SIMPLEX_TREE_H */

