#ifndef FAKE_SIMPLEX_TREE_H
#define FAKE_SIMPLEX_TREE_H

#include <gudhi/Filtered_toplex_map.h>
#include <boost/graph/adjacency_list.hpp>

namespace Gudhi {

class Fake_simplex_tree : public Filtered_toplex_map {

public:

    typedef Vertex Vertex_handle;

    typedef Simplex_ptr Simplex_handle;

    typedef void Insertion_result_type;

    /** \brief Inserts a given range `Gudhi::rips_complex::Rips_complex::OneSkeletonGraph` in the simplicial
     * complex. */
    template<class OneSkeletonGraph>
    void insert_graph(const OneSkeletonGraph& skel_graph);

    /** \brief Expands the simplicial complex containing only its one skeleton until a given maximal dimension as
     * explained in \ref ripsdefinition. */
    void expansion(int max_dim);

    /** \brief Returns the number of vertices in the simplicial complex. */
    std::size_t num_vertices();

    Simplex_ptr_set candidates() const;

    std::size_t dimension() const;

    std::size_t num_simplices() const;

    std::size_t num_vertices() const;

    Simplex simplex_vertex_range(Simplex_ptr &sptr) const;

    std::vector<Simplex_ptr> max_simplices() const;

    std::unordered_set<Simplex_ptr> filtration_simplex_range() const;

    std::unordered_set<Simplex_ptr> skeleton_simplex_range(int d=std::numeric_limits<int>::max()) const;

    std::size_t dimension(Simplex_ptr& sptr) const;

    void assign_filtration(Simplex_ptr& f_simplex, Filtration_value alpha_complex_filtration);

    void make_filtration_non_decreasing();

protected:

    /** \internal Does all the facets of the given simplex belong to the complex ?
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    bool all_facets_inside(const Input_vertex_range &vertex_range) const;

};

template<class OneSkeletonGraph>
void Fake_simplex_tree::insert_graph(const OneSkeletonGraph& skel_graph){
    typename boost::graph_traits<OneSkeletonGraph>::edge_iterator e_it,
            e_it_end;
    for (std::tie(e_it, e_it_end) = boost::edges(skel_graph); e_it != e_it_end; ++e_it) {
        auto u = source(*e_it, skel_graph);
        auto v = target(*e_it, skel_graph);
        if(u<v){
            Simplex s;
            s.emplace(u);
            s.emplace(v);
            insert_simplex_and_subfaces(s);
        }
    }
}

void Fake_simplex_tree::expansion(int max_dim){
    for(int d=3; d <= max_dim; d++){
        auto cs = candidates();
        if(cs.empty()) return;
        for(const Simplex_ptr& sptr: cs){
            Simplex sigma = *sptr;
            insert_simplex_and_subfaces(sigma);
        }
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
    std::unordered_map<Simplex_ptr, std::vector<Vertex>, Toplex_map::Sptr_hash, Toplex_map::Sptr_equal> facets_to_max;
    for(const auto& kv : filtrations){
        Simplex sigma (*(kv.first));
        for(Vertex v : sigma){
            sigma.erase(v);
            auto sptr = get_key(sigma);
            if(!facets_to_max.count(sptr)) facets_to_max.emplace(sptr, std::vector<Vertex>());
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
    return max;
}

std::size_t Fake_simplex_tree::num_simplices() const {
    return filtration_simplex_range().size();
}

std::size_t Fake_simplex_tree::num_vertices() const {
    std::unordered_set<Vertex> vertices;
    for(auto kv : filtrations)
        for (Vertex v : *(kv.first))
            vertices.emplace(v);
    return vertices.size();
}

Simplex Fake_simplex_tree::simplex_vertex_range(Simplex_ptr& sptr) const {
    return *sptr;
}

std::unordered_set<Simplex_ptr> Fake_simplex_tree::filtration_simplex_range() const{
    std::vector<Simplex_ptr> m = max_simplices();
    std::unordered_set<Simplex_ptr> seen;
    while(m.begin()!=m.end()){
        Simplex_ptr& sptr = m.back();
        m.pop_back();
        if(seen.find(sptr)!=seen.end()){
            seen.emplace(sptr);
            for(Simplex& sigma : facets(*sptr))
                m.emplace_back(get_key(sigma));
        }
    }
    return seen;
}

std::unordered_set<Simplex_ptr> Fake_simplex_tree::skeleton_simplex_range(int d) const{
    std::unordered_set<Simplex_ptr> simplices;
    for(auto sptr: filtration_simplex_range())
        if(sptr->size()<=d)
            simplices.emplace(sptr);
    return simplices;
}

std::vector<Simplex_ptr> Fake_simplex_tree::max_simplices() const{
    std::vector<Simplex_ptr> s;
    for(auto kv : filtrations)
        s.emplace_back(kv.first);
    return s;
}

std::size_t Fake_simplex_tree::dimension(Simplex_ptr& sptr) const{
    return sptr->size();
}


void Fake_simplex_tree::assign_filtration(Simplex_ptr& f_simplex, Filtration_value alpha_complex_filtration){
    filtrations.emplace(f_simplex,alpha_complex_filtration);
}

void Fake_simplex_tree::make_filtration_non_decreasing(){
    for(auto yt = filtrations.begin(); yt != filtrations.end(); ++yt)
        for (auto it = toplex_maps.begin(); it != toplex_maps.end(); ++it){
            if(it->first == yt -> second)
                break;
            if(it->second.membership(*(yt->first)))
                for(const Simplex_ptr& sptr : it->second.maximal_cofaces(*(yt->first))){
                    it->second.erase_maximal(sptr);
                    toplex_maps.at(yt->second).insert_simplex(*sptr);
                    filtrations.emplace(sptr,yt->second);
                }
        }

}



} //namespace Gudhi

#endif /* FAKE_SIMPLEX_TREE_H */

