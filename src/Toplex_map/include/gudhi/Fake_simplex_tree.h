#ifndef FAKE_SIMPLEX_TREE_H
#define FAKE_SIMPLEX_TREE_H

#include <gudhi/Simplex_tree.h>
#include <gudhi/Filtered_toplex_map.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>

#define filtration_upper_bound std::numeric_limits<Filtration_value>::max()

namespace Gudhi {

struct Visitor {
    Toplex_map* tm;

    Visitor(Toplex_map* tm)
        :tm(tm)
    {}

    template <typename Clique, typename Graph>
    void clique(const Clique& c, const Graph& g)
    {
        tm->insert_simplex(c);
    }
};

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

    Simplex_ptr_set candidates(int min_dim=-1) const;

    std::size_t dimension() const;

    std::size_t num_simplices() const;

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

template<class OneSkeletonGraph>
void Fake_simplex_tree::insert_graph(const OneSkeletonGraph& skel_graph){
    toplex_maps.emplace(filtration_upper_bound,Toplex_map());
    bron_kerbosch_all_cliques(skel_graph, Visitor(&(this->toplex_maps.at(filtration_upper_bound))));
}

void Fake_simplex_tree::expansion(int max_dim){}

template <typename Input_vertex_range>
bool Fake_simplex_tree::all_facets_inside(const Input_vertex_range &vertex_range) const{
    Simplex sigma(vertex_range);
    for(const Simplex& s : facets(sigma))
        if(!membership(s)) return false;
    return true;
}

std::size_t Fake_simplex_tree::dimension() const {
    std::size_t max = 0;
    for(const Simplex& s : max_simplices())
        max = std::max(max, s.size());
    return max-1;
}

std::size_t Fake_simplex_tree::num_simplices() const {
    //return filtration_simplex_range().size();
    return max_simplices().size();
}

std::size_t Fake_simplex_tree::num_vertices() const {
    std::unordered_set<Vertex> vertices;
    for(const Simplex& s : max_simplices())
        for (Vertex v : s)
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
    for(const Simplex& s : max_simplices())
        if(s.size()<=d)
            simplices.emplace_back(s);
    return simplices;
}

std::vector<Simplex> Fake_simplex_tree::max_simplices() const{
    std::vector<Simplex> max_s;
    for(auto kv : toplex_maps)
        for(const Simplex_ptr& sptr : kv.second.maximal_cofaces(Simplex()))
            max_s.emplace_back(*sptr);
    return max_s;
}

std::size_t Fake_simplex_tree::dimension(Simplex_ptr& sptr) const{
    return sptr->size();
}

} //namespace Gudhi

#endif /* FAKE_SIMPLEX_TREE_H */

