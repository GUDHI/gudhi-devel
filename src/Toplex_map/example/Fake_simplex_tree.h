#ifndef FAKE_SIMPLEX_TREE_H
#define FAKE_SIMPLEX_TREE_H

#include <cmath>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Filtered_toplex_map.h>
#include <gudhi/Lazy_Toplex_map.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>

namespace Gudhi {

struct Visitor {
    Lazy_Toplex_map* tm;

    Visitor(Lazy_Toplex_map* tm)
        :tm(tm)
    {}

    template <typename Clique, typename Graph>
    void clique(const Clique& c, const Graph& g)
    {
        tm->insert_simplex(c);
    }
};

/** Fake_simplex_tree is a wrapper for Filtered_toplex_map which has the interface of the Simplex_tree.
 * Mostly for retro-compatibility purpose. If you use a function that output non maximal simplices, it will be non efficient.
 * \ingroup toplex_map   */
class Fake_simplex_tree : public Filtered_toplex_map {

public:

    /** Handle type to a vertex contained in the simplicial complex.
     * \ingroup toplex_map   */
    typedef Toplex_map::Vertex Vertex_handle;

    /**  Handle type to a simplex contained in the simplicial complex.
     * \ingroup toplex_map   */
    typedef Toplex_map::Simplex Simplex_handle;

    typedef void Insertion_result_type;

    /**  Inserts the flag complex of a given range `Gudhi::rips_complex::Rips_complex::OneSkeletonGraph`
     * in the simplicial complex.
     * \ingroup toplex_map   */
    template<class OneSkeletonGraph>
    void insert_graph(const OneSkeletonGraph& skel_graph);

    /**  Do actually nothing.
     * \ingroup toplex_map   */
    void expansion(int max_dim);

    /**  Returns the number of vertices stored i.e. the number of max simplices
     *  \ingroup toplex_map   */
    std::size_t num_vertices() const;

    /**  Returns the dimension of the complex.
      * \ingroup toplex_map   */
    std::size_t dimension() const;

    /**  Returns the dimension of a given simplex in the complex.
      * \ingroup toplex_map   */
    std::size_t dimension(Simplex_ptr& sptr) const;

    /**  Returns the number of simplices stored i.e. the number of maximal simplices.
      * \ingroup toplex_map   */
    std::size_t num_simplices() const;

    /**  Returns a range over the vertices of a simplex.
      * \ingroup toplex_map   */
    Toplex_map::Simplex simplex_vertex_range(const Simplex& s) const;

    /**  Returns a set of all maximal (critical if there is filtration values) simplices.
      * \ingroup toplex_map   */
    std::vector<Toplex_map::Simplex> max_simplices() const;

    /** Returns all the simplices, of max dimension d if a parameter d is given.
      * \ingroup toplex_map   */
    std::vector<Toplex_map::Simplex> filtration_simplex_range(int d=std::numeric_limits<int>::max()) const;

    /** Returns all the simplices of max dimension d
      * \ingroup toplex_map   */
    std::vector<Toplex_map::Simplex> skeleton_simplex_range(int d) const;

    Toplex_map::Vertex contraction(const Toplex_map::Vertex x, const Toplex_map::Vertex y);


protected:

    /** \internal Does all the facets of the given simplex belong to the complex ?
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    bool all_facets_inside(const Input_vertex_range &vertex_range) const;

};

template<class OneSkeletonGraph>
void Fake_simplex_tree::insert_graph(const OneSkeletonGraph& skel_graph){
    toplex_maps.emplace(nan(""), new Lazy_Toplex_map());
    using vertex_iterator = typename boost::graph_traits<OneSkeletonGraph>::vertex_iterator;
    vertex_iterator vi, vi_end;
    for (std::tie(vi, vi_end) = boost::vertices(skel_graph); vi != vi_end; ++vi) {
        Simplex s; s.insert(*vi);
        insert_simplex_and_subfaces(s);
    }
    bron_kerbosch_all_cliques(skel_graph, Visitor(this->toplex_maps.at(nan(""))));
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

std::size_t Fake_simplex_tree::dimension(Simplex_ptr& sptr) const{
    return sptr->size();
}

std::size_t Fake_simplex_tree::num_simplices() const {
    return max_simplices().size();
}

std::size_t Fake_simplex_tree::num_vertices() const {
    std::unordered_set<Toplex_map::Vertex> vertices;
    for(const Toplex_map::Simplex& s : max_simplices())
        for (Toplex_map::Vertex v : s)
            vertices.emplace(v);
    return vertices.size();
}

Toplex_map::Simplex Fake_simplex_tree::simplex_vertex_range(const Simplex& s) const {
    return s;
}

std::vector<Toplex_map::Simplex> Fake_simplex_tree::max_simplices() const{
    std::vector<Toplex_map::Simplex> max_s;
    for(auto kv : toplex_maps)
        for(const Toplex_map::Simplex_ptr& sptr : kv.second->maximal_cofaces(Simplex()))
            max_s.emplace_back(*sptr);
    return max_s;
}

std::vector<Toplex_map::Simplex> Fake_simplex_tree::filtration_simplex_range(int d) const{
    std::vector<Toplex_map::Simplex> m = max_simplices();
    std::vector<Toplex_map::Simplex> range;
    Toplex_map::Simplex_ptr_set seen;
    while(m.begin()!=m.end()){
        Toplex_map::Simplex s(m.back());
        m.pop_back();
        if(seen.find(get_key(s))==seen.end()){
            if((int) s.size()-1 <=d)
                range.emplace_back(s);
            seen.emplace(get_key(s));
            if(s.size()>0)
                for(Simplex& sigma : facets(s))
                    m.emplace_back(sigma);
        }
    }
    return range;
}

std::vector<Toplex_map::Simplex> Fake_simplex_tree::skeleton_simplex_range(int d) const{
    return filtration_simplex_range(d);
}

Toplex_map::Vertex Fake_simplex_tree::contraction(const Toplex_map::Vertex x, const Toplex_map::Vertex y){
    for(auto kv : toplex_maps)
        kv.second->contraction(x,y);
    return y;
}

} //namespace Gudhi

#endif /* FAKE_SIMPLEX_TREE_H */

