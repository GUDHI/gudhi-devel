#include <gudhi/Skeleton_blocker.h>

#ifndef SKELETON_BLOCKER_WRAPPER_H_
#define SKELETON_BLOCKER_WRAPPER_H_

namespace Gudhi {

class Sb_wrapper{

public:

    typedef Gudhi::skeleton_blocker::Skeleton_blocker_simple_traits Traits;
    typedef Gudhi::skeleton_blocker::Skeleton_blocker_complex<Traits> Complex;

    typedef Complex::Vertex_handle Vertex_handle;

    typedef Complex::Simplex Simplex;

    typedef Toplex_map::Vertex Vertex;

    typedef Toplex_map::Simplex_ptr Simplex_ptr;

    typedef Toplex_map::Simplex_ptr_set Simplex_ptr_set;

    typedef double Filtration_value;

    template <typename Input_vertex_range>
    std::pair<bool, bool> insert_simplex_and_subfaces(const Input_vertex_range &vertex_range, double);

    template <typename Input_vertex_range>
    bool membership(const Input_vertex_range &vertex_range) const;

    typedef Toplex_map::Simplex Simplex_handle;

    typedef void Insertion_result_type;

    /**  Inserts the flag complex of a given range `Gudhi::rips_complex::Rips_complex::OneSkeletonGraph`
     * in the simplicial complex.  */
    template<class OneSkeletonGraph>
    void insert_graph(const OneSkeletonGraph& skel_graph);

    /**  Do actually nothing.  */
    void expansion(int max_dim);

    /**  Returns the number of vertices stored i.e. the number of max simplices  */
    std::size_t num_vertices() const;

    /**  Returns the dimension of the complex.  */
    std::size_t dimension() const;

    /**  Returns the dimension of a given simplex in the complex.  */
    std::size_t dimension(Simplex_ptr& sptr) const;

    /**  Returns the number of simplices stored i.e. the number of maximal simplices. */
    std::size_t num_simplices() const;

    /**  Returns a range over the vertices of a simplex. */
    Toplex_map::Simplex simplex_vertex_range(const Simplex& s) const;

    /**  Returns a set of all maximal (critical if there is filtration values) simplices.  */
    std::vector<Toplex_map::Simplex> max_simplices() const;

    /** Returns all the simplices, of max dimension d if a parameter d is given.  */
    std::vector<Toplex_map::Simplex> filtration_simplex_range(int d=std::numeric_limits<int>::max()) const;

    /** Returns all the simplices of max dimension d */
    std::vector<Toplex_map::Simplex> skeleton_simplex_range(int d) const;

private:

    Complex sb;

};


template<class OneSkeletonGraph>
void Sb_wrapper::insert_graph(const OneSkeletonGraph& skel_graph){
    using vertex_iterator = typename boost::graph_traits<OneSkeletonGraph>::vertex_iterator;
    vertex_iterator vi, vi_end;
    // for (std::tie(vi, vi_end) = boost::vertices(skel_graph); vi != vi_end; ++vi)
    //   insert_vertex(*vi);
    //edges
}

void Sb_wrapper::expansion(int max_dim){}

template <typename Input_vertex_range>
std::pair<bool, bool> Sb_wrapper::insert_simplex_and_subfaces(const Input_vertex_range &vertex_range, double){
    Complex::Simplex s;
    for (auto v : vertex_range)
        s.add_vertex(Vertex_handle(v));
    if(sb.contains(s))
        return std::make_pair(false,false);
    sb.add_simplex(s);
    return std::make_pair(true,true);
}

std::size_t Sb_wrapper::num_vertices() const{
    std::size_t num_vertices = 0;
    for(auto v : sb.vertex_range())
        ++num_vertices;
    return num_vertices;
}



}

#endif
