    #ifndef FILTERED_TOPLEX_MAP_H
#define FILTERED_TOPLEX_MAP_H

#include <gudhi/Toplex_map.h>
#include <map>
#include <limits>

namespace Gudhi {

/** A Filtered_toplex_map represents the simplicial complex with a filtration.
 * A "toplex" is a critical simplex.
 * \ingroup toplex_map   */
class Filtered_toplex_map {

public:
    /** Vertex is the type of vertices.
     * \ingroup toplex_map   */
    typedef Toplex_map::Vertex Vertex;

    /** Simplex is the type of simplices.
     * \ingroup toplex_map   */
    typedef Toplex_map::Simplex Simplex;

    /** The type of the pointers to maximal simplices.
     * \ingroup toplex_map   */
    typedef Toplex_map::Simplex_ptr Simplex_ptr;

    /** The type of the sets of Simplex_ptr.
     * \ingroup toplex_map   */
    typedef Toplex_map::Simplex_ptr_set Simplex_ptr_set;

    /** The type of the filtration values.
     * \ingroup toplex_map   */
    typedef double Filtration_value;

    /** Add a simplex and its subfaces with the given filtration value
     * in the Filtered_toplex_map.
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    std::pair<Simplex, bool> insert_simplex_and_subfaces(const Input_vertex_range &vertex_range, Filtration_value f = std::numeric_limits<double>::quiet_NaN());

    /** Gives the filtration of the input simplex.
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    Filtration_value filtration(const Input_vertex_range &vertex_range) const;

    /** Is the input simplex member of the complex ?
     * \ingroup toplex_map   */
    template <typename Input_vertex_range>
    bool membership(const Input_vertex_range &vertex_range) const;

protected:
    std::map<Filtration_value, Toplex_map*> toplex_maps;
};

template <typename Input_vertex_range>
std::pair<Toplex_map::Simplex, bool> Filtered_toplex_map::insert_simplex_and_subfaces(const Input_vertex_range &vertex_range, Filtration_value f){
    Simplex s(vertex_range.begin(),vertex_range.end());
    if(membership(s)) return make_pair(s,false);
    if(!toplex_maps.count(f)) toplex_maps.emplace(f,new Toplex_map());
    toplex_maps.at(f)->insert_simplex(vertex_range);
    return make_pair(s,true);
}


template <typename Input_vertex_range>
Filtered_toplex_map::Filtration_value Filtered_toplex_map::filtration(const Input_vertex_range &vertex_range) const{
    for(auto kv : toplex_maps)
        if(kv.second->membership(vertex_range))
            return kv.first; //min only because a map is ordered
    return std::numeric_limits<double>::quiet_NaN() ;
}

template <typename Input_vertex_range>
bool Filtered_toplex_map::membership(const Input_vertex_range &vertex_range) const{
    for(auto kv : toplex_maps)
        if(kv.second->membership(vertex_range))
            return true;
    return false;
}

} //namespace Gudhi

#endif /* FILTERED_TOPLEX_MAP_H */
