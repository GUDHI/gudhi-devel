#ifndef FILTERED_TOPLEX_MAP_H
#define FILTERED_TOPLEX_MAP_H

#include <gudhi/Toplex_map.h>
#include <map>
#include <limits>

namespace Gudhi {

class Filtered_toplex_map {

public:
    typedef double Filtration_value;

    template <typename Input_vertex_range>
    std::pair<Simplex, bool> insert_simplex_and_subfaces(const Input_vertex_range &vertex_range, Filtration_value f = nan(""));

    template <typename Input_vertex_range>
    Filtration_value filtration(const Input_vertex_range &vertex_range) const;

    template <typename Input_vertex_range>
    bool membership(const Input_vertex_range &vertex_range) const;

protected:
    std::map<Filtration_value, Toplex_map> toplex_maps;
};

template <typename Input_vertex_range>
std::pair<Simplex, bool> Filtered_toplex_map::insert_simplex_and_subfaces(const Input_vertex_range &vertex_range, Filtration_value f){
    Simplex s(vertex_range.begin(),vertex_range.end());
    if(membership(s)) return make_pair(s,false);
    if(!toplex_maps.count(f)) toplex_maps.emplace(f,Toplex_map());
    toplex_maps.at(f).insert_simplex(vertex_range);
    return make_pair(s,true);
}


template <typename Input_vertex_range>
Filtered_toplex_map::Filtration_value Filtered_toplex_map::filtration(const Input_vertex_range &vertex_range) const{
    for(auto kv : toplex_maps)
        if(kv.second.membership(vertex_range))
            return kv.first; //min only because a map is ordered
    return nan("");
}

template <typename Input_vertex_range>
bool Filtered_toplex_map::membership(const Input_vertex_range &vertex_range) const{
    for(auto kv : toplex_maps)
        if(kv.second.membership(vertex_range))
            return true;
    return false;
}

} //namespace Gudhi

#endif /* FILTERED_TOPLEX_MAP_H */
