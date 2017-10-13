#ifndef FILTERED_TOPLEX_MAP_H
#define FILTERED_TOPLEX_MAP_H

#include <gudhi/Toplex_map.h>
#include <limits>

#define filtration_upper_bound std::numeric_limits<Filtration_value>::max()

namespace Gudhi {

typedef double Filtration_value;

class Filtered_toplex_map {

public:
    template <typename Input_vertex_range>
    void insert_simplex_and_subfaces(const Input_vertex_range &vertex_range, Filtration_value f = filtration_upper_bound);

    template <typename Input_vertex_range>
    Filtration_value filtration(const Input_vertex_range &vertex_range) const;

protected:
    std::unordered_map<Filtration_value, Toplex_map> toplex_maps;
    std::unordered_map<Simplex_ptr, Filtration_value> filtrations;

};

template <typename Input_vertex_range>
void Filtered_toplex_map::insert_simplex_and_subfaces(const Input_vertex_range &vertex_range, Filtration_value f){
    if(!toplex_maps.count(f)) toplex_maps.emplace(f,Toplex_map());
    toplex_maps.at(f).insert_simplex(vertex_range);
    filtrations.emplace(get_key(vertex_range),f);
}

template <typename Input_vertex_range>
Filtration_value Filtered_toplex_map::filtration(const Input_vertex_range &vertex_range) const{
    for(auto kv : toplex_maps)
        if(kv.second.membership(vertex_range))
            return kv.first;
    return filtration_upper_bound;
}

} //namespace Gudhi

#endif /* FILTERED_TOPLEX_MAP_H */
