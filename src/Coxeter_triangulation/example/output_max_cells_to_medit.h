#ifndef OUTPUT_MAX_CELLS_TO_MEDIT_H_
#define OUTPUT_MAX_CELLS_TO_MEDIT_H_

#include <string>
#include <unordered_map>
// #include <algorithm>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>

#include "output_hasse_to_medit.h"

Hasse_cell* insert_hasse_subdiagram(const Cell_id& c_id,
				    Hasse_diagram& hd,
				    VP_map& vp_map,
				    VC_map& dictionary,
				    const Gudhi::Coxeter_triangulation_ds& ct) {
  auto res_pair = dictionary.emplace(std::make_pair(c_id, new Hasse_cell((int)(ct.dimension() - c_id.dimension()), 0.0)));
  if (res_pair.second) {
    Hasse_cell* new_cell = res_pair.first->second;
    if (new_cell->get_dimension() == 0)
      vp_map.emplace(std::make_pair(new_cell, ct.barycenter(c_id)));
    hd.emplace(new_cell);
    Hasse_boundary& boundary = new_cell->get_boundary();
    Coface_it cof_it(c_id, ct, c_id.dimension() + 1);
    Coface_it cof_end;
    for (; cof_it != cof_end; ++cof_it) {
      Hasse_cell* facet_cell = insert_hasse_subdiagram(*cof_it, hd, vp_map, dictionary, ct);
      if (facet_cell != 0)
	if (std::find(boundary.begin(), boundary.end(), std::make_pair(facet_cell, 1)) == boundary.end())
	  boundary.push_back(std::make_pair(facet_cell, 1));
    }
  }
  return res_pair.first->second;
}


template <class CellSet,
	  class Coxeter_triangulation_ds,
	  class OrderMap>
void output_max_cells_to_medit(const CellSet& max_cells,
			       const Coxeter_triangulation_ds& ct,
			       const OrderMap& order_map,
			       std::string file_name = "reconstruction")
{
  
  // Compute the Hasse diagram
  Hasse_diagram hd;
  VC_map dictionary;
  VP_map vp_map;
  for (auto c: max_cells) {
    insert_hasse_subdiagram(c, hd, vp_map, dictionary, ct);
  }  
  output_hasse_to_medit(hd, vp_map, file_name);
  hd.clear();
  
}


#endif
