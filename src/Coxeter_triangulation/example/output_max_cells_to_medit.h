#ifndef OUTPUT_MAX_CELLS_TO_MEDIT_H_
#define OUTPUT_MAX_CELLS_TO_MEDIT_H_

#include <string>
#include <unordered_map>
// #include <algorithm>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>

#include "output_hasse_to_medit.h"

void insert_hasse_subdiagram(const Cell_id& c_id,
			     Hasse_diagram& hd,
			     VP_map& vp_map,
			     VC_map& dictionary,
			     std::vector<Hasse_cell*>& cells,
			     const Gudhi::Coxeter_triangulation_ds& ct) {
  std::cout << "  Insert_hasse_subdiagram for " << c_id << "\n";
  if (dictionary.find(c_id) == dictionary.end()) {
    cells.push_back(new Hasse_cell((int)(ct.dimension() - c_id.dimension()), 0.0));
    Hasse_cell* new_cell = cells.back();
    dictionary.emplace(std::make_pair(c_id, new_cell));
    if (new_cell->get_dimension() == 0)
      vp_map.emplace(std::make_pair(new_cell, ct.barycenter(c_id)));
    if (!hd.emplace(new_cell).second)
      std::cout << "Cell " << new_cell << " Problem!\n";
    else
      std::cout << "Cell " << new_cell << " inserted. dim = " << new_cell->get_dimension() << "\n";
    Hasse_boundary& boundary = new_cell->get_boundary();
    Coface_it cof_it(c_id, ct, c_id.dimension() + 1);
    Coface_it cof_end;
    for (; cof_it != cof_end; ++cof_it) {
      insert_hasse_subdiagram(*cof_it, hd, vp_map, dictionary, cells, ct);
      boundary.push_back(std::make_pair(dictionary.at(*cof_it), 1));
    }
    std::cout << "Boundary of the cell" << c_id << " = ";
    for (auto b: boundary)
      std::cout << b.first << " ";
    std::cout  << "\n";
  }
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
  std::vector<Hasse_cell*> cells;
  VP_map vp_map;
  for (auto c: max_cells)
    insert_hasse_subdiagram(c, hd, vp_map, dictionary, cells, ct);
  output_hasse_to_medit(hd, vp_map, false, file_name+"_no_barycentric");
  output_hasse_to_medit(cells, vp_map, true, file_name+"_barycentric");
  cells.clear();
  
}


#endif
