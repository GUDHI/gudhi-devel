#ifndef OUTPUT_MAX_CELLS_TO_MEDIT_H_
#define OUTPUT_MAX_CELLS_TO_MEDIT_H_

#include <string>
#include <unordered_map>
// #include <algorithm>

#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>

#include "output_hasse_to_medit.h"

std::size_t filt_index = 1;

void insert_hasse_subdiagram(const Cell_id& c_id,
			     VP_map& vp_map,
			     VC_map& dictionary,
			     std::vector<Hasse_cell*>& cells,
			     const Gudhi::Coxeter_triangulation_ds& ct) {
  // std::cout << "  Insert_hasse_subdiagram for " << c_id << "\n";
  if (dictionary.find(c_id) == dictionary.end()) {
    cells.push_back(new Hasse_cell((int)(ct.dimension() - c_id.dimension()), filt_index++));
    Hasse_cell* new_cell = cells.back();
    dictionary.emplace(std::make_pair(c_id, new_cell));
    if (new_cell->get_dimension() == 0)
      vp_map.emplace(std::make_pair(new_cell, ct.barycenter(c_id)));
    Hasse_boundary& boundary = new_cell->get_boundary();
    Coface_it cof_it(c_id, ct, c_id.dimension() + 1);
    Coface_it cof_end;
    for (; cof_it != cof_end; ++cof_it) {
      insert_hasse_subdiagram(*cof_it, vp_map, dictionary, cells, ct);
      boundary.push_back(std::make_pair(dictionary.at(*cof_it), 1));
    }
    // std::cout << "Boundary of the cell" << c_id << " = ";
    // for (auto b: boundary)
    //   std::cout << b.first << " ";
    // std::cout  << "\n";
  }
}

Hasse_cell* insert_hasse_subdiagram_as(Hasse_cell* cell,
				       VC_map& dictionary,
				       std::map<Hasse_cell*, Hasse_cell*>& dictionary2,
				       std::vector<Hasse_cell*>& cells,
				       std::size_t cod_d,
				       std::size_t amb_d) {
  // std::cout << "  Insert_hasse_subdiagram_as for cell = " << cell << ".\n";
  auto d_it = dictionary2.find(cell);
  if (d_it == dictionary2.end()) {
    std::size_t c_dim = cell->get_dimension();
    cells.push_back(new Hasse_cell((int)(amb_d - cod_d - c_dim), cell->get_filtration()));
    Hasse_cell* new_cell = cells.back();
    dictionary2.emplace(std::make_pair(cell, new_cell));
    for (auto f_pair: cell->get_boundary()) {
      insert_hasse_subdiagram_as(f_pair.first, dictionary, dictionary2, cells, cod_d, amb_d);
      dictionary2.at(f_pair.first)->get_boundary().push_back(std::make_pair(new_cell, 1));
    }
    return new_cell;
  }
  return d_it->second;
}

template <class CellSet,
	  class Coxeter_triangulation_ds,
	  class OrderMap>
void output_max_cells_to_medit(const CellSet& max_cells,
			       const Coxeter_triangulation_ds& ct,
			       const OrderMap& order_map,
			       std::size_t cod_d,
			       std::string file_name = "reconstruction")
{
  
  // Compute the Hasse diagram
  VC_map dictionary;
  std::map<Hasse_cell*, Hasse_cell*> dictionary2;
  std::vector<Hasse_cell*> cells, cells_as;
  VP_map vp_map, vp_map2;
  for (auto c: max_cells)
    insert_hasse_subdiagram(c, vp_map, dictionary, cells, ct);
  // dictionary.clear();
  std::size_t amb_d = ct.dimension();
  for (auto c: max_cells) {
    Hasse_cell* cell = dictionary.at(c);
    auto cell2 = insert_hasse_subdiagram_as(cell, dictionary, dictionary2, cells_as, cod_d, amb_d);
    if (c.dimension() == cod_d)
      vp_map2.emplace(std::make_pair(cell2, ct.barycenter(c)));
  }
  // for (auto d_pair: dictionary2) {
  //   std::cout << d_pair.first << ": ";
  //   for (auto b_pair: d_pair.first->get_boundary())
  //     std::cout << b_pair.first << " ";
  //   std::cout << "\n";
  //   std::cout << d_pair.second << " (dim=" << d_pair.second->get_dimension() << "): "; 
  //   for (auto b_pair: d_pair.second->get_boundary())
  //     std::cout << b_pair.first << " ";
  //   std::cout << "\n\n";
  // }
  // output_hasse_to_medit(cells, vp_map, false, file_name+"_no_barycentric");
  output_hasse_to_medit(cells, vp_map, true, file_name+"_barycentric");
  allgow_switch = true;
  output_hasse_to_medit(cells_as, vp_map2, true, file_name+"_allgowerschmidt");
  for (auto c: cells)
    delete c;
  for (auto c: cells_as)
    delete c;  
}


#endif
