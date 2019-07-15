/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CELL_COMPLEX_H_
#define CELL_COMPLEX_H_

#include <Eigen/Dense>

#include <gudhi/Permutahedral_representation/Simplex_comparator.h>
#include <fstream> // for Hasse_diagram_persistence.h
#include <gudhi/Hasse_diagram_persistence.h> // for Hasse_cel

namespace Gudhi {

namespace coxeter_triangulation {


/** \class Cell_complex 
 *  \brief A class that assembles methods for manifold tracing algorithm.
 *
 *  \tparam Out_simplex_map_ The type of a map from a simplex type that is a
 *   model of SimplexInCoxeterTriangulation to Eigen::VectorXd.
 *
 *  \ingroup coxeter_triangulation
 */
template <class Out_simplex_map_>
class Cell_complex {

public:  
  typedef typename Out_simplex_map_::key_type Simplex_handle;
  using Hasse_cell = Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, double>;
  typedef std::map<Simplex_handle,
		   Hasse_cell*,
		   Simplex_comparator<Simplex_handle> > Simplex_cell_map;
  typedef std::vector<Simplex_cell_map> Simplex_cell_maps;
  
  Cell_complex() {}

private:
  Simplex_cell_maps simplex_cell_maps;
};

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
