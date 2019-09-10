/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Rips_complex.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Persistent_cohomology/Multi_field.h>
#include <gudhi/Hasse_complex.h>
#include <gudhi/Points_off_io.h>

#include <chrono>
#include <string>
#include <vector>

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = Simplex_tree::Filtration_value;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Multi_field = Gudhi::persistent_cohomology::Multi_field;
using Point = std::vector<double>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;

/* Compute the persistent homology of the complex cpx with coefficients in Z/pZ. */
template< typename FilteredComplex>
void timing_persistence(FilteredComplex & cpx
                        , int p);

/* Compute multi-field persistent homology of the complex cpx with coefficients in
 * Z/rZ for all prime number r in [p;q].*/
template< typename FilteredComplex>
void timing_persistence(FilteredComplex & cpx
                        , int p
                        , int q);

/* Timings for the computation of persistent homology with different 
 * representations of a Rips complex and different coefficient fields. The 
 * Rips complex is built on a set of 10000 points sampling a Klein bottle embedded 
 * in dimension 5.
 * We represent complexes with a simplex tree and 
 * with a Hasse diagram. The Hasse diagram represents explicitly all 
 * codimension 1 incidence relations in the complex, and hence leads to 
 * a faster computation of persistence because boundaries are precomputed. 
 * Hovewer, the simplex tree may be constructed directly from a point cloud and
 * is more compact.
 * We compute persistent homology with coefficient fields Z/2Z and Z/1223Z.
 * We present also timings for the computation of multi-field persistent 
 * homology in all fields Z/rZ for r prime between 2 and 1223.
 */
int main(int argc, char * argv[]) {
  std::chrono::time_point<std::chrono::system_clock> start, end;
  int elapsed_sec;
  {

  std::string off_file_points = "Kl.off";
  Filtration_value threshold = 0.27;
  int dim_max = 3;
  int p = 2;
  int q = 1223;

  // Extract the points from the file off_file_points
  Points_off_reader off_reader(off_file_points);

  // Compute the proximity graph of the points
  start = std::chrono::system_clock::now();
  Rips_complex rips_complex_from_file(off_reader.get_point_cloud(), threshold, Gudhi::Euclidean_distance());
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "Compute Rips graph in " << elapsed_sec << " ms.\n";

  // Construct the Rips complex in a Simplex Tree
  Simplex_tree st;
  start = std::chrono::system_clock::now();

  // insert the proximity graph in the simplex tree
  // expand the graph until dimension dim_max
  rips_complex_from_file.create_complex(st, dim_max);

  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "Compute Rips complex in " << elapsed_sec << " ms.\n";
  std::cout << "  - dimension           = " << st.dimension() << std::endl;
  std::cout << "  - number of simplices = " << st.num_simplices() << std::endl;

  // Sort the simplices in the order of the filtration
  start = std::chrono::system_clock::now();
  st.initialize_filtration();
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "Order the simplices of the filtration in " << elapsed_sec << " ms.\n";

  // Copy the keys inside the simplices
  start = std::chrono::system_clock::now();
  {
    int count = 0;
    for (auto sh : st.filtration_simplex_range())
      st.assign_key(sh, count++);
  }
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "Copied the keys inside the simplices in " << elapsed_sec << " ms.\n";

  // Convert the simplex tree into a hasse diagram
  start = std::chrono::system_clock::now();
  Gudhi::Hasse_complex<> hcpx(st);
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "Convert the simplex tree into a Hasse diagram in " << elapsed_sec << " ms.\n";


  std::cout << "Timings when using a simplex tree: \n";
  timing_persistence(st, p);
  timing_persistence(st, q);
  timing_persistence(st, p, q);

  std::cout << "Timings when using a Hasse complex: \n";
  timing_persistence(hcpx, p);
  timing_persistence(hcpx, q);
  timing_persistence(hcpx, p, q);

  start = std::chrono::system_clock::now();
  }
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "Running the complex destructors in " << elapsed_sec << " ms.\n";
  return 0;
}

template< typename FilteredComplex>
void
timing_persistence(FilteredComplex & cpx
                   , int p) {
  std::chrono::time_point<std::chrono::system_clock> start, end;
  int elapsed_sec;
  {
  start = std::chrono::system_clock::now();
  Gudhi::persistent_cohomology::Persistent_cohomology< FilteredComplex, Field_Zp > pcoh(cpx);
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "  Initialize pcoh in " << elapsed_sec << " ms.\n";
  // initializes the coefficient field for homology
  start = std::chrono::system_clock::now();
  pcoh.init_coefficients(p);
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "  Initialize the coefficient field in " << elapsed_sec << " ms.\n";

  start = std::chrono::system_clock::now();

  pcoh.compute_persistent_cohomology(INFINITY);

  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "  Compute persistent homology in Z/" << p << "Z in " << elapsed_sec << " ms.\n";
  start = std::chrono::system_clock::now();
  }
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "  Run the persistence destructors in " << elapsed_sec << " ms.\n";
}

template< typename FilteredComplex>
void
timing_persistence(FilteredComplex & cpx
                   , int p
                   , int q) {
  std::chrono::time_point<std::chrono::system_clock> start, end;
  int elapsed_sec;
  {
  start = std::chrono::system_clock::now();
  Gudhi::persistent_cohomology::Persistent_cohomology< FilteredComplex, Multi_field > pcoh(cpx);
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "  Initialize pcoh in " << elapsed_sec << " ms.\n";
  // initializes the coefficient field for homology
  start = std::chrono::system_clock::now();
  pcoh.init_coefficients(p, q);
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "  Initialize the coefficient field in " << elapsed_sec << " ms.\n";
  // compute persistent homology, disgarding persistent features of life shorter than min_persistence

  start = std::chrono::system_clock::now();

  pcoh.compute_persistent_cohomology(INFINITY);

  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "  Compute multi-field persistent homology in all coefficient fields Z/pZ "
      << "with p in [" << p << ";" << q << "] in " << elapsed_sec << " ms.\n";
  start = std::chrono::system_clock::now();
  }
  end = std::chrono::system_clock::now();
  elapsed_sec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "  Run the persistence destructors in " << elapsed_sec << " ms.\n";
}
