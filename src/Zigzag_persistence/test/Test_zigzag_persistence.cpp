

#include <chrono>

#include "io.h"
#include "graph_simplicial_complex.h"
#include "distance_functions.h"
#include "Simplex_tree.h"
#include "Hasse_complex.h"

#include "Zigzag_persistence.h"

int main() 
{
  std::chrono::time_point<std::chrono::system_clock> start, end;
  int enlapsed_sec;

  std::string      filepoints      = "../examples/Kl.txt";
  Filtration_value threshold       = 0.2;
  int              dim_max         = 3;

// Extract the points from the file filepoints
  typedef std::vector<double> Point_t;
  std::vector< Point_t > points;
  read_points( filepoints, points );

// Compute the proximity graph of the points
  start = std::chrono::system_clock::now();
  Graph_t prox_graph = compute_proximity_graph( points, threshold
                                              , euclidean_distance<Point_t> );
  end = std::chrono::system_clock::now();
  enlapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  std::cout << "Compute Rips graph in " << enlapsed_sec << " sec.\n";

// Construct the Rips complex in a Simplex Tree
  Simplex_tree<> st;        
  start = std::chrono::system_clock::now();

  st.insert_graph(prox_graph); // insert the proximity graph in the simplex tree
  st.expansion( dim_max ); // expand the graph until dimension dim_max

  end = std::chrono::system_clock::now();
  enlapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  std::cout << "Compute Rips complex in " << enlapsed_sec << " sec.\n";
  std::cout << "  - dimension           = " << st.dimension() << std::endl;
  std::cout << "  - number of simplices = " << st.num_simplices() << std::endl;

// Sort the simplices in the order of the filtration
  start = std::chrono::system_clock::now();
  st.initialize_filtration();
  end = std::chrono::system_clock::now();
  enlapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  std::cout << "Order the simplices of the filtration in " << enlapsed_sec << " sec.\n";




  Zigzag_persistence< Simplex_tree<> > zz_p (st);
  start = std::chrono::system_clock::now();
  zz_p.compute_zigzag_persistence();
  end = std::chrono::system_clock::now();
  enlapsed_sec = std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  std::cout << "  Compute zigzag persistence in " << enlapsed_sec << " sec.\n";

  return 0;
}
