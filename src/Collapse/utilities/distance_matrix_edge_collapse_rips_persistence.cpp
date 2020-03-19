#include <gudhi/FlagComplexSpMatrix.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Rips_edge_list.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>
#include <gudhi/Points_off_io.h>

#include <CGAL/Epick_d.h>

#include <boost/program_options.hpp>

// Types definition
using Point = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>::Point_d;
using Vector_of_points = std::vector<Point>;

using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = double;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Rips_edge_list = Gudhi::rips_edge_list::Rips_edge_list<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Distance_matrix = std::vector<std::vector<Filtration_value>>;

void program_options(int argc, char* const argv[], double& min_persistence, double& end_thresold,
                     int& dimension, int& dim_max, std::string& in_file_name, std::string& out_file_name) {
  namespace po = boost::program_options;
  po::options_description visible("Allowed options", 100);
      visible.add_options()
        ("help,h", "produce help message")
      	("min_persistence,m", po::value<double>(&min_persistence)->default_value(0.1),
         "Minimum persistence interval length")
        ("end_thresold,e", po::value<double>(&end_thresold)->default_value(1),
         "Final threshold for rips complex.")
        ("dimensions,D", po::value<int>(&dimension)->default_value(2),
         "Dimension of the manifold.")
        ("dim_max,k ", po::value<int>(&dim_max)->default_value(2),
         "Maximum allowed dimension of the Rips complex.")
        ("input_file_name,i", po::value<std::string>(&in_file_name),
         "The input file.")
        ("out_file_name,o", po::value<std::string>(&out_file_name),
         "The output file.");

  po::options_description all;
  all.add(visible);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << std::endl;
    std::cout << "Computes rips complexes of different threshold values, to 'end_thresold' from a n random uniform "
                 "point_vector on a selected manifold, . \n";
    std::cout << "Strongly collapses all the rips complexes and output the results in out_file. \n";
    std::cout << "The experiments are repeted 'repete' num of times for each threshold value. \n";
    std::cout << "type -m for manifold options, 's' for uni sphere, 'b' for unit ball, 'f' for file. \n";
    std::cout << "type -i 'filename' for Input file option for exported point sample. \n";
    std::cout << std::endl << std::endl;
    std::cout << "Usage: " << argv[0] << " [options]" << std::endl << std::endl;
    std::cout << visible << std::endl;
    std::abort();
  }
}

class filt_edge_to_dist_matrix {
 public:
  template <class Distance_matrix, class Filtered_sorted_edge_list>
  filt_edge_to_dist_matrix(Distance_matrix& distance_mat, Filtered_sorted_edge_list& edge_filt,
                           std::size_t number_of_points) {
    double inf = std::numeric_limits<double>::max();
    doubleVector distances;
    std::pair<std::size_t, std::size_t> e;
    for (std::size_t indx = 0; indx < number_of_points; indx++) {
      for (std::size_t j = 0; j <= indx; j++) {
        if (j == indx)
          distances.push_back(0);
        else
          distances.push_back(inf);
      }
      distance_mat.push_back(distances);
      distances.clear();
    }

    for (auto edIt = edge_filt.begin(); edIt != edge_filt.end(); edIt++) {
      e = std::minmax(std::get<1>(*edIt), std::get<2>(*edIt));
      distance_mat.at(std::get<1>(e)).at(std::get<0>(e)) = std::get<0>(*edIt);
    }
  }
};

int main(int argc, char* const argv[]) {
  auto the_begin = std::chrono::high_resolution_clock::now();
  std::string out_file_name;
  std::string in_file_name;
  std::size_t number_of_points;

  typedef size_t Vertex_handle;
  typedef std::vector<std::tuple<Filtration_value, Vertex_handle, Vertex_handle>> Filtered_sorted_edge_list;

  int dimension;
  double end_threshold;
  double min_persistence;
  int dim_max = 2;

  program_options(argc, argv, min_persistence, end_threshold, dimension, dim_max, in_file_name,
                  out_file_name);

  std::cout << "The current input values to run the program is: " << std::endl;
  std::cout << "min_persistence, end_threshold, dimension, max_complex_dimension, in_file_name, out_file_name"
            << std::endl;
  std::cout << min_persistence << ", " << end_threshold << ", " << dimension << ", " << dim_max
            << ", " << in_file_name << ", " << out_file_name << std::endl;

  Map map_empty;

  Distance_matrix distances;
  Distance_matrix sparse_distances;

  distances = Gudhi::read_lower_triangular_matrix_from_csv_file<Filtration_value>(in_file_name);
  number_of_points = distances.size();
  std::cout << "Read the distance matrix succesfully, of size: " << number_of_points << std::endl;


  std::cout << "Successfully read " << number_of_points << " point_vector.\n";

  std::cout << "Point Set Generated." << std::endl;

  Filtered_sorted_edge_list edge_t;
  std::cout << "Computing the one-skeleton for threshold: " << end_threshold << std::endl;

  Rips_edge_list Rips_edge_list_from_file(distances, end_threshold);
  Rips_edge_list_from_file.create_edges(edge_t);
  std::cout<< "Sorted edge list computed" << std::endl;
  std::cout << "Total number of edges before collapse are: " << edge_t.size() << std::endl;

  if (edge_t.size() <= 0) {
    std::cerr << "Total number of egdes are zero." << std::endl;
    exit(-1);
  }

  // Now we will perform filtered edge collapse to sparsify the edge list edge_t.
  std::cout << "Filtered edge collapse begins" << std::endl;
  FlagComplexSpMatrix mat_filt_edge_coll(number_of_points, edge_t);
  std::cout << "Matrix instansiated" << std::endl;
  Filtered_sorted_edge_list collapse_edges;
  collapse_edges = mat_filt_edge_coll.filtered_edge_collapse();
  filt_edge_to_dist_matrix(sparse_distances, collapse_edges, number_of_points);
  std::cout << "Total number of vertices after collapse in the sparse matrix are: " << mat_filt_edge_coll.num_vertices()
            << std::endl;

  Rips_complex rips_complex_after_collapse(sparse_distances, end_threshold);

  Simplex_tree stree;
  rips_complex_after_collapse.create_complex(stree, dim_max);

  std::cout << "The complex contains " << stree.num_simplices() << " simplices  after collapse. \n";
  std::cout << "   and has dimension " << stree.dimension() << " \n";

  // Sort the simplices in the order of the filtration
  stree.initialize_filtration();
  // Compute the persistence diagram of the complex
  Persistent_cohomology pcoh(stree);
  // initializes the coefficient field for homology
  pcoh.init_coefficients(3);

  pcoh.compute_persistent_cohomology(min_persistence);
  if (out_file_name.empty()) {
    pcoh.output_diagram();
  } else {
    std::ofstream out(out_file_name);
    pcoh.output_diagram(out);
    out.close();
  }

  auto the_end = std::chrono::high_resolution_clock::now();

  std::cout << "Total computation time : " << std::chrono::duration<double, std::milli>(the_end - the_begin).count()
            << " ms\n"
            << std::endl;
  return 0;
}
