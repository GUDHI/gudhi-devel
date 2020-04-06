#include <gudhi/Flag_complex_sparse_matrix.h>
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
                     int& dimension, int& dim_max, std::string& csv_matrix_file, std::string& filediag) {
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
        ("input_file_name,i", po::value<std::string>(&csv_matrix_file),
         "The input file.")
        ("filediag,o", po::value<std::string>(&filediag),
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

void program_options(int argc, char* argv[], std::string& csv_matrix_file, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence);

int main(int argc, char* argv[]) {
  auto the_begin = std::chrono::high_resolution_clock::now();

  typedef size_t Vertex_handle;
  typedef std::vector<std::tuple<Filtration_value, Vertex_handle, Vertex_handle>> Filtered_sorted_edge_list;

  std::string csv_matrix_file;
  std::string filediag;
  double threshold;
  int dim_max = 2;
  int p;
  double min_persistence;

  program_options(argc, argv, csv_matrix_file, filediag, threshold, dim_max, p, min_persistence);

  Map map_empty;

  Distance_matrix distances;
  Distance_matrix sparse_distances;

  distances = Gudhi::read_lower_triangular_matrix_from_csv_file<Filtration_value>(csv_matrix_file);
  std::size_t number_of_points = distances.size();
  std::cout << "Read the distance matrix succesfully, of size: " << number_of_points << std::endl;

  Filtered_sorted_edge_list edge_t;
  std::cout << "Computing the one-skeleton for threshold: " << threshold << std::endl;

  Rips_edge_list Rips_edge_list_from_file(distances, threshold);
  Rips_edge_list_from_file.create_edges(edge_t);
  std::cout<< "Sorted edge list computed" << std::endl;
  
  if (edge_t.size() <= 0) {
    std::cerr << "Total number of egdes are zero." << std::endl;
    exit(-1);
  }

  std::cout << "Total number of edges before collapse are: " << edge_t.size() << std::endl;

  // Now we will perform filtered edge collapse to sparsify the edge list edge_t.
  std::cout << "Filtered edge collapse begins" << std::endl;
  Flag_complex_sparse_matrix mat_filt_edge_coll(edge_t);
  std::cout << "Matrix instansiated" << std::endl;
  Filtered_sorted_edge_list collapse_edges;
  collapse_edges = mat_filt_edge_coll.filtered_edge_collapse();
  filt_edge_to_dist_matrix(sparse_distances, collapse_edges, number_of_points);
  std::cout << "Total number of vertices after collapse in the sparse matrix are: " << mat_filt_edge_coll.num_vertices()
            << std::endl;

  Rips_complex rips_complex_after_collapse(sparse_distances, threshold);

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
  if (filediag.empty()) {
    pcoh.output_diagram();
  } else {
    std::ofstream out(filediag);
    pcoh.output_diagram(out);
    out.close();
  }

  auto the_end = std::chrono::high_resolution_clock::now();

  std::cout << "Total computation time : " << std::chrono::duration<double, std::milli>(the_end - the_begin).count()
            << " ms\n"
            << std::endl;
  return 0;
}

void program_options(int argc, char* argv[], std::string& csv_matrix_file, std::string& filediag,
                     Filtration_value& threshold, int& dim_max, int& p, Filtration_value& min_persistence) {
  namespace po = boost::program_options;
  po::options_description hidden("Hidden options");
  hidden.add_options()(
      "input-file", po::value<std::string>(&csv_matrix_file),
      "Name of file containing a distance matrix. Can be square or lower triangular matrix. Separator is ';'.");

  po::options_description visible("Allowed options", 100);
  visible.add_options()("help,h", "produce help message")(
      "output-file,o", po::value<std::string>(&filediag)->default_value(std::string()),
      "Name of file in which the persistence diagram is written. Default print in std::cout")(
      "max-edge-length,r",
      po::value<Filtration_value>(&threshold)->default_value(std::numeric_limits<Filtration_value>::infinity()),
      "Maximal length of an edge for the Rips complex construction.")(
      "cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),
      "Maximal dimension of the Rips complex we want to compute.")(
      "field-charac,p", po::value<int>(&p)->default_value(11),
      "Characteristic p of the coefficient field Z/pZ for computing homology.")(
      "min-persistence,m", po::value<Filtration_value>(&min_persistence),
      "Minimal lifetime of homology feature to be recorded. Default is 0. Enter a negative value to see zero length "
      "intervals");

  po::positional_options_description pos;
  pos.add("input-file", 1);

  po::options_description all;
  all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file")) {
    std::cout << std::endl;
    std::cout << "Compute the persistent homology with coefficient field Z/pZ \n";
    std::cout << "of a Rips complex after edge collapse defined on a set of distance matrix.\n \n";
    std::cout << "The output diagram contains one bar per line, written with the convention: \n";
    std::cout << "   p   dim b d \n";
    std::cout << "where dim is the dimension of the homological feature,\n";
    std::cout << "b and d are respectively the birth and death of the feature and \n";
    std::cout << "p is the characteristic of the field Z/pZ used for homology coefficients." << std::endl << std::endl;

    std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
    std::cout << visible << std::endl;
    exit(-1);
  }
}
