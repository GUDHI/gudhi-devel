/******************************************************************************
This benchmark allows to compute the Tangential Complex from input files or 
generated point sets.

It reads the benchmark_script.txt file (located in the same folder as this 
file) and compute one or several complexes for each line. Unless TC_NO_EXPORT 
is defined, each complex is exported as an OFF file and/or as a RIB file 
(RenderMan). In addition an XML file is created at each run of the benchmark. 
It contains statistics about the complexes that were created. This XML file 
can be processed in Excel, for example.
 ******************************************************************************/

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
#define TBB_USE_THREADING_TOOL
#endif

#include <cstddef>

//#define GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
//#define TC_INPUT_STRIDES 3 // only take one point every TC_INPUT_STRIDES points
#define TC_NO_EXPORT // do not output OFF files
//#define TC_EXPORT_TO_RIB // 
//#define GUDHI_TC_EXPORT_SPARSIFIED_POINT_SET
//#define GUDHI_TC_EXPORT_ALL_COORDS_IN_OFF

const std::size_t ONLY_LOAD_THE_FIRST_N_POINTS = 20000000;

#include <gudhi/Debug_utils.h>
#include <gudhi/Clock.h>
#include <gudhi/Tangential_complex.h>
#include <gudhi/sparsify_point_set.h>
#include <gudhi/random_point_generators.h>
#include <gudhi/Tangential_complex/utilities.h>

#include <CGAL/assertions_behaviour.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/range/adaptor/strided.hpp>

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <cmath>  // for std::sqrt

#ifdef GUDHI_USE_TBB
#include <tbb/task_scheduler_init.h>
#endif
#include "XML_exporter.h"
#include "RIB_exporter.h"
#define GUDHI_TC_EXPORT_PERFORMANCE_DATA
#define GUDHI_TC_SET_PERFORMANCE_DATA(value_name, value) \
        XML_perf_data::set(value_name, value);


namespace subsampl = Gudhi::subsampling;
namespace tc = Gudhi::tangential_complex;

const char * const BENCHMARK_SCRIPT_FILENAME = "benchmark_script.txt";

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_d Point;
typedef Kernel::Vector_d Vector;
typedef tc::Tangential_complex<
Kernel, CGAL::Dynamic_dimension_tag,
CGAL::Parallel_tag> TC;
typedef TC::Simplex Simplex;
typedef TC::Simplex_set Simplex_set;

class XML_perf_data {
 public:
  typedef Streaming_XML_exporter<std::string> XML_exporter;

  XML_perf_data(const std::string &filename)
      : m_xml(filename, "ContainerPerformance", "Perf",
              construct_subelements_names()) { }

  virtual ~XML_perf_data() { }

  static XML_perf_data &get() {
    static XML_perf_data singleton(build_filename());
    return singleton;
  }

  template <typename Value_type>
  static void set(const std::string &name, Value_type value) {
    get().set_data(name, value);
  }

  static void commit() {
    get().commit_current_element();
  }

 protected:

  static std::string build_filename() {
    std::stringstream sstr;
    sstr << "perf_logs/Performance_log_" << time(0) << ".xml";
    return sstr.str();
  }

  static std::vector<std::string> construct_subelements_names() {
    std::vector<std::string> subelements;
    subelements.push_back("Input");
    subelements.push_back("Param1");
    subelements.push_back("Param2");
    subelements.push_back("Param3");
    subelements.push_back("Intrinsic_dim");
    subelements.push_back("Ambient_dim");
    subelements.push_back("Num_threads");
    subelements.push_back("Sparsity");
    subelements.push_back("Max_perturb");
    subelements.push_back("Num_points_in_input");
    subelements.push_back("Num_points");
    subelements.push_back("Perturb_technique");
    subelements.push_back("Perturb_which_points");
    subelements.push_back("Initial_num_inconsistent_local_tr");
    subelements.push_back("Best_num_inconsistent_local_tr");
    subelements.push_back("Final_num_inconsistent_local_tr");
    subelements.push_back("Init_time");
    subelements.push_back("Comput_time");
    subelements.push_back("Perturb_successful");
    subelements.push_back("Perturb_time");
    subelements.push_back("Perturb_steps");
    subelements.push_back("Result_pure_pseudomanifold");
    subelements.push_back("Result_num_wrong_dim_simplices");
    subelements.push_back("Result_num_wrong_number_of_cofaces");
    subelements.push_back("Result_num_unconnected_stars");
    subelements.push_back("Info");

    return subelements;
  }

  void set_data(const std::string &name, const std::string &value) {
    m_current_element[name] = value;
  }

  template <typename Value_type>
  void set_data(const std::string &name, Value_type value) {
    std::stringstream sstr;
    sstr << value;
    set_data(name, sstr.str());
  }

  void commit_current_element() {
    m_xml.add_element(m_current_element);
    m_current_element.clear();
  }

  XML_exporter m_xml;
  XML_exporter::Element_with_map m_current_element;
};

template<
typename Kernel, typename OutputIteratorPoints>
bool load_points_from_file(
                           const std::string &filename,
                           OutputIteratorPoints points,
                           std::size_t only_first_n_points = (std::numeric_limits<std::size_t>::max)()) {
  typedef typename Kernel::Point_d Point;

  std::ifstream in(filename);
  if (!in.is_open()) {
    std::cerr << "Could not open '" << filename << "'" << std::endl;
    return false;
  }

  Kernel k;
  Point p;
  int num_ppints;
  in >> num_ppints;

  std::size_t i = 0;
  while (i < only_first_n_points && in >> p) {
    *points++ = p;
    ++i;
  }

#ifdef DEBUG_TRACES
  std::cerr << "'" << filename << "' loaded." << std::endl;
#endif

  return true;
}

template<
typename Kernel, typename Tangent_space_basis,
typename OutputIteratorPoints, typename OutputIteratorTS>
bool load_points_and_tangent_space_basis_from_file(
                                                   const std::string &filename,
                                                   OutputIteratorPoints points,
                                                   OutputIteratorTS tangent_spaces,
                                                   int intrinsic_dim,
                                                   std::size_t only_first_n_points = (std::numeric_limits<std::size_t>::max)()) {
  typedef typename Kernel::Point_d Point;
  typedef typename Kernel::Vector_d Vector;

  std::ifstream in(filename);
  if (!in.is_open()) {
    std::cerr << "Could not open '" << filename << "'" << std::endl;
    return false;
  }

  Kernel k;
  Point p;
  int num_ppints;
  in >> num_ppints;

  std::size_t i = 0;
  while (i < only_first_n_points && in >> p) {
    *points++ = p;

    Tangent_space_basis tsb(i);
    for (int d = 0; d < intrinsic_dim; ++d) {
      Vector v;
      in >> v;
      tsb.push_back(tc::internal::normalize_vector(v, k));
    }
    *tangent_spaces++ = tsb;
    ++i;
  }

#ifdef DEBUG_TRACES
  std::cerr << "'" << filename << "' loaded." << std::endl;
#endif

  return true;
}

// color_inconsistencies: only works if p_complex = NULL
template <typename TC>
bool export_to_off(
                   TC const& tc,
                   std::string const& input_name_stripped,
                   std::string const& suffix,
                   bool color_inconsistencies = false,
                   typename TC::Simplicial_complex const* p_complex = NULL,
                   Simplex_set const *p_simpl_to_color_in_red = NULL,
                   Simplex_set const *p_simpl_to_color_in_green = NULL,
                   Simplex_set const *p_simpl_to_color_in_blue = NULL) {
#ifdef TC_NO_EXPORT
  return true;
#endif

  CGAL::Identity<Point> proj_functor;

  if (tc.intrinsic_dimension() <= 3) {
    std::stringstream output_filename;
    output_filename << "output/" << input_name_stripped << "_"
        << tc.intrinsic_dimension() << "_in_R"
        << tc.ambient_dimension() << "_"
        << tc.number_of_vertices() << "v"
        << suffix << ".off";
    std::ofstream off_stream(output_filename.str().c_str());

    if (p_complex) {
#ifndef TC_NO_EXPORT
      tc.export_to_off(
                       *p_complex, off_stream,
                       p_simpl_to_color_in_red,
                       p_simpl_to_color_in_green,
                       p_simpl_to_color_in_blue,
                       proj_functor);
#endif
    } else {
      tc.export_to_off(
                       off_stream, color_inconsistencies,
                       p_simpl_to_color_in_red,
                       p_simpl_to_color_in_green,
                       p_simpl_to_color_in_blue,
                       NULL,
                       proj_functor);
    }
    return true;
  }
  return false;
}

void make_tc(std::vector<Point> &points,
             TC::TS_container const& tangent_spaces,  // can be empty
             int intrinsic_dim,
             double sparsity = 0.01,
             double max_perturb = 0.005,
             bool perturb = true,
             bool add_high_dim_simpl = false,
             bool collapse = false,
             double time_limit_for_perturb = 0.,
             const char *input_name = "tc") {
  Kernel k;

  if (sparsity > 0. && !tangent_spaces.empty()) {
    std::cerr << "Error: cannot sparsify point set with pre-computed normals.\n";
    return;
  }

  //===========================================================================
  // Init
  //===========================================================================
  Gudhi::Clock t;

  // Get input_name_stripped
  std::string input_name_stripped(input_name);
  size_t slash_index = input_name_stripped.find_last_of('/');
  if (slash_index == std::string::npos)
    slash_index = input_name_stripped.find_last_of('\\');
  if (slash_index == std::string::npos)
    slash_index = 0;
  else
    ++slash_index;
  input_name_stripped = input_name_stripped.substr(
                                                   slash_index, input_name_stripped.find_last_of('.') - slash_index);

  GUDHI_TC_SET_PERFORMANCE_DATA("Num_points_in_input", points.size());

#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  std::vector<Point> points_not_sparse = points;
#endif

  //===========================================================================
  // Sparsify point set if requested
  //===========================================================================
  if (sparsity > 0.) {
    std::size_t num_points_before = points.size();
    std::vector<Point> sparsified_points;
    subsampl::sparsify_point_set(k, points, sparsity*sparsity,
                                 std::back_inserter(sparsified_points));
    sparsified_points.swap(points);
    std::cerr << "Number of points before/after sparsification: "
        << num_points_before << " / " << points.size() << "\n";

#ifdef GUDHI_TC_EXPORT_SPARSIFIED_POINT_SET
    std::ofstream ps_stream("output/sparsified_point_set.txt");
    tc::internal::export_point_set(k, points, ps_stream);
#endif
  }

  GUDHI_TC_SET_PERFORMANCE_DATA("Sparsity", sparsity);
  GUDHI_TC_SET_PERFORMANCE_DATA("Max_perturb", max_perturb);
  GUDHI_TC_SET_PERFORMANCE_DATA("Num_points", points.size());

  //===========================================================================
  // Compute Tangential Complex
  //===========================================================================

  TC tc(
        points,
        intrinsic_dim,
#ifdef GUDHI_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
      points_not_sparse.begin(), points_not_sparse.end(),
#endif
      k);

  if (!tangent_spaces.empty()) {
    tc.set_tangent_planes(tangent_spaces);
  }

  t.end();
  double init_time = t.num_seconds();

  t.begin();
  tc.compute_tangential_complex();
  t.end();
  double computation_time = t.num_seconds();

  //===========================================================================
  // Export to OFF
  //===========================================================================

  // Create complex
  int max_dim = -1;
  TC::Simplicial_complex complex;
  Simplex_set inconsistent_simplices;
  max_dim = tc.create_complex(complex, true, false, 2, &inconsistent_simplices);

  // TODO(CJ): TEST
  Gudhi::Simplex_tree<> stree;
  tc.create_complex(stree, true, false);
  // std::cerr << stree;

  t.begin();
  bool ret = export_to_off(
                           tc, input_name_stripped, "_INITIAL_TC", true,
                           &complex, &inconsistent_simplices);
  t.end();
  double export_before_time = (ret ? t.num_seconds() : -1);

  unsigned int num_perturb_steps = 0;
  double perturb_time = -1;
  double export_after_perturb_time = -1.;
  bool perturb_success = false;
  if (perturb) {
    //=========================================================================
    // Try to fix inconsistencies by perturbing points
    //=========================================================================
    t.begin();
    auto fix_result =
        tc.fix_inconsistencies_using_perturbation(max_perturb, time_limit_for_perturb);
    t.end();
    perturb_time = t.num_seconds();

    perturb_success = fix_result.success;
    GUDHI_TC_SET_PERFORMANCE_DATA("Initial_num_inconsistent_local_tr",
                                  fix_result.initial_num_inconsistent_stars);
    GUDHI_TC_SET_PERFORMANCE_DATA("Best_num_inconsistent_local_tr",
                                  fix_result.best_num_inconsistent_stars);
    GUDHI_TC_SET_PERFORMANCE_DATA("Final_num_inconsistent_local_tr",
                                  fix_result.final_num_inconsistent_stars);

    //=========================================================================
    // Export to OFF
    //=========================================================================

    // Re-build the complex
    Simplex_set inconsistent_simplices;
    max_dim = tc.create_complex(complex, true, false, 2, &inconsistent_simplices);

    t.begin();
    bool exported = export_to_off(
                                  tc, input_name_stripped, "_AFTER_FIX", true, &complex,
                                  &inconsistent_simplices);
    t.end();
    export_after_perturb_time = (exported ? t.num_seconds() : -1);

    //std::string fn = "output/inc_stars/";
    //fn += input_name_stripped;
    //tc.export_inconsistent_stars_to_OFF_files(fn);

#if !defined(TC_NO_EXPORT) && defined(TC_EXPORT_TO_RIB)
    std::ofstream rib(std::string("output/") + input_name_stripped + ".rib");
    RIB_exporter<TC::Points, TC::Simplicial_complex::Simplex_set> rib_exporter(
                                                                               tc.points(),
                                                                               complex.simplex_range(),
                                                                               rib,
                                                                               input_name_stripped + ".tif",
                                                                               false,  // is_preview
                                                                               std::make_tuple(2, 4, 6),
                                                                               1600, 503  // resolution
                                                                               );
    rib_exporter.write_file();

    std::ofstream rib_LQ(std::string("output/") + input_name_stripped + "_LQ.rib");
    RIB_exporter<TC::Points, TC::Simplicial_complex::Simplex_set> rib_exporter_LQ(
                                                                                  tc.points(),
                                                                                  complex.simplex_range(),
                                                                                  rib_LQ,
                                                                                  input_name_stripped + "_LQ.tif",
                                                                                  true,  // is_preview
                                                                                  std::make_tuple(0, 4, 5)
                                                                                  );
    rib_exporter_LQ.write_file();
#endif
  } else {
    GUDHI_TC_SET_PERFORMANCE_DATA("Initial_num_inconsistent_local_tr", "N/A");
    GUDHI_TC_SET_PERFORMANCE_DATA("Best_num_inconsistent_local_tr", "N/A");
    GUDHI_TC_SET_PERFORMANCE_DATA("Final_num_inconsistent_local_tr", "N/A");
  }

  max_dim = tc.create_complex(complex, true, false, 2);

  complex.display_stats();

  if (intrinsic_dim == 2)
    complex.euler_characteristic(true);

  //===========================================================================
  // Collapse
  //===========================================================================
  if (collapse) {
    complex.collapse(max_dim);
    complex.display_stats();
  }

  //===========================================================================
  // Is the result a pure pseudomanifold?
  //===========================================================================
  std::size_t num_wrong_dim_simplices,
      num_wrong_number_of_cofaces,
      num_unconnected_stars;
  Simplex_set wrong_dim_simplices;
  Simplex_set wrong_number_of_cofaces_simplices;
  Simplex_set unconnected_stars_simplices;
  bool is_pure_pseudomanifold = complex.is_pure_pseudomanifold(
                                                               intrinsic_dim, tc.number_of_vertices(),
                                                               false,  // do NOT allow borders
                                                               false, 1,
                                                               &num_wrong_dim_simplices, &num_wrong_number_of_cofaces,
                                                               &num_unconnected_stars,
                                                               &wrong_dim_simplices, &wrong_number_of_cofaces_simplices,
                                                               &unconnected_stars_simplices);

  //===========================================================================
  // Export to OFF
  //===========================================================================

  double export_after_collapse_time = -1.;
  if (collapse) {
    t.begin();
    bool exported = export_to_off(
                                  tc, input_name_stripped, "_AFTER_COLLAPSE", false, &complex,
                                  &wrong_dim_simplices, &wrong_number_of_cofaces_simplices,
                                  &unconnected_stars_simplices);
    t.end();
    std::cerr
        << " OFF colors:\n"
        << "   * Red: wrong dim simplices\n"
        << "   * Green: wrong number of cofaces simplices\n"
        << "   * Blue: not-connected stars\n";
    export_after_collapse_time = (exported ? t.num_seconds() : -1.);
  }

  //===========================================================================
  // Display info
  //===========================================================================

  std::cerr
      << "\n================================================\n"
      << "Number of vertices: " << tc.number_of_vertices() << "\n"
      << "Computation times (seconds): \n"
      << "  * Tangential complex: " << init_time + computation_time << "\n"
      << "    - Init + kd-tree = " << init_time << "\n"
      << "    - TC computation = " << computation_time << "\n"
      << "  * Export to OFF (before perturb): " << export_before_time << "\n"
      << "  * Fix inconsistencies 1: " << perturb_time
      << " (" << num_perturb_steps << " steps) ==> "
      << (perturb_success ? "FIXED" : "NOT fixed") << "\n"
      << "  * Export to OFF (after perturb): " << export_after_perturb_time << "\n"
      << "  * Export to OFF (after collapse): "
      << export_after_collapse_time << "\n"
      << "================================================\n";

  //===========================================================================
  // Export info
  //===========================================================================
  GUDHI_TC_SET_PERFORMANCE_DATA("Init_time", init_time);
  GUDHI_TC_SET_PERFORMANCE_DATA("Comput_time", computation_time);
  GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_successful",
                                (perturb_success ? 1 : 0));
  GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_time", perturb_time);
  GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_steps", num_perturb_steps);
  GUDHI_TC_SET_PERFORMANCE_DATA("Result_pure_pseudomanifold",
                                (is_pure_pseudomanifold ? 1 : 0));
  GUDHI_TC_SET_PERFORMANCE_DATA("Result_num_wrong_dim_simplices",
                                num_wrong_dim_simplices);
  GUDHI_TC_SET_PERFORMANCE_DATA("Result_num_wrong_number_of_cofaces",
                                num_wrong_number_of_cofaces);
  GUDHI_TC_SET_PERFORMANCE_DATA("Result_num_unconnected_stars",
                                num_unconnected_stars);
  GUDHI_TC_SET_PERFORMANCE_DATA("Info", "");
}

int main() {
  CGAL::set_error_behaviour(CGAL::ABORT);

#ifdef GUDHI_USE_TBB
#ifdef _DEBUG
  int num_threads = 1;
#else
  int num_threads = tbb::task_scheduler_init::default_num_threads() - 4;
#endif
#endif

  std::ifstream script_file;
  script_file.open(BENCHMARK_SCRIPT_FILENAME);
  // Script?
  // Script file format: each line gives
  //    - Filename (point set) or "generate_XXX" (point set generation)
  //    - Ambient dim
  //    - Intrinsic dim
  //    - Number of iterations with these parameters
  if (script_file.is_open()) {
    int i = 1;
#ifdef GUDHI_USE_TBB
#ifdef BENCHMARK_WITH_1_TO_MAX_THREADS
    for (num_threads = 1;
         num_threads <= tbb::task_scheduler_init::default_num_threads();
         ++num_threads)
#endif
#endif
      /*for (Concurrent_mesher_config::get().num_work_items_per_batch = 5 ;
        Concurrent_mesher_config::get().num_work_items_per_batch < 100 ;
        Concurrent_mesher_config::get().num_work_items_per_batch += 5)*/ {
#ifdef GUDHI_USE_TBB
      tbb::task_scheduler_init init(
                                    num_threads > 0 ? num_threads : tbb::task_scheduler_init::automatic);
#endif

      std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' found.\n";
      script_file.seekg(0);
      while (script_file.good()) {
        std::string line;
        std::getline(script_file, line);
        if (line.size() > 1 && line[0] != '#') {
          boost::replace_all(line, "\t", " ");
          boost::trim_all(line);
          std::cerr << "\n\n";
          std::cerr << "*****************************************\n";
          std::cerr << "******* " << line << "\n";
          std::cerr << "*****************************************\n";
          std::stringstream sstr(line);

          std::string input;
          std::string param1;
          std::string param2;
          std::string param3;
          std::size_t num_points;
          int ambient_dim;
          int intrinsic_dim;
          double sparsity;
          double max_perturb;
          char perturb, add_high_dim_simpl, collapse;
          double time_limit_for_perturb;
          int num_iteration;
          sstr >> input;
          sstr >> param1;
          sstr >> param2;
          sstr >> param3;
          sstr >> num_points;
          sstr >> ambient_dim;
          sstr >> intrinsic_dim;
          sstr >> sparsity;
          sstr >> max_perturb;
          sstr >> perturb;
          sstr >> add_high_dim_simpl;
          sstr >> collapse;
          sstr >> time_limit_for_perturb;
          sstr >> num_iteration;

          for (int j = 0; j < num_iteration; ++j) {
            std::string input_stripped = input;
            size_t slash_index = input_stripped.find_last_of('/');
            if (slash_index == std::string::npos)
              slash_index = input_stripped.find_last_of('\\');
            if (slash_index == std::string::npos)
              slash_index = 0;
            else
              ++slash_index;
            input_stripped = input_stripped.substr(
                                                   slash_index, input_stripped.find_last_of('.') - slash_index);

            GUDHI_TC_SET_PERFORMANCE_DATA("Input", input_stripped);
            GUDHI_TC_SET_PERFORMANCE_DATA("Param1", param1);
            GUDHI_TC_SET_PERFORMANCE_DATA("Param2", param2);
            GUDHI_TC_SET_PERFORMANCE_DATA("Param3", param3);
            GUDHI_TC_SET_PERFORMANCE_DATA("Ambient_dim", ambient_dim);
            GUDHI_TC_SET_PERFORMANCE_DATA("Intrinsic_dim", intrinsic_dim);
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_technique", "Tangential_translation");
            GUDHI_TC_SET_PERFORMANCE_DATA("Perturb_which_points", "Center_vertex");

#ifdef GUDHI_USE_TBB
            GUDHI_TC_SET_PERFORMANCE_DATA(
                                          "Num_threads",
                                          (num_threads == -1 ? tbb::task_scheduler_init::default_num_threads() : num_threads));
#else
            GUDHI_TC_SET_PERFORMANCE_DATA("Num_threads", "N/A");
#endif

            std::cerr << "\nTC #" << i << "...\n";

#ifdef GUDHI_TC_PROFILING
            Gudhi::Clock t_gen;
#endif

            std::vector<Point> points;
            TC::TS_container tangent_spaces;

            if (input == "generate_moment_curve") {
              points = Gudhi::generate_points_on_moment_curve<Kernel>(
                                                                      num_points, ambient_dim,
                                                                      std::atof(param1.c_str()), std::atof(param2.c_str()));
            } else if (input == "generate_plane") {
              points = Gudhi::generate_points_on_plane<Kernel>(
                                                               num_points, intrinsic_dim, ambient_dim);
            } else if (input == "generate_sphere_d") {
              points = Gudhi::generate_points_on_sphere_d<Kernel>(
                                                                  num_points, ambient_dim,
                                                                  std::atof(param1.c_str()),  // radius
                                                                  std::atof(param2.c_str()));  // radius_noise_percentage
            } else if (input == "generate_two_spheres_d") {
              points = Gudhi::generate_points_on_two_spheres_d<Kernel>(
                                                                       num_points, ambient_dim,
                                                                       std::atof(param1.c_str()),
                                                                       std::atof(param2.c_str()),
                                                                       std::atof(param3.c_str()));
            } else if (input == "generate_3sphere_and_circle_d") {
              GUDHI_CHECK(intrinsic_dim == 3,
                          std::logic_error("Intrinsic dim should be 3"));
              GUDHI_CHECK(ambient_dim == 5,
                          std::logic_error("Ambient dim should be 5"));
              points = Gudhi::generate_points_on_3sphere_and_circle<Kernel>(
                                                                            num_points,
                                                                            std::atof(param1.c_str()));
            } else if (input == "generate_torus_3D") {
              points = Gudhi::generate_points_on_torus_3D<Kernel>(
                                                                  num_points,
                                                                  std::atof(param1.c_str()),
                                                                  std::atof(param2.c_str()),
                                                                  param3 == "Y");
            } else if (input == "generate_torus_d") {
              points = Gudhi::generate_points_on_torus_d<Kernel>(
                                                                 num_points,
                                                                 intrinsic_dim,
                                                                 param1 == "Y",  // uniform
                                                                 std::atof(param2.c_str()));  // radius_noise_percentage
            } else if (input == "generate_klein_bottle_3D") {
              points = Gudhi::generate_points_on_klein_bottle_3D<Kernel>(
                                                                         num_points,
                                                                         std::atof(param1.c_str()), std::atof(param2.c_str()));
            } else if (input == "generate_klein_bottle_4D") {
              points = Gudhi::generate_points_on_klein_bottle_4D<Kernel>(
                                                                         num_points,
                                                                         std::atof(param1.c_str()), std::atof(param2.c_str()),
                                                                         std::atof(param3.c_str()));  // noise
            } else if (input == "generate_klein_bottle_variant_5D") {
              points = Gudhi::generate_points_on_klein_bottle_variant_5D<Kernel>(
                                                                                 num_points,
                                                                                 std::atof(param1.c_str()), std::atof(param2.c_str()));
            } else {
              // Contains tangent space basis
              if (input.substr(input.size() - 3) == "pwt") {
                load_points_and_tangent_space_basis_from_file
                    <Kernel, typename TC::Tangent_space_basis > (
                                                                 input, std::back_inserter(points),
                                                                 std::back_inserter(tangent_spaces),
                                                                 intrinsic_dim,
                                                                 ONLY_LOAD_THE_FIRST_N_POINTS);
              } else {
                load_points_from_file<Kernel>(
                                              input, std::back_inserter(points),
                                              ONLY_LOAD_THE_FIRST_N_POINTS);
              }
            }

#ifdef GUDHI_TC_PROFILING
            t_gen.end();
            std::cerr << "Point set generated/loaded in " << t_gen.num_seconds()
                << " seconds.\n";
#endif

            if (!points.empty()) {
#if defined(TC_INPUT_STRIDES) && TC_INPUT_STRIDES > 1
              auto p = points | boost::adaptors::strided(TC_INPUT_STRIDES);
              std::vector<Point> points(p.begin(), p.end());
              std::cerr << "****************************************\n"
                  << "WARNING: taking 1 point every " << TC_INPUT_STRIDES
                  << " points.\n"
                  << "****************************************\n";
#endif

              make_tc(points, tangent_spaces, intrinsic_dim,
                      sparsity, max_perturb,
                      perturb == 'Y', add_high_dim_simpl == 'Y', collapse == 'Y',
                      time_limit_for_perturb, input.c_str());

              std::cerr << "TC #" << i++ << " done.\n";
              std::cerr << "\n---------------------------------\n";
            } else {
              std::cerr << "TC #" << i++ << ": no points loaded.\n";
            }

            XML_perf_data::commit();
          }
        }
      }
      script_file.seekg(0);
      script_file.clear();
    }

    script_file.close();
  }    // Or not script?
  else {
    std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' NOT found.\n";
  }

  // system("pause");
  return 0;
}
