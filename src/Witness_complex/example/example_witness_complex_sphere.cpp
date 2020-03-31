#define BOOST_PARAMETER_MAX_ARITY 12

#include <gudhi/Simplex_tree.h>
#include <gudhi/Euclidean_witness_complex.h>
#include <gudhi/pick_n_random_points.h>
#include <gudhi/choose_n_farthest_points.h>
#include <gudhi/reader_utils.h>

#include <CGAL/Epick_d.h>

#include <iostream>
#include <fstream>
#include <ctime>
#include <utility>
#include <string>
#include <vector>

#include "generators.h"

/** Write a gnuplot readable file.
 *  Data range is a random access range of pairs (arg, value)
 */
template <typename Data_range>
void write_data(Data_range& data, std::string filename) {
  std::ofstream ofs(filename, std::ofstream::out);
  for (auto entry : data) ofs << entry.first << ", " << entry.second << "\n";
  ofs.close();
}

int main(int argc, char* const argv[]) {
  using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Witness_complex = Gudhi::witness_complex::Euclidean_witness_complex<Kernel>;

  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " number_of_landmarks \n";
    return 0;
  }

  int number_of_landmarks = atoi(argv[1]);

  std::vector<std::pair<int, double> > l_time;

  // Generate points
  for (int nbP = 500; nbP < 10000; nbP += 500) {
    clock_t start, end;
    // Construct the Simplex Tree
    Gudhi::Simplex_tree<> simplex_tree;
    Point_Vector point_vector, landmarks;
    generate_points_sphere(point_vector, nbP, 4);
    std::clog << "Successfully generated " << point_vector.size() << " points.\n";
    std::clog << "Ambient dimension is " << point_vector[0].size() << ".\n";

    // Choose landmarks
    start = clock();
    // Gudhi::subsampling::pick_n_random_points(point_vector, number_of_landmarks, std::back_inserter(landmarks));
    Gudhi::subsampling::choose_n_farthest_points(K(), point_vector, number_of_landmarks,
                                                 Gudhi::subsampling::random_starting_point,
                                                 std::back_inserter(landmarks));

    // Compute witness complex
    Witness_complex witness_complex(landmarks, point_vector);
    witness_complex.create_complex(simplex_tree, 0);
    end = clock();
    double time = static_cast<double>(end - start) / CLOCKS_PER_SEC;
    std::clog << "Witness complex for " << number_of_landmarks << " landmarks took " << time << " s. \n";
    std::clog << "Number of simplices is: " << simplex_tree.num_simplices() << "\n";
    l_time.push_back(std::make_pair(nbP, time));
  }
  write_data(l_time, "w_time.dat");
}
