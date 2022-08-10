/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <gudhi/Points_off_io.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Clock.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Cech_complex.h>
#include <gudhi/Simplex_tree.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>

#include "boost/filesystem.hpp"  // includes all needed Boost.Filesystem declarations

#include <string>
#include <vector>

// Types definition
using Simplex_tree = Gudhi::Simplex_tree<>;
using Filtration_value = Simplex_tree::Filtration_value;
using Point = std::vector<Filtration_value>;
using Points_off_reader = Gudhi::Points_off_reader<Point>;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;

template<typename Kernel>
Simplex_tree benchmark_cech(const std::string& off_file_points, const Filtration_value& radius, const int& dim_max, const bool exact) {
    using Point_cgal = typename Kernel::Point_d;
    using Points_off_reader_cgal = Gudhi::Points_off_reader<Point_cgal>;
    using Cech_complex = Gudhi::cech_complex::Cech_complex<Kernel, Simplex_tree>;

    // Extract the points from the file filepoints
    Points_off_reader_cgal off_reader_cgal(off_file_points);

    Gudhi::Clock cech_clock("Cech computation");
    Cech_complex cech_complex_from_points(off_reader_cgal.get_point_cloud(), radius, exact);
    Simplex_tree cech_stree;
    cech_complex_from_points.create_complex(cech_stree, dim_max);

    // ------------------------------------------
    // Display information about the Cech complex
    // ------------------------------------------
    double cech_sec = cech_clock.num_seconds();
    std::clog << cech_sec << "  ;  ";
    return cech_stree;
}

int main(int argc, char* argv[]) {
    boost::filesystem::path full_path(boost::filesystem::current_path());
    std::clog << "Current path is : " << full_path << std::endl;

    std::clog << "File name ; Radius ; Rips time ; Dim-3 Fast Cech time ; Dynamic_dim Fast Cech time ; "
                 "Dim-3 Safe Cech time ; Dynamic_dim Safe Cech time ; Dim-3 Exact Cech time ; Dynamic_dim Exact Cech time ; "
                 "Cech nb simplices ; Rips nb simplices;"
              << std::endl;
    boost::filesystem::directory_iterator end_itr;  // default construction yields past-the-end
    // For every ".off" file in the current directory, and for 3 predefined thresholds, compare Rips and various Cech constructions
    for (boost::filesystem::directory_iterator itr(boost::filesystem::current_path()); itr != end_itr; ++itr) {
        if (!boost::filesystem::is_directory(itr->status())) {
            if (itr->path().extension() == ".off") {
                Points_off_reader off_reader(itr->path().string());
                Point p0 = off_reader.get_point_cloud()[0];
                // Loop over the different thresholds
                for (Filtration_value radius = 0.1; radius < 0.35; radius += 0.1) {
                    std::clog << itr->path().stem() << "  ;  ";
                    std::clog << radius << "  ;  ";

                    Gudhi::Clock rips_clock("Rips computation");
                    Rips_complex rips_complex_from_points(off_reader.get_point_cloud(), radius, Gudhi::Euclidean_distance());
                    Simplex_tree rips_stree;
                    int dim_max = p0.size() - 1;
                    rips_complex_from_points.create_complex(rips_stree, dim_max);
                    // ------------------------------------------
                    // Display information about the Rips complex
                    // ------------------------------------------
                    double rips_sec = rips_clock.num_seconds();
                    std::clog << rips_sec << "  ;  ";

                    // --------------
                    // Cech complex
                    // --------------
                    // Fast
                    benchmark_cech<CGAL::Epick_d<CGAL::Dimension_tag<3>>>(itr->path().string(), radius, dim_max, false);
                    benchmark_cech<CGAL::Epick_d<CGAL::Dynamic_dimension_tag>>(itr->path().string(), radius, dim_max, false);
                    // Safe
                    benchmark_cech<CGAL::Epeck_d<CGAL::Dimension_tag<3>>>(itr->path().string(), radius, dim_max, false);
                    benchmark_cech<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>>(itr->path().string(), radius, dim_max, false);
                    // Exact
                    benchmark_cech<CGAL::Epeck_d<CGAL::Dimension_tag<3>>>(itr->path().string(), radius, dim_max, true);
                    auto cech_stree = benchmark_cech<CGAL::Epeck_d<CGAL::Dynamic_dimension_tag>>(itr->path().string(), radius, dim_max, true);

                    std::clog << cech_stree.num_simplices() << "  ;  ";
                    std::clog << rips_stree.num_simplices() << ";" << std::endl;
                }
            }
        }
    }

    return 0;
}
