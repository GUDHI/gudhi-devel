#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
//#include <gudhi/Fake_simplex_tree.h>
//#include <gudhi/Sb_wrapper.h>
#include <gudhi/distance_functions.h>
#include <gudhi/Points_fvecs_reader.h>


#include <CGAL/Epick_d.h>

#include <iostream>
#include <string>
#include <vector>
#include <ctime>

void usage(int nbArgs, char * const progName) {
    std::cerr << "Error: Number of arguments (" << nbArgs << ") is not correct\n";
    std::cerr << "Usage: " << progName << " filename.fvecs threshold dim_max [ouput_file.txt]\n";
    std::cerr << "       i.e.: " << progName << " ../../data/points/alphacomplexdoc.fvecs 60.0\n";
    exit(-1);  // ----- >>
}

int main(int argc, char **argv) {
    if (argc != 4) usage(argc, (argv[0] - 1));

    std::string file_name(argv[1]);
    double threshold = atof(argv[2]);
    int dim_max = atoi(argv[3]);

    // Type definitions
    using K = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
    using Point = typename K::Point_d;
    using Simplex_tree = Gudhi::Simplex_tree<>;
    //using Simplex_tree = Gudhi::Fake_simplex_tree;
    //using Simplex_tree = Gudhi::Sb_wrapper;
    using Filtration_value = Simplex_tree::Filtration_value;
    using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
    using Point_vector = std::vector<Point>;

    // ----------------------------------------------------------------------------
    // Init of a Rips complex from an fvecs file
    // ----------------------------------------------------------------------------
    Point_vector point_vector;
    Gudhi::load_points_from_fvecs_file<K>(file_name, std::back_insert_iterator< Point_vector >(point_vector));

    Rips_complex rips_complex_from_file(point_vector, threshold, Gudhi::Euclidean_distance());


    Simplex_tree stree;

    clock_t start, end;
    start = clock();
    rips_complex_from_file.create_complex(stree, dim_max);
    end = clock();

    std::cout << "Strong witness complex took "<< static_cast<double>(end - start) / CLOCKS_PER_SEC << " s." << std::endl;
    //std::cout << "Rips complex is of dimension " << stree.dimension() <<  " - " << stree.num_simplices() << " simplices." << std::endl;

    return 0;
}
