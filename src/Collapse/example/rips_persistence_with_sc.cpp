#include <gudhi/FlagComplexSpMatrix.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Rips_edge_list.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>
#include <gudhi/PointSetGen.h>

// Types definition
using Vector_of_points         = std::vector<Point>;

using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = double;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Rips_edge_list      = Gudhi::rips_edge_list::Rips_edge_list<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Distance_matrix = std::vector<std::vector<Filtration_value>>;


class filt_edge_to_dist_matrix
{
public:
    template<class Distance_matrix, class Filtered_sorted_edge_list>
    filt_edge_to_dist_matrix(Distance_matrix & distance_mat, Filtered_sorted_edge_list & edge_filt, std::size_t number_of_points)
    {
        double inf = std::numeric_limits<double>::max();
        doubleVector distances ; 
        std::pair<std::size_t, std::size_t> e;
        for(std::size_t indx = 0; indx < number_of_points; indx++) { 
            for (std::size_t j = 0; j <= indx; j++) {
               if( j == indx)
                    distances.push_back(0);
                
                else
                     distances.push_back(inf);
            }
            distance_mat.push_back(distances);
            distances.clear();
        }

        for(auto edIt = edge_filt.begin(); edIt != edge_filt.end(); edIt++) {
            e=std::minmax(std::get<1>(*edIt),std::get<2>(*edIt));
            distance_mat.at(std::get<1>(e)).at(std::get<0>(e)) = std::get<0>(*edIt);
        } 
    }
};


int main(int argc, char * const argv[]) {
	
    auto the_begin = std::chrono::high_resolution_clock::now();
    PointSetGen point_generator;
    std::string out_file_name   = "default";
    std::string in_file_name    = "default";
    std::size_t number_of_points;
    
    typedef size_t Vertex_handle;
    typedef std::vector< std::tuple<Filtration_value, Vertex_handle, Vertex_handle > > Filtered_sorted_edge_list;

    int     dimension;
    double  end_threshold;
    double  steps;
    char    manifold;

    Vector_of_points * point_vector;
    
    int dim_max  = 2;

    point_generator.program_options(argc, argv, steps, end_threshold, manifold, dimension, dim_max, in_file_name, out_file_name);
    
    std::cout << "The current input values to run the program is: "<< std::endl;
    std::cout << "steps, end_threshold, manifold, dimension, max_complex_dimension, in_file_name, out_file_name" << std::endl;
    std::cout << steps << ", " << end_threshold << ", " << manifold << ", " << dimension << ", " << dim_max << ", " << in_file_name << ", " << out_file_name << std::endl;
    
    Map map_empty;
    
    std::string filediag_aft ("./PersistenceOutput/collapsed_persistence_diags") ;
    
    filediag_aft = filediag_aft+"_"+ out_file_name+ ".txt";

    Distance_matrix distances;
    Distance_matrix *sparse_distances = new Distance_matrix();


    if(manifold == 'f') {
        Gudhi::Points_off_reader<Point> off_reader(in_file_name);
        if (!off_reader.is_valid()) {
            std::cerr << "Unable to read file " << in_file_name << "\n";
            exit(-1);  // ----- >>
        }

        point_vector = new Vector_of_points(off_reader.get_point_cloud().begin(), off_reader.get_point_cloud().end());
        dimension = point_vector->at(0).dimension() ;
        number_of_points = point_vector->size();
        std::cout << "Successfully read " << number_of_points << " point_vector.\n";
        std::cout << "Ambient dimension is " << dimension << ".\n";
    }
    else if (manifold == 'm'){
        std::string csv_file_name(in_file_name);
        distances = Gudhi::read_lower_triangular_matrix_from_csv_file<Filtration_value>(csv_file_name);
        number_of_points = distances.size();
        std::cout << "Read the distance matrix succesfully, of size: " << number_of_points << std::endl;
    }
    else {
        std::cerr << "Wrong parameters for input manifold..." <<std::endl;  
        exit(-1); 
    }
    
    std::cout << "Point Set Generated."  <<std::endl;
 
    Filtered_sorted_edge_list * edge_t = new Filtered_sorted_edge_list();
    std::cout << "Computing the one-skeleton for threshold: " << end_threshold << std::endl; 
    
    if(manifold == 'm') {                                           //Input is a distance 'm'atrix //Creating the edge list
        Rips_edge_list Rips_edge_list_from_file(distances, end_threshold);
        Rips_edge_list_from_file.create_edges(*edge_t);
        std::cout<< "Sorted edge list computed" << std::endl;
    }
    else{ //Point cloud input  //Creating the edge list
        Rips_edge_list Rips_edge_list_from_points(*point_vector, end_threshold, Gudhi::Euclidean_distance());
        Rips_edge_list_from_points.create_edges(*edge_t);
        std::cout<< "Sorted edge list computed" << std::endl;
        std::cout << "Total number of edges before collapse are: " << edge_t->size() << std::endl;
    }
    
    //Now we will perform filtered edge collapse to sparsify the edge list edge_t.
    std::cout<< "Filtered edge collapse begins" << std::endl;
    FlagComplexSpMatrix * mat_filt_edge_coll  = new FlagComplexSpMatrix(number_of_points,*edge_t);
    std::cout<< "Matrix instansiated" << std::endl;
    if(edge_t->size() >0){
        delete edge_t;
        edge_t = new Filtered_sorted_edge_list();
        *edge_t = mat_filt_edge_coll->filtered_edge_collapse(); 
        filt_edge_to_dist_matrix(*sparse_distances, *edge_t, number_of_points); 
        std::cout << "Total number of vertices after collapse in the sparse matrix are: " << mat_filt_edge_coll->num_vertices() << std::endl;
    }
    else
    {
        std::cerr << "Total number of egdes are zero." <<std::endl;
        exit(-1) ;
    }
 
    // Rips_complex rips_complex_before_collapse(distances, end_threshold);
    Rips_complex rips_complex_after_collapse(*sparse_distances, end_threshold);
    
    Simplex_tree   simplex_tree_aft;
    rips_complex_after_collapse.create_complex(simplex_tree_aft, dim_max);

    std::cout << "The complex contains " << simplex_tree_aft.num_simplices() << " simplices  after collapse. \n";
    std::cout << "   and has dimension " << simplex_tree_aft.dimension() << " \n";

    // Sort the simplices in the order of the filtration
    simplex_tree_aft.initialize_filtration();
    // Compute the persistence diagram of the complex
    Persistent_cohomology pcoh_aft(simplex_tree_aft);
    // initializes the coefficient field for homology
    pcoh_aft.init_coefficients(2);

    // pcoh_bfr.compute_persistent_cohomology(steps);
    pcoh_aft.compute_persistent_cohomology(steps);
    if (filediag_aft.empty()) {
        pcoh_aft.output_diagram();
    } 
    else {
        std::ofstream out(filediag_aft);
        pcoh_aft.output_diagram(out);
        out.close();
      }
    
    auto the_end = std::chrono::high_resolution_clock::now();   

    std::cout << "Total computation time : " <<  std::chrono::duration<double, std::milli>(the_end- the_begin).count()
                  << " ms\n" << std::endl;         
    return 0;

}
  
