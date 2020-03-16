#include <gudhi/TowerAssembler_FlagComplex.h>
#include <gudhi/Rips_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
#include <gudhi/Rips_edge_list.h>
#include <gudhi/distance_functions.h>
#include <gudhi/reader_utils.h>
#include <gudhi/PointSetGen.h>

// Types definition
using Vector_of_points         = std::vector<Point>;
using Vector_of_SM_pointers    = std::vector<FlagComplexSpMatrix*>;

using Simplex_tree = Gudhi::Simplex_tree<Gudhi::Simplex_tree_options_fast_persistence>;
using Filtration_value = double;
using Rips_complex = Gudhi::rips_complex::Rips_complex<Filtration_value>;
using Rips_edge_list      = Gudhi::rips_edge_list::Rips_edge_list<Filtration_value>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp>;
using Distance_matrix = std::vector<std::vector<Filtration_value>>;
class extract_sub_one_skeleton
{
public:
    template<class Filtered_sorted_edge_list,  class Fil_vector >
    extract_sub_one_skeleton(double threshold, Filtered_sorted_edge_list & current_edge_t, Filtered_sorted_edge_list & edge_t, Fil_vector & edge_filt ) {
 
        auto end_it = std::upper_bound(edge_filt.begin(), edge_filt.end(), threshold); // find_index(edge_t, threshold, 0, end_idx);
        size_t end_idx = std::distance(edge_filt.begin(), end_it);

        for( size_t idx = 0; idx < end_idx ; idx++) {
           current_edge_t.push_back(*edge_t.begin()); 
           edge_filt.erase(edge_filt.begin());
           edge_t.erase(edge_t.begin());
        }
        
    }
};
class extract_one_new_edge
{
public:
    template<class Filtered_sorted_edge_list,  class Fil_vector >
    extract_one_new_edge(Filtered_sorted_edge_list & current_edge_t, Filtered_sorted_edge_list & edge_t, Fil_vector & edge_filt ) {
        current_edge_t.push_back(*edge_t.begin()); 
        edge_filt.erase(edge_filt.begin());
        edge_t.erase(edge_t.begin());
      
    }
};
class strong_filtered_vertex_collapse
{
public:
    template<class Filtered_sorted_edge_list, class Distance_matrix>
    strong_filtered_vertex_collapse(std::size_t number_of_points, Filtered_sorted_edge_list & edge_t, double steps, Distance_matrix & sparse_distances)
    {
         auto the_collapse_begin = std::chrono::high_resolution_clock::now();
        //An additional vector <edge_filt> to perform binary search to find the index of given threshold
        std::vector<Filtration_value> * edge_filt = new std::vector<Filtration_value>();
        edge_filt->clear();
        for(auto edIt = edge_t.begin(); edIt != edge_t.end(); edIt++) {
            edge_filt->push_back(std::get<0>(*edIt));
        }
        double begin_thresold = edge_filt->at(0);
        double end_threshold = edge_filt->at(edge_filt->size()-1);
        double totAssembleTime = 0.0;
        Filtered_sorted_edge_list * sub_skeleton  = new Filtered_sorted_edge_list();  
        TowerAssembler_FlagComplex twr_assembler(number_of_points) ;

        std::cout<< "Begin and end thresholds (filteration values) are , " << begin_thresold <<", " << end_threshold << "." <<std::endl;

        // int iterations = (end_threshold - begin_thresold)/steps;
        // std::cout << "Total number of iterations to be run are: " << iterations << std::endl;
       
        auto threshold =  begin_thresold;  

        FlagComplexSpMatrix * mat_coll       = new FlagComplexSpMatrix(); 
        FlagComplexSpMatrix * mat_prev_coll  = new FlagComplexSpMatrix(number_of_points); 

        std::cout << "Going for vertex collapse and tower assembly" << std::endl;

        int i = 1;
        Map * redmap;
        while(edge_filt->size() > 0) {
           
            if(steps > 0)
                extract_sub_one_skeleton(threshold, *sub_skeleton, edge_t, *edge_filt);

            else
                extract_one_new_edge(*sub_skeleton, edge_t, *edge_filt);


            mat_coll = new FlagComplexSpMatrix(number_of_points, *sub_skeleton);
            // mat_coll->strong_vertex_edge_collapse();
            mat_coll->strong_vertex_collapse();
            // std::cout<< " Total number of current edges are " << sub_skeleton->size() << std::endl;
                       
            redmap = new Map();
            *redmap = mat_coll->reduction_map(); 
            
            // std::cout << "Subcomplex #" << i << " of threshold "<< threshold << " Collapsed" << std::endl;
            totAssembleTime += twr_assembler.build_tower_for_two_cmplxs(*mat_prev_coll, *mat_coll, *redmap, threshold, "./PersistenceOutput/CollapsedTowerRips.txt");
            // std::cout << "Tower updated for subcomplex #" << i << std::endl; 
            
            delete mat_prev_coll;
            mat_prev_coll = new FlagComplexSpMatrix();
            mat_prev_coll = mat_coll;
            mat_coll  = new FlagComplexSpMatrix();
            i++;
            delete redmap;
            if(steps > 0) {
                threshold = threshold+steps;
                if(threshold > end_threshold)
                    threshold = end_threshold; 
            }
            else {
                if(edge_filt->size() > 0)
                    threshold = edge_filt->at(0);  
            }    
        }
        std::cout << "Tower updated for subcomplex #" << i << std::endl; 
        auto the_collapse = std::chrono::high_resolution_clock::now();    
        sparse_distances = twr_assembler.distance_matrix();
        std::cout << "Collapse And assembly time : " <<  std::chrono::duration<double, std::milli>(the_collapse- the_collapse_begin).count() << " ms\n" << std::endl;
    }    

};
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
            for (int j = 0; j <= indx; j++) {
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
    double  begin_thresold;
    double  end_threshold;
    double  steps;
    int     repetetions = 1;
    char    manifold;

    Vector_of_points * point_vector;
    Vector_of_points file_all_points;

    std::string manifold_full = "sphere";
    
    double radius  = 1;
    double r_min  = 0.6;
    double r_max = 1;
    int dim_max  = 2;

    point_generator.program_options(argc, argv, number_of_points, begin_thresold, steps, end_threshold, repetetions, manifold, dimension, dim_max, in_file_name, out_file_name);
    
    std::cout << "The current input values to run the program is: "<< std::endl;
    std::cout << "number_of_points, begin_thresold, steps, end_threshold, repetetions, manifold, dimension, max_complex_dimension, in_file_name, out_file_name" << std::endl;
    std::cout << number_of_points << ", " << begin_thresold << ", " << steps << ", " << end_threshold << ", " << repetetions << ", " << manifold << ", " << dimension << ", " << dim_max << ", " << in_file_name << ", " << out_file_name << std::endl;
    
    if(manifold == 'f' || manifold =='F') {
        Gudhi::Points_off_reader<Point> off_reader(in_file_name);
        if (!off_reader.is_valid()) {
            std::cerr << "Unable to read file " << in_file_name << "\n";
            exit(-1);  // ----- >>
        }

        file_all_points = Vector_of_points(off_reader.get_point_cloud());
        dimension = file_all_points[0].dimension() ;
        std::cout << "Successfully read " << file_all_points.size() << " point_vector.\n";
        std::cout << "Ambient dimension is " << dimension << ".\n";
    }
   
    Map map_empty;

    std::string origFile ("./PersistenceOutput/original_tower_rips" );
    std::string collFile  ("./PersistenceOutput/collapsed_tower_rips") ;
    
    std::string filediag_bfr ("./PersistenceOutput/uncollapsed_persistence_diags") ;
    std::string filediag_aft ("./PersistenceOutput/collapsed_persistence_diags") ;
    
     std::string origPoints ("./PersistenceOutput/pointsamaple.off");
    // std::string otherStats ("./PersistenceOutput/maximal_simplx_cnt");
    // otherStats = otherStats+"_"+ out_file_name+ ".txt";
    filediag_bfr = filediag_bfr+"_"+ out_file_name+ ".txt";
    filediag_aft = filediag_aft+"_"+ out_file_name+ ".txt";

    double currentCreationTime = 0.0;
    double maxCreationTime     = 0.0;

   
    point_vector = new Vector_of_points();
    Distance_matrix distances;
    Distance_matrix *sparse_distances = new Distance_matrix();


    if(manifold == 's' || manifold == 'S'){
        // point_generator.generate_points_sphere(*point_vector, number_of_points, dimension, radius); 
        point_generator.generate_grid_2sphere(*point_vector, number_of_points, radius); 
        origFile = origFile+"_sphere_"+out_file_name+".txt";
        collFile = collFile+"_sphere_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< dimension <<"-sphere of radius " << radius << std::endl;
    }
    else if(manifold == 'b' || manifold == 'B'){
        point_generator.generate_points_ball(*point_vector, number_of_points, dimension, radius); 
        origFile = origFile+"_ball_"+out_file_name+".txt";
        collFile = collFile+"_ball_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< dimension <<"-ball of radius " << radius << std::endl;
    
    }
    else if( (manifold == 'a' || manifold == 'A')&& dimension == 2){
        point_generator.generate_points_2annulus(*point_vector, number_of_points, r_min, r_max); 
        origFile = origFile+"_annulus_"+out_file_name+".txt";
        collFile = collFile+"_annulus_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from "<< 2 <<"-annulus of radii (" << r_min << ',' << r_max << ") " << std::endl;
    }
    else if( (manifold == 'a' || manifold == 'A') && dimension == 3){
        point_generator.generate_points_spherical_shell(*point_vector, number_of_points, r_min, r_max); 
        origFile = origFile+"_annulus_"+out_file_name+".txt";
        collFile = collFile+"_annulus_"+out_file_name+".txt";
        std::cout << number_of_points << " points successfully chosen randomly from spherical shell of radii (" << r_min << ',' << r_max << ") " << std::endl;
    }
    
    else if(manifold == 'f' || manifold =='f') {
        // Subsampling from all points for each iterations
        Gudhi::subsampling::pick_n_random_points(file_all_points, number_of_points, std::back_inserter(*point_vector));
        origFile = origFile+"_"+ out_file_name+ ".txt";
        collFile = collFile+"_"+ out_file_name+ ".txt";
        std::cout << number_of_points << " points succesfully chosen randomly of dimension "<< dimension << " ." << std::endl;
    }
    else if (manifold == 'm'){
        std::string csv_file_name(in_file_name);
        distances = Gudhi::read_lower_triangular_matrix_from_csv_file<Filtration_value>(csv_file_name);
        number_of_points = distances.size();
        std::cout << "Read the distance matrix succesfully, of size: " << number_of_points << std::endl;
        origFile = origFile+"_"+ out_file_name+ ".txt";
        collFile = collFile+"_"+ out_file_name+ ".txt";
    }
    else {
        std::cerr << "Wrong parameters for input manifold..." <<std::endl;  
        exit(-1); 
    }
    
    if( manifold != 'm')
        number_of_points = point_vector->size();
    std::cout << "Point Set Generated."  <<std::endl;
 
    // for(int i = 0; i < number_of_points; i++ )
    //     point_generator.print_point(point_vector->at(i));
    // point_generator.output_points(*point_vector, origPoints);
   
    Filtered_sorted_edge_list * edge_t = new Filtered_sorted_edge_list();
    std::cout << "Computing the one-skeleton for threshold: " << end_threshold << std::endl; 
    
    auto begin_full_cmplx = std::chrono::high_resolution_clock::now();
    if(manifold == 'm') {                                           //Input is a distance 'm'atrix //Creating the edge list
        Rips_edge_list Rips_edge_list_from_file(distances, end_threshold);
        Rips_edge_list_from_file.create_edges(*edge_t);
        std::cout<< "Sorted edge list computed" << std::endl;

        //Creating the Rips Complex
        //Rips_complex rips_complex_from_file(distances, end_threshold);
        //rips_complex_from_file.create_complex(*subComplex, dim_max);
        //std::cout<< "Rips complex computed" << std::endl;
    }
    else{ //Point cloud input  //Creating the edge list
        Rips_edge_list Rips_edge_list_from_points(*point_vector, end_threshold, Gudhi::Euclidean_distance());
        Rips_edge_list_from_points.create_edges(*edge_t);
        std::cout<< "Sorted edge list computed" << std::endl;
        std::cout << "Total number of edges before collapse are: " << edge_t->size() << std::endl;

        // Creating the Rips Complex
        // Rips_complex rips_complex_after_collapse(*point_vector, end_threshold, Gudhi::Euclidean_distance());
        // rips_complex_from_points.create_complex(*subComplex, dim_max);
        // std::cout<< "Rips complex computed" << std::endl;
    }
    
    //Now we will perform filtered edge collapse to sparsify the edge list edge_t.
    std::cout<< "Filtered edge collapse begins" << std::endl;
    FlagComplexSpMatrix * mat_filt_edge_coll  = new FlagComplexSpMatrix(number_of_points,*edge_t,true);
    std::cout<< "Matrix instansiated" << std::endl;
    if(edge_t->size() >0){
        delete edge_t;
        edge_t = new Filtered_sorted_edge_list();
        *edge_t = mat_filt_edge_coll->filtered_edge_collapse(0); 
        filt_edge_to_dist_matrix(*sparse_distances, *edge_t, number_of_points); 
        std::cout << "Total number of vertices after collapse in the sparse matrix are: " << mat_filt_edge_coll->num_vertices() << std::endl;
    }
    else
    {
        std::cerr << "Total number of egdes are zero." <<std::endl;
        exit(-1) ;
    }
    

    // std::cout<< "Vertex strong-collapse begins" << std::endl;
    // strong_filtered_vertex_collapse(number_of_points, *edge_t, steps, *sparse_distances); 

    auto end_full_cmplx = std::chrono::high_resolution_clock::now();
    currentCreationTime = std::chrono::duration<double, std::milli>(end_full_cmplx - begin_full_cmplx).count();
    maxCreationTime = currentCreationTime;
   
 
    // Rips_complex rips_complex_before_collapse(distances, end_threshold);
    Rips_complex rips_complex_after_collapse(*sparse_distances, end_threshold);
     // Rips_complex rips_complex_after_collapse(*point_vector, end_threshold, Gudhi::Euclidean_distance());
    // Construct the Rips complex in a Simplex Tree
    
    Simplex_tree   simplex_tree_aft;
    // Simplex_tree   simplex_tree_bfr;
    // rips_complex_before_collapse.create_complex(simplex_tree_bfr, dim_max);
    rips_complex_after_collapse.create_complex(simplex_tree_aft, dim_max);

    // std::cout << "The complex contains " << simplex_tree_bfr.num_simplices() << " simplices before collapse. \n";
    // std::cout << "   and has dimension " << simplex_tree_bfr.dimension() << " \n";

    std::cout << "The complex contains " << simplex_tree_aft.num_simplices() << " simplices  after collapse. \n";
    std::cout << "   and has dimension " << simplex_tree_aft.dimension() << " \n";

    // Sort the simplices in the order of the filtration
    // simplex_tree_bfr.initialize_filtration();
    simplex_tree_aft.initialize_filtration();
    // Compute the persistence diagram of the complex
    // Persistent_cohomology pcoh_bfr(simplex_tree_bfr);
    Persistent_cohomology pcoh_aft(simplex_tree_aft);
    // initializes the coefficient field for homology
    // pcoh_bfr.init_coefficients(2);
    pcoh_aft.init_coefficients(2);

    // pcoh_bfr.compute_persistent_cohomology(steps);
    pcoh_aft.compute_persistent_cohomology(steps);
    // Output the diagram in filediag
    // if (filediag_bfr.empty()) {
    //     pcoh_bfr.output_diagram();
    // } 
    // else {
    //     std::ofstream out(filediag_bfr);
    //     pcoh_bfr.output_diagram(out);
    //     out.close();
    //   }

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
  