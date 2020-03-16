#pragma once
#include <gudhi/FlagComplexSpMatrix.h>
#include <gudhi/Rips_edge_list.h>

#include <set>
#include <fstream>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>


typedef std::size_t Vertex;
using Edge                  = std::pair<Vertex,Vertex>;
using edge_list             = std::vector<Edge>;
using Simplex               = std::vector<Vertex>;

using vectorVertex          = std::vector<Vertex>;
using vert_unSet            = std::unordered_set<Vertex>;
using Map                   = std::unordered_map<Vertex,Vertex>;
using Distance_matrix       = std::vector<std::vector<double>>;
// using infinity              = std::numeric_limits<double>::max();
// assumptions : (1) K1 and K2 have the same vertex set
//               (2) The set of simplices of K1 is a subset of set of simplices of K2
// K1  ->  K2    [Original Simplicial Complexes]
// |       |
// |       |
// K1c ->  K2c   [Strongly Collapsed Flag Complexes]

class TowerAssembler_FlagComplex
{
  private:
	Map renamedVertices; 
    Map representative_map;
	std::size_t current_rename_counter;
    Distance_matrix * distance_mat;
    size_t total_vertices;

    // Filtered_sorted_edge_list * edge_t = new Filtered_sorted_edge_list();
    FlagComplexSpMatrix * flag_Filtration; 
    typedef std::vector< std::tuple< double, Vertex, Vertex > > Filtered_sorted_edge_list;

	
  public:
    
    TowerAssembler_FlagComplex(std::size_t numVert)
    {
    	for (std::size_t i = 0; i <= numVert; ++i){
    		renamedVertices[i] = i;
        }
    	total_vertices = 0;
        current_rename_counter = numVert+1;
        
        flag_Filtration = new FlagComplexSpMatrix((numVert*log2(numVert)) +1);
        distance_mat    = new Distance_matrix();
        // distance_mat->push_back({0});
    }
    
    ~TowerAssembler_FlagComplex(){};
    double build_tower_for_two_cmplxs(FlagComplexSpMatrix mat_1, const FlagComplexSpMatrix & mat_2,  Map redmap_2,  double filtration_value, std::string outFile) // mat_1 and mat_2 are simplex_trees of K1c and K2c (the collapsed ones), redmap_2 is the map of K2 -> K2c
    {
        auto begin_print  = std::chrono::high_resolution_clock::now();
        std::ofstream myfile (outFile, std::ios::app);
        if (myfile.is_open())
        {   
            for (auto & v : mat_1.vertex_set()) {
                auto collapsed_to = redmap_2.find(v);  // If v collapsed to something?
                if(collapsed_to != redmap_2.end()) {  // Collapse happened, because there is a vertex in the map
                    if(mat_1.membership(collapsed_to->second)) { // Collapsed to an existing vertex in mat_1.

                    	// myfile << filtration_value  << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
                    	std::cout << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl;
                        flag_Filtration->active_strong_expansion(renamedVertices.at(v), renamedVertices.at(collapsed_to->second), filtration_value);
                        renamedVertices.at(v) = current_rename_counter;
                        current_rename_counter++;                        
                    }
                    else {
	                    // myfile << filtration_value << " i " << renamedVertices.at(collapsed_to->second) << std::endl;
	                    // myfile  << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
	                    flag_Filtration->active_strong_expansion(renamedVertices.at(v), renamedVertices.at(collapsed_to->second),filtration_value);
                        std::cout << filtration_value << " i " << renamedVertices.at(collapsed_to->second) << std::endl;
                        std::cout  << filtration_value << " c " << renamedVertices.at(v) << " " << renamedVertices.at(collapsed_to->second) << std::endl; 
                        renamedVertices.at(v) = current_rename_counter;
	                    current_rename_counter++;
                    }
                     // If the vertex "collapsed_to->second" is not a member of mat_1, the contraction function will simply add and then collapse
                     mat_1.contraction(v, collapsed_to->second);
                     // std::cout << " Contraction Done " << std::endl;

                }
            }

            //The core K1c (mat_1) has gone through the transformation(re-labeling)/collapse and it is now a subcomplex of K2c, the remaining simplices need to be included
            // Writing the inclusion of all remaining simplices...
            // std::cout << "Begining the inclusion of edges " << std::endl;
            for( const Edge & e  : mat_2.all_edges()) {
                auto u = std::get<0>(e);
                auto v = std::get<1>(e);
                // std::cout << "Going to insert the edge :" << renamedVertices.at(u) << ", " << renamedVertices.at(v) << std::endl;
                // std::cout << "Going to insert the vertex :" << u << std::endl;
                if(!mat_1.membership(u)) {
                    flag_Filtration->insert_vertex(renamedVertices.at(u),filtration_value);
                    // std::cout << "Inserted the vertex :" << renamedVertices.at(u) <<  " in the new distance matrix"<< std::endl;
                    // myfile << filtration_value << " i";
                    // myfile  << " " << renamedVertices.at(u);
                    // myfile  << std::endl;
                    mat_1.insert_vertex(u,1);
                    // std::cout << "Inserted the vertex :" << renamedVertices.at(u) <<  " in the old skeleton matrix"<< std::endl;


                }

                // std::cout << "Going to insert the vertex :" << v << std::endl ;
                if(!mat_1.membership(v)) {
                    // std::cout << "Begining the insertion the vertex :" << renamedVertices.at(v) <<  " in the new distance matrix"<< std::endl;

                    flag_Filtration->insert_vertex(renamedVertices.at(v),filtration_value);
                    // std::cout << "Inserted the vertex :" << renamedVertices.at(v) <<  " in the new distance matrix"<< std::endl;

                    // myfile << filtration_value << " i";
                    // myfile  << " " << renamedVertices.at(v);
                    // myfile  << std::endl;
                    mat_1.insert_vertex(v,1);
                    // std::cout << "Inserted the vertex :" << v <<  " in the old skeleton matrix"<< std::endl;

                }
                // std::cout << "Going to insert the edge :" << u << ", " << v << std::endl;
                if(!mat_1.membership(e)){
                    flag_Filtration->insert_new_edges(renamedVertices.at(u),renamedVertices.at(v), filtration_value);
                    // std::cout << "Inserted the edge :" << renamedVertices.at(u) << ","<< renamedVertices.at(v) <<  " in the new distance matrix"<< std::endl;

                    // myfile << filtration_value << " i";
                    // myfile  << " " <<  renamedVertices.at(u) << ", " << renamedVertices.at(v);
                    // myfile  << std::endl;
                    // std::cout << "Going to insert the edge :" << u << ","<< v <<  " in the old skeleton"<< std::endl;
                    mat_1.insert_new_edges(u,v,1);
                    // std::cout << "Inserted the edge :" << u << ","<< v <<  " in the old skeleton"<< std::endl;
                } 
                 // std::cout << " Insertion Done " << std::endl;        
            }

            // myfile << "# Tower updated for the additional subcomplex.\n";
            myfile.close();
        }
        else {
            std::cerr << "Unable to open file";
            exit(-1) ;
        }
        
        auto end_print  = std::chrono::high_resolution_clock::now();
        auto printTime = std::chrono::duration<double, std::milli>(end_print- begin_print).count();
        // std::cout << " Time to print the tower : " << printTime << " ms\n" << std::endl;
        return printTime;
    } 
    Distance_matrix distance_matrix()
    {
        size_t non_zero_rw = flag_Filtration->num_vertices();
        double inf = std::numeric_limits<double>::max();
        sparseRowMatrix mat = flag_Filtration->uncollapsed_matrix();
        doubleVector distances ;  
        for(size_t indx = 0; indx < non_zero_rw; indx++){ 
            rowInnerIterator it(mat, indx);   
            for(size_t j = 0; j <= indx; j++) {             // Iterate over the non-zero columns
                //std::cout << "j = " << j << " - it.index = " << it.index() << " - indx = " << indx << std::endl;
                if(it.index() == j && j != indx){
                    distances.push_back(it.value());  // inner index, here it is equal to it.columns()
                    ++it;
                }
                else if( j == indx)
                    distances.push_back(0);
                else
                     distances.push_back(inf);

            }
            distance_mat->push_back(distances);
            distances.clear();
        }        
        return *distance_mat;
    }      
   void print_sparse_matrix(){
        flag_Filtration->print_sparse_skeleton();
   }
    
};
