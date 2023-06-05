/* This file is part of the Gudhi Library - https://gudhi.inria.fr/ -
*  which is released under MIT.
*  See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full
*  license details.
*  Author(s):       Hugo Passe
*
*  Copyright (C) Year Inria
*
*  Modification(s):
*    - YYYY/MM Author: Description of the modification
*/

#include <gudhi/Bitmap_cubical_complex.h>

#include <future>
#include <thread>
#include <chrono>

#ifndef EMBEDDED_CUBICAL_COMPLEX_H_
#define EMBEDDED_CUBICAL_COMPLEX_H_

template <typename T>
void print_vector(std::vector<T> vect){
    if(vect.size() == 0){
        std::cout << "[]\n";
    }else{
        std::cout << "[" << vect[0];
        for(std::size_t i = 1; i < vect.size(); i++){
            std::cout << ", " << vect[i];
        }
        std::cout << "]\n";
    }
}

template <typename T>
class Embedded_cubical_complex : public Gudhi::cubical_complex::Bitmap_cubical_complex<T>
{

    private:

        typedef std::size_t Simplex_key;
        typedef typename T::filtration_type Filtration_value;
        typedef Simplex_key Simplex_handle;

    protected:
        std::vector<std::vector<double>> embedding;         //Array of the cubical complexe's points' coordinates
        std::vector<int> embedding_index;                   //Array to link vertices index in the cubical complex to their index in the embedding

        std::vector<int> sizes_pdt;                         //Products of the sizes from s_0 up to s_0*s_1*...*s_(n-2)

        std::vector<std::vector<int>> critical_vertices;
        std::vector<std::vector<int>> critical_multiplicity;

        //If we have a vector e = (e_1,...,e_n), its index is the sum of the 2^i were i are the indexes such that e_i >= 0
        //If index_v = index_w they have the same critical points. The critical points are stored in critical_vertices[index]

    public:
        //*********************************************//
        //Constructors
        //*********************************************//
        Embedded_cubical_complex(const std::vector<unsigned>& dimensions,
            const std::vector<Filtration_value>& top_dimensional_cells, unsigned int num_jobs = 0):Gudhi::cubical_complex::Bitmap_cubical_complex<T>(dimensions, top_dimensional_cells)
            {
                sizes_pdt.push_back(2*this->sizes[0]+1);
                //In this loop we compute the product of the number of cubes in each direction, it optimizes the complexity of the get_coordinates_in_complex function
                for(Simplex_handle i = 1; i < this->dimension()-1; i++){
                    sizes_pdt.push_back(sizes_pdt[i-1]*(2*this->sizes[i]+1));
                }

                std::cout << "Init embedding ...\n";
                initalize_embedding();
                initalize_embedding_index();

                std::cout << "Computing critical vertices\n";
                if(num_jobs > std::thread::hardware_concurrency() || num_jobs <= 0){
                    num_jobs = std::thread::hardware_concurrency();
                }
                compute_critical_vertices(num_jobs);
                std::cout << "Cubical complex successfully created\n\n";
            }   

        //*********************************************//
        //Functions for pretratement
        //*********************************************//

        void initalize_embedding(){
            for(Simplex_handle i = 0; i < this->num_simplices();i++){
                if(this->dimension(i) == 0){
                    std::vector<int> coords = get_coordinates_in_complex(i);
                    std::vector<double> embedded_coords(coords.size());
                    for(Simplex_handle j = 0; j < coords.size(); j++){
                        embedded_coords[j] = coords[j] / 2./ this->sizes[j] - 0.5;
                    }
                    embedding.push_back(embedded_coords);
                }
            }
        }

        void initalize_embedding_index(){
            int index = 0;
            for(Simplex_handle handle = 0; handle < this->num_simplices(); handle++){
                if(this->dimension(handle) == 0){
                    embedding_index.push_back(index);
                     index++;
                }else{
                    embedding_index.push_back(-1);
                }
            }
        }

        //*********************************************//
        //Functions for critical points
        //*********************************************//
        void compute_critical_vertices(int num_jobs){
            int dim = this->dimension();
            Simplex_key n = this->num_simplices();
            std::chrono::seconds zero_sec{0};

            int num_vertex = 1;
            for(std::size_t  i = 0; i < this->sizes.size(); i++){
                num_vertex *= this->sizes[i] + 1;
            }

            std::vector<int> direction(dim,-1);
            long unsigned int index = 0;

            for(int i = 0; i < (1 << dim); i++){    //Loop on every possible critical direction 

                std::vector<int> tmp;
                critical_vertices.push_back(tmp);
                std::vector<int> tmp1;
                critical_multiplicity.push_back(tmp1);

                std::vector<std::promise<std::vector<std::vector<int>>>> promise_vect;
                std::vector<std::future<std::vector<std::vector<int>>>> future_vect;
                std::vector<std::thread> thread_vector;

                for(int i = 0; i < num_jobs; i++){
                    std::promise<int> promiseObj;
            
                    promise_vect.push_back(std::promise<std::vector<std::vector<int>>>());
                    future_vect.push_back(promise_vect[i].get_future());
                    //Objects to get return values from the thread
                    
                    thread_vector.push_back(std::thread(&Embedded_cubical_complex::compute_critical_vertices_subroutine, this, std::move(promise_vect[i]), direction, dim, n, i, num_jobs));
                    thread_vector[i].detach();  //Thread is now running concurently
                }

                int job = 0;
                while(job < num_jobs){
                    if(future_vect[job].wait_for(zero_sec) == std::future_status::ready){
                        std::vector<std::vector<int>> thread_res = future_vect[job].get();
                        critical_vertices[index].insert(critical_vertices[index].end(), thread_res[0].begin(), thread_res[0].end());
                        critical_multiplicity[index].insert(critical_multiplicity[index].end(), thread_res[1].begin(), thread_res[1].end());
                        job++;
                    }
                }

                //Finding next direction vector
                int k=0;
                while(direction[k] == 1 && k < dim){
                    direction[k] = -1;
                    index -= (1u << k);
                    k++;
                }
                direction[k] = 1;
                index += (1u << k);
            }
        }

        //This function is not yet finished, it will be used as a parallel subroutine of the previous one.
        //Subroutine, we parallelize on the first coordinate (potentially not a good choice on a n*m complex where n < num_threads and m >> n)
        void compute_critical_vertices_subroutine(std::promise<std::vector<std::vector<int>>> promiseObj, std::vector<int> direction, int dim, Simplex_key n_simplices, int job_index, int num_jobs){
            std::vector<int> crit_subroutine;
            std::vector<int> mult_subroutine;
            
            std::vector<int> coords(dim,0);     //On each pass throught the loop, it will store the coordinates of the simplex
            coords[0] = 2*job_index;
            Simplex_key vertex = get_key_from_coordinates(coords);   //Will store the index of the vertex

            while(vertex < n_simplices){  //Loop on all vertices
                int multiplicity = is_vertex_critical(vertex,direction);    //Computing critical point multiplicity in direction
                if(multiplicity != 0){                                      //If relevant we add it to the critical points
                    crit_subroutine.push_back(vertex);     
                    mult_subroutine.push_back(multiplicity);
                }

                //Finding next vertex
                coords[0] = coords[0] + 2*num_jobs;
                for(int j = 0; j < dim-1; j++){  
                    if((unsigned)coords[j] > 2*this->sizes[j]+1){
                        if(j == 0){
                            coords[0] = 2*job_index;
                        }else{
                            coords[j] = 0;
                        }
                        coords[j+1] = coords[j+1] + 2;
                    }else{
                        break;
                    }
                }
                vertex = get_key_from_coordinates(coords);
            }
            std::vector<std::vector<int>> res;
            res.push_back(crit_subroutine);
            res.push_back(mult_subroutine);
            promiseObj.set_value(res);
        }
        
        //Returns 0 if the vertex isn't critial, its multiplicity if it is, works on 0-1 complexes.
        int is_vertex_critical(Simplex_handle vertex, std::vector<int> direction){
            int euler_1 = compute_euler_car_in_direction(vertex,direction,1);
            int euler_2 = compute_euler_car_in_direction(vertex,direction,-1);

            return euler_1 - euler_2;
        }

        //Return euler car with multiplicity of the intersection between the cells in the neigbourhood of the vertex and the hyperplane orthogonal to direction in the neighbourhood of the vertex
        int compute_euler_car_in_direction(Simplex_handle vertex, std::vector<int> direction, int reverse_vector){
            int euler_car = 0;
            int dim = direction.size();

            std::vector<int> coordinates = get_coordinates_in_complex(vertex); //This vector will successively take the coordinates of adjacent cells involved in the Euler's caracteristic calculation
            
            std::vector<int> tmp(dim);  //This vector will help us to find all the adjacent cells involved in calculations
            int simplex_dim_sign = 1;
            
            for(int i = 0; i < (1 << dim)-1; i++){  //Looping on all of the adjacent cell in direction
                int k = 0;
                while(tmp[k] == 1){             //Finding the next adjacent cell 
                    tmp[k] = 0;
                    coordinates[k] -= reverse_vector * direction[k];
                    simplex_dim_sign *= -1;     //Gives us the sign to put in front (e.g : edges must be taken positively were faces must be taken negatively (it is the oppostie to the formulae due to the intersection with the hyperplane that transforms n dimensional cells in n-1 dimensional cells))
                    k++;
                }
                coordinates[k] += reverse_vector * direction[k];

                if(k < dim){
                    tmp[k] = 1;
                    simplex_dim_sign *= -1;
                }

                if(are_coordinates_in_complex(coordinates) == 1){   //If the cell exists, adding a term to the caracteristic
                    Simplex_key key = get_key_from_coordinates(coordinates);
                    euler_car += this->filtration(key) * simplex_dim_sign;
                }
            }

            return euler_car;
        }


        //*********************************************//
        //Printing functions
        //*********************************************//

        void print_embedding(){
            std::cout << "[";
            for(Simplex_handle i = 0; i < embedding_index.size(); i++){
                if(embedding_index[i] != -1){
                    print_vector(embedding[embedding_index[i]]);
                }
            }
            std::cout << "]\n";
        }

        void print_critical_vertices(){
            std::cout << "Critical vertices test : \n";
            for(std::size_t i = 0; i < critical_vertices.size(); i++){
                print_vector(critical_vertices[(int)i]);
            }
            std::cout << "\n";
        }

        void print_critical_multiplicity(){
            std::cout << "Critical multiplicity : \n";
            for(std::size_t i = 0; i < critical_multiplicity.size(); i++){
                print_vector(critical_multiplicity[(int)i]);
            }
            std::cout << "\n";
        }

        void print_filtration(){
            std::cout << "Filtration : \n[";
            int n = this->num_simplices()-1;
            for(int i = 0; i < 2*this->sizes[1]+1; i++){
                for(int j = 0; j < 2*this->sizes[0]+1; j++){
                    std::cout << this->filtration(j*(2*this->sizes[0]+1)+2*this->sizes[1]-i) << ", ";
                }
                std::cout << "\n";
            }
            std::cout << "]\n";
        }

        //*********************************************//
        //Those are arithmetic functions to find the index of the vertex within the vertices
        //*********************************************//

        //This function gives the coordinates of the cell given by key
        //The coordinates are the cartesians' ones where the complex is seen as a cartesian coordinate system
        std::vector<int> get_coordinates_in_complex(Simplex_key key){

            int n = (int)this->dimension();
            std::vector<int> coordinates;      //Coordinates of the face indexed by key to return
            
            for(int i = n-1; i > 0; i--){
                coordinates.insert(coordinates.begin(),key / sizes_pdt[i-1]);
                key = key - coordinates[0]*sizes_pdt[i-1];
            }

            coordinates.insert(coordinates.begin(),key);

            return coordinates;
        }

        //This functions gives you the key of the vertex given the coordinates of it
        //The opposite operation than the previous function
        Simplex_key get_key_from_coordinates(std::vector<int> coordinates){
            Simplex_key key = 0;
            for(int i = this->dimension()-1; i >= 0; i--){
                key = key*(2*this->sizes[i] + 1) + coordinates[i];
            }
            return key;
        }

        //This function returns a vector with the keys of the vertices of the cell given by key
        std::vector<int> get_cell_vertices(Simplex_key key){
            std::vector<int> cell_coordinates = get_coordinates_in_complex(key);    //Get the coordinates of cel indexed by key
            int n = cell_coordinates.size();

            std::vector<int> odd_coordinates;
            std::vector<int> cell_vertex(n);                //The first identified vertex 
            int n_odd_coords = 0;                           //Gives us the dimension of the cell, as well as the number of vertices in it
            
            for(int i = 0; i < n; i++){                         //Computing the number and indexes of odd coordinates in vector cell_coordinates and finding the first vertex
                if(cell_coordinates[i]%2 == 1){
                    odd_coordinates.push_back(i);
                    cell_vertex[i] = cell_coordinates[i]-1;
                    n_odd_coords++;
                }else{
                    cell_vertex[i] = cell_coordinates[i];
                }
            }

            std::vector<int> cell_vertices;

            if(n_odd_coords == 0){  //If key is a vertex we return key
                cell_vertices.push_back(key);
                return cell_vertices;  
            }

            std::vector<int> tmp_vect(n_odd_coords);

            cell_vertices.push_back(get_key_from_coordinates(cell_vertex));

            for(int i = 0; i < (1 << n_odd_coords)-1;i++){
                for(int j = 0; j < n_odd_coords;j++){                   //We use a binary counter on the odd cordinates of the simplex to compute all of the vertices coordinates
                    if(tmp_vect[j] == 1){                           //Cell is vertex if and only if it all of its coordinates are even
                        tmp_vect[j] = 0;
                        cell_vertex[odd_coordinates[j]] = cell_vertex[odd_coordinates[j]] - 2;
                    }else{
                        tmp_vect[j] = 1;
                        cell_vertex[odd_coordinates[j]] = cell_vertex[odd_coordinates[j]] + 2;
                        break;
                    }
                }
                
                cell_vertices.push_back(get_key_from_coordinates(cell_vertex));     //Adding the key of the vertex we found in cell_vertices
            }

            return cell_vertices;
        }

        //Given a direction vector e, return the index of the subvector that contains critical points in direction e
        int get_vector_index(std::vector<double> e){
            int index = 0;
            int muliplier = 1;

            for(int i = 0; i < (int)e.size(); i++){
                if(e[i] >= 0){
                    index += muliplier;
                }
                muliplier = muliplier*2;
            }

            return index;
        }

        //Check if given coordinates are in complex, useful when computing the adjacent cells indexes to compute the euler caracteristic
        int are_coordinates_in_complex(std::vector<int> coordinates){
            std::size_t coord_dim = coordinates.size();
            if(coord_dim != this->sizes.size()){
                return 0;
            }

            for(std::size_t i = 0; i < coord_dim; i++){
                if(coordinates[i] < 0 || coordinates[i] > 2*((int)this->sizes[i])){
                    return 0;
                }
            }
            return 1;
        }

        //*********************************************//
        //Functions to compute hybrid transform, needs a functions to compute kernel's antiderivative (abusevely called kernel) and a direction vector e
        //*********************************************//
        
        double compute_hybrid_transform(double (*kernel)(double), std::vector<double> e){
            int index = get_vector_index(e);
            int reverse_vector = 1;

            //As multiplicity values are stored assuming that the last coordinate of direction is > 0, potentialy change index and set reverse_vector to -1
            if(index >= (int)critical_vertices.size()){
                reverse_vector = -1;
                index = (1 << e.size()) - 1 - index;
            }
            
            double sum = 0.0;

            for(std::size_t i = 0; i < critical_vertices[index].size(); i++){   //Looping on critical vertices
                sum +=  critical_multiplicity[index][i] * kernel(std::inner_product(e.begin(),e.end(),embedding[embedding_index[critical_vertices[index][i]]].begin(),0.0));
            }
            return reverse_vector * sum;
        }

        //An overload of previous function to support multithreading
        std::vector<double> compute_hybrid_transform(double (*kernel)(double), std::vector<std::vector<double>> vect_list, unsigned int num_threads = -1){
            std::vector<double> results;
            std::chrono::seconds zero_sec{0};

            std::size_t num_vectors = vect_list.size();
            
            if(num_threads > std::thread::hardware_concurrency() || num_threads <= 0){
                num_threads = std::thread::hardware_concurrency();
            }

            int step = (int)num_vectors / (int)num_threads;
            
            std::vector<std::promise<std::vector<double>>> promise_vect;
            std::vector<std::future<std::vector<double>>> future_vect;
            std::vector<std::thread> thread_vector;
            
            for(std::size_t i = 0; i < num_threads; i++){   //We create threads
            
                promise_vect.push_back(std::promise<std::vector<double>>());
                future_vect.push_back(promise_vect[i].get_future());
                //Objects to get return values from the thread

                int begin = i*step;
                int end = (i+1)*step;

                if(i == num_threads-1){
                    end = (int)num_vectors;
                }
                
                thread_vector.push_back(std::thread(&Embedded_cubical_complex::compute_hybrid_transform_subvector, this, std::move(promise_vect[i]), kernel, vect_list, begin, end));
                thread_vector[i].detach();  //Thread is now running concurently
            }

            int b = 1;
            while(b){
                b = 0;
                for(std::size_t i = 0; i < num_threads; i++){       //Waiting for all the threads to finish their job
                    if(future_vect[i].wait_for(zero_sec) != std::future_status::ready){
                        b = 1;
                        break;
                    }
                }
            }
            
            for(std::size_t i = 0; i < num_threads; i++){   //Merging answers in one vector
                std::vector<double> thread_res = future_vect[i].get();
                results.insert(results.end(),thread_res.begin(),thread_res.end());
            }
            
            return results;
        }

        //Computing multiple transforms on one kernel, used by previous function
        void compute_hybrid_transform_subvector(std::promise<std::vector<double>> promiseObj, double (*kernel)(double), std::vector<std::vector<double>> vect_list, std::size_t begin_index, std::size_t end_index){
            std::vector<double> results;
            for(std::size_t i = begin_index; i < end_index; i++){
                results.push_back(compute_hybrid_transform(kernel,vect_list[i]));
            }
            promiseObj.set_value(results);
        }
};

#endif