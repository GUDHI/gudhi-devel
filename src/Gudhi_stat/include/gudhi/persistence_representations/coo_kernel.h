//this is an implementation of a kernel by Mathieu Carrière,
//Steve Y. Oudot, and Maksim Ovsjanikov presended in Stable Topological
//Signatures for Points on 3D Shapes. Proc. Sympos. on Geometry Processing, 2015

#pragma once
#include <cmath>
#include <climits>

#include "persistence_intervals.h"
#include "kernel.h"
#include "Distances_of_points_in_diagram.h"

template <typename T , typename F>
class coo_kernel : public Persistence_kernel<T>
{
public:
    //Distances_of_points_in_diagram( std::vector< std::pair< T , T > > intervals , F f , size_t where_to_cut );
    //Distances_of_points_in_diagram( char* filename , F f , size_t where_to_cut );

    coo_kernel( std::vector< Persistence_intervals<T> >& intervals , F f , size_t where_to_cut = 100 );
    coo_kernel( char* filename , F f , size_t where_to_cut = 100 );

    //i tutaj pojawia sie problem, bo reprezentacja jest liczona wielokrotnie, strata czasu.
    double compute_scalar_product( size_t number_of_first_barcode , size_t number_of_second_barcode );
private:
    //in this case, we will construct the representation of diagrams in the vector space. That is why, we store here some more stuff compared to
    //the Persistence_kernel class.
    std::vector< Distances_of_points_in_diagram<T,F> > representation;
    void constuct_representation( std::vector< Persistence_intervals<T> >& intervals , F f , size_t where_to_cut);
};

template <typename T , typename F>
coo_kernel<T,F>::coo_kernel( std::vector< Persistence_intervals<T> >& intervals , F f , size_t where_to_cut ):Persistence_kernel<T>(intervals)
{
    this->constuct_representation(intervals, f , where_to_cut);
}

template <typename T , typename F>
coo_kernel<T,F>::coo_kernel( char* filename , F f , size_t where_to_cut )
{
    cerr << "Constuct representation \n";
    std::vector< Persistence_intervals<T> > intervals;
    std::vector<std::string> names = readFileNames( filename );
    for ( size_t file_no = 0 ; file_no != names.size() ; ++file_no )
    {
        cout << "Reading file : " << names[file_no] << endl;
        Persistence_intervals<T> interval( (char*)names[file_no].c_str() );
        intervals.push_back( interval );
    }
    this->constuct_representation(intervals, f , where_to_cut);
}

template <typename T , typename F>
void coo_kernel<T,F>::constuct_representation( std::vector< Persistence_intervals<T> >& intervals , F f , size_t where_to_cut)
{
    std::vector< Distances_of_points_in_diagram<T,F> > representation( intervals.size() );
    this->representation = representation;
    #pragma omp parallel for
    for ( size_t i = 0 ; i < intervals.size() ; ++i )
    {
        this->representation[i] = Distances_of_points_in_diagram<T,F>( intervals[i].intervals , f , where_to_cut );

    }
    this->number_of_intervals = intervals.size();
}

template <typename T , typename F>
double coo_kernel<T,F>::compute_scalar_product( size_t number_of_first_barcode , size_t number_of_second_barcode )
{
    double result = 0;
    size_t range = std::min( this->representation[number_of_first_barcode].size() , this->representation[number_of_second_barcode].size() );
    for ( size_t i = 0 ; i != range ; ++i )
    {
        result += this->representation[number_of_first_barcode].vector_in_position(i)*this->representation[number_of_second_barcode].vector_in_position(i);
    }
    return result;
}


