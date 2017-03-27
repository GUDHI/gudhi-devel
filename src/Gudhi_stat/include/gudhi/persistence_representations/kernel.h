#pragma once
#include "persistence_intervals.h"
#include "read_files_names.h"

template <typename T>
class Persistence_kernel
{
public:
    Persistence_kernel(){}
    //Persistence_kernel( std::vector< Persistence_intervals<T> >& inter );
    //Persistence_kernel( const char* filename );
    Persistence_kernel(size_t number_of_intervals_):number_of_intervals(number_of_intervals_){}
	std::vector< std::vector<double> > compute_Gram_matrix();
	virtual double compute_scalar_product( size_t number_of_first_barcode , size_t number_of_second_barcode )=0;
	size_t& number(){return this->number_of_intervals;}
protected:
    //std::vector< Persistence_intervals<T> > intervals;
    size_t number_of_intervals;
};
/*
template <typename T>
Persistence_kernel<T>::Persistence_kernel( std::vector< Persistence_intervals<T> >& inter )
{
    this->intervals(inter);
}

template <typename T>
Persistence_kernel<T>::Persistence_kernel( const char* filename )
{
    std::vector<std::string> names = readFileNames( filename );
    for ( size_t file_no = 0 ; file_no != names.size() ; ++file_no )
    {
        cout << "Reading file : " << names[file_no] << endl;
        Persistence_intervals<double> interval( (char*)names[file_no].c_str() );
        this->intervals.push_back( interval );
    }
}
*/

template <typename T>
std::vector< std::vector<double> > Persistence_kernel<T>::compute_Gram_matrix()
{
    cerr << "compute_Gram_matrix \n";
    bool dbg = true;
    std::vector< std::vector<double> >  result( this->number_of_intervals );
    for ( size_t i = 0 ; i != this->number_of_intervals ; ++i )
    {
        result[i] = std::vector<double>( this->number_of_intervals );
    }
    for ( size_t i = 0 ; i != this->number_of_intervals ; ++i )
    {
        cerr << "Computing row number : " << i << endl;
        #pragma omp parallel for
        for ( size_t j = i ; j < this->number_of_intervals ; ++j )
        {
            cerr << j << " ";
            result[i][j] = result [j][i] = this->compute_scalar_product( i,j );
        }
    }
    return result;
}//compute_Gram_matrix

void pint_Gram_matrix_to_file( std::vector< std::vector<double> >& matrix , char* filename )
{
    ofstream out;
    out.open(filename);
    for ( size_t i = 0 ; i != matrix.size() ; ++i )
    {
        for ( size_t j = 0 ; j != matrix[i].size() ; ++j )
        {
            out << matrix[i][j] << " ";
        }
        out << std::endl;
    }
    out.close();
}
