#pragma once

#include <string>

#include "kernel.h"
#include "one_d_gaussians.h"


template <typename T>
class one_d_gaussians_kernel : public Persistence_kernel<T>
{
public:
    one_d_gaussians_kernel( std::vector< Persistence_intervals<T> >&  intervals_ , size_t number_of_points , double min_ = -1 , double max_ = -1 , size_t range = -1 );
    one_d_gaussians_kernel( char* filename , size_t number_of_points , double min_ = -1 , double max_ = -1 , size_t range = -1 );

    double compute_scalar_product( size_t number_of_first_barcode , size_t number_of_second_barcode );
private:
    std::vector< one_d_gaussians<T> > gaussians;
    void construct( std::vector< Persistence_intervals<T> >&  intervals , size_t number_of_points , double min_ = -1 , double max_ = -1 , size_t range = -1 );
};


template <typename T>
void one_d_gaussians_kernel<T>::construct( std::vector< Persistence_intervals<T> >&  intervals_ , size_t number_of_points , double min_ , double max_ , size_t range )
{
    for ( size_t nr = 0 ; nr != intervals_.size() ; ++nr )
    {
        std::vector< Persistence_intervals<T> > aa;
        aa.push_back( intervals_[nr] );
        cerr << "intervals_[nr].size() : " << intervals_[nr].size() << endl;
        cerr << "A \n";
        this->gaussians.push_back( one_d_gaussians<T>( aa , number_of_points , min_ , max_ , range ) );
        cerr << "B \n";
    }
}

template <typename T>
one_d_gaussians_kernel<T>::one_d_gaussians_kernel( std::vector< Persistence_intervals<T> >& intervals , size_t number_of_points , double min_ , double max_ , size_t range ):Persistence_kernel<T>(intervals.size())
{
    this->construct( intervals , number_of_points , min_ , max_ , range );
}

template <typename T>
one_d_gaussians_kernel<T>::one_d_gaussians_kernel( char* filename , size_t number_of_points , double min_ , double max_ , size_t range )
{
    std::vector<std::string> names = readFileNames( filename );
    std::vector< Persistence_intervals<double> > intervals;
    for ( size_t file_no = 0 ; file_no != names.size() ; ++file_no )
    {
        cout << "Reading file : " << names[file_no] << endl;
        Persistence_intervals<double> interval( (char*)names[file_no].c_str() );
        intervals.push_back( interval );
    }

    //in this case, if  min_ and max_ are not specified by the user, we should compute it, so that they are the same for all the considered cases:
    if ( min_ == max_ == -1 )
    {
        min_ = INT_MAX;
        max_ = INT_MIN;
        //in this case, the values of min_ and max_ are not specified by the user, and we should specify them here.
        for ( size_t diagNo = 0 ; diagNo != intervals.size() ; ++diagNo )
        {
            std::pair<T,T> mM = intervals[diagNo].min_max();
            if ( min_ > mM.first )min_ = mM.first;
            if ( max_ < mM.second )max_ = mM.second;
        }
    }
    this->construct( intervals , number_of_points , min_ , max_ , range );
    this->number_of_intervals = names.size();
}

template <typename T>
double one_d_gaussians_kernel<T>::compute_scalar_product( size_t number_of_first_barcode , size_t number_of_second_barcode )
{
    double result = 0;
    for ( size_t i = 0 ; i != this->gaussians[number_of_first_barcode].size() ; ++i )
    {
        result += this->gaussians[number_of_first_barcode][i] * this->gaussians[number_of_second_barcode][i];
    }
    //cerr << "Scalar product : " << result << endl; getchar();
    return result;
}
