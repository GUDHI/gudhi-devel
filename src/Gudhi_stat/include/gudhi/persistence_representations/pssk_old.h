//This is an implementation of a kernel described in //http://arxiv.org/pdf/1412.6821v1.pdf
//by Reininghaus, Huber, Bauer and Kwitt.

#pragma once
#include <cmath>
#include <climits>

#include "persistence_intervals.h"
#include "kernel.h"

inline double norm_squared( double x , double y )
{
    return x*x + y*y;
};

template <typename T>
class pssk_kernel : public Persistence_kernel<T>
{
public:
    pssk_kernel( std::vector< Persistence_intervals<T> >& inter , double stdev_ ):Persistence_kernel<T>(inter.size()),stdev(stdev_)
    {
        this->intervals(inter);
    }
    pssk_kernel( const char* filename , double stdev_ = -1 )
    {
        //first, read the intervals:
        std::vector<std::string> names = readFileNames( filename );
        this->number_of_intervals = names.size();
        std::vector< Persistence_intervals<T> > intervals;
        this->intervals = intervals;
        for ( size_t file_no = 0 ; file_no != names.size() ; ++file_no )
        {
            cout << "Reading file : " << names[file_no] << endl;
            Persistence_intervals<T> interval( (char*)names[file_no].c_str() );
            this->intervals.push_back(interval);
            cerr << "Number of intervals : " << interval.size() << endl;

        }
        if ( stdev_ != -1 )
        {
            this->stdev = stdev_;
        }
        else
        {
            T min_ = INT_MAX;
            T max_ = INT_MIN;
            //compute range of the considered diagrams:
            for ( size_t i = 0 ; i != this->intervals.size() ; ++i )
            {
                std::pair< T , T > min_max = this->intervals[i].min_max();
                if ( min_ > min_max.first )min_ = min_max.first;
                if ( max_ < min_max.second )max_ = min_max.second;
            }
            //and arbitrary choice based on min and max:
            this->stdev = (max_ - min_)/100;
        }

    }
    double compute_scalar_product( size_t number_of_first_barcode , size_t number_of_second_barcode );
private:
    //since we do not construct explicit representation of the diagrams in the feature space, we do not store it here.
    double stdev;
    std::vector< Persistence_intervals<T> > intervals;
};

//First version, a naive implementation:
template <typename T>
double pssk_kernel<T>::compute_scalar_product( size_t number_of_first_barcode , size_t number_of_second_barcode )
{
    bool dbg = false;
    //cerr << "Entering \n";getchar();
    double result = 0;
    for ( size_t i = 0 ; i != this->intervals[number_of_first_barcode].intervals.size() ; ++i )
    {
        //cerr << "i :" << i <<endl;
        #pragma omp parallel for
        for ( size_t j = 0 ; j != this->intervals[number_of_second_barcode].intervals.size() ; ++j )
        {
            if ( dbg )
            {
                cerr << "( " << this->intervals[number_of_first_barcode].intervals[i].first << "," << this->intervals[number_of_first_barcode].intervals[i].second << ") \n";
                cerr << "( " << this->intervals[number_of_second_barcode].intervals[j].first << "," << this->intervals[number_of_second_barcode].intervals[j].second << ") \n";
                getchar();
            }
            T toAdd =
            exp( -( norm_squared( this->intervals[number_of_first_barcode].intervals[i].first-this->intervals[number_of_second_barcode].intervals[j].first, this->intervals[number_of_first_barcode].intervals[i].second- this->intervals[number_of_second_barcode].intervals[j].second ) ) / (8*this->stdev) )
            +
            exp( -( norm_squared( this->intervals[number_of_first_barcode].intervals[i].first-this->intervals[number_of_second_barcode].intervals[j].second, this->intervals[number_of_first_barcode].intervals[i].second- this->intervals[number_of_second_barcode].intervals[j].first ) ) / (8*this->stdev) );
            result += toAdd;
            if ( dbg )
            {
                cerr << "toAdd : " << toAdd << endl;
                cerr << result << " : " << result << endl;
                getchar();
            }
        }
    }
    result *= 1/( 8*this->stdev*3.141592 );
    //cerr << "result : " << result << endl; getchar();
    return result;
}
