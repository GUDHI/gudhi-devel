#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>



#include "persistence_intervals.h"
#include "Distances_of_points_in_diagram.h"
#include "read_files_names.h"




using namespace std;


class normal_dist
{
public:
    normal_dist(double mean , double stdev):stdev(stdev),mean(mean){}
    double operator() (double x)
    {
        //compute value of the normal distribution at x:
        return 1.0/(this->stdev*2.506628)*exp(-(x-this->mean)*(x-this->mean)/(2*this->stdev*this->stdev));
    }
private:
    double mean;
    double stdev;
};


template <typename T>
class one_d_gaussians
{
public:
    one_d_gaussians( std::vector< Persistence_intervals<T> > intervals , size_t number_of_points , double min_ = -1 , double max_ = -1 , size_t range = -1 );
    void print_to_file( char* filename  );
    size_t size()const{return this->gaussians.size();}
    inline T operator [](size_t i)
    {
        if ( i >= this->gaussians.size() )throw("Index out of range! Operator [], one_d_gaussians class\n");
        return this->gaussians[i].second;
    }
protected:
    std::vector< std::pair<T,T> > gaussians;
    size_t number_of_points;
    T min_;
    T max_;
    size_t range;
};


template <typename T>
one_d_gaussians<T>::one_d_gaussians( std::vector< Persistence_intervals<T> > intervals , size_t number_of_points , double min_ , double max_ , size_t range )
{
    bool dbg = true;

    this->number_of_points = number_of_points;
    //setting up min max and range:
    if ( range == -1 )
    {
        this->range = number_of_points/10.00;
    }
    else
    {
        this->range = range;
    }

    if ( (min_ == -1) && (max_ == -1) )
    {
        if ( dbg )cerr << "Setting up min and max \n";
        //in this case, we need to find min and max:
        min_ = INT_MAX;
        max_ = -INT_MAX;
        for ( size_t i = 0 ; i != intervals.size() ; ++i )
        {
            for ( size_t j = 0 ; j != intervals[i].intervals.size() ; ++j )
            {
                if ( intervals[i].intervals[j].first < min_ )min_ = intervals[i].intervals[j].first;
                if (intervals[i].intervals[j].first > max_ )max_ = intervals[i].intervals[j].first;
            }
        }
    }

    this->min_= min_;
    this->max_ = max_;

    if ( this->min_ == this->max_ )return;

    if ( dbg )
    {
        cerr << "min_ : " << this->min_ << endl;
        cerr << "max_ : " << this->max_ << endl;
        cerr << "number_of_points : " << this->number_of_points << endl;
        cerr << "range : " << this->range << endl;

    }

    //initialization of the output variable:
    std::vector< std::pair<double,double> > result( this->number_of_points+1 );
    T dx = (this->max_-this->min_)/(T)this->number_of_points;
    T current = this->min_;
    for ( size_t i = 0 ; i != this->number_of_points ; ++i )
    {
        result[i] = std::make_pair( current , 0 );
        current += dx;
    }

    //for every diagram
    for ( size_t diag_no = 0 ; diag_no != intervals.size() ; ++diag_no )
    {
        if ( dbg )
        {
            cerr << "We are considering diagram number : " << diag_no << endl;
        }
        //for every point in this diagram:
        for ( size_t point_no = 0 ; point_no != intervals[diag_no].intervals.size() ; ++point_no )
        {
            if (dbg){cerr << "(intervals[diag_no].intervals[point_no].first-min_) : " << (intervals[diag_no].intervals[point_no].first-min_) << "\n";}

            size_t position = ((intervals[diag_no].intervals[point_no].first-this->min_)/( this->max_ - this->min_ )) * this->number_of_points;
            //now we need to draw a gaussian distribution 'range' value to left, and to right from this point

            if ( dbg )
            {
                cerr << "point : (" << intervals[diag_no].intervals[point_no].first << "," << intervals[diag_no].intervals[point_no].second << "), its position : " << position << endl;
            }
            normal_dist d(intervals[diag_no].intervals[point_no].first,dx*this->range);
            double length = ( intervals[diag_no].intervals[point_no].second - intervals[diag_no].intervals[point_no].first );

            if ( position > this->range )
            {
                cerr << "Position of this element extends the current range. This point will not be taken into consideration \n";
                continue;
            }

            for ( size_t i = 0 ; i != this->range ; ++i )
            {
                //cerr << "i : " << i << " range : " << this->range << endl;
                if ( position+i < number_of_points )
                {
                    result[position+i].second += length*d( min_+(position+i)*dx );
                    //cerr << "position+i : " << position+i << endl;
                }
                if ( position >= i )
                {
                    result[position-i].second += length*d( min_+(position-i)*dx );
                    //cerr << "position-i : " << position-i << endl;
                }
            }

            //if ( dbg )cerr << "Done \n";
        }
    }

    //weighting, normalization. This is needed, othervise the result depends on the number of diagram involved, which should not be the case for averages.
    for ( size_t i = 0 ; i != result.size() ; ++i )
    {
        //cerr << "result[i].second  : " << result[i].second  << endl;
        std::pair<double,double> newResult = std::make_pair( result[i].first , result[i].second / (double)intervals.size() );
        result[i] = newResult;
        this->gaussians.push_back(newResult);
        //cerr << "result[i].second  : " << result[i].second  << endl;getchar();
    }

}//plot_gaussians_centered_at_birth_times



template <typename T>
void one_d_gaussians<T>::print_to_file( char* filename  )
{
    cerr << "print_gaussians_centered_at_birth_times_to_file\n";
    ofstream out;
    out.open( filename );
    for ( size_t i = 0 ; i != gaussians.size() ; ++i )
    {
        out << gaussians[i].first << " " << gaussians[i].second << endl;
    }
    out.close();
}//print_to_file
