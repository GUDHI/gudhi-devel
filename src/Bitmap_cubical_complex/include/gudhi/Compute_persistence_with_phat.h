/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once


#include "phat/compute_persistence_pairs.h"
#include "phat/representations/vector_vector.h"
#include "phat/algorithms/standard_reduction.h"
#include "phat/algorithms/chunk_reduction.h"
#include "phat/algorithms/row_reduction.h"
#include "phat/algorithms/twist_reduction.h"


namespace Gudhi
{


//the only aim of this class is to have a ability to compute persistence with phat.
template <typename K>
void writeBettiNumbersAndPersistenceIntervalsToFile( char* prefix , std::pair< std::vector<std::vector< K > > , std::vector< std::vector< std::pair<K,K> > > > resutsFromPhat )
{
    std::ostringstream filenameStr;
    filenameStr << prefix << "_bettiNumbers";
    std::string str = filenameStr.str();
    const char* filename = str.c_str();
    ofstream out;
    out.open( filename );
    for ( size_t dim = 0 ; dim != resutsFromPhat.first.size() ; ++dim )
    {
        out << "Dimension : " << dim << endl;
        for ( size_t i = 0 ; i != resutsFromPhat.first[dim].size() ; ++i )
        {
            out << resutsFromPhat.first[dim][i] << endl;
        }
        out << endl;
    }
    out.close();


    cerr << "Write persistence to file \n";
    for ( size_t dim = 0 ; dim != resutsFromPhat.second.size() ; ++dim )
    {
        cerr << "resutsFromPhat.second[dim].size() : " << resutsFromPhat.second[dim].size()  << endl;
        if ( resutsFromPhat.second[dim].size() == 0 )continue;
        std::ostringstream filenameStr;
        filenameStr << prefix << "_persistence_" << dim;
        std::string str = filenameStr.str();
        const char* filename = str.c_str();
        ofstream out1;
        out1.open( filename );
        for ( size_t i = 0 ; i != resutsFromPhat.second[dim].size() ; ++i )
        {
            out1 << resutsFromPhat.second[dim][i].first << " " << resutsFromPhat.second[dim][i].second << endl;
        }
        out1.close();
    }
}//writeBettiNumbersAndPersistenceIntervalsToFile


template <typename T , typename K>
class Compute_persistence_with_phat
{
public:
    Compute_persistence_with_phat( T* data_structure_ );
    std::pair< std::vector< std::vector<K> > , std::vector< std::vector< std::pair<K,K> > > >  get_the_intervals( phat::persistence_pairs pairs );

    phat::persistence_pairs compute_persistence_pairs_dualized_chunk_reduction();
    phat::persistence_pairs compute_persistence_pairs_twist_reduction();
    phat::persistence_pairs compute_persistence_pairs_standard_reduction();
    //phat::persistence_pairs compute_persistence_pairs_spectral_sequence_reduction();
private:
    void print_bd_matrix();
    phat::boundary_matrix< phat::vector_vector > boundary_matrix;
    T* data_structure;
};

template <typename T , typename K>
void Compute_persistence_with_phat<T,K>::print_bd_matrix()
{
    std::cout << "The boundary matrix has " << this->boundary_matrix.get_num_cols() << " columns: " << std::endl;
    for( phat::index col_idx = 0; col_idx < this->boundary_matrix.get_num_cols(); col_idx++ ) {
        std::cout << "Colum " << col_idx << " represents a cell of dimension " << (int)this->boundary_matrix.get_dim( col_idx ) << ". ";
        if( !this->boundary_matrix.is_empty( col_idx ) ) {
            std::vector< phat::index > temp_col;
            this->boundary_matrix.get_col( col_idx, temp_col );
            std::cout << "Its boundary consists of the cells";
            for( phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++ )
                std::cout << " " << temp_col[ idx ];
        }
        std::cout << std::endl;
    }
}

template <typename T , typename K>
phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_dualized_chunk_reduction()
{
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs_dualized< phat::chunk_reduction >( pairs, this->boundary_matrix );
    return pairs;
}

template <typename T , typename K>
phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_twist_reduction()
{
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs< phat::twist_reduction >( pairs, this->boundary_matrix );
    return pairs;
}

template <typename T , typename K>
phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_standard_reduction()
{
    phat::persistence_pairs pairs;
    phat::compute_persistence_pairs< phat::standard_reduction >( pairs, this->boundary_matrix );
    return pairs;
}

//template <typename T , typename K>
//phat::persistence_pairs Compute_persistence_with_phat<T,K>::compute_persistence_pairs_spectral_sequence_reduction()
//{
//    phat::persistence_pairs pairs;
//    phat::compute_persistence_pairs< phat::spectral_sequence_reduction >( pairs, this->boundary_matrix );
//    return pairs;
//}

template <typename T , typename K>
Compute_persistence_with_phat<T,K>::Compute_persistence_with_phat( T* data_structure_ ):data_structure( data_structure_ )
{
    bool dbg = false;
    this->boundary_matrix.set_num_cols( this->data_structure->num_simplices() );

    //setting up the dimensions of cells:
    for ( size_t i = 0 ; i != this->data_structure->num_simplices() ; ++i )
    {
        this->boundary_matrix.set_dim( i, this->data_structure->dimension( this->data_structure->simplex(i) ) );
    }


    //now it is time to set up the boundary matrix:
    typename T::Filtration_simplex_range range = this->data_structure->filtration_simplex_range();
    std::vector< phat::index > temp_col;
    for ( typename T::Filtration_simplex_iterator it = range.begin() ; it != range.end() ; ++it )
    {
        typename T::Boundary_simplex_range boundary_range = this->data_structure->boundary_simplex_range( *it );
        for ( typename T::Boundary_simplex_iterator bd = boundary_range.begin() ; bd != boundary_range.end() ; ++bd )
        {
            temp_col.push_back( this->data_structure->key( *bd ) );
        }
        //we do not know if the boundary elements are sorted according to filtration, that is why I am enforcing it here:
        this->boundary_matrix.set_col( this->data_structure->key( *it ) , temp_col );
        temp_col.clear();
    }
}

template <typename T , typename K>
std::pair< std::vector< std::vector<K> > , std::vector< std::vector< std::pair<K,K> > > >  Compute_persistence_with_phat<T,K>::get_the_intervals( phat::persistence_pairs pairs )
{
    bool dbg = false;
    //in order to find the birth times of the infinite homology classes, we need to know which elements are not paired. To search for them, we will use this vector:
    std::vector<bool> isTheElementPaired( this->data_structure->num_simplices() , false );

    //now it is time to recover the finite persistence pairs and the Betti numbers:
    std::vector< std::vector< std::pair<K,K> > > finitePersistencePairs( this->data_structure->dimension() );
    for( phat::index idx = 0; idx < pairs.get_num_pairs(); idx++ )
    {
        typename T::Simplex_key positionOfBeginOfInterval = pairs.get_pair( idx ).first;
        typename T::Simplex_key positionOfEndOfInterval = pairs.get_pair( idx ).second;

        typename T::Simplex_handle first_simplex = this->data_structure->simplex(positionOfBeginOfInterval);
        typename T::Simplex_handle second_simplex = this->data_structure->simplex(positionOfEndOfInterval);

        typename T::Filtration_value valueFirst = this->data_structure->filtration( first_simplex );
        typename T::Filtration_value valueSecond = this->data_structure->filtration( second_simplex );

        if ( valueFirst > valueSecond ){std::swap( valueFirst , valueSecond );}

        unsigned dimFirst = this->data_structure->dimension(first_simplex);
        unsigned dimSecond = this->data_structure->dimension(second_simplex);
        unsigned dim = std::min( dimFirst , dimSecond );


        //we are ignoring trivial barcodes
        if ( valueFirst != valueSecond )
        {
            finitePersistencePairs[ dim ].push_back( std::make_pair(valueFirst , valueSecond) );
            if ( dbg ){cerr << "Adding barcode : " << valueFirst << "," << valueSecond << endl;}
        }

        //isTheElementPaired[ positionOfBeginOfIntervalInBitmap ] = true;
        //isTheElementPaired[ positionOfEndOfIntervalInBitmap ] = true;
        isTheElementPaired[ pairs.get_pair( idx ).first ] = true;
        isTheElementPaired[ pairs.get_pair( idx ).second ] = true;
    }


    std::vector< std::vector<K> > birthTimesOfInfinitePersistnceClasses(this->data_structure->dimension()+1 );
    for ( size_t i = 0 ; i != this->data_structure->dimension()+1 ; ++i )
    {
        std::vector<K> v;
        birthTimesOfInfinitePersistnceClasses[i] = v;
    }
    for ( size_t i = 0 ; i != isTheElementPaired.size() ; ++i )
    {
        if ( isTheElementPaired[i] == false )
        {
            //i-th element is not paired, therefore it gives an infinite class
            typename T::Simplex_handle simplex = this->data_structure->simplex(i);
            birthTimesOfInfinitePersistnceClasses[ this->data_structure->dimension( simplex ) ].push_back( this->data_structure->filtration(simplex) );
        }
    }

    //sorting finite persistence pairs:
    for ( size_t dim = 0 ; dim != finitePersistencePairs.size() ; ++dim )
    {
        std::sort( finitePersistencePairs[dim].begin() , finitePersistencePairs[dim].end() );
    }
    return std::make_pair( birthTimesOfInfinitePersistnceClasses , finitePersistencePairs );
}//Compute_persistence_with_phat



}//namespace Gudhi
