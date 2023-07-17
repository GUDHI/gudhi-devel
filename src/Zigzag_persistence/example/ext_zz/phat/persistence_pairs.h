/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Michael Kerber, Jan Reininghaus

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */

#pragma once

#include <phat/helpers/misc.h>

namespace phat {
    class persistence_pairs {

    protected:
        std::vector< std::pair< index, index > > pairs;

    public:
        index get_num_pairs() const { 
            return (index)pairs.size();
        }

        void append_pair( index birth, index death ) {
            pairs.push_back( std::make_pair( birth, death ) );
        }

        std::pair< index, index > get_pair( index idx ) const {
            return pairs[ idx ];
        }

        void set_pair( index idx, index birth, index death )  {
            pairs[ idx ] = std::make_pair( birth, death );
        }
        
        void clear() {
            pairs.clear();
        }

        void sort() {
            std::sort( pairs.begin(), pairs.end() );
        }

        // Loads the persistence pairs from given file in asci format 
        // Format: nr_pairs % newline % birth1 % death1 % newline % birth2 % death2 % newline ...
        bool load_ascii( std::string filename ) {
            std::ifstream input_stream( filename.c_str() );
            if( input_stream.fail() )
                return false;

            int64_t nr_pairs;
            input_stream >> nr_pairs;
            pairs.clear();
            for( index idx = 0; idx < nr_pairs; idx++ ) {
                int64_t birth;
                input_stream >> birth;
                int64_t death;
                input_stream >> death;
                append_pair(  (index)birth,  (index)death );
            }

            input_stream.close();
            return true;
        }

        // Saves the persistence pairs to given file in binary format 
        // Format: nr_pairs % newline % birth1 % death1 % newline % birth2 % death2 % newline ...
        bool save_ascii( std::string filename ) {
            std::ofstream output_stream( filename.c_str() );
            if( output_stream.fail() )
                return false;

            this->sort();
            output_stream << get_num_pairs() << std::endl;
            for( std::size_t idx = 0; idx < pairs.size(); idx++ ) {
                output_stream << pairs[idx].first << " " << pairs[idx].second << std::endl;
            }

            output_stream.close();
            return true;
        }

        // Loads the persistence pairs from given file in binary format 
        // Format: nr_pairs % birth1 % death1 % birth2 % death2 ...
        bool load_binary( std::string filename ) {
            std::ifstream input_stream( filename.c_str(), std::ios_base::binary | std::ios_base::in );
            if( input_stream.fail() )
                return false;

            int64_t nr_pairs;
            input_stream.read( (char*)&nr_pairs, sizeof( int64_t ) );
            for( index idx = 0; idx < nr_pairs; idx++ ) {
                int64_t birth;
                input_stream.read( (char*)&birth, sizeof( int64_t ) );
                int64_t death;
                input_stream.read( (char*)&death, sizeof( int64_t ) );
                append_pair(  (index)birth,  (index)death );
            }

            input_stream.close();
            return true;
        }

        // Saves the persistence pairs to given file in binary format 
        // Format: nr_pairs % birth1 % death1 % birth2 % death2 ...
        bool save_binary( std::string filename ) {
            std::ofstream output_stream( filename.c_str(), std::ios_base::binary | std::ios_base::out );
            if( output_stream.fail() )
                return false;

            this->sort();
            int64_t nr_pairs = get_num_pairs();
            output_stream.write( (char*)&nr_pairs, sizeof( int64_t ) );
            for( std::size_t idx = 0; idx < pairs.size(); idx++ ) {
                int64_t birth = pairs[ idx ].first;
                output_stream.write( (char*)&birth, sizeof( int64_t ) );
                int64_t death = pairs[ idx ].second;
                output_stream.write( (char*)&death, sizeof( int64_t ) );
            }

            output_stream.close();
            return true;
        }

        bool operator==( persistence_pairs& other_pairs ) {
            this->sort();
            other_pairs.sort();
            if( pairs.size() != (std::size_t)other_pairs.get_num_pairs() )
                return false;

            for( index idx = 0; idx < (index)pairs.size(); idx++ )
                if( get_pair( idx ) != other_pairs.get_pair( idx ) )
                    return false;

            return true;
        }

        bool operator!=( persistence_pairs& other_pairs ) {
            return !( *this == other_pairs );
        }
    };
    

    
}
