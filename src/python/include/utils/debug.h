/*    This file is part of the MMA Library - https://gitlab.inria.fr/dloiseau/multipers - which is released under MIT.
 *    See file LICENSE for full license details.
 *    Author(s):       David Loiseaux
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2022/03 Hannah Schreiber: Integration of the new Vineyard_persistence class, renaming and cleanup.
 */
/**
 * @file debug.h
 * @author David Loiseaux,
 * @brief Display functions for debug purposes
 */


#ifndef DEBUG_H_INCLUDED
#define DEBUG_H_INCLUDED

#include <chrono>
#include <iostream>
#include <vector>
#include <list>

#include <functional>

namespace utils {

using clk = std::chrono::high_resolution_clock;
using tp = clk::time_point;

constexpr bool debug = false;


class Timer
{
public:
    Timer() : activated_(false) {}
    Timer(const std::string &string, bool verbose)
        : timer_(clk::now()), activated_(verbose)
    {
        if(verbose){
            std::cout << string << std::flush;
        }
    }
    ~Timer(){
        if (activated_)
        {
            std::chrono::duration<double> elapsed =
                    std::chrono::duration_cast<std::chrono::duration<double>>(
                            clk::now() - timer_);
            std::cout << " Done ! It took "<< elapsed.count()
                      << " seconds." << std::endl;
        }
    }

private:
    tp timer_;
    bool activated_;
};




template<typename T>
void disp_vect(std::vector<T> v){
    for(uint i=0; i< v.size(); i++){
        std::cout << v[i] << " ";
    }
    std::cout << std::endl;
}

template<typename T>
void disp_vect(std::list<T> v){
    while(!v.empty()){
        std::cout << v.front() << " ";
        v.pop_front();
    }
    std::cout << std::endl;
}

template<typename T>
void disp_vect(std::vector<std::pair<T,T> > v){
    for(unsigned int i=0; i< v.size(); i++){
        std::cout << "(" << v[i].first << " " << v[i].second <<")  ";
    }
}

template<typename T>
void disp_vect(std::vector<std::vector<T>> v, bool show_small = true){
    for(uint i = 0; i < v.size(); i++){
        if(v[i].size() <= 1 && !show_small) continue;
        std::cout << "(";
        for (uint j = 0; j < v[i].size(); j++){
            std::cout << v[i][j];
            if(j < v[i].size() - 1) std::cout << " ";
        }
        std::cout << ") ";
    }
    std::cout << std::endl;
}


} //namespace Debug
namespace std{
    template<typename T>
    std::ostream& operator<<(std::ostream& stream, const std::vector<T> truc){
        stream << "[";
        for(unsigned int i = 0; i < truc.size()-1; i++){
            stream << truc[i] << ", ";
        }
        if(!truc.empty()) stream << truc.back();
        stream << "]";
        return stream;
    }
}






#endif // DEBUG_H_INCLUDED


