/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DISCRETE_MORSE_THEORY_H_
#define DISCRETE_MORSE_THEORY_H_

#include <set>

template<typename ComplexWithMatching>
class Discrete_morse_theory {
public:
  typedef ComplexWithMatching Complex;
  typedef typename Complex::Simplex_handle Cell_handle;

/** Restrict the matching the a range of cells. value_type in range must match 
  * Complex::Simplex_handle. The available cells belong to the Complex cpx.
  *
  * The range must be ordered in filtration order (i.e., faces before cofaces). We 
  * read it from right to left.
  *
  * Assume all faces are marked critical by default.
  * Benedetti-Lutz heuristic.
  */
  template<typename RangeCells>
  void compute_matching(const RangeCells &available_cells, Complex *cpx) {
    compute_matching(available_cells.rbegin(), available_cells.rend(), cpx);
  }
  template<typename IteratorCells>
  void compute_matching(IteratorCells av_begin, IteratorCells av_end, Complex *cpx)
  {
    //sort by address, which is a O(1) comparison function
    struct Address_cmp {
      //compare Nodes pointed by Simplex_handle by address
      bool operator()(const Cell_handle& lhs, const Cell_handle& rhs) const { 
          return &(lhs->second) < &(rhs->second); 
      }
    };
    std::set<Cell_handle, Address_cmp> available;
    
#ifdef GUDHI_USE_TBB
    std::vector<Cell_handle> tmp_av; tmp_av.reserve(av_end-av_begin);
    for(auto it_sh = av_begin; it_sh != av_end; ++it_sh) {
      tmp_av.push_back(*it_sh);
    }
    tbb::parallel_sort(tmp_av.begin(),tmp_av.end(), Address_cmp());
    available = std::set<Cell_handle, Address_cmp>(tmp_av.begin(),tmp_av.end());
#else
    for(auto it_sh = av_begin; it_sh != av_end; ++it_sh) 
    {  available.insert(*it_sh);  }
#endif

    auto it_sh = av_begin;//in reverse inclusion order.
    
    bool free_pair;
    while(!available.empty()) {
      //max cell to be paired
      typename std::set<Cell_handle>::iterator it_av_sh = available.find(*it_sh);
      typename std::set<Cell_handle>::iterator it_av_b_sh;//to find available facets
      free_pair = false;
      if(it_av_sh != available.end()) {//if it has not been paired yet
        //traverse all facets to potentially find an available one
        for(auto b_sh : cpx->boundary_simplex_range(*it_sh)) {
          it_av_b_sh = available.find(b_sh); 
          //is the facet b_sh available for match, i.e. in range and not yet paired?
          if(it_av_b_sh != available.end()) {//yes, pair sh with b_sh
            cpx->assign_pairing(b_sh,*it_sh);
            free_pair = true;
            break;
          }
        }
        if(free_pair) {//we have found a free pair, remove both cells from the set
          available.erase(it_av_b_sh);
        }
        else { cpx->make_critical(*it_sh); }
        available.erase(it_av_sh);//always remove the max face we have considered
      }
      ++it_sh;
    }
    for(auto remain_sh : available) {
      //finish the range -> all remaining cells are critical
      cpx->make_critical(remain_sh);
    }
  }



  void clear_matching(Complex *cpx) {
    for(auto sh : cpx->complex_simplex_range()) { cpx->make_critical(sh); }
  }
};

#endif //DISCRETE_MORSE_THEORY_H_