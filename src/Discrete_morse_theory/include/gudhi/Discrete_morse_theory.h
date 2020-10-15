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

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

template<typename ComplexWithMatching>
class Discrete_morse_theory {
public:
  typedef ComplexWithMatching Complex;
  typedef typename Complex::Simplex_handle Cell_handle;

/** Restrict the matching the a range of cells. value_type in range must match 
  * Complex::Cell_handle. The available cells belong to the Complex cpx.
  *
  * The range presented by pointers av_begin and av_end must satisfy the following 
  * property: an iterator 'it' in this range always points to a Cell_handle *it 
  * which is maximal w.r.t. to inclusion among all simplices in the range 
  * [it, av_end).
  *
  * Similarly, but in reverse, the range available_cells must satisfy: any iterator 
  * 'it' in the range points to a simplex *it which is maximal w.r.t. inclusion 
  * among all simplices in the range [available_cells.begin(), it). The reversed 
  * convention ensures that compute_matching can be called on a range representing 
  * a filtration, as used in persistent homology. 
  *
  * There is no assumption on the status of cells (critical or not) in the input 
  * range. The algorithm implement the Benedetti-Lutz heuristic to compute a 
  * (not necessarily optimal) Morse matching.
  *
  * Assign a new ordering of the Cell_handles in available_cells satisfying:
  * - the inclusion order, i.e., t \subset s => t comes before s in available_cells
  *   (which is already satisfied by the input available_cells), and
  * - two paired simplices (t,s), t \subset s, must be consecutive.
  */
  template<typename RangeCells>
  void compute_matching(RangeCells &available_cells, Complex *cpx) {
    std::vector<Cell_handle> new_order = //in reverse inclusion order
            compute_matching(available_cells.rbegin(), available_cells.rend(), cpx);
    auto it2 = new_order.rbegin();//read it in reverse
    for(auto it1=available_cells.begin(); it1 != available_cells.end(); 
             ++it1, ++it2) {  *it1 = *it2;  }//copy new_order into available_cells
  }
  template<typename IteratorCells>
  std::vector<Cell_handle> 
        compute_matching(IteratorCells av_begin, IteratorCells av_end, Complex *cpx)
  {
    //use the key to mark simplices: 
    // Any face NOT in the input range must have a Simplex_key != -1 and != -2.
    //-1 = free to pair with another cell in the range, 
    //-2 already treated (paired with something), postpone insertion. 
    //Return a different vector where paired simplices are consecutive 
    //-> ensure that they get consecutive keys later.
    // When a simplex s is paired with a facet t, s is maximal in [s ... end()) but
    // there may be other cofaces of t in (s ... t].
    //We postpone the insertion of s in new_filt until it_sh falls on t, we then 
    //insert consecutively s,t in new_filt.
    for(auto it_sh = av_begin; it_sh != av_end; ++it_sh) {
      cpx->assign_key(*it_sh,-1);
    }
    //the output filtration in reverse inclusion order (cofaces first) + paired 
    //simplices are consecutive
    std::vector<Cell_handle> new_filtration; 
    new_filtration.reserve(av_end-av_begin);
    
    for(auto it_sh = av_begin; it_sh != av_end; ++it_sh) {
      switch(cpx->key(*it_sh)) {
        //
        case -1: {//not yet paired, try to pair with a facet in range
          bool critical = true;//is the face critical?
          Cell_handle paired_sh;
          for(auto b_sh : cpx->boundary_simplex_range(*it_sh)) {//all facets
            if(cpx->key(b_sh) == -1) {//free facet in range, pair *it_sh with b_sh
              cpx->assign_pairing(b_sh,*it_sh);
              cpx->assign_key(b_sh,-2);//do not try to pair b_sh in the future
              critical = false;//not critical
              paired_sh = b_sh;
              break;            
            }
          }
          if(critical) {//*it_sh is definitely critical
            cpx->make_critical(*it_sh); 
            new_filtration.push_back(*it_sh);
          }//else postpone the insertion into new_filt
          break;
        }
        //       
        case -2: {//it_sh has been paired earlier with a coface, insert both in new_or
          new_filtration.push_back(cpx->paired_with(*it_sh));//the coface s
          new_filtration.push_back(*it_sh);//the facet t, in pair (s,t)
          break;
        }
      }
    }
    return new_filtration;
  }

  /** Remove the matching by turning all cells in cpx critical.*/
  void clear_matching(Complex *cpx) {
    for(auto sh : cpx->complex_simplex_range()) { cpx->make_critical(sh); }
  }
};

#endif //DISCRETE_MORSE_THEORY_H_