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
#include <map>
#include <list>
#include <gudhi/Simplex_tree.h>

#ifdef GUDHI_USE_TBB
#include <tbb/parallel_sort.h>
#endif

namespace Gudhi {
namespace dmt {
/** \brief Compute a Morse matching on a range of cells in a cell complex. 
 * 
 * \details The iterators of the range must have value_type  
 * Complex::Simplex_handle. All such simplex handles point to valid simplices in 
 * the complex cpx. The cells of the input range must all be critical in cpx. For a 
 * given cell c in the range, all of its cofaces must be in the range as well as, 
 * otherwise, we may create a partial matching with cycles. All simplices of the complex cpx, not in the input range, must have a key >= 0.
 *
 * The output vector of simplex handles contains all simplex handles of the input 
 * range, in reverse inclusion order, i.e., any iterator 'it' in the output vector 
 * points to a simplex *it which is maximal w.r.t. inclusion among all simplices in 
 * the range (it, ...]. Additionally, if a simplex t has been paired with a simplex 
 * s, then s and t are consecutive in the output range.
 * 
 * The algorithm implements the Benedetti-Lutz heuristic to compute a 
 * (not necessarily optimal) Morse matching.
 *
 * The complexity is at most quadratic in the size of the input range, times the 
 * complexity of computing the boundary of a simplex. 
 */
  template<typename IteratorCells, typename Complex>
  std::vector<typename Complex::Simplex_handle> 
    compute_matching(IteratorCells av_begin, IteratorCells av_end, Complex *cpx)
  { //use the key to mark simplices: 
    // Any face NOT in the input range must have a Simplex_key >= 0.
    // We use the key to count the number of cofacets in the range, i.e., we assign 
    // the key -(x+1) to sh (in the range) iff sh has x cofacets in the range.
    //
    //if a simplex in the range has already been considered (set critical or paired 
    //in the range, we give it key 0 to laxily remove it later).

    for(auto it_sh = av_begin; it_sh != av_end; ++it_sh) {
      cpx->assign_key(*it_sh,-1);//mark all simplices in the range
      if(!cpx->critical(*it_sh)) { std::cout << "Simplex not critical ; undefined behavior.\n"; }
    }
    //count the number of cofacets of each simplex in the range
    for(auto it_sh = av_begin; it_sh != av_end; ++it_sh) {
      for(auto b_sh : cpx->boundary_simplex_range(*it_sh)) {
        auto curr_key = cpx->key(b_sh);
        //negative => in the range, by assumption
        if(curr_key <0) { cpx->assign_key(b_sh, curr_key-1); }
      }
    }
    //now, a simplex with key -1 is maximal in the range, if one of its facets has 
    //key -2, they form a free pair

    //the output filtration in reverse inclusion order (cofaces first) + paired 
    //simplices are consecutive -> todo
    std::vector<typename Complex::Simplex_handle> new_filtration; 
    new_filtration.reserve(av_end-av_begin);
    //all available simplices
    std::list< typename Complex::Simplex_handle > available_list(av_begin,av_end);

    typename std::list< typename Complex::Simplex_handle >::iterator some_max_simplex;//record a max simplex
    while(!available_list.empty()) {
      bool found_free_pair = false;

      for(auto it = available_list.begin(); it != available_list.end(); ++it) {
        //when a simplex is not available, we mark it with a key==0 for a lazy remove
        if(cpx->key(*it) == 0) { //simplex mark -> effectively remove it form av_list
          auto tmp_it = it;  ++it;
          available_list.erase(tmp_it);
        }
        else {//simplex really available
          if(cpx->key(*it) == -1) {//maximal simplex
            //look for a facet with a single coface, i.e. key==-2
            for(auto b_sh : cpx->boundary_simplex_range(*it)) {
              if(cpx->key(b_sh) == -2) {//we have a free pair, remove both
                //one less coface for the simplices of the boundary of *it
                for(auto bb_sh : cpx->boundary_simplex_range(*it)) {
                  auto curr_key = cpx->key(bb_sh);
                  if(curr_key < 0) { cpx->assign_key(bb_sh, curr_key+1); }
                }
                //one less coface for the simplices of the boundary of b_sh
                for(auto bb_sh : cpx->boundary_simplex_range(b_sh)) {
                  auto curr_key = cpx->key(bb_sh);
                  if(curr_key < 0) { cpx->assign_key(bb_sh, curr_key+1); }
                }
                //pair the two simplices
                cpx->assign_pairing(b_sh,*it);
                //remove effectively *it from available simplices
                auto tmp_it = it; ++it;
                available_list.erase(tmp_it);
                //mark b_sh for lazy removal
                cpx->assign_key(b_sh,0);
                //we've found a free pair
                found_free_pair = true;
                //record them consecutively in new_filtration, coface first
                new_filtration.push_back(*it);
                new_filtration.push_back(b_sh);
                break; //done with *it
              }
            }
            if(!found_free_pair) { some_max_simplex = it; }//record an available max simplex
          }
        }
      }
      //if we haven't found a free pair, we remove the recorded max simplex
      if(!found_free_pair) {
        //one less coface for the simplices of the boundary of *some_max_simplex
        for(auto bb_sh : cpx->boundary_simplex_range(*some_max_simplex)) {
          auto curr_key = cpx->key(bb_sh);
          if(curr_key < 0) { cpx->assign_key(bb_sh, curr_key+1); }
        } 
        new_filtration.push_back(*some_max_simplex);//record in new_filtration
        available_list.erase(some_max_simplex);//remove from available simplices
      }
    }
    return new_filtration;
  }




  //   //-1 = free to pair with another cell in the range, 
  //   //-2 already treated (paired with something), postpone insertion. 
  //   //Return a different vector where paired simplices are consecutive 
  //   //-> ensure that they get consecutive keys later.
  //   // When a simplex s is paired with a facet t, s is maximal in [s ... end()) but
  //   // there may be other cofaces of t in (s ... t].
  //   //We postpone the insertion of s in new_filt until it_sh falls on t, we then 
  //   //insert consecutively s,t in new_filt.
  //   for(auto it_sh = av_begin; it_sh != av_end; ++it_sh) {
  //     cpx->assign_key(*it_sh,-1);
  //   }

  //   for(auto it_sh = av_begin; it_sh != av_end; ++it_sh) {
  //     switch(cpx->key(*it_sh)) {
  //       //
  //       case -1: {//not yet paired, try to pair with a facet in range
  //         bool critical = true;//is the face critical?
  //         typename Complex::Simplex_handle paired_sh;
  //         for(auto b_sh : cpx->boundary_simplex_range(*it_sh)) {//all facets
  //           if(cpx->key(b_sh) == -1) {//free facet in range, pair *it_sh with b_sh
  //             cpx->assign_pairing(b_sh,*it_sh);
  //             cpx->assign_key(b_sh,-2);//do not try to pair b_sh in the future
  //             critical = false;//not critical
  //             paired_sh = b_sh;
  //             break;            
  //           }
  //         }
  //         if(critical) {//*it_sh is definitely critical
  //           cpx->make_critical(*it_sh); 
  //           new_filtration.push_back(*it_sh);
  //         }//else postpone the insertion into new_filt
  //         break;
  //       }
  //       //       
  //       case -2: {//it_sh has been paired earlier with a coface, insert both in new_or
  //         new_filtration.push_back(cpx->paired_with(*it_sh));//the coface s
  //         new_filtration.push_back(*it_sh);//the facet t, in pair (s,t)
  //         break;
  //       }
  //     }
  //   }
  //   return new_filtration;
  // }

/** \brief Compute a Morse matching on a range of cells in a cell complex, and 
 * re-order cells in the input range.
 * 
 * \details The new order of cells in the input range satisfies the following 
 * properties:
 * - a subface comes before any of its cofaces,
 * - two simplices that have been paired together in the Morse matching are 
 * consecutive in the new ordering. 
 *
 * Any cell in the input range must have all its cofaces in the range, otherwise we 
 * may compute a partial matching with cycles.
 */ 


  template<typename RangeCells, typename Complex>
  void compute_matching(RangeCells &available_cells, Complex *cpx) 
  { //compute the Morse matching, the output is in reverse inclusion order
    std::vector<typename Complex::Simplex_handle> new_order = 
            compute_matching(available_cells.begin(), available_cells.end(), cpx);
    auto it2 = new_order.rbegin();//read it in reverse
    for(auto it1=available_cells.begin(); it1 != available_cells.end(); ++it1, ++it2)
    {  *it1 = *it2;  }//copy new_order into available_cells
  }

  /** Remove the matching by turning all cells in cpx critical.*/
  template< typename Complex >
  void clear_matching(Complex *cpx) {
    for(auto sh : cpx->complex_simplex_range()) { cpx->make_critical(sh); }
  }

/** \brief Type storing a chain in a chain complex, as a map of 
 * Simplex_handle -> coefficient.
 * 
 * \details Supports coefficients as int.
 */
template< typename Complex >
struct Chain {
  //coefficient in a ring
  typedef int                                   Coefficient;
  typedef typename Complex::Simplex_handle      Simplex_handle;
  typedef std::map<Simplex_handle, Coefficient> Chain_t;

  Chain() {}

  void add(Simplex_handle sh, Coefficient w) {
    auto res_insert = chain_.emplace(sh,w);
    if(!res_insert.second) {//sh already in map
      res_insert.first.second += w;
      if(res_insert.first.second == 0) { chain_.erase(res_insert.first); }
    } 
  }
  void add(Chain &other) {
    for(auto &sh_w : other.chain_) { add(sh_w.first, sh_w.second); }
  }

  typename Chain_t::iterator begin() { return chain_.begin(); }
  typename Chain_t::iterator end ()  { return chain_.end();   }

  Chain_t chain() { return chain_; }

private:
  Chain_t chain_;
};



/* Computation of boundaries in a Morse complex. */

  template< typename Complex, typename PartialBoundariesMap >
  typename std::map< typename Complex::Simplex_handle, Chain<Complex> >::iterator 
  rec_boundary_morse_complex( 
    Complex *cpx, typename Complex::Simplex_handle t_sh, int dim_sh, 
    PartialBoundariesMap &partial_boundaries);

/** \brief Computes the boundary of a cell in a Morse complex.
 * 
 * \details Uses a depth first search approach, with a recursive implementation. 
 * Every time the search finishes with a non-critical simplex, it stores a chain in
 * this complex. The chain corresponds to a set of critical simplices and their 
 * coefficients (int), that can be reached in the Hasse diagram starting from that 
 * simplex. This allows to not traverse several times the same arrow.
 * 
 * In a simplicial complex, with a Morse maching admitting k critical cells in 
 * dimension d-1, and m non-critical cells of dimension d,  
 * and where computing the simplicial complex boundary of a simplex of dimension d 
 * costs O(B) operations, 
 * the complexity of computing the Morse boundary of a critical cell c of dimension 
 * d is: \f$O( m \cdot (B + d k + \log (m) + d)\f$ as,
 * - each edge of the Hasse diagram is traversed at most once (m d),
 * - we compute the complex boundary of a non-critical cell at most once (mB),
 * - we store a partial boundary in each non-critical cell (agglomeration of the 
 * recursion in \f$m \cdot d \cdot k\f$ opeartions),
 * - we maintain the non-critical cells already visited in a map (log m per call).
 * 
 * Complex must be a Simplex_tree (the computation of coefficient only works for a 
 * Simplex_tree.)
 *  
 * todo: 1/ Complex::Boundary_simplex_iterator gives a coefficient .coefficient() of 
 * type in a Ring.
 * 
 * 2/ Store intermediate data in Node to avoid passing several times through the same 
 * non-critical cell.
 *
 * 3/ Strict weak ordering of Simplex_handles in the Simplex_tree
 * 
 * 4/ Class Chain for insert map, etc.
 */ 
  template< typename Complex >
  std::map< typename Complex::Simplex_handle, int > 
    boundary_morse_complex( Complex *cpx, typename Complex::Simplex_handle sh ) 
  {
    //For any non-critical cell t of dimension dim(sh)-1, records the chain of 
    //critical simplices obtained by a depth first search from t. 
    std::map< typename Complex::Simplex_handle, Chain<Complex>
            , typename Complex::cmp_simplices >          partial_boundaries;

    //The multiplicity of a gradient path from s (dim d) to t (dim d-1) is:
    //the product of -1 / (incidence rel up arrow) for arrows going up (Morse pair),
    //times (incidence rel down arrow) for arrows going down
    Chain<Complex> boundary_sh;//record the final boundary
    int dim_sh = cpx->dimension(sh);//dimension d of the simplex sh
    int incidence_coefficient = 1 - 2 * (dim_sh % 2);//+/-1 for incidence relation
    for(auto b_sh : cpx->boundary_simplex_range(sh)) {
      if(cpx->critical(b_sh)) { boundary_sh.add(b_sh,incidence_coefficient); }
      else { 
        //iterator in partial_boundaries corresponding to the facet b_sh.
        auto it = rec_boundary_morse_complex<Complex>(cpx, b_sh, dim_sh, partial_boundaries);
        boundary_sh.add(it->first, incidence_coefficient);
      }
      incidence_coefficient = -1 * incidence_coefficient;
    }

// //forget +/-1 copefficients
//     std::vector< typename Complex::Simplex_handle > morse_boundary;
//     morse_boundary.reserve(boundary_sh.size());
//     for(auto &sh_w : boundary_sh) { 
//       morse_boundary.push_back(sh_w.first); 
//     }
//     return morse_boundary;
    return boundary_sh.chain();
  }

  /* Computes recursively the partial boundary of t_sh, and insert it in 
   * partial_boundaries. Returns an iterator to 
   * partial_boundaries[t_sh]=partial_boundary_t_sh
   * 
   * @param[in] t_sh non-critical simplex of dimension dim(sh)-1
   */
  template< typename Complex, typename PartialBoundariesMap >
  typename std::map< typename Complex::Simplex_handle, Chain<Complex> >::iterator 
  rec_boundary_morse_complex( 
    Complex *cpx, typename Complex::Simplex_handle t_sh, int dim_sh, 
    PartialBoundariesMap &partial_boundaries)
  {
    //test wether the partial boundary of t_sh has already been computed
    auto res_insert = partial_boundaries.emplace(t_sh, Chain<Complex>());
    if(!res_insert.second) { return res_insert.first; }//already computed
    //otherwise, compute it recursively
    auto c_sh = cpx->paired_with(t_sh);//Morse pair (t,c)
    int incidence_c_t;//records weight of edge t->c                     //iterator
    //else compute the partial boundary at tsh
    auto it_tsh = res_insert.first;
    //pre-compute the chain boundary of c_sh, minus t_sh:
    std::vector< std::pair<typename Complex::Simplex_handle, int> > tmp_boundary;
    tmp_boundary.reserve(dim_sh);

    int incidence_coefficient = 1 - 2 * (dim_sh % 2);//+/-1 for incidence relation

    for(auto b_sh : cpx->boundary_simplex_range(c_sh)) {
      //true iff b_sh == t_sh, record the weight of (t_sh->c_sh)
      if(cpx->same_simplex(b_sh, t_sh)) { incidence_c_t = incidence_coefficient; }
      else { tmp_boundary.emplace_back(b_sh, incidence_coefficient); }
      incidence_coefficient = -1 * incidence_coefficient;
    }

    typename Complex::Simplex_handle sh;   int weight;
    for(auto &sh_w : tmp_boundary) {
      std::tie(sh,weight) = sh_w;
      //if critical, add the simplex to partial boundary of t, weight is the 
      //product: -1 * 1/[t_sh:c_sh] *[sh:c_sh] (length 2 path from t to this simplex)
      if(cpx->critical(sh)) { 
        it_tsh->second.add(sh, weight * (1/incidence_c_t) * (-1));
      }
      else { //recursively compute/output the partial boundary of sh
          auto it_sh =rec_boundary_morse_complex(cpx, sh, dim_sh, partial_boundaries);
          //add to partial boundary of t       
          it_tsh->second.add(it_sh->second, weight * (1/incidence_c_t) * (-1));
      }
    }
    return it_tsh;
  }

} //namespace dmt
} //namespace Gudhi

#endif //DISCRETE_MORSE_THEORY_H_