
#ifndef DISCRETE_MORSE_THEORY_H_
#define DISCRETE_MORSE_THEORY_H_

#include <set>

template<typename ComplexWithMatching>
class Discrete_morse_theory {
  typedef ComplexWithMatching Complex;
  typedef Complex::Simplex_handle Cell_handle;

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
  void compute_matching(RangeCells available_cells, Complex *cpx) {
    //sort by address, which is a O(1) comparison function
    struct Address_cmp {
      //compare Nodes pointed by Simplex_handle by address
      bool operator()(const Cell_handle& lhs, const Cell_handle& rhs) const { 
          return &(lhs->second) < &(rhs->second); 
      }
    };
    std::set<Cell_handle, Address_cmp> available;
    for(auto sh : available_cells) { available.insert(sh); }

    auto it_sh = available_cells_.rbegin();//in filtration order already, read from the end
    
    bool free_pair;
    while(!available.empty()) {
      //max cell to be paired
      std::set<Cell_handle>::iterator it_av_sh = available.find(*it_sh);
      std::set<Cell_handle>::iterator it_av_b_sh;//to find available facets
      free_pair = false;
      if(av_sh != available.end()) {//if it has not been paired yet
        //traverse all facets to potentially find an available one
        for(auto b_sh : cpx->boundary_simplex_range(*it)) {
          it_av_b_sh = available.find(b_sh); 
          //is the facet b_sh available for match, i.e. in range and not yet paired?
          if(it_av_b_sh != available.end()) {//yes, pair sh with b_sh
            cpx->assign_morse_pairing(b_sh,*it_sh);
            free_pair = true;
            break;
          }
        }
        if(free_pair) {//we have found a free pair, remove both cells from the set
          available.erase(it_av_b_sh);
        }
        available.erase(it_av_sh);//always remove the max face we have considered
      }
      ++it_sh;
    }
    // while(it_sh != available_cells.rend()) {
    //   //finish the range -> not necessary if cells are critical by default
    // }
  }
};

#endif //DISCRETE_MORSE_THEORY_H_