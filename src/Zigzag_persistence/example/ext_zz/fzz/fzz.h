#ifndef _FZZ_H_
#define _FZZ_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <memory>

namespace FZZ { 

typedef int Integer;
typedef std::vector<Integer> Simplex;

class FastZigzag {
public:
    /*
      'filt_simp' and 'filt_op' should have the same length which altogether
      specify the input zigzag filtration. 'filt_simp' specifies the simplices
      being added or deleted (following the order of the filtration) and 
      'filt_op' specifies whether it's an addition (true) or deletion (false).
      'persistence' returns the barcode, with the first element of the tuple
      being the birth, the second element being the death, and the third
      being the dimension.
    */
    void compute(
        const std::vector<Simplex> &filt_simp, 
        const std::vector<bool> &filt_op,
        std::vector< std::tuple<Integer, Integer, Integer> > *persistence);

private:
    void mapOrdIntv(Integer &b, Integer &d) {
        // assert(b-1 > 0);
        // assert(d < orig_f_add_id.size());

        // Up-down interval is same, 
        // so directly map to interval of input filtration
        b = orig_f_add_id[b-1] + 1;
        d = orig_f_add_id[d];
    }

    void mapRelExtIntv(Integer &p, Integer &b, Integer &d) {
        // assert(d >= simp_num);

        if (b > simp_num) { // Open-closed
            // Map to up-down interval
            std::swap(b, d);
            b = 3*simp_num - b;
            d = 3*simp_num - d;
            p --;

            // Map to interval of input filtration
            b = orig_f_del_id[b-1-simp_num] + 1;
            d = orig_f_del_id[d-simp_num];
        } else { // Closed-closed
            // Map to up-down interval
            d = 3*simp_num - d-1;
            
            // Map to interval of input filtration
            b = orig_f_add_id[b-1];
            d = orig_f_del_id[d-simp_num];

            if (b < d) {
                b = b+1;
            } else {
                std::swap(b, d);
                b = b+1;
                p = p-1;
            }
        }
    }

private:
    // 'orig_f_add_id' and 'orig_f_del_id' form a mapping 
    // from the up-down filtration to the original filtration
    std::vector<Integer> orig_f_add_id;
    std::vector<Integer> orig_f_del_id;

    Integer simp_num;
};

}

#endif
