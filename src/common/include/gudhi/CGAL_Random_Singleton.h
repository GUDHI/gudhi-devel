/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2025 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CGAL_RANDOM_SINGLETON_H_
#define CGAL_RANDOM_SINGLETON_H_

#include <CGAL/Random.h>
#include <CGAL/version.h>  // for CGAL_VERSION_NR

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error random_point_generators is only available for CGAL >= 4.11
#endif

namespace Gudhi {

  class CGAL_Random_Singleton {
   public:
    static CGAL::Random* get() {
      if(singleton_ == nullptr){
        singleton_ = new CGAL::Random();
      }
      return singleton_;
    }
    static CGAL::Random* get(unsigned int seed) {
      set_seed(seed);
      return singleton_;
    }
    static void set_seed(unsigned int seed) {
      if(singleton_ != nullptr){
        delete singleton_;
      }
      singleton_ = new CGAL::Random(seed);
    }
   private:
    CGAL_Random_Singleton() {};
    static CGAL::Random* singleton_;
  };

  CGAL::Random* CGAL_Random_Singleton::singleton_= nullptr;

}  // namespace Gudhi

#endif // CGAL_RANDOM_SINGLETON_H_
