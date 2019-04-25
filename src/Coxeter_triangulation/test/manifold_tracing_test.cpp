#include <gudhi/Manifold_tracing.h>
#include <gudhi/Implicit_function_intersection_oracle.h>

#include <gudhi/Coxeter_triangulation/Freudenthal_representation.h>
#include <gudhi/Coxeter_triangulation.h>

#include "../example/functions/sphere_S1_in_R2.h"

using Vertex = std::vector<int>;
using Part = std::vector<uint>;
using Ordered_partition = std::vector<Part>;
using FR = Freudenthal_representation<Vertex, Ordered_partition>;

/* a hash function for the FR */
namespace std {
  template<>
  struct hash<FR> {
    typedef FR argument_type;
    typedef std::size_t result_type;
    result_type operator()(const argument_type& c) const noexcept {
      return c.vertex[0];
    }
  };
}

bool compare_parts(const Part& lhs, const Part& rhs) {
  Part lhs_sorted(lhs);
  Part rhs_sorted(rhs);
  std::sort(lhs_sorted.begin(), lhs_sorted.end());
  std::sort(rhs_sorted.begin(), rhs_sorted.end());
  return lhs_sorted < rhs_sorted;
}

/* It is supposed that the two part(itions have the same size */
bool compare_partitions(const Ordered_partition& lhs, const Ordered_partition& rhs) {
  assert(lhs.size() == rhs.size());
  for (std::size_t k = 0; k < lhs.size(); k++) {
    if (compare_parts(lhs[k], rhs[k]))
      return true;
    else if (compare_parts(rhs[k], lhs[k]))
      return false;
  }
  return false;
  // return std::sort(lhs.begin(), lhs.end()) < std::sort(rhs.begin(), rhs.end());
}


/** \brief The default comparison operator between two tuples of cell coordinates.
    \detail If the sizes of the two tuples are different, the smallest tuple is returned as the smallest. If the sizes are the same, then order for the comparison is the lexicographical order, with a true coordinate smaller than a false coordinate, in the case of difference.
*/
namespace std {
  bool operator< (const FR& lhs, const FR& rhs) {
    if (lhs.dimension() < rhs.dimension())
      return true;
    else if (lhs.dimension() > rhs.dimension())
      return false;
    else
      for (std::size_t k = 0; k < lhs.vertex.size(); ++k) {
	if (lhs.vertex[k] < rhs.vertex[k])
	  return true;
	else if (lhs.vertex[k] > rhs.vertex[k])
	  return false;
      }
    if (compare_partitions(lhs.partition, rhs.partition))
      return true;
    else
      return false;
  }
}


int main() {
  using uint = unsigned;
  using Kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Point_d = Kernel::Point_d;
  
  Gudhi::Manifold_tracing mt;

  Kernel k;
  /* function definition */ 
  double r = 5;
  Function_S1_in_R2 fun(r);
  Gudhi::Implicit_function_intersection_oracle<Function_S1_in_R2> oracle(fun);

  /* triangulation definition */
  Gudhi::Coxeter_triangulation triangulation(2);
  FR s = {Vertex(2,0), Ordered_partition({Part(1,1), Part(1,2)})};
  std::cout << "Intersects = " << oracle.intersects(s, triangulation) << "\n";
  
  std::vector<Point_d> seed_points = {Gudhi::construct_point(k, r+fun.off_[0], fun.off_[1])};
  std::unordered_set<FR> output;

  mt.compute_complex(seed_points, triangulation, 2, oracle, output);
  std::cout << "Number of simplices in the output = " << output.size() << "\n";
  return 0;
}

