// SPDX-License-Identifier: GPL-3.0-or-later
// Author(s)    : Marc Glisse

#include <CGAL/Delaunay_triangulation.h>
#include <iterator>
#include <vector>

namespace Gudhi {

/** \ingroup alpha_complex
 *
 * This function constructs the incremental Delaunay complex, as defined in \cite incremental-delaunay .
 * Given a list of points (p1, ..., pn), the incremental Delaunay complex contains the Delaunay triangulation
 * of (p1, ..., pk) for all k≤n. For any simplex s of the triangulation of (p1, ..., pk) that is not in the
 * triangulation of (p1, ..., pn), the complex also contains a simplex formed with the same points as s plus the
 * first pj that conflicts with s (is inside its circumsphere). For points in R^d, the complex usually has
 * dimension d+1.
 *
 * This function calls `complex.insert_simplex_and_subfaces` on a set of simplices large enough that, with their
 * faces, they define the incremental Delaunay complex. However, this set is not supposed to be minimal, there
 * may be redundancies.
 *
 * No filtration values are computed here, those can be computed separately afterwards, for instance
 * with `Gudhi::cech_complex::assign_MEB_filtration`.
 *
 * @tparam K CGAL kernel: either Epick_d or Epeck_d.
 * @tparam SimplicialComplex Simplicial complex. It must provide a type `Vertex_handle` constructible from `int`
 *   to name vertices, and a member function `void insert_simplex_and_subfaces(vrange)` which accepts a range
 *   of `Vertex_handle` (representing a simplex).
 * @tparam PointRange Forward range of `K::Point_d`.
 *
 * @param[in] k The geometric kernel.
 * @param[out] complex The simplicial complex.
 * @param[in] points Embedding of the points.
 */
template<typename K, typename SimplicialComplex, typename PointRange>
void construct_incremental_delaunay(K const&k, SimplicialComplex& complex, PointRange const& points) {
  using TDS = CGAL::Triangulation_data_structure<
    typename K::Dimension,
    CGAL::Triangulation_vertex<K, typename SimplicialComplex::Vertex_handle>,
    CGAL::Triangulation_ds_full_cell<void, CGAL::TDS_full_cell_mirror_storage_policy> >;
  using Triangulation = CGAL::Delaunay_triangulation<K, TDS>;
  auto pt_it = std::begin(points);
  auto pt_end = std::end(points);
  if(pt_it == pt_end) return;
  int d = k.point_dimension_d_object()(*pt_it);
  Triangulation tri(d);

  // local variables pulled out for recycling
  typedef std::vector<typename Triangulation::Full_cell_handle> Full_cell_h_vector;
  Full_cell_h_vector cs;
  std::vector<typename SimplicialComplex::Vertex_handle> splx;

  for(int idx = 0; pt_it != pt_end; ++pt_it, ++idx) {
    // Adapted from Delaunay_triangulation::insert(Point)
    auto&& p = *pt_it;
    typename Triangulation::Locate_type lt;
    typename Triangulation::Face f(d);
    typename Triangulation::Facet ft;
    typename Triangulation::Full_cell_handle s = tri.locate(p, lt, f, ft);
    // Adapted from Delaunay_triangulation::insert(Point,Locate_type,Face,Facet,Full_cell_handle)
    if (lt == Triangulation::OUTSIDE_AFFINE_HULL) {
      // nothing disappears
      tri.insert_outside_affine_hull(p)->data() = idx;
    }
    else if (lt == Triangulation::ON_VERTEX) {
      // nothing disappears
    }
    else {
      if(tri.current_dimension() == 1) {
	if(lt == Triangulation::OUTSIDE_CONVEX_HULL) {
	  // only an infinite face disappears
	  tri.insert_outside_convex_hull_1(p, s)->data() = idx;
	} else {
	  complex.insert_simplex_and_subfaces({s->vertex(0)->data(), s->vertex(1)->data(), idx});
	  typename Triangulation::Vertex_handle v = tri.tds().insert_in_full_cell(s);
	  v->set_point(p);
	  v->data() = idx;
	}
      } else { // main case
        // Adapted from Delaunay_triangulation::insert_in_conflicting_cell(Point,Full_cell_handle)
	/* Full_cell_h_vector cs; */ cs.clear();
	std::back_insert_iterator<Full_cell_h_vector> out(cs);
	typename Triangulation::Facet ft = tri.compute_conflict_zone(p, s, out);
	for(auto c : cs) {
	  /* std::vector<int> splx; */ splx.clear();
	  for(int i = 0; i <= tri.current_dimension(); ++i) {
	    // We could skip the infinite vertex and output a (new) simplex, but that's not so useful
	    if(tri.is_infinite(c->vertex(i))) goto next_conflicting_cell;
	    splx.push_back(c->vertex(i)->data());
	  }
          // We could put idx at the beginning of the vector once instead of reinserting it each time
	  splx.push_back(idx);
	  complex.insert_simplex_and_subfaces(splx);
next_conflicting_cell:;
	}
	tri.insert_in_hole(p, cs.begin(), cs.end(), ft)->data() = idx;
      }
    }
  }
  // We have added all the conflict cells, but we are missing some current cells.
  // There will be a lot of redundancy though.
  for(auto it = tri.finite_full_cells_begin(); it != tri.finite_full_cells_end(); ++it) {
    /* std::vector<int> splx; */ splx.clear();
    for(int i = 0; i <= tri.current_dimension(); i++) {
      splx.push_back(it->vertex(i)->data());
    }
    complex.insert_simplex_and_subfaces(splx);
    // Or we could call insert_simplex_and_subfaces on a transform_iterator.
  }
}

/* Strategies
 1. Current
 - every time a finite simplex disappears, output the conflict simplex
 - output the final triangulation at the end (a bit costly)
 2. Alternative from Kerber's function_delaunay
 - every time a full-dimensional simplex disappears, output the finite faces of the conflict simplex
 - this only works with some assumptions, for instance that the first d+1 points are affinely independent and
   that the simplex they define is not in the final Delaunay triangulation.
 3. Mix?
 - we could use 1. until the triangulation reaches the ambient dimension then switch to 2.
 - we could pretend the current dimension is going to be the final one, use 2, but fix things up every time
   we insert a point outside the affine hull. Probably not very interesting, since triangulations of non-full
   dimension should be rare.
*/

/* Possibilities:
 - output conflict cells, new boundary cells, and final cells separately in case a user wants to handle them differently?
 - we could remove duplicates in `points` so the numbering becomes contiguous, but then we have to output the new `points`.
 - use less undocumented functions. Calling `insert` shouldn't cost that much more.
*/
} // namespace Gudhi
