/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONTRACTION_POLICIES_CONTRACTION_VISITOR_H_
#define CONTRACTION_POLICIES_CONTRACTION_VISITOR_H_

#include <gudhi/Contraction/Edge_profile.h>
#include <boost/optional.hpp>

namespace Gudhi {

namespace contraction {

/**
 *@class Contraction_visitor
 *@brief Interface for a visitor of the edge contraction process.
 *@ingroup contr
 */
template <typename EdgeProfile>
class Contraction_visitor {  // : public Dummy_complex_visitor<typename EdgeProfile::Vertex_handle>  {
 public:
  // typedef typename ComplexType::GeometryTrait GT;
  typedef EdgeProfile Profile;
  typedef double FT;
  typedef typename Profile::Complex Complex;
  typedef Complex ComplexType;
  typedef typename ComplexType::Point Point;

  virtual ~Contraction_visitor() { }

  /**
   * @brief Called before the edge contraction process starts.
   */
  virtual void on_started(ComplexType & complex) { }

  /**
   * @brief Called when the algorithm stops.
   */
  virtual void on_stop_condition_reached() { }

  /**
   * @brief Called during the collecting phase (when a cost is assigned to the edges), for each edge collected.
   */
  virtual void on_collected(const Profile &profile, boost::optional< FT > cost) { }

  /**
   * @brief Called during the processing phase (when edges are contracted), for each edge that is selected.
   */
  virtual void on_selected(const Profile &profile, boost::optional< FT > cost, int initial_count, int current_count) { }

  /**
   * @brief Called when an edge is about to be contracted and replaced by a vertex whose position is *placement.
   */
  virtual void on_contracting(const Profile &profile, boost::optional< Point > placement) { }

  /**
   * @brief Called when after an edge has been contracted onto a new point placement.
   * A possibility would to remove popable blockers at this point for instance.
   */
  virtual void on_contracted(const Profile &profile, boost::optional< Point > placement) { }

  /**
   * @brief Called for each selected edge which cannot be contracted because the ValidContractionPredicate is false
   */
  virtual void on_non_valid(const Profile &profile) { }
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // CONTRACTION_POLICIES_CONTRACTION_VISITOR_H_
