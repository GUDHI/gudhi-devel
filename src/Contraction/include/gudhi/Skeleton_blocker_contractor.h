/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - 2019/08 Vincent Rouvreau: Fix issue #10 for CGAL
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef SKELETON_BLOCKER_CONTRACTOR_H_
#define SKELETON_BLOCKER_CONTRACTOR_H_

#include <gudhi/Contraction/Edge_profile.h>
#include <gudhi/Contraction/policies/Cost_policy.h>
#include <gudhi/Contraction/policies/Edge_length_cost.h>
#include <gudhi/Contraction/policies/Placement_policy.h>
#include <gudhi/Contraction/policies/First_vertex_placement.h>
#include <gudhi/Contraction/policies/Valid_contraction_policy.h>
#include <gudhi/Contraction/policies/Dummy_valid_contraction.h>  // xxx remove
#include <gudhi/Contraction/policies/Link_condition_valid_contraction.h>  // xxx remove
#include <gudhi/Contraction/policies/Contraction_visitor.h>

#include <gudhi/Skeleton_blocker/Skeleton_blocker_complex_visitor.h>
#include <gudhi/Debug_utils.h>

// todo remove the queue to be independent from cgal
#include <CGAL/Modifiable_priority_queue.h>
#include <CGAL/version.h>  // for CGAL_VERSION_NR

#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>

#include <memory>
#include <cassert>
#include <list>
#include <utility>  // for pair
#include <vector>

// Make compilation fail - required for external projects - https://github.com/GUDHI/gudhi-devel/issues/10
#if CGAL_VERSION_NR < 1041101000
# error Skeleton_blocker_contractor is only available for CGAL >= 4.11
#endif

namespace Gudhi {

namespace contraction {

template <class Profile>
Placement_policy<Profile>* make_first_vertex_placement() {
  return new First_vertex_placement<Profile>();
}

template <class Profile>
Valid_contraction_policy<Profile>* make_link_valid_contraction() {
  return new Link_condition_valid_contraction<Profile>();
}

/**
 *@brief Visitor to remove popable blockers after an edge contraction.
 */
template <class Profile>
class Contraction_visitor_remove_popable : public Contraction_visitor<Profile> {
 public:
  typedef typename Profile::Point Point;
  typedef typename Profile::Complex Complex;
  typedef typename Complex::Vertex_handle Vertex_handle;

  void on_contracted(const Profile &profile, boost::optional< Point > placement) override {
    profile.complex().remove_all_popable_blockers(profile.v0_handle());
  }
};

template <class Profile>
Contraction_visitor<Profile>* make_remove_popable_blockers_visitor() {
  return new Contraction_visitor_remove_popable<Profile>();
}

/**
 *@class Skeleton_blocker_contractor
 *@brief Class that allows to contract iteratively edges of a simplicial complex.
 *@ingroup contr
 *
 * @details The simplification algorithm consists in iteratively picking the
 * edge with lowest cost and performing an edge contraction if the contraction is valid.
 * This class is policy based (and much inspired from the edge collapse package of CGAL http://doc.cgal.org/latest/Surface_mesh_simplification/index.html).
 *
 * Policies that can be changed are :
 *  - the cost policy : how much cost an edge contraction
 *  - the placement policy : where will be placed the contraction point
 *  - the valid contraction policy : is the contraction valid. For instance, it can be
 *  a topological condition (link condition) or a geometrical condition (normals messed up).
 *
 */
template<class GeometricSimplifiableComplex, class EdgeProfile = Edge_profile<GeometricSimplifiableComplex>>
class Skeleton_blocker_contractor : private skeleton_blocker::Dummy_complex_visitor<
typename GeometricSimplifiableComplex::Vertex_handle> {
  GeometricSimplifiableComplex& complex_;

 public:
  typedef typename GeometricSimplifiableComplex::Graph_vertex Graph_vertex;
  typedef typename GeometricSimplifiableComplex::Vertex_handle Vertex_handle;
  typedef typename GeometricSimplifiableComplex::Simplex Simplex;
  typedef typename GeometricSimplifiableComplex::Root_vertex_handle Root_vertex_handle;
  typedef typename GeometricSimplifiableComplex::Graph_edge Graph_edge;
  typedef typename GeometricSimplifiableComplex::Edge_handle Edge_handle;
  typedef typename GeometricSimplifiableComplex::Point Point;

  typedef EdgeProfile Profile;


  typedef Cost_policy<Profile> Cost_policy_;
  typedef Placement_policy<Profile> Placement_policy_;
  typedef Valid_contraction_policy<Profile> Valid_contraction_policy_;
  typedef Contraction_visitor<EdgeProfile> Contraction_visitor_;
  typedef Edge_profile_factory<EdgeProfile> Edge_profile_factory_;



  typedef boost::optional<double> Cost_type;
  typedef boost::optional<Point> Placement_type;

  typedef size_t size_type;

  typedef Skeleton_blocker_contractor Self;

 private:
  struct Compare_id {
    Compare_id() : algorithm_(0) { }

    Compare_id(Self const* aAlgorithm) : algorithm_(aAlgorithm) { }

    bool operator()(Edge_handle a, Edge_handle b) const {
      return algorithm_->get_undirected_edge_id(a) < algorithm_->get_undirected_edge_id(b);
    }

    Self const* algorithm_;
  };

  struct Compare_cost {
    Compare_cost() : algorithm_(0) { }

    Compare_cost(Self const* aAlgorithm) : algorithm_(aAlgorithm) { }

    bool operator()(Edge_handle a, Edge_handle b) const {
      // NOTE: A cost is an optional<> value.
      // Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
      // In consequence, edges with undefined costs will be promoted to the top of the priority queue and popped out
      // first.
      return algorithm_->get_data(a).cost() < algorithm_->get_data(b).cost();
    }

    Self const* algorithm_;
  };

  struct Undirected_edge_id : boost::put_get_helper<size_type, Undirected_edge_id> {
    typedef boost::readable_property_map_tag category;
    typedef size_type value_type;
    typedef size_type reference;
    typedef Edge_handle key_type;

    Undirected_edge_id() : algorithm_(0) { }

    Undirected_edge_id(Self const* aAlgorithm) : algorithm_(aAlgorithm) { }

    size_type operator[](Edge_handle e) const {
      return algorithm_->get_undirected_edge_id(e);
    }

    Self const* algorithm_;
  };

  typedef CGAL::Modifiable_priority_queue<Edge_handle, Compare_cost, Undirected_edge_id> PQ;
  typedef typename PQ::handle pq_handle;


  // An Edge_data is associated with EVERY edge in the complex (collapsible or not).
  // It relates the edge with the PQ-handle needed to update the priority queue
  // It also relates the edge with a policy-based cache

  class Edge_data {
   public:
    Edge_data() : PQHandle_(), cost_() { }

    Cost_type const& cost() const {
      return cost_;
    }

    Cost_type & cost() {
      return cost_;
    }

    pq_handle PQ_handle() const {
      return PQHandle_;
    }

    bool is_in_PQ() const {
      return PQHandle_ != PQ::null_handle();
    }

    void set_PQ_handle(pq_handle h) {
      PQHandle_ = h;
    }

    void reset_PQ_handle() {
      PQHandle_ = PQ::null_handle();
    }

   private:
    pq_handle PQHandle_;
    Cost_type cost_;
  };
  typedef Edge_data* Edge_data_ptr;
  typedef boost::scoped_array<Edge_data> Edge_data_array;

  int get_undirected_edge_id(Edge_handle edge) const {
    return complex_[edge].index();
  }

  const Edge_data& get_data(Edge_handle edge) const {
    return edge_data_array_[get_undirected_edge_id(edge)];
  }

  Edge_data& get_data(Edge_handle edge) {
    return edge_data_array_[get_undirected_edge_id(edge)];
  }

  Cost_type get_cost(const Profile & profile) const {
    return (*cost_policy_)(profile, get_placement(profile));
  }

  Profile create_profile(Edge_handle edge) const {
    if (edge_profile_factory_)
      return edge_profile_factory_->make_profile(complex_, edge);
    else
      return Profile(complex_, edge);
  }

  void insert_in_PQ(Edge_handle edge, Edge_data& data) {
    data.set_PQ_handle(heap_PQ_->push(edge));
    ++current_num_edges_heap_;
  }

  void update_in_PQ(Edge_handle edge, Edge_data& data) {
    data.set_PQ_handle(heap_PQ_->update(edge, data.PQ_handle()));
  }

  void remove_from_PQ(Edge_handle edge, Edge_data& data) {
    data.set_PQ_handle(heap_PQ_->erase(edge, data.PQ_handle()));
    --current_num_edges_heap_;
  }

  boost::optional<Edge_handle> pop_from_PQ() {
    boost::optional<Edge_handle> edge = heap_PQ_->extract_top();
    if (edge)
      get_data(*edge).reset_PQ_handle();
    return edge;
  }

 private:
  /**
   * @brief Collect edges.
   *
   * Iterates over all edges of the simplicial complex and
   * 1) inserts them in the priority queue sorted according to the Cost policy.
   * 2) set the id() field of each edge
   */
  void collect_edges() {
    //
    // Loop over all the edges in the complex in the heap
    //
    size_type size = complex_.num_edges();
    DBG("Collecting edges ...");
    DBGMSG("num edges :", size);

    edge_data_array_.reset(new Edge_data[size]);

    heap_PQ_.reset(new PQ(size, Compare_cost(this), Undirected_edge_id(this)));

    std::size_t id = 0;

    // xxx do a parralel for
    for (auto edge : complex_.edge_range()) {
      complex_[edge].index() = id++;
      Profile const& profile = create_profile(edge);
      Edge_data& data = get_data(edge);
      data.cost() = get_cost(profile);
      ++initial_num_edges_heap_;
      insert_in_PQ(edge, data);
      if (contraction_visitor_) contraction_visitor_->on_collected(profile, data.cost());
    }

    DBG("Edges collected.");
  }

  bool should_stop(double lCost, const Profile &profile) const {
    return false;
  }

  boost::optional<Point> get_placement(const Profile& profile) const {
    return (*placement_policy_)(profile);
  }

  bool is_contraction_valid(Profile const& profile, Placement_type placement) const {
    if (!valid_contraction_policy_) return true;
    return (*valid_contraction_policy_)(profile, placement);
  }


 public:
  /**
   * \brief Contract edges.
   *
   * While the heap is not empty, it extracts the edge with the minimum
   * cost in the heap then try to contract it.
   * It stops when the Stop policy says so or when the number of contractions
   * given by 'num_max_contractions' is reached (if this number is positive).
   */
  void contract_edges(int num_max_contractions = -1) {
    DBG("\n\nContract edges");
    int num_contraction = 0;

    bool unspecified_num_contractions = (num_max_contractions == -1);
    //
    // Pops and processes each edge from the PQ
    //
    boost::optional<Edge_handle> edge;
    while ((edge = pop_from_PQ()) && ((num_contraction < num_max_contractions) || (unspecified_num_contractions))) {
      Profile const& profile = create_profile(*edge);
      Cost_type cost(get_data(*edge).cost());
      if (contraction_visitor_) contraction_visitor_->on_selected(profile, cost, 0, 0);

      DBGMSG("\n\n---- Pop edge - num vertices :", complex_.num_vertices());

      if (cost) {
        DBGMSG("sqrt(cost):", std::sqrt(*cost));
        if (should_stop(*cost, profile)) {
          if (contraction_visitor_) contraction_visitor_->on_stop_condition_reached();
          DBG("should_stop");
          break;
        }
        Placement_type placement = get_placement(profile);
        if (is_contraction_valid(profile, placement) && placement) {
          DBG("contraction_valid");
          contract_edge(profile, placement);
          ++num_contraction;
        } else {
          DBG("contraction not valid");
          if (contraction_visitor_) contraction_visitor_->on_non_valid(profile);
        }
      } else {
        DBG("uncomputable cost");
      }
    }
    if (contraction_visitor_) contraction_visitor_->on_stop_condition_reached();
  }

  bool is_in_heap(Edge_handle edge) const {
    if (heap_PQ_->empty()) {
      return false;
    } else {
      return edge_data_array_[get_undirected_edge_id(edge)].is_in_PQ();
    }
  }

  bool is_heap_empty() const {
    return heap_PQ_->empty();
  }

  /**
   * @brief Returns an Edge_handle and a Placement_type. This pair consists in
   * the edge with the lowest cost in the heap together with its placement.
   * The returned value is initialized iff the heap is non-empty.
   */
  boost::optional<std::pair<Edge_handle, Placement_type > > top_edge() {
    boost::optional<std::pair<Edge_handle, Placement_type > > res;

    if (!heap_PQ_->empty()) {
      auto edge = heap_PQ_->top();
      Profile const& profile = create_profile(edge);
      Placement_type placement = get_placement(profile);
      res = std::make_pair(edge, placement);
      DBGMSG("top edge:", complex_[edge]);
    }
    return res;
  }

  /**
   * @brief Constructor with default policies.
   *
   * @details The default cost, placement, valid and visitor policies
   * are respectively : the edge length, the first point, the link condition
   */
  Skeleton_blocker_contractor(GeometricSimplifiableComplex& complex)
      : complex_(complex),
      cost_policy_(new Edge_length_cost<Profile>),
      placement_policy_(new First_vertex_placement<Profile>),
      valid_contraction_policy_(new Link_condition_valid_contraction<Profile>),
      contraction_visitor_(new Contraction_visitor_()),
      edge_profile_factory_(0),
      initial_num_edges_heap_(0),
      current_num_edges_heap_(0) {
    complex_.set_visitor(this);
    if (contraction_visitor_) contraction_visitor_->on_started(complex);
    collect_edges();
  }

  /**
   * @brief Constructor with customed policies.
   * @remark Policies destruction is handle by the class with smart pointers.
   */
  Skeleton_blocker_contractor(GeometricSimplifiableComplex& complex,
                              Cost_policy_ *cost_policy,
                              Placement_policy_ * placement_policy = new First_vertex_placement<Profile>,
                              Valid_contraction_policy_ * valid_contraction_policy =
                              new Link_condition_valid_contraction<Profile>,
                              Contraction_visitor_* contraction_visitor = new Contraction_visitor_(),
                              Edge_profile_factory_* edge_profile_factory = NULL) :
      complex_(complex),
      cost_policy_(cost_policy),
      placement_policy_(placement_policy),
      valid_contraction_policy_(valid_contraction_policy),
      contraction_visitor_(contraction_visitor),
      edge_profile_factory_(edge_profile_factory),
      initial_num_edges_heap_(0),
      current_num_edges_heap_(0) {
    complex_.set_visitor(this);
    if (contraction_visitor) contraction_visitor->on_started(complex);
    collect_edges();
  }

  ~Skeleton_blocker_contractor() {
    complex_.set_visitor(0);
  }

 private:
  void contract_edge(const Profile& profile, Placement_type placement) {
    if (contraction_visitor_) contraction_visitor_->on_contracting(profile, placement);

    assert(complex_.contains_vertex(profile.v0_handle()));
    assert(complex_.contains_vertex(profile.v1_handle()));
    assert(placement);

    profile.complex().point(profile.v0_handle()) = *placement;

    // remark : this is not necessary since v1 will be deactivated
    // profile.complex().point(profile.v1_handle()) = *placement;

    complex_.contract_edge(profile.v0_handle(), profile.v1_handle());

    assert(complex_.contains_vertex(profile.v0_handle()));
    assert(!complex_.contains_vertex(profile.v1_handle()));

    update_changed_edges();

    // the visitor could do something as complex_.remove_popable_blockers();
    if (contraction_visitor_) contraction_visitor_->on_contracted(profile, placement);
  }

 private:
  // every time the visitor's method on_changed_edge is called, it adds an
  // edge to changed_edges_
  std::vector< Edge_handle > changed_edges_;

  /**
   * @brief we update the cost and the position in the heap of an edge that has
   * been changed
   */
  inline void on_changed_edge(Vertex_handle a, Vertex_handle b) override {
    boost::optional<Edge_handle> ab(complex_[std::make_pair(a, b)]);
    assert(ab);
    changed_edges_.push_back(*ab);
  }

  void update_changed_edges() {
    // xxx do a parralel for
    DBG("update edges");

    // sequential loop
    for (auto ab : changed_edges_) {
      // 1-get the Edge_handle corresponding to ab
      // 2-change the data in mEdgeArray[ab.id()]
      // 3-update the heap
      Edge_data& data = get_data(ab);
      Profile const& profile = create_profile(ab);
      data.cost() = get_cost(profile);
      if (data.is_in_PQ()) {
        update_in_PQ(ab, data);
      } else {
        insert_in_PQ(ab, data);
      }
    }
    changed_edges_.clear();
  }


 private:
  void on_remove_edge(Vertex_handle a, Vertex_handle b) override {
    boost::optional<Edge_handle> ab((complex_[std::make_pair(a, b)]));
    assert(ab);
    Edge_data& lData = get_data(*ab);
    if (lData.is_in_PQ()) {
      remove_from_PQ(*ab, lData);
    }
  }

 private:
  /**
   * @brief Called when the edge 'ax' has been added while the edge 'bx'
   * is still there but will be removed on next instruction.
   * We assign the index of 'bx' to the edge index of 'ax'
   */
  void on_swaped_edge(Vertex_handle a, Vertex_handle b, Vertex_handle x) override {
    boost::optional<Edge_handle> ax(complex_[std::make_pair(a, x)]);
    boost::optional<Edge_handle> bx(complex_[std::make_pair(b, x)]);
    assert(ax && bx);
    complex_[*ax].index() = complex_[*bx].index();
  }

 private:
  /**
   * @brief Called when a blocker is removed.
   * All the edges that passes through the blocker may be edge-contractible
   * again and are thus reinserted in the heap.
   */
  void on_delete_blocker(const Simplex * blocker) override {
    // we go for all pairs xy that belongs to the blocker
    // note that such pairs xy are necessarily edges of the complex
    // by definition of a blocker

    // todo uniqument utile pour la link condition
    // laisser a l'utilisateur ? booleen update_heap_on_removed_blocker?
    Simplex blocker_copy(*blocker);
    for (auto x = blocker_copy.begin(); x != blocker_copy.end(); ++x) {
      for (auto y = x; ++y != blocker_copy.end();) {
        auto edge_descr(complex_[std::make_pair(*x, *y)]);
        assert(edge_descr);
        Edge_data& data = get_data(*edge_descr);
        Profile const& profile = create_profile(*edge_descr);
        data.cost() = get_cost(profile);

        // If the edge is already in the heap
        // its priority has not changed.
        // If the edge is not present, we reinsert it
        // remark : we could also reinsert the edge
        // only if it is valid
        if (!data.is_in_PQ()) {
          insert_in_PQ(*edge_descr, data);
        }
      }
    }
  }


 private:
  std::shared_ptr<Cost_policy_> cost_policy_;
  std::shared_ptr<Placement_policy_> placement_policy_;
  std::shared_ptr<Valid_contraction_policy_> valid_contraction_policy_;
  std::shared_ptr<Contraction_visitor_> contraction_visitor_;

  // in case the user wants to do something special when the edge profile
  // are created (for instance add some info)
  std::shared_ptr<Edge_profile_factory_> edge_profile_factory_;
  Edge_data_array edge_data_array_;

  boost::scoped_ptr<PQ> heap_PQ_;
  int initial_num_edges_heap_;
  int current_num_edges_heap_;
};

}  // namespace contraction

}  // namespace Gudhi

#endif  // SKELETON_BLOCKER_CONTRACTOR_H_
