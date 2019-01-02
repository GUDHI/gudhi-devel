#ifndef COXETER_TRIANGULATION_DS_H_
#define COXETER_TRIANGULATION_DS_H_

#include <stack>
#include <map>

#include <gudhi/Coxeter_triangulation/Cell_id.h>
#include <gudhi/Coxeter_triangulation/Subset_chooser.h>

#include <boost/range/iterator_range.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/property_map/property_map.hpp>

#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/SVD>

#include "comb/mixedradix-init.cc"
#include "comb/print-mixedradix.cc"
#include "comb/print-setpart.cc"
#include "comb/print-perm.cc"
#include "comb/mixedradix-sod-lex.h"
#include "comb/setpart-p-rgs-lex.h"
#include "comb/perm-lex.h"
#include "fxtalloca.h"

namespace Gudhi {
  
class Coxeter_triangulation_ds {

  typedef Eigen::MatrixXd Matrix;
  typedef Eigen::SparseMatrix<double> SparseMatrix;
  typedef Eigen::Triplet<double> Triplet;
  
public:

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Constructors
  //////////////////////////////////////////////////////////////////////////////////////////////////

  
  Coxeter_triangulation_ds() {
  }

  Coxeter_triangulation_ds(unsigned dimension)
    : root_t_(dimension, dimension), dimension_(dimension){
    unsigned short d = dimension;
    Matrix cartan(d,d);
    for (unsigned i = 0; i < d; i++) {
      cartan(i,i) = 1.0;
    }
    for (unsigned i = 1; i < d; i++) {
      cartan(i-1,i) = -0.5;
      cartan(i,i-1) = -0.5;
    }
    for (int i = 0; i < d; i++)
      for (int j = 0; j < d; j++)
	if (j < i-1 || j > i+1) 
	  cartan(i,j) = 0;    
    // std::cout << "cartan =" << std::endl << cartan << std::endl;
    Eigen::SelfAdjointEigenSolver<Matrix> saes(cartan);
    Eigen::VectorXd sqrt_diag(d);
    for (int i = 0; i < d; ++i)
      sqrt_diag(i) = std::sqrt(saes.eigenvalues()[i]);
    root_t_ = saes.eigenvectors()*sqrt_diag.asDiagonal();
    colpivhouseholderqr_ = root_t_.colPivHouseholderQr();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Access functions
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  unsigned dimension() const {
    return dimension_;
  }

  unsigned pos_root_count() const {
    return dimension_*(dimension_ + 1)/2;
  }

  Eigen::MatrixXd const& simple_root_matrix() const {
    return root_t_;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Alcove dimension
  //////////////////////////////////////////////////////////////////////////////////////////////////
  unsigned alcove_dimension(const Cell_id& a_id) const {
    std::size_t i = 0, j = 1, k = 0;
    unsigned return_value = dimension_;
    while (k < a_id.size()) {
      if (a_id.mask(k)) {
	std::size_t l = i + 1;
	bool lin_independent = true;
	while (l < j && lin_independent) {
	  std::size_t k1 = (l*l+l-2)/2 - i;
	  std::size_t k2 = (j*j+j-2)/2 - l;
	  lin_independent = (!a_id.mask(k1)) || (!a_id.mask(k2));
	  l++;
	}
	if (lin_independent) {
	  return_value--;
	  if (!return_value)
	    return 0;
	}
      }
      k++;
      if (i == 0) {
	j++;
	i = j - 1;
      }
      else
	i--;
    }
    return return_value;
  }

  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Point location
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  template <class Point_d>
  Cell_id locate_point(const Point_d& point, double level=1) {
    unsigned short d = point.size();
    assert(d == dimension_);
    double error = 1e-9;
    Cell_id c_id(level, d);
    Eigen::VectorXd p_vect(d);
    for (short i = 0; i < d; i++)
      p_vect(i) = point[i];
    Eigen::VectorXd scalprod_vect = root_t_ * p_vect;
    for (short i = 0; i < d; i++) {
      double root_scalprod = 0;
      for (short j = i; j >= 0; j--) {
	root_scalprod += scalprod_vect(j);
	double value = level * root_scalprod;
	if (std::abs(value - std::round(value)) >= error)
	  c_id.push_back(std::floor(level * root_scalprod), false);
	else
	  c_id.push_back(std::round(level * root_scalprod), true);
      }
    }
    c_id.set_dimension(alcove_dimension(c_id));
    return c_id;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Cartesian coordinates and barycenter
  //////////////////////////////////////////////////////////////////////////////////////////////////
  
  Eigen::VectorXd cartesian_coordinates(const Cell_id& v_id) const {
    assert(v_id.dimension() == 0);
    Eigen::VectorXd val_vector(dimension_);
    std::size_t k = 0, j = 1;
    for (; j < (unsigned)dimension_ + 1; k += j, j++)
      val_vector(j-1) = v_id.value(k) / v_id.level();
    // return root_t_.colPivHouseholderQr().solve(val_vector);
    return colpivhouseholderqr_.solve(val_vector);
  }

  Eigen::VectorXd barycenter(const Cell_id& c_id) const {
    Eigen::VectorXd res_vector(dimension_);
    for (size_t i = 0; i < dimension_; ++i)
      res_vector(i) = 0;
    for (auto v: vertex_range(c_id)) {
      res_vector += cartesian_coordinates(v);
    }
    return (1/(c_id.dimension()+1)) * res_vector;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Vertex computation
  //////////////////////////////////////////////////////////////////////////////////////////////////

  class Vertex_iterator : public boost::iterator_facade< Vertex_iterator,
							 Cell_id const,
							 boost::forward_traversal_tag> {
  protected:
    friend class boost::iterator_core_access;

    bool equal(Vertex_iterator const& other) const {
      return (is_end_ && other.is_end_);
    }

    Cell_id const& dereference() const {
      return value_;
    }

    void update_value() {
      while (!is_end_) {
	std::size_t k = value_.size();
	std::size_t j = std::floor(0.5*(1 + std::sqrt(1+8*k)));
	if (state_stack_.top())
	  value_.push_back(c_id_.value(k)+1, true);
	else
	  value_.push_back(c_id_.value(k), true);
	std::size_t i = j-2;
	for (; i < j; --i) { // i is unsigned: if goes past 0, then becomes >j
	  std::size_t
	    k_s = i*(i+1)/2, // index of the simple root s_i
	    k_r = value_.size()-1; // index of the previous positive root r_{i+1,j}
	  value_.push_back(value_.value(k_r) + value_.value(k_s), true);
	  if (c_id_.mask(k_r+1)) {
	    if (value_.value(k_r+1) != c_id_.value(k_r+1))
	      break;
	  }
	  else // !c_id_.mask(k_r+1)
	    if (value_.value(k_r+1) < c_id_.value(k_r+1) ||
		value_.value(k_r+1) > c_id_.value(k_r+1) + 1)
	    break;
	}
	if (i >= j) { // success
	  if (value_.size() == c_id_.size()) // finished
	    return;
	  state_stack_.push(false);
	}
	if (i < j) { // fail
	  value_.resize(k);
	  elementary_increment();
	}
      }
    }

    void elementary_increment() {
      if (is_end_)
        return;
      while (!state_stack_.empty()) {
	bool p = state_stack_.top();
	state_stack_.pop();
	std::size_t j = state_stack_.size();
	std::size_t k_s = j*(j+1)/2;
	value_.resize(k_s);
	if (c_id_.mask(k_s))
	  continue;
	else if (p)
	  continue;
	else {
	  state_stack_.push(true);
	  return;
	}
      }
      if (state_stack_.empty())
	is_end_ = true;
    }
    
    void increment() {
      is_end_ = is_itself_;
      elementary_increment();
      update_value();
    }

  public:
    Vertex_iterator(const Cell_id& c_id,
		    const Coxeter_triangulation_ds& scs)
      :
      value_(c_id.level(), 0),
      is_end_(false),
      is_itself_(c_id.dimension() == 0)
    {
      if (is_itself_) {
	value_ = c_id;
	return;
      }
      c_id_ = c_id;
      state_stack_.push(false);
      update_value();
    }

    Vertex_iterator() : is_end_(true) {}
  
  
  protected:
    Cell_id c_id_; 
    Cell_id value_;
    bool is_end_;
    bool is_itself_;
    std::stack<bool> state_stack_; // true means +1, false means +0
  };
  
  typedef boost::iterator_range<Vertex_iterator> Vertex_range;
  Vertex_range vertex_range(const Cell_id& a_id) const {
    return Vertex_range(Vertex_iterator(a_id, *this),
			Vertex_iterator());
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Face computation
  //////////////////////////////////////////////////////////////////////////////////////////////////

  class Face_iterator : public boost::iterator_facade< Face_iterator,
						       Cell_id const,
						       boost::forward_traversal_tag> {
  protected:
    friend class boost::iterator_core_access;
    
    void update_value() {
      if (is_end_)
        return;
      value_.clear();
      for (unsigned k = 0; k < vertices_[state_[0]].size(); ++k) {
        int val = vertices_[state_[0]].value(k);
        bool is_true = true;
        for (auto p: state_)
          if (vertices_[p].value(k) < val) {
            val = vertices_[p].value(k);
            is_true = false;
            break;
          }
          else if (vertices_[p].value(k) > val) {
            is_true = false;
            break;
          }
        value_.push_back(val, is_true);
      }
    }
    
    bool equal(Face_iterator const& other) const {
      return (is_end_ && other.is_end_); 
    }

    Cell_id const& dereference() const {
      return value_;
    }

    void increment() {
      is_end_ = is_itself_ || !increment_state();
      update_value();
    }

    bool increment_state() {
      do {
        if (state_.empty()) {
          state_.push_back(0);
        } else {
          if (state_.size() < state_.capacity()) {
            state_.push_back(state_.back());
          }
          
          // Roll over when the remaining elements wouldn't fill the subset.
          while (++state_.back() == vertices_.size()) {
            state_.pop_back();
            if (state_.empty())
              // we have run out of possibilities
              return false;
          }
        }
      } while (state_.size() < state_.capacity());
      return true;
    }
    
  public:
    Face_iterator(const Cell_id& c_id,
		  const Coxeter_triangulation_ds& scs,
		  std::size_t value_dim)
      : value_(c_id.level(), value_dim),
	is_end_(false),
	is_itself_(value_dim == c_id.dimension())
    {
      if (is_itself_) {
	value_ = c_id;
	return;
      }
      for (const Cell_id& v_id: scs.vertex_range(c_id))
        vertices_.push_back(v_id);
      state_.reserve(value_dim + 1);
      if (value_dim > c_id.dimension())
        is_end_ = true;
      is_end_ = !increment_state();
      update_value();
    }

    Face_iterator() : is_end_(true) {}
    
    
  protected:
    Cell_id value_;
    std::vector<Cell_id> vertices_;
    SubsetChooser<std::vector<Cell_id>::const_iterator> chooser_;
    std::vector<std::size_t> state_;
    bool is_end_;
    bool is_itself_;
  };

  
  typedef boost::iterator_range<Face_iterator> Face_range;
  Face_range face_range(const Cell_id& c_id,
			std::size_t value_dim) const {
    return Face_range(Face_iterator(c_id, *this, value_dim),
		      Face_iterator());
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Coface computation
  //////////////////////////////////////////////////////////////////////////////////////////////////

  class Coface_iterator : public boost::iterator_facade< Coface_iterator,
							 Cell_id const,
							 boost::forward_traversal_tag> {
  protected:
    friend class boost::iterator_core_access;
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Graph_node;
    typedef boost::graph_traits<Graph>::edge_descriptor Graph_edge;
    typedef boost::graph_traits<Graph>::adjacency_iterator Adj_it;
    typedef std::pair<Adj_it, Adj_it> Out_edge_it;
    typedef std::vector<std::size_t> Rank;
    typedef std::vector<std::size_t> Parent;
    typedef boost::disjoint_sets<std::size_t*, std::size_t*> Disjoint_sets;
    typedef std::map<int, Graph_node> IN_map;
    typedef typename IN_map::iterator IN_map_iterator;
    typedef std::map<Graph_node, int> NI_map;
    typedef std::vector<std::vector<std::size_t> > Partition;
    
    void update_value() {
      if (is_end_)
        return;
      value_.clear();
      std::vector<bool> mask_vect(v_id_.size(), false);
      // std::vector<long int> val_vect(v_id_.size(), 0);
      std::vector<long int> val_vect(v_id_.begin(), v_id_.end());
      std::vector<std::vector<std::size_t> > refined_partition;
      std::size_t offset = 0;
      for (std::size_t n = 0; n < nonrefined_partition_.size(); ++n) {
	for (std::size_t m = 0; m < msl_->data()[n] + 1; ++m)
	  refined_partition.push_back(std::vector<std::size_t>());
	if (msl_->data()[n] == 0)
	  for (std::size_t l = 0; l < nonrefined_partition_[n].size(); ++l)
	    refined_partition[offset + sprl_vector_[n]->data()[l]].push_back(nonrefined_partition_[n][l]);
	else
	  for (std::size_t l = 0; l < nonrefined_partition_[n].size(); ++l)
	    refined_partition[offset + pl_vector_[n]->data()[sprl_vector_[n]->data()[l]]].push_back(nonrefined_partition_[n][l]);
	  // std::cout << "sprl_vector[" << n << "].data()[" << l << "] = " << sprl_vector_[n]->data()[l] << ", nonrefined_partition_[" << n << "][" << l << "] = " << nonrefined_partition_[n][l] << "\n";
	offset += msl_->data()[n] + 1;
      }
      // std::cout << "refined_partition = " << refined_partition << "\n";
      std::size_t i, j, k;
      for (std::size_t n1 = 0; n1 < refined_partition.size(); ++n1)
	for (std::size_t m1 = 0; m1 < refined_partition[n1].size(); ++m1) 
	  for (std::size_t n2 = n1; n2 < refined_partition.size(); ++n2) 
	    for (std::size_t m2 = (n2 == n1 ? m1+1 : 0); m2 < refined_partition[n2].size(); ++m2) {
	      std::tie(i, j) = std::minmax(refined_partition[n1][m1], refined_partition[n2][m2]);
	      k = (j*j+j-2)/2 - i;
	      if (n1 == n2)
		mask_vect[k] = true;
	      else if (refined_partition[n1][m1] > refined_partition[n2][m2])
		val_vect[k]--;
	    }
      for (k = 0; k < v_id_.size(); ++k)
	value_.push_back(val_vect[k], mask_vect[k]);
      // std::cout << "v_id_ = " << v_id_ << ", value_ = " << value_ << "\n\n";
    }
    
    bool equal(Coface_iterator const& other) const {
      return (is_end_ && other.is_end_); 
    }

    Cell_id const& dereference() const {
      return value_;
    }

    void increment() {
      is_end_ = is_itself_;
      while (true) {
	if (pl_vector_.size() == 0) {
	  if (!msl_->next())
	    is_end_ = true;
	  break;
	}
	std::size_t n = pl_vector_.size()-1;
	if (!msl_->data()[n] || !pl_vector_[n]->next()) {
	  delete pl_vector_[n];
	  pl_vector_.pop_back();
	  if (!sprl_vector_[n]->next()) {
	    delete sprl_vector_[n];
	    sprl_vector_.pop_back();
	    continue;
	  }
	  else
	    pl_vector_.push_back(new perm_lex(msl_->data()[n]+1));
	}
	else
	  break;
      }
      if (!is_end_)
	while (pl_vector_.size() != nonrefined_partition_.size()) {
	  std::size_t n = pl_vector_.size();
	  sprl_vector_.push_back(new setpart_p_rgs_lex(nonrefined_partition_[n].size(), msl_->data()[n]+1));
	  // sprl_vector_.back()->print_sets("    "); std::cout << "\n";
	  pl_vector_.push_back(new perm_lex(msl_->data()[n]+1));
	  // pl_vector_.back()->print("    "); std::cout << "\n";
	}      
      update_value();
    }

    Partition nonrefined_partition(const Cell_id& c_id, const Cell_id& v_id) {
      Graph graph;
      IN_map in_map;
      NI_map ni_map;
      std::size_t k = 0;
      
      for (std::size_t j = 0; k < c_id.size(); ++j) {
      	Graph_node node = boost::add_vertex(graph);
      	in_map.emplace(std::make_pair(j, node));
      	ni_map.emplace(std::make_pair(node, j));
      	for (std::size_t i = j-1; i < j; --i, ++k)
      	  if (!c_id.mask(k)) {
      	    if (c_id.value(k) == v_id.value(k))
	      boost::add_edge(j, i, graph);
	    else
	      boost::add_edge(i, j, graph);
	  }
      }
      std::vector<IN_map_iterator> iterators;
      for (auto it = in_map.begin(); it != in_map.end(); ++it)
	iterators.push_back(it);
      auto compare = [graph](const IN_map_iterator& lhs, const IN_map_iterator& rhs)
		     {return boost::out_degree(lhs->second, graph) < boost::out_degree(rhs->second, graph);};
      std::sort(iterators.begin(), iterators.end(), compare);
      std::size_t curr_degree = 0;
      Partition output;
      if (!iterators.empty())
	output.push_back(std::vector<std::size_t>());
      for (auto it: iterators) {
	std::size_t node_degree = boost::out_degree(it->second, graph);
	if (node_degree == curr_degree)
	  output.back().push_back(it->first);
	else {
	  curr_degree = node_degree;
	  output.push_back(std::vector<std::size_t>(1, it->first));
	}
      }
      // std::cout << output << "\n";
      return output;
    }
    
  public:
    Coface_iterator(const Cell_id& c_id,
		    const Coxeter_triangulation_ds& scs,
		    std::size_t value_dim)
      : value_(c_id.level(), value_dim),
	v_id_(*Coxeter_triangulation_ds::Vertex_iterator(c_id, scs)),
	nonrefined_partition_(nonrefined_partition(c_id, v_id_)),
	is_end_(false),
	is_itself_(value_dim == c_id.dimension())
    {
      if (is_itself_) {
	value_ = c_id;
	msl_ = new mixedradix_sod_lex(1, 0);
	return;
      }
      // std::size_t d = scs.dimension();
      std::size_t n = nonrefined_partition_.size();
      ALLOCA(std::size_t, radix, n); // as in 'radix' in the mixed radix
      for (std::size_t i = 0; i < nonrefined_partition_.size(); ++i) {
	radix[i] = nonrefined_partition_[i].size();
        // std::cout << i << ": " << radix[i] << "\n";
      }
      msl_ = new mixedradix_sod_lex(n, 0, radix);
      msl_->first(value_dim+1-n);
      // msl_->print_nines("Nines: "); std::cout << "\n";
      // std::cout << "n = " << n << ", value_dim = " << value_dim << "\n";
      // msl_->print("  ", true); std::cout << "\n";
      
      for (std::size_t i = 0; i < n; ++i) {
	sprl_vector_.push_back(new setpart_p_rgs_lex(radix[i], msl_->data()[i]+1));
	// sprl_vector_.back()->print_sets("    "); std::cout << "\n";
	pl_vector_.push_back(new perm_lex(msl_->data()[i]+1)); 
	// pl_vector_.back()->print("    "); std::cout << "\n";
      }
      update_value();
    }

    Coface_iterator() : is_end_(true) {
      msl_ = new mixedradix_sod_lex(1, 0);
    }

    ~Coface_iterator() {
      delete msl_;
      for (setpart_p_rgs_lex* sprl_: sprl_vector_)
	delete sprl_;
      for (perm_lex* pl_: pl_vector_)
	delete pl_;
    }
    
  protected:
    Cell_id value_;
    Cell_id v_id_;
    Partition nonrefined_partition_;
    mixedradix_sod_lex* msl_;
    std::vector<setpart_p_rgs_lex*> sprl_vector_;
    std::vector<perm_lex*> pl_vector_;
    bool is_end_;
    bool is_itself_;
  };

  
  typedef boost::iterator_range<Coface_iterator> Coface_range;
  Coface_range coface_range(const Cell_id& c_id,
			    std::size_t value_dim) const {
    return Coface_range(Coface_iterator(c_id, *this, value_dim),
			Coface_iterator());
  }  

  
protected:
  Matrix root_t_;
  unsigned dimension_;
  Eigen::ColPivHouseholderQR<Matrix> colpivhouseholderqr_;
};

    
}

#endif
