/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Clément Maria
 *
 *    Copyright (C) 2021 Inria
 *
 *    Modification(s):
 *      - 2023/05 Hannah Schreiber: Rework of the interface, reorganization and debug
 *      - 2023/05 Hannah Schreiber: Addition of infinit bars
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ZIGZAG_PERSISTENCE_H_
#define ZIGZAG_PERSISTENCE_H_

#include <boost/tuple/tuple.hpp>
#include <boost/intrusive/list.hpp>
#include <boost/intrusive/set.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/timer/progress_display.hpp>

#include <cmath>
#include <limits>
#include <map>
#include <list>
#include <ostream>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#include <gudhi/Debug_utils.h>
#include <gudhi/Simple_object_pool.h>

namespace Gudhi {
namespace zigzag_persistence {
//represent matrix columns with sets.
struct Zigzag_persistence_colset;
//----------------------------------------------------------------------------------
/** \class Zigzag_persistence Zigzag_persistence.h gudhi/Zigzag_persistence.h
  * \brief Computation of the zigzag persistent homology of a zigzag
  * filtered complex.
  *
  * \details The type ZigzagFilteredComplex::Simplex_key counts the number of
  * insertions and
  * deletions of simplices, which may be large in zigzag persistence and require
  * more than 32 bits of storage. The type used (int, long, etc) should be chosen in
  * consequence. Simplex_key must be signed.
  *
  * Over all insertions, the Simplex_key must be positive and strictly increasing
  * when forward iterating along the zigzag filtration.
  */
template < typename ZigzagFilteredComplex
		   , typename ZigzagPersistenceOptions = Zigzag_persistence_colset >
class Zigzag_persistence {
public:
	typedef ZigzagFilteredComplex                 Complex;
	typedef ZigzagPersistenceOptions              Options;
	/*** Types defined in the complex ***/
	// Data attached to each simplex to interface with a Property Map.
	typedef typename Complex::Simplex_key         Simplex_key;//must be signed
	typedef typename Complex::Simplex_handle      Simplex_handle;
	typedef typename Complex::Vertex_handle       Vertex_handle;
	typedef typename Complex::Filtration_value    Filtration_value;
	//
private:
	/*** Matrix cells and columns types ***/
	struct matrix_row_tag; // for horizontal traversal in the persistence matrix
	struct matrix_column_tag; // for vertical traversal in the persistence matrix
	typedef boost::intrusive::list_base_hook<
				boost::intrusive::tag < matrix_row_tag >              //allows .unlink()
			  , boost::intrusive::link_mode < boost::intrusive::auto_unlink >
	>														base_hook_matrix_row_list;
	//hook for a column represented by an intrusive list
	typedef boost::intrusive::list_base_hook <  //faster hook, less safe
				boost::intrusive::tag < matrix_column_tag >
			  //, boost::intrusive::link_mode < boost::intrusive::auto_unlink >
			  , boost::intrusive::link_mode < boost::intrusive::safe_link >
	>														base_hook_matrix_column_list;
	//hook for a column represented by an intrusive set
	typedef boost::intrusive::set_base_hook <  //faster hook, less safe
				boost::intrusive::tag < matrix_column_tag >
			  , boost::intrusive::optimize_size<true>
			  //, boost::intrusive::link_mode < boost::intrusive::auto_unlink >
			  , boost::intrusive::link_mode < boost::intrusive::safe_link >
	>														base_hook_matrix_column_set;
	//the data structure for columns is selected in Options::searchable_column
	typedef typename std::conditional<Options::searchable_column,
				base_hook_matrix_column_set,
				base_hook_matrix_column_list
	>::type													base_hook_matrix_column;
	//the only option for rows is the intrusive list
	typedef base_hook_matrix_row_list                       base_hook_matrix_row;

	/* Cell for the persistence matrix. Contains a key for the simplex index, and
	 * horizontal and vertical hooks for connections within sparse rows and columns.
	 */
	struct matrix_chain;//defined below, a chain contains a row, a column, and more
	/** Type of cell in the sparse homology matrix.
	 *  For now, only coefficients in Z/2Z, so only the row index (called key) is
	 * stored in the cell.
	 */
	struct Zigzag_persistence_cell
			: public base_hook_matrix_row, public base_hook_matrix_column
	{
		Zigzag_persistence_cell(Simplex_key key, matrix_chain *self_chain)
			: key_(key)
			, self_chain_(self_chain)
		{}

		Simplex_key key() const { return key_; }
		//compare by increasing key value
		friend bool operator<( const Zigzag_persistence_cell& c1
							   , const Zigzag_persistence_cell& c2) {
			return c1.key() < c2.key();
		}
		/* In a matrix M, if M[i][j] == x not 0, we represent a cell with key_=i
		 * (the row index),
		 * self_chain_ points to the chain corresponding the j-th column, and x_=x.
		 * Currently, only Z/2Z coefficients are implemented, so x=1.
		 *
		 * A cell is connected to all cells of the same row, and all cells of the same
		 * column, via the two boost::intrusive hooks (row and column).
		 */
		Simplex_key       key_;
		matrix_chain    * self_chain_;
	};

	//Homology matrix cell
	typedef Zigzag_persistence_cell												Cell;
	// Remark: constant_time_size must be false because base_hook_matrix_row and
	// base_hook_matrix_column have auto_unlink link_mode
	//vertical list of cells, forming a matrix column stored as an intrusive list
	typedef boost::intrusive::list <
				Cell
			  , boost::intrusive::constant_time_size<false>
			  , boost::intrusive::base_hook< base_hook_matrix_column_list >  >	Column_list;
	//vertical list of cells, forming a matrix column stored as an intrusive set
	typedef boost::intrusive::set <
				Cell
			  , boost::intrusive::constant_time_size<false>
			  , boost::intrusive::base_hook< base_hook_matrix_column_set >  >	Column_set;
	//choice encoded in Options::searchable_column. a column can be
	//iterated through, and keys are read in strictly increasing natural order.
	typedef typename std::conditional<
				Options::searchable_column,
				Column_set,
				Column_list >::type												Column;
	//horizontal list of cells, forming a matrix row, no particular order on keys.
	typedef boost::intrusive::list <
				Cell
			  , boost::intrusive::constant_time_size<false>
			  , boost::intrusive::base_hook< base_hook_matrix_row >  >			Row_list;
	//rows are encoded by lists. need to be sorted and traversed
	typedef Row_list															Row;

	/* Chain for zigzag persistence. A chain stores:
	 * - a matrix column (col_i) that represents the chain as a sum of simplices
	 * (represented by their unique key, stored in the cells of the sparse column),
	 * - a matrix row of all elements of index the lowest index of the column (row_i),
	 * - is paired with another chain, indicating its type F, G, or H,
	 * - has a direct access to its lowest index.
	 */
	struct matrix_chain {
		/* Trivial constructor, birth == -3 */
		matrix_chain() : column_(nullptr), row_(nullptr), paired_col_(nullptr),
			birth_(-3), lowest_idx_(-1) {}

		/* Creates a matrix chain of type F with one cell of index 'key'. */
		matrix_chain(Simplex_key key)
			: paired_col_(nullptr), birth_(key), lowest_idx_(key)
		{
			// Cell *new_cell = new Cell(key, this);
			Cell *new_cell = cellPool_.construct(key, this);
			if constexpr(Options::searchable_column) { column_.insert(*new_cell); }
			else                                     { column_.push_back(*new_cell); }
			row_.push_back(*new_cell);
		}
		/* Creates a matrix chain of type F with new cells of key indices given by a
		 * range. Birth and lowest indices are given by 'key'.
		 * The range [beg,end) must be sorted by increasing key values, the same
		 * order as the column_ when read from column_.begin() to column_.end().
		 *
		 * SimplexKeyIterator value_type must be Simplex_key.
		 * KeyToMatrixChain must be of type
		 *         std::map< Simplex_key, typename std::list<matrix_chain>::iterator >
		 */
		template< typename SimplexKeyIterator, typename KeyToMatrixChain >
		matrix_chain(Simplex_key key, SimplexKeyIterator beg, SimplexKeyIterator end, KeyToMatrixChain &lowidx_to_matidx)
			: paired_col_(nullptr), birth_(key), lowest_idx_(key)
		{
			for(SimplexKeyIterator it = beg; it != end; ++it)
			{
				// Cell *new_cell = new Cell(*it, this);//create a new cell
				Cell *new_cell = cellPool_.construct(*it, this);
				//insertion in the column
				if constexpr(Options::searchable_column) {
					column_.insert(column_.end(), *new_cell); //ordered range
				}
				else { column_.push_back(*new_cell); }
				//insertion in a row, not corresponding to the row stored in this->row_.
				lowidx_to_matidx[*it]->row_.push_back( *new_cell );
			}
			//Add the bottom coefficient for the chain
			// Cell *new_cell = new Cell(key, this);
			Cell *new_cell = cellPool_.construct(key, this);
			//insertion in the column, key is larger than any *it in [beg, end) above.
			if constexpr(Options::searchable_column) {
				column_.insert(column_.end(), *new_cell);
			}
			else { column_.push_back(*new_cell); }
			//insertion of row_, that stores no particular order.
			row_.push_back( *new_cell );
		}

		/* Creates a matrix chain of type H with new cells of key indices given by a
		 * range. Birth and lowest indices are given by 'key'.
		 * The range [beg,end) must be sorted by increasing key values, the same
		 * order as the column_ when read from column_.begin() to column_.end().
		 *
		 * SimplexKeyIterator value_type must be Simplex_key.
		 * KeyToMatrixChain must be of type
		 *             std::map< Simplex_key, typename std::list<matrix_chain>::iterator >
		 */
		template< typename SimplexKeyIterator, typename KeyToMatrixChain >
		matrix_chain(Simplex_key key, matrix_chain *paired_col, SimplexKeyIterator beg, SimplexKeyIterator end, KeyToMatrixChain &lowidx_to_matidx)
			: paired_col_(paired_col), birth_(-2), lowest_idx_(key)
		{
			for(SimplexKeyIterator it = beg; it != end; ++it)
			{
				// Cell * new_cell = new Cell(*it, this);//create a new cell
				Cell *new_cell = cellPool_.construct(*it, this);
				//insertion in the column
				if constexpr(Options::searchable_column) {
					column_.insert(column_.end(), *new_cell); //ordered range
				}
				else { column_.push_back(*new_cell); }
				//insertion in a row, not corresponding to the row stored in this->row_.
				lowidx_to_matidx[*it]->row_.push_back( *new_cell );
			}
			//Add the bottom coefficient for the chain
			// Cell * new_cell = new Cell(key, this);
			Cell *new_cell = cellPool_.construct(key, this);
			//insertion in the column, key is larger than any *it in [beg, end) above.
			if constexpr(Options::searchable_column) {
				column_.insert(column_.end(), *new_cell);
			}
			else { column_.push_back(*new_cell); }
			//insertion of row_, that stores no particular order.
			row_.push_back( *new_cell );
		}

		/* Erase the chain, all cells were allocated with operator new. */
		~matrix_chain()
		{ //empty the column, call delete on all cells
			for(typename Column::iterator c_it = column_.begin(); c_it != column_.end(); )
			{
				auto tmp_it = c_it; ++c_it;
				Cell * tmp_cell = &(*tmp_it);
				tmp_it->base_hook_matrix_row::unlink(); //rm from row
				column_.erase(tmp_it);
				cellPool_.destroy(tmp_cell);
				// delete tmp_cell;
			}
		}

		/* Returns the chain with which *this is paired in the F,G,H classification.
		 * If in F (i.e., paired with no other column), return nullptr.*/
		matrix_chain * paired_chain() const { return paired_col_; }
		/* Assign a paired chain. */
		void assign_paired_chain(matrix_chain *other_col) { paired_col_ = other_col; }
		/* Access the column. */
		Column & column() { return column_; }
		/* Returns the birth index (b >= 0) of the chain if the column is in F.
		 * Returns -2 if the chain is in H, and -1 if the chain is in G. */
		Simplex_key birth() const                       { return birth_; }
		/* Assign a birth index to the chain. */
		void assign_birth(Simplex_key b)          { birth_ = b; }
		void assign_birth(matrix_chain *other) { birth_ = other->birth_; }
		/* Returns true iff the chain is indexed in F. */
		bool inF() const { return birth_ >  -1; }
		/* Returns true iff the chain is indexed in G. */
		bool inG() const { return birth_ == -1; }
		/* Returns true iff the chain is indexed in H. */
		bool inH() const { return birth_ == -2; }
		Simplex_key lowest_idx() const { return lowest_idx_; }

		Column            column_      ; //col at index i, with lowest index i
		Row               row_         ; //row at index i
		matrix_chain    * paired_col_  ; //\in F -> nullptr, \in H -> g, \in G -> h
		Simplex_key       birth_       ; //\in F -> b, \in H -> -2 \in G -> -1
		Simplex_key       lowest_idx_  ; //lowest_idx_ = i (upper triangular matrix)
		inline static Simple_object_pool<Cell> cellPool_;
	};

public:
	/** \brief Structure to store persistence intervals by their filtration values.
	 *
	 * \details By convention, interval \f$[b;d]\f$ are
	 * closed for finite indices b and d, and open for left-infinite and/or
	 * right-infinite endpoints.*/
	struct interval_filtration {
		interval_filtration() {}
		interval_filtration(int dim, Filtration_value b, Filtration_value d) : dim_(dim), b_(b), d_(d) {}
		/** Returns the absolute length of the interval \f$|d-b|\f$. */
		Filtration_value length() {
			if(b_ == d_) { return 0; } //otherwise inf - inf would return nan.
			return std::abs(b_ - d_);
		}
		/** Returns the absolute length of the log values of birth and death, i.e.  \f$|\log d - \log b|\f$.. */
		Filtration_value log_length() {//return the log-length
			if(b_ == d_) { return 0; } //otherwise inf - inf would return nan.
			return std::abs(log2((double)b_) - log2((double)d_));
		}
		/** Returns the dimension of the homological feature corresponding to the
		 * interval. */
		int dim() const { return dim_; }//return the homological dimension of the interval
		/** Returns the birth of the interval.*/
		Filtration_value birth() const { return b_; }//return the birth value
		/** Returns the death of the interval.*/
		Filtration_value death() const { return d_; }//return the death value
		/** Swaps the values of birth and death.*/
		void swap_birth_death() { std::swap(b_,d_); }

	private://note that we don't assume b_ <= d_
		int              dim_; //homological dimension
		Filtration_value b_; //filtration value associated to birth index
		Filtration_value d_; //filtration value associated to death index
	};

	/** \brief Structure to store persistence intervals by their index values.
	 *
	 * \details By convention, interval [b;d] are
	 * closed for finite indices b and d, and open for left-infinite and/or
	 * right-infinite endpoints.
	 */
	struct interval_index {
		interval_index() {}
		interval_index(int dim, Simplex_key b, Simplex_key d) : dim_(dim), b_(b), d_(d) {}
		/** Returns the dimension of the homological feature corresponding to the
		 * interval. */
		int dim() const { return dim_; }//return the homological dimension of the interval
		/** Returns the birth index of the interval.*/
		Filtration_value birth() const { return b_; }//return the birth value
		/** Returns the death index of the interval.*/
		Filtration_value death() const { return d_; }//return the death value

	private://note that we don't assume b_ <= d_
		int              dim_; //homological dimension
		Simplex_key b_; //filtration value associated to birth index
		Simplex_key d_; //filtration value associated to death index
	};

private:
	/* Comparison function to sort intervals by decreasing log-length in the
	 * output persistence diagram, i.e.,
	 * [f(b),f(d)]<[f(b'),f(d')] iff |log2(f(b))-log2(f(d))|> |log2(f(b'))-log2(f(d'))|
	 */
	struct cmp_intervals_by_log_length {
		cmp_intervals_by_log_length(){}
		bool operator()( interval_filtration p, interval_filtration q)
		{
			if(p.dim() != q.dim()) {return p.dim() < q.dim();}//lower dimension first
			if(p.log_length() != q.log_length()) {return p.log_length() > q.log_length();}
			if(p.birth() != q.birth()) {return p.birth() < q.birth();}//lex order
			return p.death() < q.death();
		}
	};
	/* Comparison function to sort intervals by decreasing length in the
	 * output persistence diagram, i.e.,
	 * [f(b),f(d)]<[f(b'),f(d')] iff  |f(b)-f(d)| > |f(b')-f(d')|
	 */
	struct cmp_intervals_by_length {
		cmp_intervals_by_length(){}
		bool operator()( interval_filtration p, interval_filtration q)
		{
			if(p.length() != q.length()) { return p.length() > q.length(); }//longest 1st
			if(p.dim() != q.dim()) {return p.dim() < q.dim();}//lower dimension first
			if(p.birth() != q.birth()) {return p.birth() < q.birth();}//lex order
			return p.death() < q.death();
		}
	};

public:
	/** \brief Initialization of the Zigzag_persistence class.
	 *
	 * \param[in] cpx   A model of ZigzagFilteredComplex.
	 * */
	Zigzag_persistence(int ignore_cycles_above_dim = -1)
		: cpx_()
		, dim_max_(ignore_cycles_above_dim)
		, lowidx_to_matidx_()
		, matrix_()
		, birth_ordering_()
		, persistence_diagram_()
		, num_arrow_(-1)
		, previous_filtration_value_(std::numeric_limits<Filtration_value>::infinity())
		, filtration_values_() {}

private:
	/* Set c1 <- c1 + c2, assuming canonical order of indices induced by the order in
	 * the vertical lists. self1 is the matrix_chain whose column is c1, for self
	 * reference of the new cells.
	 */
	void plus_equal_column(matrix_chain * self1, Column & c1, Column & c2)
	{
		//insert all elements of c2 in c1, in O(|c2| * log(|c1|+|c2|))
		if constexpr (Options::searchable_column) {
			for(auto &cell : c2) {
				auto it1 = c1.find(cell);
				if(it1 != c1.end()) {//already there => remove as 1+1=0
					Cell * tmp_ptr = &(*it1);
					it1->base_hook_matrix_row::unlink(); //unlink from row
					c1.erase(it1); //remove from col
					matrix_chain::cellPool_.destroy(tmp_ptr);
					// delete tmp_ptr;
				}
				else {//not there, insert new cell
					// Cell *new_cell = new Cell(cell.key(), self1);
					Cell *new_cell = matrix_chain::cellPool_.construct(cell.key(), self1);
					c1.insert(*new_cell);
					lowidx_to_matidx_[cell.key()]->row_.push_back(*new_cell);//row link,no order
				}
			}
		}
		else {//traverse both columns doing a standard column addition, in O(|c1|+|c2|)
			auto it1 = c1.begin();   auto it2 = c2.begin();
			while(it1 != c1.end() && it2 != c2.end())
			{
				if(it1->key() < it2->key()) { ++it1; }
				else {
					if(it1->key() > it2->key()) {
						// Cell * new_cell = new Cell(it2->key(), self1);
						Cell *new_cell = matrix_chain::cellPool_.construct(it2->key(), self1);
						c1.insert(it1, *new_cell); //col link, in order
						lowidx_to_matidx_[it2->key()]->row_.push_back(*new_cell);//row link,no order
						++it2;
					}
					else { //it1->key() == it2->key()
						auto tmp_it = it1;    ++it1; ++it2;
						Cell * tmp_ptr = &(*tmp_it);
						tmp_it->base_hook_matrix_row::unlink(); //unlink from row
						c1.erase(tmp_it); //remove from col
						matrix_chain::cellPool_.destroy(tmp_ptr);
						// delete tmp_ptr;
					}
				}
			}
			while(it2 != c2.end()) {//if it1 reached the end of its column, but not it2
				// Cell * new_cell = new Cell(it2->key(),self1);
				Cell *new_cell = matrix_chain::cellPool_.construct(it2->key(), self1);
				lowidx_to_matidx_[it2->key()]->row_.push_back(*new_cell); //row links
				c1.push_back(*new_cell);
				++it2;
			}
		}
	}

	/** Maintains the birth ordering <=b. Contains an std::map of size the number of
	 * non-zero rows of the homology matrix, at any time during the computation of
	 * zigzag persistence.
	 *
	 * By construction, we maintain the map satisfying
	 * 'birth_to_pos_[i] < birth_to_pos_[j]',
	 * with 0 <= i,j <= k indices in the quiver '0 \leftrightarrow ... \leftrightarrow i \leftrightarrow .. \leftrightarrow k'
	 * visited at time k of the algorithm (prefix of length k of the full zigzag
	 * filtration '0 \leftrightarrow ... \leftrightarrow i \leftrightarrow .. \leftrightarrow k \leftrightarrow ... \leftrightarrow n' that is studied),
	 * iff i <b j for the birth ordering.
	 *
	 * By construction, when adding index k+1 to '0 \leftrightarrow ... \leftrightarrow i \leftrightarrow .. \leftrightarrow k \leftrightarrow k+1',
	 * we have:
	 * - if k -> k+1 forward, then j <b k+1 for all indices j < k+1, otherwise
	 * - if k <- k+1 backward, then k+1 <b j for all indices j < k+1.
	 */
	struct birth_ordering {
		//example quiver indices    empty_cpx -> 0 -> 1 -> 2 <- 3 <- 4 -> 5 <- 6 etc
		birth_ordering() : birth_to_pos_(), max_birth_pos_(0), min_birth_pos_(-1) {}

		//when the arrow key-1 -> key is forward, key is larger than any other index
		//i < key in the birth ordering <b. We give key the largest value max_birth_pos_
		void add_birth_forward(Simplex_key key) { //amortized constant time
			birth_to_pos_.emplace_hint(birth_to_pos_.end(), key, max_birth_pos_);
			++max_birth_pos_;
		}
		//when the arrow key-1 <- key is backward, key is smaller than any other index
		//i < key in the birth ordering <b. We give key the smallest value min_birth_pos_
		void add_birth_backward(Simplex_key key) { //amortized constant time
			birth_to_pos_.emplace_hint(birth_to_pos_.end(), key, min_birth_pos_);
			--min_birth_pos_;
		}
		//when the row at index key is removed from the homology matrix, we do not need
		//to maintain its position in <b anymore
		void remove_birth(Simplex_key key) { birth_to_pos_.erase(key); }
		//increasing birth order <=b, true iff k1 <b k2
		bool birth_order(Simplex_key k1, Simplex_key k2) {
			return birth_to_pos_[k1] < birth_to_pos_[k2];
		}
		//decreasing birth order <=b, true iff k1 >b k2
		bool reverse_birth_order(Simplex_key k1, Simplex_key k2) {
			return birth_to_pos_[k1] > birth_to_pos_[k2];
		}

	private:
		//birth_to_pos_[i] < birth_to_pos_[j] iff i <b j
		std::map< Simplex_key, Simplex_key > birth_to_pos_;
		//by construction, max_birth_pos_ (resp. min_birth_pos_) is strictly larger
		//(resp. strictly smaller) than any value assigned to a key so far.
		Simplex_key                          max_birth_pos_;
		Simplex_key                          min_birth_pos_;
	};

public:
	/** \brief Computes the zigzag persistent homology of a zigzag filtered complex,
	 * using the reflection and transposition algorithm of \cite zigzag_reflection.
	 *
	 * \details After computation, the persistence diagram can be accessed via
	 * member method <CODE>persistence_diagram</CODE>, for the diagram with filtration
	 * values, or member method <CODE>index_persistence_diagram</CODE>, for the
	 * diagram with
	 * indices of paired simplices.
	 *
	 *
	 * <CODE>matrix_</CODE>, originally empty, maintains the set of chains, with a
	 * partition \f$ F \sqcup G \sqcup H\f$
	 * representing a compatible homology basis as in \cite zigzag_reflection.
	 *
	 * Each simplex in the complex stores a key field that stores the index of
	 * its insertion in the zigzag filtration.
	 *
	 * The algorithm maintains a compatible homology basis for the zigzag filtration.
	 *
	 * \f$$\emptyset = K_0 \leftrightarrow (...) \leftrightarrow K_i \leftarrow ... \leftarrow \emptyset\f$$
	 *
	 * where the prefix from \f$K_0\f$ to \f$K_i\f$ is equal to the i-th prefix of
	 * the input zigzag
	 * filtration given by <CODE>cpx_.filtration_simplex_range()</CODE>, and
	 * the suffix
	 * (from \f$K_i\f$
	 * to the right) is a sequence of simplex removals. Due to the structure of
	 * reflection diamonds, the removals are in reverse order of the insertions, to
	 * reduce the amount of transposition diamonds.
	 *
	 * Consequently, using <CODE>cpx_.key(zzsh)</CODE> as indexing for the matrix
	 * rows/cells,
	 * with the natural order on integers, makes our homology matrix <CODE>matrix_</CODE> upper
	 * triangular for the suffix \f$K_i \leftarrow ... \leftarrow 0\f$, seen as a
	 * standard persistence
	 * filtration. At \f$K_i\f$, the natural order on integers is also equivalent to the
	 * death-order \f$\leq_d\f$ (because all arrows in the suffix are backward).
	 *
	 * Insertion: <CODE>cpx_.key(*zzit)</CODE> is a strictly increasing sequence
	 * for <CODE>zzit</CODE>
	 * insertion of cells (does not need to be contiguous). However, for every forward
	 * arrow, we have <CODE>cpx_.key(*zzit) == num_arrows_</CODE>.
	 * Removal: <CODE>cpx_.key(*zzit)</CODE> gives the assigned key (during past
	 * insertion) of a
	 * <CODE>cell == *zzit</CODE> during a removal. We use <CODE>num_arrows_</CODE>
	 * to record the deaths in the
	 * persistence diagram.
	 * Insertion and Removal: <CODE>zzit.filtration()</CODE> is totally monotone.
	 * Note that the
	 * iterator encodes the filtration, and not the cells within the complex structure.
	 */
//	void zigzag_persistent_homology()
//	{ //compute index persistence, interval are closed, i.e., [b,d) is stored as
//		//[b,d-1]. The filtration values are maintained in field filtration_values_
//		Filtration_value prev_fil_, curr_fil_;

//		assert(num_arrow_ == 0);
//		auto zzrg = cpx_.filtration_simplex_range();
//		auto zzit = zzrg.begin();
//		dim_max_ = zzit.dim_max();

//		num_arrow_ = cpx_.key(*zzit);//should be 0

//		prev_fil_ = zzit.filtration();
//		filtration_values_.emplace_back(num_arrow_, prev_fil_);

//		while( zzit != zzrg.end() )
//		{ //insertion of a simplex
//			if(zzit.arrow_direction()) { num_arrow_ = cpx_.key(*zzit); }
//			else { ++num_arrow_; } //removal of a simplex, a simplex key corresponds to the index of its INSERTION
//			curr_fil_ = zzit.filtration();//cpx_.filtration(*zzit) is invalid for (<-);
//			if(curr_fil_ != prev_fil_) //check whether the filt value has changed
//			{ //consecutive pairs (i,f), (j,f') mean simplices of index k in [i,j-1] have
//				prev_fil_ = curr_fil_;                                //filtration value f
//				filtration_values_.emplace_back(num_arrow_, prev_fil_);
//			}
//			if(zzit.arrow_direction()) { //forward arrow, only consider critical cells
//				forward_arrow(*zzit);
//			}
//			else { //backward arrow
//				backward_arrow(*zzit);
//			}
//			++zzit;
//		}

////		if(!matrix_.empty()) {
////			std::cout << "There remain " << matrix_.size() << " columns in the matrix.\n";
////		}
//	}

	template<class VertexRange = std::initializer_list<Vertex_handle>>
	void insert_simplex(const VertexRange& simplex, Filtration_value filtration_value)
	{
		if (dim_max_ != -1 && simplex.size() > static_cast<unsigned int>(dim_max_) + 1) return;

		++num_arrow_;

		if (filtration_value != previous_filtration_value_) //check whether the filt value has changed
		{ //consecutive pairs (i,f), (j,f') mean simplices of index k in [i,j-1] have
			previous_filtration_value_ = filtration_value;                                //filtration value f
			filtration_values_.emplace_back(num_arrow_, previous_filtration_value_);
		}

		std::pair<Simplex_handle, bool> res = cpx_.insert_simplex(simplex, filtration_value);
		GUDHI_CHECK(res.second, "Zigzag_persistence::insert_simplex - insertion of a simplex already in the complex");
		cpx_.assign_key(res.first, num_arrow_);
		forward_arrow(res.first);
	}

	template<class VertexRange = std::initializer_list<Vertex_handle>>
	void remove_simplex(const VertexRange& simplex, Filtration_value filtration_value)
	{
		if (dim_max_ != -1 && simplex.size() > static_cast<unsigned int>(dim_max_) + 1) return;

		++num_arrow_;

		Simplex_handle sh = cpx_.find(simplex);
		GUDHI_CHECK(sh != cpx_.null_simplex(), "Zigzag_persistence::remove_simplex - removal of a simplex not in the complex");

		if (filtration_value != previous_filtration_value_) //check whether the filt value has changed
		{ //consecutive pairs (i,f), (j,f') mean simplices of index k in [i,j-1] have
			previous_filtration_value_ = filtration_value;                                //filtration value f
			filtration_values_.emplace_back(num_arrow_, previous_filtration_value_);
		}

		backward_arrow(sh);
		cpx_.remove_maximal_simplex(sh);
	}

	template<class SimplexRange = std::initializer_list<std::initializer_list<Vertex_handle>>,
			 class FiltrationRange = std::initializer_list<Filtration_value>>
	void insert_simplices_contiguously(const SimplexRange& simplices, const FiltrationRange& filtration_values)
	{
		auto simplexIt = simplices.begin();
		auto filIt = filtration_values.begin();
		for (; simplexIt != simplices.end(); ++simplexIt, ++filIt) {
			insert_simplex(*simplexIt, *filIt);
		}
	}

	template<class SimplexRange = std::initializer_list<std::initializer_list<Vertex_handle>>,
			 class FiltrationRange = std::initializer_list<Filtration_value>>
	void remove_simplices_contiguously(const SimplexRange& simplices, const FiltrationRange& filtration_values)
	{
		auto simplexIt = simplices.begin();
		auto filIt = filtration_values.begin();
		for (; simplexIt != simplices.end(); ++simplexIt, ++filIt) {
			remove_simplex(*simplexIt, *filIt);
		}
	}

	template<class SimplexRangeIterators, class FiltrationRangeIterators>
	void insert_simplices_contiguously(SimplexRangeIterators simplex_range_start,
									   SimplexRangeIterators simplex_range_end,
									   FiltrationRangeIterators filtration_range_start)
	{
		for (; simplex_range_start != simplex_range_end; ++simplex_range_start, ++filtration_range_start) {
			insert_simplex(*simplex_range_start, *filtration_range_start);
		}
	}

	template<class SimplexRangeIterators, class FiltrationRangeIterators>
	void remove_simplices_contiguously(SimplexRangeIterators simplex_range_start,
									   SimplexRangeIterators simplex_range_end,
									   FiltrationRangeIterators filtration_range_start)
	{
		for (; simplex_range_start != simplex_range_end; ++simplex_range_start, ++filtration_range_start) {
			remove_simplex(*simplex_range_start, *filtration_range_start);
		}
	}

	void print_current_complex(){
		for (auto& sh : cpx_.complex_simplex_range()){
			for (auto v : cpx_.simplex_vertex_range(sh)){
				std::cout << v << " ";
			}
			std::cout << " - " << cpx_.filtration(sh) << "\n";
		}
	}

private: 
	/** \brief Computes the boundary cycle of the new simplex zzsh, and express it as a
	 * sum of cycles. If all cycles are boundary cycles, i.e., columns with G-index
	 * in the matrix, then [\partial zzsh] = 0 and we apply an injective diamond to
	 * the zigzag module. Otherwise, we keep reducing with boundary- and live- cycles,
	 * i.e., columns with (F \cup G)-indices, and then apply a surjective diamond to
	 * the zigzag module.
	 */
	void forward_arrow( Simplex_handle zzsh )
	{ //maintain the <=b order
		birth_ordering_.add_birth_forward(num_arrow_);

		//Reduce the boundary of zzsh in the basis of cycles.
		//Compute the simplex keys of the simplices of the boundary of zzsh.
		std::set< Simplex_key > col_bsh; //set maintains the natural order on indices
		for( auto b_sh : cpx_.boundary_simplex_range(zzsh) )
		{ col_bsh.insert(cpx_.key(b_sh)); }

		//If empty boundary (e.g., when zzsh is a vertex in a simplicial complex)
		//Add a non-trivial cycle [c = zzsh] to the matrix, lowidx_to_matidx_make it a creator in F.
		if(col_bsh.empty()) // -> creator
		{ //New row and column with a bottom-right non-zero element, at index key(zzsh)
			//i.e., create a new cycle in F, equal to    *zzsh alone.
			matrix_.emplace_front(num_arrow_);
			auto new_chain_it = matrix_.begin();//the new chain
			//Update the map [index idx -> chain with lowest index idx] in matrix_
			lowidx_to_matidx_[num_arrow_] = new_chain_it;
			return;
		}

		// col_bsh.rbegin()) is idx of lowest element in col_bsh, because it is a set.
		matrix_chain *col_low = &(*lowidx_to_matidx_[*(col_bsh.rbegin())]);
		auto paired_idx = col_low->paired_col_; //col with which col_low is paired
		std::vector< matrix_chain * > chains_in_H; //for corresponding indices in H
		std::vector< matrix_chain * > chains_in_G;

		//Reduce col_bsh with boundary cycles, i.e., indices in G.
		std::pair< typename std::set< Simplex_key >::iterator, bool > res_insert;
		while( paired_idx != nullptr )
		{
			chains_in_H.push_back(paired_idx);//keep the col_h with which col_g is paired
			chains_in_G.push_back(col_low);   //keep the col_g
			for(auto &cell : (col_low->column())) { //Reduce with the column col_g
				res_insert = col_bsh.insert(cell.key());
				if( !res_insert.second ) { col_bsh.erase(res_insert.first); } //1+1 = 0
				//o.w. insertion has succeeded.
			}
			//If col_bsh is entirely reduced, \partial zzsh is a boundary cycle.
			if(col_bsh.empty()) {
				// if(cpx_.dimension(zzsh) >= max_dim_) {return;} we need max_dim creators
				injective_reflection_diamond(zzsh, chains_in_H);
				return;
			}
			//Continue the reduction
			col_low     =  &(*lowidx_to_matidx_[*(col_bsh.rbegin())]);//curr low index col
			paired_idx  =  col_low->paired_col_;//col with which col_low is paired
		}

		//Continue reducing with boundary and 'live' cycles, i.e., indices in G U F.
		std::vector< matrix_chain * > chains_in_F;
		while(true)
		{
			if(paired_idx == nullptr) { chains_in_F.push_back(col_low); }//col_low is in F
			else { chains_in_H.push_back(paired_idx); } //col_low in G, paired_idx is in H
			//Reduce with the column col_g or col_f
			for(auto &cell : (col_low->column())) {
				res_insert = col_bsh.insert(cell.key());
				if( !res_insert.second ) { col_bsh.erase(res_insert.first); } //1+1 = 0
				//o.w. insertion has succeeded.
			}
			//If col_bsh is entirely reduced, i.e. col_bsh == \emptyset.
			if(col_bsh.empty())
			{
				surjective_reflection_diamond(zzsh, chains_in_F, chains_in_H);
				return;
			}
			//Else, keep reducing.
			col_low = &(*lowidx_to_matidx_[*(col_bsh.rbegin())]); //curr low index col
			paired_idx = col_low->paired_col_;//col with which col_low is paired
		}
	}

	/** \brief Computes an injective diamond in the zigzag module, by inserting a new
	 * column for the chain zzsh - \sum col_h, for all col_h in chains_in_H, and a
	 * new row for the simplex zzsh.
	 */
	void injective_reflection_diamond ( Simplex_handle zzsh
										, std::vector< matrix_chain * > & chains_in_H )
	{ //Compute the chain   zzsh + \sum col_h, for col_h \in chains_in_H
		std::set< Simplex_key > col_bsh;
		std::pair< typename std::set< Simplex_key >::iterator, bool > res_insert;
		//produce the sum of all col_h in chains_in_H
		for( matrix_chain *idx_h : chains_in_H ) {
			for(auto &cell : (idx_h->column()) ) {
				res_insert = col_bsh.insert(cell.key());
				if( !res_insert.second ) { col_bsh.erase(res_insert.first); }
			}
		}
		//create a new cycle (in F) sigma - \sum col_h
		matrix_.emplace_front(num_arrow_, col_bsh.begin(), col_bsh.end(),
							  lowidx_to_matidx_);
		//Update the map 'index idx -> chain with lowest index idx' in matrix_
		auto chain_it = matrix_.begin();
		lowidx_to_matidx_[num_arrow_] = chain_it;
	}

	/** The vector chains_in_F is sorted by decreasing lowest index values in the
	 * columns corresponding to the chains, due to its computation in the reduction of
	 * \partial zzsh in forward_arrow(...). It is equivalent to decreasing death index
	 * order w.r.t. the <d ordering.
	 */
	void surjective_reflection_diamond( Simplex_handle zzsh
										, std::vector< matrix_chain * > & chains_in_F
										, std::vector< matrix_chain * > & chains_in_H )
	{ //fp is the largest death index for <=d
		//Set col_fp: col_fp <- col_f1+...+col_fp (now in G); preserves lowest idx
		auto chain_fp = *(chains_in_F.begin()); //col_fp, with largest death <d index.

		for(auto other_col_it = chains_in_F.begin()+1;
			other_col_it != chains_in_F.end(); ++other_col_it)
		{ plus_equal_column(chain_fp, chain_fp->column(), (*other_col_it)->column()); }
		//doesn't change the lowest idx as chain_fp has maximal lowest idx of all

		//chains_in_F is ordered, from .begin() to end(), by decreasing lowest_idx_. The
		//lowest_idx_ is also the death of the chain in the right suffix of the
		//filtration (all backward arrows). Consequently, the chains in F are ordered by
		//decreasing death for <d.
		//Pair the col_fi, i = 1 ... p-1, according to the reflection diamond principle
		//Order the fi by reverse birth ordering <=_b
		auto cmp_birth = [this](Simplex_key k1, Simplex_key k2)->bool
		{ return birth_ordering_.reverse_birth_order(k1,k2); };//true iff b(k1) >b b(k2)

		//available_birth: for all i by >d value of the d_i,
		//contains at step i all b_j, j > i, and maybe b_i if not stolen
		std::set< Simplex_key, decltype(cmp_birth) > available_birth(cmp_birth);
		//for f1 to f_{p} (i by <=d), insertion in available_birth_to_fidx sorts by >=b
		for(auto &chain_f : chains_in_F) { available_birth.insert(chain_f->birth()); }

		auto maxb_it = available_birth.begin();//max birth cycle
		auto maxb = *maxb_it; //max birth value, for persistence diagram
		available_birth.erase(maxb_it); //remove max birth cycle (stolen)

		auto last_modified_chain_it = chains_in_F.rbegin();

		//consider all death indices by increasing <d order i.e., increasing lowest_idx_
		for(auto chain_f_it  = chains_in_F.rbegin(); //by increasing death order <d
			*chain_f_it != chain_fp; ++chain_f_it )//chain_fp=*begin() has max death
		{ //find which reduced col has this birth
			auto birth_it = available_birth.find((*chain_f_it)->birth());
			if(birth_it == available_birth.end()) //birth is not available. *chain_f_it
			{ //must become the sum of all chains in F with smaller death index.
				//this gives as birth the maximal birth of all chains with strictly larger
				//death <=> the maximal availabe death.
				//Let c_1 ... c_f be the chains s.t. <[c_1+...+c_f]> is the kernel and
				// death(c_i) >d death(c_i-1). If the birth of c_i is not available, we set
				//c_i <- c_i + c_i-1 + ... + c_1, which is [c_i + c_i-1 + ... + c_1] on
				//the right (of death the maximal<d death(c_i)), and is [c_i + c_i-1 + ... +
				//c_1] + kernel = [c_f + c_f-1 + ... + c_i+1] on the left (of birth the max<b
				//of the birth of the c_j, j>i  <=> the max<b available birth).
				//N.B. some of the c_k, k<i, ahve already been modified to be equal to
				//c_k + c_k-1 + ... + c_1. The largest k with this property is maintained in
				//last_modified_chain_it (no need to compute from scratch the full sum).

				//last_modified is equal to c_k+...+c_1, all c_j, i>j>k, are indeed c_j
				//set c_i <- c_i + (c_i-1) + ... + (c_k+1) + (c_k + ... + c_1)
				for(auto chain_passed_it = last_modified_chain_it;//all with smaller <d death
					chain_passed_it != chain_f_it;  ++chain_passed_it)
				{
					plus_equal_column( (*chain_f_it), (*chain_f_it)->column()
									   , (*chain_passed_it)->column() );
				}
				last_modified_chain_it = chain_f_it;//new cumulated c_i+...+c_1
				//remove the max available death
				auto max_avail_b_it = available_birth.begin();//max because order by deacr <b
				Simplex_key max_avail_b = *max_avail_b_it;//max available birth

				(*chain_f_it)->assign_birth(max_avail_b); //give new birth
				available_birth.erase(max_avail_b_it); //remove birth from availability
			}
			else { available_birth.erase(birth_it); } //birth not available anymore, do not
		}                                          //modify *chain_f_it.
		//Compute the new column zzsh + \sum col_h, for col_h in chains_in_H
		std::set< Simplex_key > col_bsh;
		std::pair< typename std::set< Simplex_key >::iterator, bool > res_insert;
		for(auto other_col : chains_in_H)
		{ //Compute (\sum col_h) in a set
			for(auto &cell : (other_col->column()))
			{
				res_insert = col_bsh.insert(cell.key());
				if( !res_insert.second ) { col_bsh.erase(res_insert.first); } //1+1=0
			}
		}
		//Create and insert (\sum col_h) + sigma (in H, paired with chain_fp) in matrix_
		matrix_.emplace_front(cpx_.key(zzsh), chain_fp, col_bsh.begin(), col_bsh.end(), lowidx_to_matidx_);
		//record that the chain with lowest index key(zzsh) is the one just created
		auto chain_it = matrix_.begin();
		lowidx_to_matidx_[cpx_.key(zzsh)] = chain_it;//new row

		chain_fp->assign_paired_chain( &(*chain_it) );//pair chain_fp with the new chain
		chain_fp->assign_birth(-1); //now belongs to G now -> right interval [m-1,g]

		//Update persistence diagram with left interval [fil(b_max) ; fil(m))
		persistence_diagram_.emplace_back( cpx_.dimension(zzsh)-1
										   , maxb
										   , cpx_.key(zzsh));//-1);//
	}

	//cpx_.key(zzsh) is the key of the simplex we remove, not a new one
	void backward_arrow( Simplex_handle zzsh )
	{
		//maintain the <=b order
		birth_ordering_.add_birth_backward(num_arrow_);
		//column whose key is the one of the removed simplex
		auto curr_col_it = lowidx_to_matidx_.find(cpx_.key(zzsh));
		//corresponding chain
		matrix_chain * curr_col    = &(*(curr_col_it->second));
		//Record all columns that get affected by the transpositions, i.e., have a coeff
		std::vector< matrix_chain * > modified_columns;//in the row of idx key(zzsh)
		for(auto & hcell : (curr_col->row_)) {
			modified_columns.push_back(hcell.self_chain_);
		}
		//Sort by left-to-right order in the matrix_ (no order maintained in rows)
		std::stable_sort( modified_columns.begin(),modified_columns.end()
						  , [](matrix_chain *mc1, matrix_chain *mc2)
		{ return mc1->lowest_idx_ < mc2->lowest_idx_;} );

		//Modifies the pointer curr_col, not the other one.
		for(auto other_col_it = modified_columns.begin()+1;
			other_col_it != modified_columns.end(); ++other_col_it) {
			curr_col = arrow_transposition_case_study(curr_col, *other_col_it);
		}

		//curr_col points to the column to remove by restriction of K to K-{\sigma}
		if( curr_col->paired_col_ == nullptr ) { // in F
			int dim_zzsh = cpx_.dimension(zzsh);
			if(dim_max_ == -1 || (dim_max_ != -1 && dim_zzsh < dim_max_)) { //don't record intervals of max dim
				persistence_diagram_.emplace_back( dim_zzsh
												   , curr_col->birth()
												   , num_arrow_);// -1);
			}
		}
		else { //in H    -> paired with c_g, that now belongs to F now
			curr_col->paired_col_->assign_paired_chain(nullptr);
			curr_col->paired_col_->assign_birth(num_arrow_); //closed interval
		}

		//cannot be in G as the removed simplex is maximal
		matrix_.erase(curr_col_it->second);
		lowidx_to_matidx_.erase(curr_col_it);
	}

	/* Exchanges members of matrix_chains, except the column_ pointer. Modify
	 * also the lowidx_to_matidx_ data structure, considering that the matrix chains
	 * also exchange their lowest_idx_. Specifically, it is called by
	 * arrow_transposition_case_study when:
	 * c_s has originally birth b_s and low idx s, and c_t has birth b_t and low idx t
	 * however, c_s becomes c_s <- c_s+c_t with b_t <b b_s, and has SAME birth b_s but
	 * NEW lowest index t. c_t remains the same but, because the positions of s and t
	 * are transposed in the filtration, c_t keeps its birth b_t but has NEW lowest
	 * idx s.
	 *
	 * Note that Cells in the matrix store a pointer to their matrix_chain_.
	 * Consequently, exchanging columns would require to update all such pointers. That
	 * is why we avoid doing it, and prefer exchanging all other attributes.
	 */
	void exchange_lowest_indices_chains( matrix_chain * curr_col//c_s+c_t (prev. c_s)
										 , matrix_chain * other_col )//c_t
	{ //lowidx_to_matidx_[i]==matrix_chain* whose lowest index is i
		auto it_s = lowidx_to_matidx_.find(curr_col->lowest_idx_);
		auto it_t = lowidx_to_matidx_.find(other_col->lowest_idx_);

		std::swap(it_s->second, it_t->second);//swap matrix_chain* in lowidx_to_matidx_
		std::swap(curr_col->row_, other_col->row_);//swap associated row of lowest idx
		std::swap(curr_col->lowest_idx_, other_col->lowest_idx_);//swap lowest idx.
	}

	/**
	 * Permutes s and t, s goes up, whose insertions are adjacent, i.e., the following
	 * transformation (from (a) to (b)) in the filtration:
	 * from (a) ... \leftrightarrow K \leftarrow ... \leftarrow K' U {s,t} \leftarrow K' U {s} \leftarrow K' \leftarrow ...,
	 * where K' is at matrix index i, K' U {s} at matrix index i+1 and K' U {s,t} at
	 * matrix index i+2,
	 *
	 * to (b)   ... \leftrightarrow K \leftarrow ... \leftarrow K' U {s,t} \leftarrow K' U {t} \leftarrow K' \leftarrow ...,
	 *
	 * and the chain c_t has a non-trivial coefficient for s, i.e.,
	 * the bloc matrix gives (and becomes):
	 *                                 c_t                  c_t
	 *                                  +                    +
	 *      c_s c_t                    c_s c_s              c_s c_t
	 *   s   1   1                  t   1   0            t   1   1
	 *   t   0   1      --> either  s   0   1      or    s   0   1
	 *
	 * By construction, s is a simplex that we want to remove in the complex K. It is
	 * consequently maximal in K, and all complexes between K and K' U {s} in filtration
	 * (a).
	 *
	 * If c_s and c_t are both cycles (in F)that, before the permutation, are carried by
	 * respectively the closed intervals [b_s, i+1] and [b_t, i+2], then the sum
	 * c_s + c_t is a cycle carried by the interval
	 *                                          [max<b {b_s,b_t} , max<d {i,i+1} = i+1].
	 *
	 * If c_s and c_t are both chains in H, paired with c_gs and c_gt (in G) resp.,
	 * where c_gs is carried by [i;gs] and c_gt is carried by [i+1;gt] before
	 * permutation, then c_gt+c_gs is carried by
	 * [max<b {i+1,i}=i+1 ; max<d {gs,gt} = max {gs,gt}] because, all arrows being
	 * backward the max<b birth is the leftmost, and the max<d death is the leftmost.
	 *
	 * We have \partial(c_s+c_t) = c_gs+c_gt. Because both c_s and c_t contain s in their
	 * sum, c_s+c_t (in Z/2Z) contains t but not s in its sum, and all other simplices
	 * are in K. Consequently, the chain exists in any complex containing K U {t} as
	 * subcomplex.
	 *
	 * Note that because all arrows are backward on the right side of the quiver
	 * ... \leftrightarrow K \leftarrow ... \leftarrow K U {s,t} \leftarrow K U {t} \leftarrow K \leftarrow ..., we always have
	 * i+1 >d i (i+1 \leftarrow i backward).
	 * If j \leftarrow  ...   \leftarrow k are both birth indices on the right part of the quiver (all
	 * backward arrows) then systematically k <b j.
	 * If b \leftrightarrow ... \leftrightarrow K \leftarrow ... \leftarrow j \leftarrow ... are both birth indices, with j on the right
	 * part of the quiver (all backward arrows), and b in the prefix, then
	 * systematically j <b b.
	 *
	 * Because s is maximal in K, none of c_s or c_t, that have a non-trivial
	 * coefficient for s, can belong to G because cycles in G are the boundary of some
	 * chain in H.
	 *
	 * We get the following cases:
	 * c_s | c_t in:
	 *  F  |  F   keep (c_s+c_t), c_x: such that c_x has the min birth of c_s and c_t.
	 *                                 Both chains are cycles in F.
	 *  F  |  H   keep (c_s+c_t), c_s: c_s+c_t still in H (boundary is unchanged) and
	 *                                 c_s still in F
	 *  H  |  F   keep (c_s+c_t), c_t: c_s+c_t still in H (boundary is unchanged) and
	 *                                 c_t still in F
	 *  H  |  H   keep (c_s+c_t), c_x: let c_s be paired with c_gs in G and c_t be
	 *            paired with c_gt in G. Then (c_s+c_t) is paired with c_gs+c_gt in G,
	 *            of death index max<d d_gs,d_gt. With only backward arrow, the maximal
	 *            death for <d is the leftmost (i.e., == to the highest key of c_gs or
	 *            c_gt in the matrix), and the maximal birth for <b is the leftmost.
	 *            Because c_s+c_t exists in any complex containing K U {t}, c_gs+c_gt is
	 *            killed by the insertion of t before and after permutation of arrows.
	 *            In particular, [c_gs+c_gt] has the same lifespan as the cycle c_gs or
	 *            c_gt with the maximal death for <d.
	 *
	 *            c_x is the chain c_s or c_t that is paired with the chain in G with
	 *            smaller death index in <d. We also update c_gy <- c_gs+c_gt, where c_gy
	 *            is the cycle c_gs or c_gt with larger death index for <d (now paired
	 *            with c_gs+c_gt in H).
	 *
	 * Returns the new value of curr_col we continue with, i.e., the chain whose lowest
	 * index is the one of s, the simplex we are percolating left in the filtration,
	 * after permutation of arrows.
	 */
	matrix_chain * arrow_transposition_case_study( matrix_chain * curr_col//c_s
												   , matrix_chain * other_col )//c_t
	{ //c_s has low idx s and c_t low idx t
		if(curr_col->inF())
		{//case F x *
			if(other_col->inH()) { //                                    case F x H
				plus_equal_column( other_col, other_col->column()//c_t <- c_s+c_t still in H
								   , curr_col->column() );//(birth -2) and low idx t
				return curr_col;
			}//end case F x H
			else //                                                      case F x F
			{ //in F x F:           c_s+c_t has max<=b birth between b_t and b_s:
				if(birth_ordering_.birth_order(curr_col->birth(), other_col->birth()))
				{ //max<=b is the birth of other_col i.e., b_s <b b_t
					plus_equal_column( other_col, other_col->column()//c_t <- c_s+c_t of birth
									   , curr_col->column() );//b_t and lowest idx t. (same)
					//c_s still has birth b_s (minimal) and lowest idx s
					return curr_col;//continue with c_s of smaller birth b_s and lowest idx s
				}//endif
				else
				{ //max<=b is the birth of curr_col, i.e., b_t <b b_s
					plus_equal_column( curr_col, curr_col->column()//c_s <- c_s+c_t of birth
									   , other_col->column() );//b_s and of NEW lowest idx t
					//now c_t has (same) birth b_t (minimal) but NEW lowest idx s, so
					//exchange lowest_idx, the rows, and update lowidx_to_matidx structure
					exchange_lowest_indices_chains(curr_col, other_col);
					return other_col;//continue with c_t of (smaller) birth b_t and low idx s
				}//end else
			}//end case F x F
		}//end case F x *
		else {//curr_col->inH() == true,                                case H x *
			if(other_col->inH()) {//                                      case H x H
				//Case H x H, c_s+c_t paired w/ c_gs+c_gt, of death
				//max<d { key(c_gs),key(c_gt) } == usual max { key(c_gs),key(c_gt) }:
				auto curr_p_col  = curr_col->paired_col_; //c_s paired with c_gs, death d_gs
				auto other_p_col = other_col->paired_col_;//c_t paired with c_gt, death d_gt
				if( curr_p_col->lowest_idx_ < other_p_col->lowest_idx_)//<=> d_gs <d d_gt
				{
					plus_equal_column( other_p_col, other_p_col->column()//c_gt <- c_gs+c_gt,
									   , curr_p_col->column() );//of death d_gt, low idx d_gt
					//(same because bigger), paired with c_s+c_t (now &c_t, updated below)
					plus_equal_column( other_col, other_col->column()//c_t <- c_t+c_s, still
									   , curr_col->column() );//in H, low idx t (same)
					return curr_col;//continue with c_s, paired with c_gs of min death d_gs
				}
				else
				{// d_gt <d d_gs
					plus_equal_column( curr_p_col, curr_p_col->column()//c_gs <- c_gs+c_gt,
									   , other_p_col->column() );//of death d_gs, low idx d_gs
					//(same because bigger), paired with c_s+c_t (now &c_s, updated below)
					plus_equal_column( curr_col, curr_col->column()//c_s <- c_s+c_t, of NEW
									   , other_col->column());//low idx t (still in H)
					//now c_s is still in H (birth -2) but has NEW lowest idx t, and c_t has
					//low idx s after transposition.
					//exchange lowest_idx, the rows, and update lowidx_to_matidx structure
					exchange_lowest_indices_chains(curr_col, other_col);
					return other_col; //continue with c_t, paired w. c_g' of min death g'
				}
			}//end case H x H
			else {//other_col->inF() == true,                             case H x F
				plus_equal_column( curr_col, curr_col->column() //c_s <- c_s+c_t still in H,
								   , other_col->column());       //(birth -2) and NEW low idx t
				//now c_t, still in F, has (same) birth b_t but NEW lowest idx s, so
				//exchange lowest_idx, the rows, and update lowidx_to_matidx structure
				exchange_lowest_indices_chains(curr_col, other_col);
				return other_col; //continue with c_t, still in F, of birth b_t and low idx s
			}
		}
	}


public:
	/** \brief Returns the index persistence diagram as an std::list of intervals.*/
	const std::list< interval_index > & index_persistence_diagram() const
	{ 
		return persistence_diagram_; 
	}

	/** \brief Returns the filtration values \f$[f(b),f(d)]\f$ (generally real-valued)
	 * associated to the indices \f$[b,d]\f$ (integer valued) of the insertion or
	 * deletion of a simplex in the zigzag filtration.
	 *
	 * \details Used to convert a persistent interval \f$[b,d]\f$, computed by the
	 * persistent homology algorithm, into its filtration valued version
	 * \f$[f(b),f(d)]\f$ used in the persistence barcode. The information
	 * index->filtration is stored in the member <CODE>filtration_values_</CODE> of
	 * the class <CODE>Zigzag_persistence</CODE>.
	 *
	 * @param[in] b_key, d_key The indices of birth and death of a persistent
	 *                         interval.
	 *
	 * @param[out] std::pair A pair of real values \f$(f(b),f(d))\f$.
	 */
	std::pair<Filtration_value,Filtration_value> index_to_filtration(
			Simplex_key b_key, Simplex_key d_key) {
		// filtration_values_ must be sorted by increasing keys.
		auto it_b = //lower_bound(x) returns leftmost y s.t. x <= y
				std::lower_bound( filtration_values_.begin(), filtration_values_.end()
								  , std::pair<Simplex_key, Filtration_value>(b_key
																			 , std::numeric_limits<Filtration_value>::infinity() )
								  , []( std::pair<Simplex_key, Filtration_value> p1
								  , std::pair<Simplex_key, Filtration_value> p2)
		{ return p1.first < p2.first; }
				);
		if(it_b == filtration_values_.end() || it_b->first > b_key) { --it_b; }
		//it points to the rightmost z such that z <= x

		auto it_d = //
				std::lower_bound( filtration_values_.begin(), filtration_values_.end()
								  , std::pair<Simplex_key, Filtration_value>(d_key
																			 , std::numeric_limits<Filtration_value>::infinity() )
								  , []( std::pair<Simplex_key, Filtration_value> p1
								  , std::pair<Simplex_key, Filtration_value> p2)
		{ return p1.first < p2.first; }
				);
		if(it_d == filtration_values_.end() || it_d->first > d_key) { --it_d; }

		return std::make_pair(it_b->second, it_d->second);
	}

	/** \brief Writes the persistence diagram in a file.
	 *
	 * \details The filtration values are given by the zigzag persistence iterator, that assigns
	 * to any insertion or deletion of a simplex a filtration value ; we say that an
	 * arrow has an index \f$i\f$, and a corresponding filtration value \f$f(i)\f$.
	 * Reading a zigzag filtration from left to right, indices are strictly
	 * monotonically increasing, and the associated filtration values are monotonous
	 * (not necessarily
	 * strictly, either increasing or decreasing).
	 *
	 * Consider two consecutive arrows (insertion or deletion):
	 *
	 * \f$$K_1 \leftrightarrow K_2 \leftrightarrow K_3\f$$
	 *
	 * with respectively indices \f$i\f$ (left) and \f$i+1\f$ (right), and associated
	 * filtration values \f$f(i)\f$ and \f$f(i+1)\f$ respectively.
	 *
	 * If, the arrow \f$K_2 \leftrightarrow K_3\f$ leads to the creation of a new
	 * homology feature in \f$K_3\f$, it creates an (indexed) persistent interval
	 * \f$[\f$i+1; \cdot\f$, and a corresponding (filtration) persistent interval
	 * \f$[f(i+1); \cdot]\f$ in the persistence diagram.
	 *
	 * If a homology feature in \f$K_2\f$ is destroyed by the arrow \f$K_2 \leftrightarrow K_3\f$, it closes an (indexed)
	 * interval \f$[\cdot ; i+1]\f$, and a corresponding (filtration) persistent
	 * interval \f$[\cdot ; f(i+1)]\f$ in the persistence diagram.
	 *
	 * For example, in an oscillating Rips zigzag filtration, if, in the following
	 * chunk of filtration:
	 *
	 * \f$R_{\eta \varepsilon_i}(P_i) \rightarrow \cdots \leftarrow R_{\eta \varepsilon_{i+1}}(P_{i+1}),\f$
	 *
	 * if anything is created by any of the arrows above, it leads to an interval
	 * \f$[\varepsilon_{i+1}; \cdot]\f$. If anything is destroyed by any of the arrows
	 * above, if leads to an interval \f$[\cdot;\varepsilon_i]\f$. Note that we may
	 * have \f$\varepsilon_i > \varepsilon_{i+1}\f$.
	 *
	 * The bars are ordered by decreasing length.
	 *
	 * @param[in] os   the output stream in which the diagram is written.
	 * @param[in] shortest_interval   all intervals of lenght smaller or equal to
	 *                                this value are ignore. Default is 0.
	 */
	void persistence_diagram( std::ostream& os
							  , Filtration_value shortest_interval = 0.) {

		std::stable_sort(filtration_values_.begin(), filtration_values_.end(),
						 []( std::pair< Simplex_key, Filtration_value > p1
						 , std::pair< Simplex_key, Filtration_value > p2 )
		{ return p1.first < p2.first; }
		);

		std::vector< interval_filtration > tmp_diag;
		tmp_diag.reserve(persistence_diagram_.size());
		for(auto bar : persistence_diagram_)
		{
			Filtration_value birth,death;
			std::tie(birth,death) = index_to_filtration(bar.birth(), bar.death());

			if( std::abs(birth - death) > shortest_interval ) {
				tmp_diag.emplace_back(bar.dim(), birth, death );
			}
		}
		// cmp_intervals_by_length cmp;
		std::stable_sort(tmp_diag.begin(), tmp_diag.end(), cmp_intervals_by_length());

		os << "# dim  birth  death  [length]\n";
		for(auto bar : tmp_diag) {
			if(bar.birth() > bar.death()) { bar.swap_birth_death(); }
			os << bar.dim() << " " << bar.birth() << " " << bar.death() <<
				  " - [" << bar.length() << "] \n";
		}
	}

	/** \brief Returns the persistence diagram as a vector of real-valued intervals. */
	std::vector< interval_filtration >
	persistence_diagram(Filtration_value shortest_interval = 0., bool include_infinit_bars = false)
	{
		std::stable_sort(filtration_values_.begin(), filtration_values_.end(),
						 []( std::pair< Simplex_key, Filtration_value > p1
						 , std::pair< Simplex_key, Filtration_value > p2 )
		{ return p1.first < p2.first; }
		);

		std::vector< interval_filtration > diag;
		diag.reserve(persistence_diagram_.size());
		for(auto bar : persistence_diagram_)
		{
			Filtration_value birth,death;
			std::tie(birth,death) = index_to_filtration(bar.birth(), bar.death());

			if( std::abs(birth - death) > shortest_interval ) {
				diag.emplace_back(bar.dim(), birth, death );
			}
		}
		//put lower value as birth
		for(auto &bar : diag) {
			if( bar.birth() > bar.death() ) { bar.swap_birth_death(); }
		}
		std::stable_sort(diag.begin(), diag.end(), cmp_intervals_by_length());

		auto birth =
			[this](Simplex_key b_key) {
				auto it_b =	 // lower_bound(x) returns leftmost y s.t. x <= y
					std::lower_bound(filtration_values_.begin(), filtration_values_.end(),
									 std::pair<Simplex_key, Filtration_value>(
										 b_key, std::numeric_limits<Filtration_value>::infinity()),
									 [](std::pair<Simplex_key, Filtration_value> p1,
										std::pair<Simplex_key, Filtration_value> p2) { return p1.first < p2.first; });
				if (it_b == filtration_values_.end() || it_b->first > b_key) {
					--it_b;
				}
				return it_b->second;
			};

		//TODO: dimension value
		if (include_infinit_bars) {
			for (const matrix_chain &col : matrix_) {
				if (col.inF())
					diag.emplace_back(-1, birth(col.birth()), std::numeric_limits<Filtration_value>::infinity());
			}
		}

		return diag;
	}

private:
	Complex                                            cpx_; // complex
	int                                                 dim_max_;//max dim complex
	//idx -> chain with lowest element at index idx in matrix_
	std::map< Simplex_key, typename std::list<matrix_chain>::iterator >
	lowidx_to_matidx_;
	//arbitrary order for the matrix chains
	std::list< matrix_chain >                            matrix_; // 0 ... m-1
	// birth_vector                                           birth_vector_; //<=b order
	birth_ordering                                      birth_ordering_;
	std::list< interval_index >                             persistence_diagram_;
	Simplex_key                                         num_arrow_; //current index
	Filtration_value                                         previous_filtration_value_;
	// filtration_values stores consecutive pairs (i,f) , (j,f') with f != f',
	// meaning that all inserted simplices with key in [i;j-1] have filtration value f
	//i is the smallest simplex index whose simplex has filtration value f.
	std::vector< std::pair< Simplex_key, Filtration_value > > filtration_values_;
};//end class Zigzag_persistence


/** ZigzagPersistenceOptions, represents matrix columns by intrusive lists.*/
struct Zigzag_persistence_collist {
	static const bool searchable_column = false;
};
/** ZigzagPersistenceOptions, represents matrix columns by intrusive sets.*/
struct Zigzag_persistence_colset {
	static const bool searchable_column = true;
};

} //namespace zigzag_persistence

} //namespace Gudhi

#endif //ZIGZAG_PERSISTENCE_H_

