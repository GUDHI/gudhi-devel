/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siddharth Pritam
 *
 *    Copyright (C) 2019 Inria
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef FLAG_COMPLEX_SPARSE_MATRIX_H_
#define FLAG_COMPLEX_SPARSE_MATRIX_H_

#include <iostream>
#include <utility>
#include <vector>
#include <queue>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <unordered_set>
#include <ctime>
#include <fstream>

#include <Eigen/Sparse>

namespace Gudhi {

namespace strong_collapse {

using Vertex = int;
using Edge = std::pair<Vertex, Vertex>;
using Edge_list = std::vector<Edge>;

using Reduction_map = std::unordered_map<Vertex, Vertex>;

using Sparse_row_matrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using Sparse_row_iterator = Sparse_row_matrix::InnerIterator;

using Filtered_sorted_edge_list = std::vector<std::tuple<double, Vertex, Vertex> >;

//!  Class Flag_complex_sparse_matrix
/*!
  The class for storing the Vertices v/s MaxSimplices Sparse Matrix and performing collapses operations using the N^2()
  Algorithm.
*/
class Flag_complex_sparse_matrix {
 private:
  /** \brief Stores the vertices (or unordered set of vertex) of the original simplicial complex. */
  std::unordered_set<Vertex> vertices_;

  /** \brief Stores the 1-simplices (or list of edges) of the original simplicial complex. The list is updated in
   * `after_collapse()` method. */
  Edge_list one_simplices_;

  /** \brief The row in the sparse matrix does not correspond to the vertex number. This map helps to switch from row
   * to vertex 
   * */
  std::unordered_map<int, Vertex> row_to_vertex_;

  /** \brief The row in the sparse matrix does not correspond to the vertex number. This map helps to switch from
   * vertex to row
   * */
  std::unordered_map<Vertex, int> vertex_to_row_;

  //! Stores the number of vertices in the original Simplicial Complex.
  /*!
    This stores the count of vertices (which is also the number of rows in the Matrix).
  */
  std::size_t rows_;

  /** \brief Stores the collapsed sparse matrix representation. Initialized by the `after_collapse()` method.
   * */
  Sparse_row_matrix sparse_collapsed_matrix_;
  /** \brief Stores the sparse matrix representation of the 1-simplices. Initialized by the constructor.
   * */
  Sparse_row_matrix sparse_matrix_;

  /** \brief Stores <I>true</I> for dominated rows and  <I>false</I> for undominated rows.
   * Initialised to a vector of length equal to the value of the variable <B>rows</B> with all <I>false</I> values.
   * Subsequent removal of dominated vertices is reflected by concerned entries changing to <I>true</I> in this vector.
  */
  std::vector<bool> domination_indicator_;

  /** \brief Stores the indices of the rows to be checked for domination in the current iteration.
   * Initialised to a queue with all row-indices inserted.
   * Subsequently once the row is checked for dominated the row-index is poped out from the queue. A row-index is
   * inserted once again if it is a non-zero element of a dominated column.
   */
  std::vector<bool> contraction_indicator_;

  /** \brief Stores the indices of the row to-be checked for domination in the current iteration.
   * Initialised to an empty queue.
   * Subsequently once a dominated column is found, its non-zero row indices are inserted.
   */
  std::queue<std::size_t> row_iterator_;

  /** \brief Stores <I>true</I> if the current row is inserted in the queue <B>row_iterator_<B> otherwise its value is
   * <I>false<I>.
   *
   * Initialised to a boolean vector of length equal to the value of the variable <B>rows</B> with all <I>true</I>
   * values. Subsequent removal/addition of a row from <B>row_iterator_<B> is reflected by concerned entries changing
   * to <I>false</I>/<I>true</I> in this vector.
   */
  std::vector<bool> row_insert_indicator_;

  /** \brief Stores the Reduction / Collapse of vertices.
   * This is empty to begin with. As and when collapses are done (let's say from dominated vertex <I>v</I> to
   * dominating vertex <I>v'</I>) : <br> <B>reduction_map_</B>[<I>v</I>] = <I>v'</I> is entered into the map. <br>
   * This does not store uncollapsed vertices. What it means is that say vertex <I>x</I> was never collapsed onto
   * any other vertex.
   * Then, this map <B>WILL NOT</B> have any entry like <I>x</I> -> <I>x</I>. Basically, it will have no entry
   * corresponding to vertex <I>x</I> at all.
   */
  Reduction_map reduction_map_;

  int expansion_limit_;

  void init() {
    row_to_vertex_.clear();
    vertex_to_row_.clear();
    one_simplices_.clear();
    reduction_map_.clear();

    domination_indicator_.clear();
    row_insert_indicator_.clear();
    // VR: row_iterator_.push(0);
    // VR: row_iterator_.pop();

    rows_ = 0;

    expansion_limit_ = 2;
  }

  //!	Function for computing the Fake Simplex_tree corresponding to the core of the complex.
  /*!
    First calls strong_collapse(), and then computes the Fake Simplex_tree of the core using the Sparse matrix that we
    have. How does it compute the Fake simplex tree ? <br> Goes over all the columns (remaining MaximalSimplices) and
    for each of them, inserts that simplex <br>
    ['that simplex' means the maximal simplex with all the (remaining) vertices] with all subfaces using the <br>
    <I>insert_new_edges()</I> function from Gudhi's Fake_simplex_tree.
  */
  void after_collapse() {
    sparse_collapsed_matrix_ = Sparse_row_matrix(rows_, rows_);  // Just for debugging purpose.
    one_simplices_.clear();
    for (std::size_t rw = 0; rw < rows_; ++rw) {
      if (not domination_indicator_[rw])  // If the current column is not dominated
      {
        auto nbhrs_to_insert = read_row_index(rw);  // returns row indices of the non-dominated vertices.
        for (auto& v : nbhrs_to_insert) {
          sparse_collapsed_matrix_.insert(rw, v) = 1;
          if (rw < v) {
            one_simplices_.push_back({row_to_vertex_[rw], row_to_vertex_[v]});
          }
        }
      }
    }
    return;
  }
  //! Function to fully compact a particular vertex of the reduction_map_.
  /*!
    It takes as argument the iterator corresponding to a particular vertex pair (key-value) stored in the reduction_map_.
    <br> It then checks if the second element of this particular vertex pair is present as a first element of some other
    key-value pair in the map. If no, then the first element of the vertex pair in consideration is fully compact. If
    yes, then recursively call fully_compact_this_vertex() on the second element of the original pair in consideration
    and assign its resultant image as the image of the first element of the original pair in consideration as well.
  */
  void fully_compact_this_vertex(Reduction_map::iterator iter) {
    Reduction_map::iterator found = reduction_map_.find(iter->second);
    if (found == reduction_map_.end()) return;

    fully_compact_this_vertex(found);
    iter->second = reduction_map_[iter->second];
  }

  //! Function to fully compact the Reduction Map.
  /*!
    While doing strong collapses, we store only the immediate collapse of a vertex. Which means that in one round,
    vertex <I>x</I> may collapse to vertex <I>y</I>. And in some later round it may be possible that vertex <I>y</I>
    collapses to <I>z</I>. In which case our map stores : <br> <I>x</I> -> <I>y</I> and also <I>y</I> -> <I>z</I>. But
    it really should store : <I>x</I> -> <I>z</I> and <I>y</I> -> <I>z</I>. This function achieves the same. <br> It
    basically calls fully_compact_this_vertex() for each entry in the map.
  */
  void fully_compact() {
    Reduction_map::iterator it = reduction_map_.begin();
    while (it != reduction_map_.end()) {
      fully_compact_this_vertex(it);
      it++;
    }
  }

  void sparse_strong_collapse() {
    complete_domination_check();
    return;
  }

  // Complete check for rows in row_iterator_, row_insert_indicator_ is a list of boolean
  // indicator if a vertex is already inserted in the working row_queue (row_iterator_)
  void complete_domination_check() {
    // row_iterator_ contains list (FIFO) of rows to be considered for domination check
    while (not row_iterator_.empty())
    {
      double k = row_iterator_.front();
      row_iterator_.pop();
      row_insert_indicator_[k] = false;
      if (not domination_indicator_[k])  // Check if is  already dominated
      {
        std::vector<double> non_zero_inner_indices = read_row_index(k);
        for (auto index : non_zero_inner_indices) {
          // "true" for row domination comparison
          int check_domination = pair_domination_check(k, index);
          if (check_domination == 1)
          {
            // row k is dominated by index, k <= index
            set_zero(k, index);
            break;
          } else if (check_domination == -1) {
            // row index is dominated by k, index <= k;
            set_zero(index, k);
          }
        }
      }
    }
  }

  // True for row comparison, false for column comparison
  int pair_domination_check(double i, double j)
  {
    if (i != j) {
      std::vector<double> list_i = read_row_index(i);
      std::vector<double> list_j = read_row_index(j);
      if (list_j.size() <= list_i.size()) {
        if (std::includes(list_i.begin(), list_i.end(), list_j.begin(), list_j.end()))
          // list_j is a subset of list_i
          return -1;
      }

      else if (std::includes(list_j.begin(), list_j.end(), list_i.begin(), list_i.end()))
        // list_i is a subset of list_j
        return 1;
    }
    return 0;
  }

  // Returns list of non-zero columns of the particular indx.
  std::vector<double> read_row_index(double indx)
  {
    std::vector<double> non_zero_indices;
    if (not domination_indicator_[indx])
      // Iterate over the non-zero columns
      for (Sparse_row_iterator it(sparse_matrix_, indx); it; ++it) {
        if (not domination_indicator_[it.index()]) {
          // inner index, here it is equal to it.columns()
          non_zero_indices.push_back(it.index());
        }
      }
    return non_zero_indices;
  }

  void set_zero(double dominated, double dominating) {
    domination_indicator_[dominated] = true;
    reduction_map_[row_to_vertex_[dominated]] = row_to_vertex_[dominating];

    vertex_to_row_.erase(row_to_vertex_[dominated]);
    vertices_.erase(row_to_vertex_[dominated]);
    row_to_vertex_.erase(dominated);

    // Iterate over the non-zero rows
    for (Sparse_row_iterator it(sparse_matrix_, dominated); it; ++it)
      // Checking if the row is already dominated(set zero) or inserted
      if (not domination_indicator_[it.index()] &&
          not row_insert_indicator_[it.index()])
      {
        row_iterator_.push(it.index());
        row_insert_indicator_[it.index()] = true;
      }
  }

  // Returns list of non-zero "vertices" of the particular colIndx. the difference
  // is in the return type
  std::vector<Vertex> read_row(double row_index)
  {
    std::vector<Vertex> colmns;
    // Iterate over the non-zero columns
    for (Sparse_row_iterator non_zero_column_it(sparse_matrix_, row_index); non_zero_column_it; ++non_zero_column_it)
      // Check if the row corresponds to a dominated vertex
      if (not domination_indicator_[non_zero_column_it.index()])
        // inner index, here it is equal to it.col()
        colmns.push_back(row_to_vertex_[non_zero_column_it.index()]);
    std::sort(colmns.begin(), colmns.end());
    return colmns;
  }

  // Returns list of all non-zero "vertices" of the particular colIndx which
  // are currently active. the difference is in the return type.
  std::vector<Vertex> read_active_row(double row_index)
  {
    std::vector<Vertex> colmns;
    // Iterate over the non-zero columns
    for (Sparse_row_iterator non_zero_column_it(sparse_matrix_, row_index); non_zero_column_it; ++non_zero_column_it)
      // Check if the row corresponds to a contracted vertex
      if (not contraction_indicator_[non_zero_column_it.index()])
        // inner index, here it is equal to it.col()
        colmns.push_back(row_to_vertex_[non_zero_column_it.index()]);

    std::sort(colmns.begin(), colmns.end());
    return colmns;
  }

  // Returns list of all non-zero "vertices" of the particular colIndx whether
  // dominated or not. the difference is in the return type.
  std::vector<Vertex> read_all_row(double row_index)
  {
    std::vector<Vertex> colmns;
    // Iterate over the non-zero columns
    for (Sparse_row_iterator non_zero_column_it(sparse_matrix_, row_index); non_zero_column_it; ++non_zero_column_it)
      // inner index, here it is equal to it.row()
      colmns.push_back(row_to_vertex_[non_zero_column_it.index()]);
    std::sort(colmns.begin(), colmns.end());
    return colmns;
  }

  // swap the rows of v and w. Both should be members of the skeleton
  void swap_rows(const Vertex& v, const Vertex& w) {
    if (membership(v) && membership(w)) {
      auto rw_v = vertex_to_row_[v];
      auto rw_w = vertex_to_row_[w];
      vertex_to_row_[v] = rw_w;
      vertex_to_row_[w] = rw_v;
      row_to_vertex_[rw_v] = w;
      row_to_vertex_[rw_w] = v;
    }
  }

 public:
  //! Default Constructor
  /*!
    Only initialises all Data Members of the class to empty/Null values as appropriate.
    One <I>WILL</I> have to create the matrix using the Constructor that has an object of the Simplex_tree class as
    argument.
  */

  Flag_complex_sparse_matrix() { init(); }

  Flag_complex_sparse_matrix(std::size_t expRows) {
    init();
    // Initializing sparse_matrix_, This is a row-major sparse matrix.
    sparse_matrix_ = Sparse_row_matrix(
        expansion_limit_ * expRows,
        expansion_limit_ * expRows);
  }

  //! Main Constructor
  /*!
    Argument is an instance of Fake_simplex_tree. <br>
    This is THE function that initialises all data members to appropriate values. <br>
    <B>row_to_vertex_</B>, <B>vertex_to_row_</B>, <B>rows</B>, <B>cols</B>, <B>sparseMxSimplices</B> are initialised here.
    <B>domination_indicator_</B>, <B>row_insert_indicator_</B>
    ,<B>row_iterator_<B>,<B>simpDomnIndicator<B>,<B>colInsertIndicator<B> and <B>columnIterator<B> are initialised by
    init_lists() function which is called at the end of this. <br> What this does:
      1. Populate <B>row_to_vertex_</B> and <B>vertex_to_row_</B> by going over through the vertices of the Fake_simplex_tree
    and assign the variable <B>rows</B> = no. of vertices
      2. Initialise the variable <B>cols</B> to zero and allocate memory from the heap to <B>MxSimplices</B> by doing
    <br> <I>MxSimplices = new std::vector<bool>[rows];</I>
      3. Iterate over all maximal simplices of the Fake_simplex_tree and populates the column of the sparseMatrix.
      4. Initialise <B>active_rows</B> to an array of length equal to the value of the variable <B>rows</B> and all
    values assigned true. [All vertices are there in the simplex to begin with]
      5. Initialise <B>active_cols</B> to an array of length equal to the value of the variable <B>cols</B> and all
    values assigned true. [All maximal simplices are maximal to begin with]
      6. Calls the private function init_lists().
  */
  Flag_complex_sparse_matrix(const std::size_t num_vertices, const Filtered_sorted_edge_list& edge_t) {
    init();
    // Initializing sparse_matrix_, This is a row-major sparse matrix.
    sparse_matrix_ = Sparse_row_matrix(
        expansion_limit_ * num_vertices,
        expansion_limit_ * num_vertices);

    for (std::size_t bgn_idx = 0; bgn_idx < edge_t.size(); bgn_idx++) {
      insert_new_edges(std::get<1>(edge_t.at(bgn_idx)), std::get<2>(edge_t.at(bgn_idx)), 1);
    }
    sparse_matrix_.makeCompressed();
  }

  //!	Function for performing strong collapse.
  /*!
    calls sparse_strong_collapse(), and
    Then, it compacts the reduction_map_ by calling the function fully_compact().
  */
  void strong_collapse() {
    sparse_strong_collapse();
    // Now we complete the Reduction Map
    fully_compact();
    // Post processing...
    after_collapse();
    return;
  }

  bool membership(const Vertex& v) {
    auto rw = vertex_to_row_.find(v);
    if (rw != vertex_to_row_.end())
      return true;
    else
      return false;
  }

  bool membership(const Edge& e) {
    auto u = std::get<0>(e);
    auto v = std::get<1>(e);
    if (membership(u) && membership(v)) {
      auto rw_u = vertex_to_row_[u];
      auto rw_v = vertex_to_row_[v];
      if (rw_u <= rw_v)
        // Taking advantage of sorted lists.
        for (auto x : read_row_index(rw_v)) {
          if (rw_u == x)
            return true;
          else if (rw_u < x)
            return false;
        }
      else
        // Taking advantage of sorted lists.
        for (auto x : read_row_index(rw_u)) {
          if (rw_v == x)
            return true;
          else if (rw_v < x)
            return false;
        }
    }
    return false;
  }
  void insert_vertex(const Vertex& vertex, double filt_val) {
    auto rw = vertex_to_row_.find(vertex);
    if (rw == vertex_to_row_.end()) {
      // Initializing the diagonal element of the adjency matrix corresponding to rw_b.
      sparse_matrix_.insert(rows_, rows_) = filt_val;
      domination_indicator_.push_back(false);
      row_insert_indicator_.push_back(true);
      contraction_indicator_.push_back(false);
      row_iterator_.push(rows_);
      vertex_to_row_.insert(std::make_pair(vertex, rows_));
      row_to_vertex_.insert(std::make_pair(rows_, vertex));
      vertices_.emplace(vertex);
      rows_++;
    }
  }

  // The edge must not be added before, it should be a new edge.
  void insert_new_edges(const Vertex& u, const Vertex& v, double filt_val) {
    insert_vertex(u, filt_val);
    if (u != v) {
      insert_vertex(v, filt_val);

      auto rw_u = vertex_to_row_.find(u);
      auto rw_v = vertex_to_row_.find(v);

      sparse_matrix_.insert(rw_u->second, rw_v->second) = filt_val;
      sparse_matrix_.insert(rw_v->second, rw_u->second) = filt_val;
      one_simplices_.emplace_back(u, v);
    }
  }

  std::size_t num_vertices() const { return vertices_.size(); }

  //!	Function for returning the reduction_map_.
  /*!
    This is the (stl's unordered) map that stores all the collapses of vertices. <br>
    It is simply returned.
  */

  Reduction_map reduction_map() const { return reduction_map_; }

  std::unordered_set<Vertex> vertex_set() const { return vertices_; }
  
  Sparse_row_matrix collapsed_matrix() const { return sparse_collapsed_matrix_; }

  Sparse_row_matrix uncollapsed_matrix() const { return sparse_matrix_; }

  Edge_list all_edges() const { return one_simplices_; }

  void contraction(const Vertex& del, const Vertex& keep) {
    if (del != keep) {
      bool del_mem = membership(del);
      bool keep_mem = membership(keep);
      if (del_mem && keep_mem) {
        std::vector<double> del_indcs, keep_indcs, diff;
        auto row_del = vertex_to_row_[del];
        auto row_keep = vertex_to_row_[keep];
        del_indcs = read_row_index(row_del);
        keep_indcs = read_row_index(row_keep);
        std::set_difference(del_indcs.begin(), del_indcs.end(), keep_indcs.begin(), keep_indcs.end(),
                            std::inserter(diff, diff.begin()));
        for (auto& v : diff) {
          if (v != row_del) {
            sparse_matrix_.insert(row_keep, v) = 1;
            sparse_matrix_.insert(v, row_keep) = 1;
          }
        }
        vertex_to_row_.erase(del);
        vertices_.erase(del);
        row_to_vertex_.erase(row_del);
      } else if (del_mem && not keep_mem) {
        vertex_to_row_.insert(std::make_pair(keep, vertex_to_row_.find(del)->second));
        row_to_vertex_[vertex_to_row_.find(del)->second] = keep;
        vertices_.emplace(keep);
        vertices_.erase(del);
        vertex_to_row_.erase(del);

      } else {
        std::cerr << "The first vertex entered in the method contraction() doesn't exist in the skeleton." << std::endl;
        exit(-1);
      }
    }
  }

  void print_sparse_skeleton() { std::cout << sparse_matrix_ << std::endl; }

  // Returns the contracted edge. along with the contracted vertex in the begining of the list at {u,u} or {v,v}
  void active_strong_expansion(const Vertex& v, const Vertex& w, double filt_val) {
    if (membership(v) && membership(w) && v != w) {
      auto active_list_v_w = active_relative_neighbors(v, w);
      auto active_list_w_v = active_relative_neighbors(w, v);
      // simulate the contraction of w by expanding the star of v
      if (active_list_w_v.size() <
          active_list_v_w.size()) {
        // simulate the contraction of w by expanding the star of v
        for (auto& x : active_list_w_v) {
          active_edge_insertion(v, x, filt_val);
        }
        swap_rows(v, w);
      } else {
        for (auto& y : active_list_v_w) {
          active_edge_insertion(w, y, filt_val);
        }
      }
      auto rw_v = vertex_to_row_.find(v);
      contraction_indicator_[rw_v->second] = true;
    }
    if (membership(v) && !membership(w)) {
      relable(v, w);
    }
  }

 private:
  std::vector<Vertex> active_relative_neighbors(const Vertex& v, const Vertex& w) {
    std::vector<Vertex> diff;
    if (membership(v) && membership(w)) {
      auto nbhrs_v = active_neighbors(v);
      auto nbhrs_w = active_neighbors(w);
      std::set_difference(nbhrs_v.begin(), nbhrs_v.end(), nbhrs_w.begin(), nbhrs_w.end(),
                          std::inserter(diff, diff.begin()));
    }
    return diff;
  }

  std::vector<Vertex> active_neighbors(const Vertex& v) {
    std::vector<Vertex> nb;
    auto rw_v = vertex_to_row_.find(v);
    if (rw_v != vertex_to_row_.end()) nb = read_active_row(rw_v->second);

    return nb;
  }

  // relable v as w.
  void relable(const Vertex& v, const Vertex& w) {
    if (membership(v) and v != w) {
      auto rw_v = vertex_to_row_[v];
      row_to_vertex_[rw_v] = w;
      vertex_to_row_.insert(std::make_pair(w, rw_v));
      vertices_.emplace(w);
      vertex_to_row_.erase(v);
      vertices_.erase(v);
    }
  }

  void active_edge_insertion(const Vertex& v, const Vertex& w, double filt_val) {
    insert_new_edges(v, w, filt_val);
  }

};

}  // namespace strong_collapse

}  // namespace Gudhi

#endif  // FLAG_COMPLEX_SPARSE_MATRIX_H_
