/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Hannah Schreiber
 *
 *    Copyright (C) 2024 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#include <vector>

#include <gudhi/matrix.h>
#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Simplex_tree.h>

using Gudhi::persistence_matrix::Default_options;
using Gudhi::persistence_matrix::Column_types;
using Gudhi::Simplex_tree;

struct RU_options : Default_options<Column_types::INTRUSIVE_LIST, true> 
{
  static const bool has_column_pairings = true;
  static const bool can_retrieve_representative_cycles = true;
};

struct Chain_options : Default_options<Column_types::INTRUSIVE_LIST, true> 
{
  static const bool has_column_pairings = true;
  static const bool is_of_boundary_type = false;
};

using Base_matrix = Gudhi::persistence_matrix::Matrix<>;
using RU_matrix = Gudhi::persistence_matrix::Matrix<RU_options>;
using Chain_matrix = Gudhi::persistence_matrix::Matrix<Chain_options>;

void build_simplex_tree(Simplex_tree<>& st){
  std::vector<std::vector<int> > simplices = {
    {0},
    {1},
    {2},
    {3},
    {4},
    {5},
    {6},
    {2, 3},
    {4, 5},
    {0, 2},
    {0, 1},
    {1, 3},
    {1, 2},
    {1, 2, 3},
    {0, 1, 2},
    {5, 6},
    {2, 4},
    {4, 6},
    {4, 5, 6},
    {3, 6}
  };

  //insertion in the simplex tree
  for (unsigned int i = 0; i < simplices.size(); ++i){
    st.insert_simplex(simplices[i], i);
  }
}

template <class Column>
void print_column(Column& col, unsigned int num, unsigned int size) {
  std::cout << "col " << num << ":  ";
  if (num < 10) std::cout << " ";
  for (const auto& e : col.get_content(size)) {
    if (e == 0u)
      std::cout << "- ";
    else
      std::cout << e << " ";
  }
  std::cout << "\n";
}

void print_matrix(Base_matrix& bm, unsigned int size){
  std::cout << "Base_matrix:\n";
  for (unsigned int i = 0; i < size; ++i) {
    print_column(bm.get_column(i), i, size);
  }
  std::cout << "\n";
}

void print_matrix(const Chain_matrix& cm, unsigned int size){
  std::cout << "Chain_matrix:\n";
  // just note that if some vine swaps or removals occured, this would
  // not give us the columns in the order of filtration anymore, but just
  // in the order they are stored in the matrix
  for (unsigned int i = 0; i < size; ++i) {
    print_column(cm.get_column(i), i, size);
  }
  std::cout << "\n";
}

void print_matrix(RU_matrix& rum, unsigned int size){
  std::cout << "RU_matrix:\n";
  std::cout << "R:\n";
  for (unsigned int i = 0; i < size; ++i) {
    print_column(rum.get_column(i), i, size);
  }
  std::cout << "\n";
  std::cout << "U:\n";
  for (unsigned int i = 0; i < size; ++i) {
    print_column(rum.get_column(i, false), i, size);
  }
  std::cout << "\n";
}

int main() {
  Simplex_tree<> st;
  build_simplex_tree(st); //could be any other way to build a simplex tree
  auto size = st.num_simplices();

  //reserving space with `size` is not mandatory but recommended for better performances.
  Base_matrix bm(size);
  RU_matrix rum(size);
  Chain_matrix cm(size);

  //filling the matrices with the boundary matrix computed from the simplex tree
  unsigned int id = 0;
  for (auto sh : st.filtration_simplex_range()){
    //identifying the simplex such that the IDs are strictly increasing in the order of filtration.
    st.assign_key(sh, id++);

    //creating boundary
    std::vector<unsigned int> boundary;
    for (auto b : st.boundary_simplex_range(sh)){
      boundary.push_back(st.key(b));
    }
    std::sort(boundary.begin(), boundary.end());  //boundaries have to ordered

    //insertion in the matrices, the id is continuously increasing from 0 so no need to give it as an argument.
    bm.insert_boundary(boundary);
    rum.insert_boundary(boundary);
    cm.insert_boundary(boundary);
  }

  //content of the matrices
  print_matrix(bm, size);
  print_matrix(rum, size);
  print_matrix(cm, size);
}
