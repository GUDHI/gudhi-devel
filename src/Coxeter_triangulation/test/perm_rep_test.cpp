#include <gudhi/Permutahedral_representation.h>

int main() {
  typedef std::vector<int> Vertex;
  typedef std::vector<std::size_t> Part;
  typedef std::vector<Part> Partition;
  typedef Gudhi::coxeter_triangulation::Permutahedral_representation<Vertex, Partition> Simplex_handle;
  Vertex v = {0,0,0};
  Partition p = {Part({1}), Part({0,2,3})};
  Simplex_handle s(v, p);
  std::cout << s << "\n";
  return 0;
}
