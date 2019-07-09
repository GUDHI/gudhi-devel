#include <gudhi/Permutahedral_representation.h>

int main() {
  typedef std::vector<int> Vertex;
  typedef std::vector<std::size_t> Part;
  typedef std::vector<Part> Partition;
  typedef Gudhi::coxeter_triangulation::Permutahedral_representation<Vertex, Partition> Simplex_handle;
  Vertex v0 = {0,0,0};
  Partition omega = {Part({1}), Part({2}), Part({0,3})};
  Simplex_handle s(v0, omega);
  std::cout << "Simplex" << s << ", its dimension is " << s.dimension() << "\n";
  auto print_vertex =
    [](const Vertex& v) {
      std::cout << " (";
      if (v.empty())
	std::cout << ")";
      auto v_it = v.begin();
      std::cout << *v_it++;
      for (; v_it != v.end(); ++v_it)
	std::cout << ", " << *v_it;  
      std::cout << ")\n";
    };
  std::cout << "Vertices of " << s << ":\n";
  for (auto& v: s.vertex_range())
    print_vertex(v);  
  std::cout << "Facets of " << s << ":\n";
  for (auto& f: s.facet_range())
    std::cout << f << "\n";
  std::cout << "Faces of " << s << " of dimension 0:\n";
  for (auto& f: s.face_range(0))
    std::cout << f << "\n";
  std::cout << "Cofacets of " << s << ":\n";
  for (auto& f: s.cofacet_range())
    std::cout << f << "\n";
  std::cout << "Cofaces of " << s << " of dimension 2:\n";
  for (auto& f: s.coface_range(2))
    std::cout << f << "\n";

  Vertex v1 = {0,1,1};
  Partition omega1 = {Part({0,1,2,3})};
  Simplex_handle s1(v1, omega1);
  Vertex v2 = {0,1,0};
  Partition omega2 = {Part({0}), Part({1,2,3})};
  Simplex_handle s2(v2, omega2);
  std::cout << "Simplex " << s1 << " is " << (s1.is_face_of(s)? "": "not ") << "a face of " << s << "\n";
  std::cout << "Simplex " << s2 << " is " << (s2.is_face_of(s)? "": "not ") << "a face of " << s << "\n";

  return 0;
}
