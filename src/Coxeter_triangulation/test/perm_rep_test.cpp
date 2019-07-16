#include <gudhi/Permutahedral_representation.h>
#include <gudhi/Permutahedral_representation/join.h>

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

  Vertex v1 = {0,0,0};
  Partition omega1 = {Part({0,1,2,3})};
  Simplex_handle s1(v1, omega1);
  Vertex v2 = {0,1,0};
  Partition omega2 = {Part({0}), Part({1,2,3})};
  Simplex_handle s2(v2, omega2);
  std::cout << "Simplex " << s1 << " is " << (s1.is_face_of(s)? "": "not ") << "a face of " << s << "\n";
  std::cout << "Simplex " << s2 << " is " << (s2.is_face_of(s)? "": "not ") << "a face of " << s << "\n";

  v1 = {1,1,0};
  omega1 = {Part({2}), Part({0,1,3})};
  s1 = Simplex_handle(v1, omega1);
  std::cout << "The smallest coface of " << s1 << " is " << s1.smallest_coface() << "\n";
  std::cout << "The greatest face of " << s1 << " is " << s1.greatest_face() << "\n";

  std::cout << "The smallest coface of " << s2 << " is " << s2.smallest_coface() << "\n";
  std::cout << "The greatest face of " << s2 << " is " << s2.greatest_face() << "\n";

  Vertex v3 = {0,0,-1};
  Partition omega3 = {Part({0,1,2}), Part({3})};
  Simplex_handle s3(v3, omega3);
  std::vector<Simplex_handle> simplices = {s2,s3};
  std::cout << "The join of " << s2 << " and " << s3 << " is " << join(simplices) << "\n";
  
  return 0;
}
