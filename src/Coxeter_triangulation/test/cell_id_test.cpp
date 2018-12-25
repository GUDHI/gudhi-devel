#include <gudhi/Coxeter_triangulation/Cell_id.h>


int main() {
  Gudhi::Cell_id c_id_0;
  std::cout << "c_id_0: " << c_id_0 << " Level: " << c_id_0.level() << std::endl;
  Gudhi::Cell_id c_id_1(1);
  std::cout << "c_id_1: " << c_id_1 << " Level: " << c_id_1.level() << std::endl;
  Gudhi::Cell_id c_id_2(1, 2);
  std::cout << "c_id_2: " << c_id_2 << " Level: " << c_id_2.level() << std::endl;

  c_id_2.push_back(0);
  std::cout << "c_id_2: " << c_id_2 << " Level: " << c_id_2.level() << std::endl;
  c_id_2.push_back(-1, true);
  std::cout << "c_id_2: " << c_id_2 << " Level: " << c_id_2.level() << std::endl;
  c_id_2.pop_back();
  std::cout << "c_id_2: " << c_id_2 << " Level: " << c_id_2.level() << std::endl;
  c_id_2.reserve(6);
  std::cout << "c_id_2: " << c_id_2 << " Level: " << c_id_2.level() << std::endl;
  c_id_2.resize(6);
  std::cout << "c_id_2: " << c_id_2 << " Level: " << c_id_2.level() << std::endl;
  c_id_2.resize(4);
  std::cout << "c_id_2: " << c_id_2 << " Level: " << c_id_2.level() << std::endl;
  c_id_2.set_dimension(3);
  c_id_2.set_level(1.5);
  std::cout << "c_id_2: " << c_id_2 << " Level: " << c_id_2.level() << std::endl;
  std::cout << "c_id_0 < c_id_2: " << (c_id_0 < c_id_2) << std::endl;
  std::cout << "*c_id_2.begin(): " << *c_id_2.begin() << std::endl;
  std::cout << "*c_id_2.mask_begin(): " << *c_id_2.mask_begin() << std::endl;
  std::cout << "c_id_0 == c_id_1: " << (c_id_0 == c_id_1) << std::endl;
  return 0;
}

