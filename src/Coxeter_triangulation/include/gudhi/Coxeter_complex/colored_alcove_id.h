#ifndef COLORED_ALCOVE_ID_OUT_
#define COLORED_ALCOVE_ID_OUT_

#include <gudhi/Simple_coxeter_system.h>

namespace Gudhi {

std::ostream& operator<<(std::ostream & os, const typename Simple_coxeter_system::Alcove_id& a_id) {
  int i = 0, j = 0;
  os << "[";
  if (a_id.empty())
    return os;
  auto a_it = a_id.begin();
  if (a_it == a_id.end() - 1)
    os << "\033[1;32m" << *a_it++ << "\033[0m";
  else if (j == i)
    os << "\033[1;31m" << *a_it++ << "\033[0m";
  else
    os << *a_it++;
  if (j == 0) {
    i++;
    j = i;
  }
  else
    j--;
  while (a_it != a_id.end()) {
    if (a_it == a_id.end() - 1)
      os << ", \033[1;32m" << *a_it++ << "\033[0m";
    else if (j == i)
      os << ", \033[1;31m" << *a_it++ << "\033[0m";
    else
      os << ", " << *a_it++;
    if (j == 0) {
      i++;
      j = i;
    }
    else
      j--;
  }
  std::cout << "]";
  return os;
}

}

#endif
