#include <gudhi/Persistence_on_a_line.h>
#include <vector>
#include <iostream>

int main() {
  std::vector<float> data{ 0., 1.5, .7, 2.8, 3.1, -1., 2. };
  Gudhi::persistent_cohomology::compute_persistence_of_function_on_line(data,
      [](float b, float d){ std::cout << b << ' ' << d << '\n'; });
}
