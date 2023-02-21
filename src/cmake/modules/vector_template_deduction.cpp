#include <vector>

// cf. https://gitlab.kitware.com/cmake/cmake/-/blob/master/Source/Checks/cm_cxx17_check.cpp
// CMake only tests that std::invoke_result and std::optional works to tell if c++17 is available or not
// This is a workaround to https://github.com/GUDHI/gudhi-devel/issues/816
template<class Range>
void test_cpp17(Range diag1) {
  // std::vector template deduction
  std::vector v1(std::begin(diag1), std::end(diag1));
}

int main()
{
  std::vector<int> range {1,2,3};
  test_cpp17(range);
}