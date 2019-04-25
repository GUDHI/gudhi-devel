#include <gudhi/Manifold_tracing.h>
#include <gudhi/Implicit_function_intersection_oracle.h>

#include "../example/functions/sphere_S1_in_R2.h"

int main() {
  Gudhi::Manifold_tracing mt;

  /* function definition */ 
  double r = 5;
  Function_S1_in_R2 fun(r);
  Gudhi::Implicit_function_intersection_oracle<Function_S1_in_R2> oracle(fun);
  
  
  return 0;
}

