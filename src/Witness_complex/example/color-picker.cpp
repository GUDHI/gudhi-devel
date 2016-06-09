#include <iostream>
#include <assert.h>

struct Color {
  double r, g, b;
  Color() : r(0), g(0), b(0) { }
  Color(double r_, double g_, double b_)
    : r(r_), g(g_), b(b_) { } 
};

// Convert a double in [0,1] to a color
template< typename ColorType >
void double_to_color(double t, ColorType & color) {
  assert(t >= 0 && t <= 1);
  double s = 1.0/6;
  if (t >= 0 && t <= s)
    color = ColorType(1, 6*t, 0);
  else if (t > s && t <= 2*s)
    color = ColorType(1-6*(t-s), 1, 0);
  else if (t > 2*s && t <= 3*s)
    color = ColorType(0, 1, 6*(t-2*s));
  else if (t > 3*s && t <= 4*s)
    color = ColorType(0, 1-6*(t-3*s), 1);
  else if (t > 4*s && t <= 5*s)
    color = ColorType(6*(t-4*s), 0, 1);
  else
    color = ColorType(1, 0, 1-6*(t-5*s));
}

int main (int argc, char * const argv[]) {
  double number_to_convert = 0;
  std::cout << "Enter a real number in [0,1]: ";
  std::cin >> number_to_convert;
  Color color;
  double_to_color(number_to_convert, color);
  std::cout << "The corresponding color in rgb is: (" << color.r << ", " << color.g << ", " << color.b << ")\n";
}
