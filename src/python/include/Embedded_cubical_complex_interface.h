#ifndef INCLUDE_EMBEDDED_CUBICAL_COMPLEX_INTERFACE_H_
#define INCLUDE_EMBEDDED_CUBICAL_COMPLEX_INTERFACE_H_

#include <gudhi/embedded_complex_release.hpp>

#include <iostream>
#include <vector>
#include <string>

namespace Gudhi {

namespace cubical_complex {

    template<typename Embedded_cubical_complex_options = Bitmap_cubical_complex_base<double>>;
    class Embedded_cubical_complex_interface : public Embedded_cubical_complex<Embedded_cubical_complex_options>{
        Embedded_cubical_complex_interface(const std::vector<unsigned>& dimensions,
                            const std::vector<double>& top_dimensional_cells)
        : Embedded_cubical_complex<CubicalComplexOptions>(dimensions, top_dimensional_cells) {
        }
    };

}   //namespace cubical_complex

}   //namespace Gudhi

#endif