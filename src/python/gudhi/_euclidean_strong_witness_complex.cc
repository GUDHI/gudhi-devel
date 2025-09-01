/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - 2025/03 Thibaud Kloczko: Use nanobind instead of Cython for python bindings.
 *      - 2025/04 Hannah Schreiber: Re-add possibility of tensors (numpy, torch etc.) as input.
 *      - YYYY/MM Author: Description of the modification
 */

#include <optional>
#include <vector>
#include <cstddef>
#include <limits>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>

#include <CGAL/Epick_d.h>

#include <gudhi/Simplex_tree.h>
#include <gudhi/Euclidean_strong_witness_complex.h>
#include <python_interfaces/Simplex_tree_interface.h>
#include <python_interfaces/points_utils.h>

namespace Gudhi {
namespace witness_complex {

class Euclidean_strong_witness_complex_interface
{
  using Dynamic_kernel = CGAL::Epick_d<CGAL::Dynamic_dimension_tag>;
  using Point_d = Dynamic_kernel::Point_d;
  using Simplex_key = typename Simplex_tree<>::Simplex_key;

 public:
  Euclidean_strong_witness_complex_interface() : number_of_landmarks_(0) {}

  Euclidean_strong_witness_complex_interface(const Sequence2D& landmarks, const Sequence2D& witnesses)
      : number_of_landmarks_(landmarks.size())
  {
    // Cannot build an Euclidean_strong_witness_complex with empty landmarks, as it segfaults
    if (number_of_landmarks_ != 0) {
      // TODO: copy could be avoided here with custom iterator or some transform (boost::adaptors::transform?)
      // witness complex just needs a begin/end method deferrenceating in a Point_d
      std::vector<Point_d> lm;
      lm.reserve(landmarks.size());
      for (auto& landmark : landmarks) lm.emplace_back(landmark.begin(), landmark.end());
      witness_complex_.emplace(lm, witnesses);
    }
  }

  // TODO: remove one copy by directly constructing std::vector<Point_d> instead of std::vector<std::vector<double>>
  Euclidean_strong_witness_complex_interface(const Tensor2D& landmarks, const Tensor2D& witnesses)
      : Euclidean_strong_witness_complex_interface(_get_sequence_from_tensor(landmarks),
                                                   _get_sequence_from_tensor(witnesses))
  {}

  void create_simplex_tree(Simplex_tree_interface* simplex_tree,
                           double max_alpha_square,
                           std::size_t limit_dimension = std::numeric_limits<std::size_t>::max())
  {
    if (witness_complex_) witness_complex_->create_complex(*simplex_tree, max_alpha_square, limit_dimension);
  }

  std::vector<double> get_point(unsigned int vh)
  {
    std::vector<double> vd;
    if (vh < number_of_landmarks_) {
      Point_d ph = witness_complex_->get_point(vh);
      for (auto coord = ph.cartesian_begin(); coord < ph.cartesian_end(); ++coord) {
        vd.push_back(*coord);
      }
    }
    return vd;
  }

 private:
  unsigned int number_of_landmarks_;
  std::optional<Euclidean_strong_witness_complex<Dynamic_kernel> > witness_complex_;
};

}  // namespace witness_complex
}  // namespace Gudhi

namespace nb = nanobind;
namespace gwc = Gudhi::witness_complex;
using gwci = gwc::Euclidean_strong_witness_complex_interface;

NB_MODULE(_euclidean_strong_witness_complex_ext, m)
{
  m.attr("__license__") = "GPL v3";

  nb::class_<gwci>(m, "Euclidean_strong_witness_complex_interface")
      .def(nb::init<>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Sequence2D&, const Sequence2D&>(), nb::call_guard<nb::gil_scoped_release>())
      .def(nb::init<const Tensor2D&, const Tensor2D&>(), nb::call_guard<nb::gil_scoped_release>())
      .def("create_simplex_tree",
           &gwci::create_simplex_tree,
           nb::arg("simplex_tree"),
           nb::arg("max_alpha_square"),
           nb::arg("limit_dimension") = std::numeric_limits<std::size_t>::max(),
           nb::call_guard<nb::gil_scoped_release>())
      .def("get_point",
           &gwci::get_point,
           R"doc(
This function returns the point corresponding to a given vertex.
Args:
    arg (int): The vertex id.
    
Returns:
    list of float: The coordinates of the vertex.
           )doc");
}

//
// _euclidean_strong_witness_complex.cc ends here
