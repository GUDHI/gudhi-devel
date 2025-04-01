/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2020 Inria
 *
 *    Modification(s):
 *      - 2024/03 Vincent Rouvreau: Renamed Alpha_complex_factory as Delaunay_complex_factory for DelaunayCechComplex.
 *                                  Factorize create_complex
 *      - 2024/10 Vincent Rouvreau: Add square root filtration values interface
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef INCLUDE_DELAUNAY_COMPLEX_FACTORY_H_
#define INCLUDE_DELAUNAY_COMPLEX_FACTORY_H_

#include <gudhi/Simplex_tree.h>
#include <gudhi/Alpha_complex.h>
#include <gudhi/Alpha_complex_3d.h>
#include <gudhi/Alpha_complex_options.h>
#include <gudhi/MEB_filtration.h>

#include <boost/range/adaptor/transformed.hpp>

#include "Simplex_tree_interface.h"

#include <vector>
#include <memory>  // for std::unique_ptr
#include <stdexcept>
#include <cmath>  // for std::sqrt

namespace Gudhi {

namespace delaunay_complex {

enum class Delaunay_filtration : char {
    NONE  = 'n',   ///< Delaunay complex
    CECH  = 'c',   ///< Delaunay Cech complex
    ALPHA = 'a',   ///< Alpha complex
};

// template Functor that transforms a CGAL point to a vector of double as expected by cython
template<typename CgalPointType, bool Weighted>
struct Point_cgal_to_cython;

// Specialized Unweighted Functor
template<typename CgalPointType>
struct Point_cgal_to_cython<CgalPointType, false>
{
    std::vector<double> operator()(CgalPointType const& point) const
    {
        std::vector<double> vd;
        vd.reserve(point.dimension());
        for (auto coord = point.cartesian_begin(); coord != point.cartesian_end(); coord++)
            vd.push_back(CGAL::to_double(*coord));
        return vd;
    }
};

// Specialized Weighted Functor
template<typename CgalPointType>
struct Point_cgal_to_cython<CgalPointType, true>
{
    std::vector<double> operator()(CgalPointType const& weighted_point) const
    {
        const auto& point = weighted_point.point();
        return Point_cgal_to_cython<decltype(point), false>()(point);
    }
};

// Function that transforms a cython point (aka. a vector of double) to a CGAL point
template <typename CgalPointType>
static CgalPointType pt_cython_to_cgal(std::vector<double> const& vec)
{
    return CgalPointType(vec.size(), vec.begin(), vec.end());
}

template <typename Delaunay_complex, typename Kernel, bool Weighted, typename Point_cloud>
bool create_complex(Delaunay_complex& delaunay_complex, Simplex_tree_interface* simplex_tree,
        const Point_cloud& points, double max_alpha_square,
        bool exact_version, Delaunay_filtration filtration, bool output_squared_values)
{
    if (filtration == Delaunay_filtration::CECH) {
        if (Weighted)
            throw std::invalid_argument("Weighted Delaunay-Cech complex is not available");
        // Construct the Delaunay complex
        bool result = delaunay_complex.create_complex(*simplex_tree,
                std::numeric_limits<Simplex_tree_interface::Filtration_value>::infinity(),
                exact_version,
                true);
        if (result == true) {
            // Construct the Delaunay-Cech complex by assigning filtration values with MEB
            if (!output_squared_values) {
                Gudhi::cech_complex::assign_MEB_filtration<false>(Kernel(), *simplex_tree, points);
                simplex_tree->prune_above_filtration(std::sqrt(max_alpha_square));
            } else {
                Gudhi::cech_complex::assign_MEB_filtration<true>(Kernel(), *simplex_tree, points);
                simplex_tree->prune_above_filtration(max_alpha_square);
            }
        }
        return result;
    } else {
        if (output_squared_values)
            return delaunay_complex.template create_complex<true>(*simplex_tree, max_alpha_square,
                    exact_version, filtration == Delaunay_filtration::NONE);
        else
            return delaunay_complex.template create_complex<false>(*simplex_tree, max_alpha_square,
                    exact_version, filtration == Delaunay_filtration::NONE);
    }

}

class Abstract_delaunay_complex
{
public:
    virtual std::vector<double> get_point(int vh) = 0;

    virtual bool create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square,
            Delaunay_filtration filtration, bool output_squared_values = false) = 0;

    virtual std::size_t num_vertices() const = 0;

    virtual ~Abstract_delaunay_complex() = default;

};

template <typename Kernel, bool Weighted = false>
class Delaunay_complex_t final : public Abstract_delaunay_complex
{
private:
    using Bare_point = typename Kernel::Point_d;
    using Point = std::conditional_t<Weighted, typename Kernel::Weighted_point_d,
            typename Kernel::Point_d>;
    using Delaunay_complex = Gudhi::alpha_complex::Alpha_complex<Kernel, Weighted>;

public:
    Delaunay_complex_t(const std::vector<std::vector<double>>& points, bool exact_version)
        : exact_version_(exact_version),
          points_(boost::begin(boost::adaptors::transform(points, pt_cython_to_cgal<Bare_point>)),
                  boost::end(boost::adaptors::transform(points, pt_cython_to_cgal<Bare_point>))),
          delaunay_complex_(points_) {
    }

    Delaunay_complex_t(const std::vector<std::vector<double>>& points,
            const std::vector<double>& weights, bool exact_version)
        : exact_version_(exact_version),
          points_(boost::begin(boost::adaptors::transform(points, pt_cython_to_cgal<Bare_point>)),
                  boost::end(boost::adaptors::transform(points, pt_cython_to_cgal<Bare_point>))),
          delaunay_complex_(points_, weights) {
    }

    virtual std::vector<double> get_point(int vh) override
    {
        // Can be a Weighted or a Bare point in function of Weighted
        return Point_cgal_to_cython<Point, Weighted>()(delaunay_complex_.get_point(vh));
    }

    virtual bool create_simplex_tree(Simplex_tree_interface* simplex_tree, double max_alpha_square,
            Delaunay_filtration filtration, bool output_squared_values) override
    {
        return create_complex<Delaunay_complex, Kernel, Weighted, std::vector<Bare_point>>(delaunay_complex_, simplex_tree,
                points_, max_alpha_square,
                exact_version_, filtration,
                output_squared_values);
    }

    virtual std::size_t num_vertices() const override {
        return delaunay_complex_.num_vertices();
    }

private:
    bool exact_version_;
    std::vector<Bare_point> points_;
    Delaunay_complex delaunay_complex_;
};

}  // namespace delaunay_complex

}  // namespace Gudhi

#endif  // INCLUDE_DELAUNAY_COMPLEX_FACTORY_H_
