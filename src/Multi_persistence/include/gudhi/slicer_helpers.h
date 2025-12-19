/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Loiseaux, Hannah Schreiber
 *
 *    Copyright (C) 2023-25 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

/**
 * @file slicer_helpers.h
 * @author David Loiseaux, Hannah Schreiber
 * @brief Contains the helper methods @ref Gudhi::multi_persistence::build_complex_from_scc_file,
 * @ref Gudhi::multi_persistence::write_complex_to_scc_file, @ref Gudhi::multi_persistence::build_slicer_from_scc_file,
 * @ref Gudhi::multi_persistence::build_complex_from_bitmap and @ref Gudhi::multi_persistence::build_slicer_from_bitmap.
 */

#ifndef MP_SLICER_HELPERS_H_
#define MP_SLICER_HELPERS_H_

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <ostream>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <set>
#include <limits>
#include <iomanip>
#include <cmath>

#ifdef GUDHI_USE_TBB
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/parallel_for.h>
#endif

#include <boost/range/iterator_range_core.hpp>

#include <gudhi/Debug_utils.h>
#include <gudhi/simple_mdspan.h>
#include <gudhi/Multi_parameter_filtered_complex.h>
#include <gudhi/Bitmap_cubical_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Multi_filtration/multi_filtration_utils.h>
#include <gudhi/Multi_filtration/multi_filtration_conversions.h>
#include <gudhi/Multi_persistence/Line.h>

namespace Gudhi {
namespace multi_persistence {

/**
 * @ingroup multi_persistence
 *
 * @brief Builds a complex for the scc format file given. Assumes that every index appearing in a boundary in the file
 * corresponds to a real line in the file (for example, the lowest dimension has always empty boundaries).
 *
 * @tparam MultiFiltrationValue Filtration value class respecting the @ref MultiFiltrationValue concept. It will be
 * used as filtration value type of the new complex.
 * @param inFilePath Path to scc file.
 * @param isRivetCompatible Set to true if the file is written such that Rivet can read it. See TODO ref.
 * Default value: false.
 * @param isReversed Set to true if the cells in the file are written in increasing dimension order instead of
 * the standard decreasing order. Default value: false.
 * @param shiftDimensions Indicates if there is a shift in the dimension written in the file: if the value is 0, it
 * means that the smallest dimension is 0, if the value is positive, the smallest dimension is assumed to be
 * `shiftDimensions` instead of 0, and if the value is negative, the `abs(shiftDimensions)` smallest dimensions in
 * the file are ignored and the smallest remaining dimension is interpreted as 0. Default value: 0.
 */
template <class MultiFiltrationValue>
inline Multi_parameter_filtered_complex<MultiFiltrationValue> build_complex_from_scc_file(
    const std::string& inFilePath,
    bool isRivetCompatible = false,
    bool isReversed = false,
    int shiftDimensions = 0)
{
  using Fil = MultiFiltrationValue;
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using Index = typename Complex::Index;

  std::string line;
  std::ifstream file(inFilePath);
  unsigned int numberOfParameters;

  if (!file.is_open()) {
    // TODO: throw instead?
    std::cerr << "Unable to open input file: " << inFilePath << '\n';
    file.setstate(std::ios::failbit);
    return Complex();
  }

  auto error = [&file](const std::string& msg) {
    file.close();
    throw std::invalid_argument(msg);
  };
  auto is_comment_or_empty_line = [](const std::string& line) -> bool {
    size_t current = line.find_first_not_of(' ', 0);
    if (current == std::string::npos) return true;  // is empty line
    if (line[current] == '#') return true;          // is comment
    return false;
  };

  while (getline(file, line, '\n') && is_comment_or_empty_line(line));
  if (!file) error("Empty file!");

  if (isRivetCompatible && line != "firep") error("Wrong file format. Should start with 'firep'.");
  if (!isRivetCompatible && line != "scc2020") error("Wrong file format. Should start with 'scc2020'.");

  while (getline(file, line, '\n') && is_comment_or_empty_line(line));
  if (!file) error("Premature ending of the file. Stops before numbers of parameters.");

  if (isRivetCompatible) {
    numberOfParameters = 2;
    getline(file, line, '\n');  // second rivet label
  } else {
    std::size_t current = line.find_first_not_of(' ', 0);
    std::size_t next = line.find_first_of(' ', current);
    numberOfParameters = std::stoi(line.substr(current, next - current));
  }

  while (getline(file, line, '\n') && is_comment_or_empty_line(line));
  if (!file) error("Premature ending of the file. Not a single cell was specified.");

  std::vector<unsigned int> counts;
  Index numberOfCells = 0;
  counts.reserve(line.size() + shiftDimensions);
  std::size_t current = line.find_first_not_of(' ', 0);
  if (shiftDimensions != 0 && isReversed && current != std::string::npos) {
    if (shiftDimensions > 0) {
      counts.resize(shiftDimensions, 0);
    } else {
      for (int i = shiftDimensions; i < 0 && current != std::string::npos; ++i) {
        std::size_t next = line.find_first_of(' ', current);
        current = line.find_first_not_of(' ', next);
      }
    }
  }
  while (current != std::string::npos) {
    std::size_t next = line.find_first_of(' ', current);
    counts.push_back(std::stoi(line.substr(current, next - current)));
    numberOfCells += counts.back();
    current = line.find_first_not_of(' ', next);
  }
  if (shiftDimensions != 0 && !isReversed) {
    counts.resize(counts.size() + shiftDimensions, 0);
  }

  std::size_t dimIt = 0;
  while (dimIt < counts.size() && counts[dimIt] == 0) ++dimIt;

  if (dimIt == counts.size()) return Complex();

  std::size_t shift = isReversed ? 0 : counts[dimIt];
  unsigned int nextShift = isReversed ? 0 : counts.size() == 1 ? 0 : counts[dimIt + 1];
  unsigned int tmpNextShift = counts[dimIt];

  auto get_boundary = [&isReversed, &numberOfCells](
                          const std::string& line, std::size_t start, std::size_t shift) -> std::vector<Index> {
    std::vector<Index> res;
    res.reserve(line.size() - start);
    std::size_t current = line.find_first_not_of(' ', start);
    while (current != std::string::npos) {
      std::size_t next = line.find_first_of(' ', current);
      Index idx = std::stoi(line.substr(current, next - current)) + shift;
      res.push_back(isReversed ? idx : numberOfCells - 1 - idx);
      current = line.find_first_not_of(' ', next);
    }
    std::sort(res.begin(), res.end());
    return res;
  };
  auto get_filtration_value = [numberOfParameters, &error](const std::string& line, std::size_t end) -> Fil {
    std::vector<typename Fil::value_type> res;
    res.reserve(end);
    bool isPlusInf = true;
    bool isMinusInf = true;
    std::size_t current = line.find_first_not_of(' ', 0);
    while (current < end) {
      std::size_t next = line.find_first_of(' ', current);
      res.push_back(std::stod(line.substr(current, next - current)));
      if (isPlusInf && res.back() != Fil::T_inf) isPlusInf = false;
      if (isMinusInf && res.back() != Fil::T_m_inf) isMinusInf = false;
      current = line.find_first_not_of(' ', next);
    }
    if (isPlusInf) return Fil::inf(numberOfParameters);
    if (isMinusInf) return Fil::minus_inf(numberOfParameters);
    if (res.size() % numberOfParameters != 0) error("Wrong format. The number of parameters does not match.");
    return Fil(res.begin(), res.end(), numberOfParameters);
  };

  typename Complex::Boundary_container boundaries(numberOfCells);
  typename Complex::Dimension_container dimensions(numberOfCells);
  typename Complex::Filtration_value_container filtrationValues(numberOfCells);
  std::size_t i = 0;
  // because of possible negative dimension shifts, the document should not always be read to the end
  // therefore `dimIt < counts.size()` is also a stop condition
  while (getline(file, line, '\n') && dimIt < counts.size()) {
    if (!is_comment_or_empty_line(line)) {
      std::size_t sep = line.find_first_of(';', 0);
      filtrationValues[i] = get_filtration_value(line, sep);
      boundaries[i] = get_boundary(line, sep + 1, shift);
      dimensions[i] = isReversed ? dimIt : counts.size() - 1 - dimIt;

      --counts[dimIt];
      while (dimIt < counts.size() && counts[dimIt] == 0) {
        ++dimIt;
        if (dimIt != counts.size()) {
          shift += nextShift;
          nextShift = isReversed ? tmpNextShift : dimIt < counts.size() - 1 ? counts[dimIt + 1] : 0;
          tmpNextShift = counts[dimIt];
        }
      }
      ++i;
    }
  }

  if (!isReversed) {  // to order by dimension
    std::reverse(dimensions.begin(), dimensions.end());
    std::reverse(boundaries.begin(), boundaries.end());
    std::reverse(filtrationValues.begin(), filtrationValues.end());
  }

  file.close();

  return Complex(std::move(boundaries), std::move(dimensions), std::move(filtrationValues));
}

/**
 * @ingroup multi_persistence
 *
 * @brief Writes the given complex into a file with scc format. Assumes that every index appearing in a boundary of
 * the complex corresponds to an existing index in the complex (for example, the lowest dimension has always empty
 * boundaries).
 *
 * @tparam MultiFiltrationValue Filtration value of the given complex.
 * @param outFilePath Path with file name into which to write.
 * @param complex Complex to write. Every index appearing in a boundary of the complex has to correspond to an existing
 * index in the complex
 * @param degree TODO Default value: -1.
 * @param rivetCompatible Set to true if the written file has to be Rivet compatible. Note that Rivet only accepts
 * bi-filtrations. Default value: false.
 * @param ignoreLastGenerators Set to true, if the generators with last dimension in the list should be ignored
 * (maximal dimension by default, minimal dimension if `reverse` is true). Default value: false.
 * @param stripComments Set to true, if no comment should be written in the file (comments are lines starting with `#`
 * and which are ignored when read). Default value: false.
 * @param reverse Set to true if the generators should be written in increasing order of dimension instead of
 * decreasing. Default value: false.
 */
template <class MultiFiltrationValue>
inline void write_complex_to_scc_file(const std::string& outFilePath,
                                      const Multi_parameter_filtered_complex<MultiFiltrationValue>& complex,
                                      int degree = -1,
                                      bool rivetCompatible = false,
                                      bool ignoreLastGenerators = false,
                                      bool stripComments = false,
                                      bool reverse = false)
{
  if (!complex.is_ordered_by_dimension()) {
    // other solution would be to call build_permuted_complex ourself, but this is a good way to make the
    // user aware of it.
    throw std::invalid_argument(
        "The given complex has to be ordered by dimension. If it is not the case, call this method with "
        "`build_permuted_complex(complex).first` or `build_permuted_complex(complex, permutation_by_dim)` instead.");
    return;
  }

  unsigned int numberOfParameters = complex.get_number_of_parameters();

  std::ofstream file(outFilePath);

  if (rivetCompatible)
    file << "firep\n";
  else
    file << "scc2020\n";

  // TODO: change line for gudhi
  if (!stripComments && !rivetCompatible)
    file << "# This file was generated by multipers (https://github.com/DavidLapous/multipers).\n";

  if (!stripComments && !rivetCompatible) file << "# Number of parameters\n";

  if (rivetCompatible) {
    GUDHI_CHECK(numberOfParameters == 2, "Rivet only handles bifiltrations.");
    file << "Filtration 1\n";
    file << "Filtration 2\n";
  } else {
    file << std::to_string(numberOfParameters) << "\n";
  }

  if (!stripComments) file << "# Sizes of generating sets\n";

  using Fil = MultiFiltrationValue;

  int maxDim = complex.get_max_dimension();
  int minDim = maxDim;
  const auto& dimensions = complex.get_dimensions();

  std::vector<std::vector<std::size_t> > indicesByDim(maxDim + 1);
  std::vector<std::size_t> shiftedIndices(complex.get_number_of_cycle_generators());
  for (std::size_t i = 0; i < complex.get_number_of_cycle_generators(); ++i) {
    auto dim = dimensions[i];
    minDim = dim < minDim ? dim : minDim;
    auto& atDim = indicesByDim[reverse ? dim : maxDim - dim];
    shiftedIndices[i] = atDim.size();
    atDim.push_back(i);
  }
  if (degree < 0) degree = minDim;
  int minIndex = reverse ? degree - 1 : 0;
  int maxIndex = reverse ? maxDim : maxDim - degree + 1;
  maxIndex = std::max(maxIndex, -1);
  if (ignoreLastGenerators) maxIndex--;
  if (rivetCompatible) minIndex = maxIndex - 2;

#ifdef DEBUG_TRACES
  std::cout << "minDim = " << minDim << " maxDim = " << maxDim << " minIndex = " << minIndex
            << " maxIndex = " << maxIndex << " degree = " << degree << std::endl;
#endif

  auto print_fil_values = [&](const Fil& fil) {
    GUDHI_CHECK(fil.num_parameters() == numberOfParameters, "Filtration value has wrong number of parameters.");
    for (unsigned int g = 0; g < fil.num_generators(); ++g) {
      for (unsigned int p = 0; p < fil.num_parameters(); ++p) {
        file << fil(g, p) << " ";
      }
    }
  };

  if (minIndex < 0) file << 0 << " ";
  for (int i = 0; i < minIndex; ++i) file << 0 << " ";
  for (int i = std::max(minIndex, 0); i <= std::min(maxDim, maxIndex); ++i) {
    file << indicesByDim[i].size() << " ";
  }
  if (!rivetCompatible)
    for (int i = maxIndex + 1; i <= maxDim; ++i) file << 0 << " ";
  if (maxIndex > maxDim) file << 0;
  file << "\n";

  file << std::setprecision(std::numeric_limits<typename Fil::value_type>::digits);

  std::size_t startIndex = reverse ? minIndex + 1 : minIndex;
  std::size_t endIndex = reverse ? maxIndex : maxIndex - 1;
  const auto& filtValues = complex.get_filtration_values();
  const auto& boundaries = complex.get_boundaries();
  int currDim;
  if (reverse)
    currDim = minIndex == -1 ? 0 : minIndex;
  else
    currDim = maxIndex == maxDim + 1 ? maxDim + 1 : maxDim;

  if (reverse) {
    if (!stripComments) file << "# Block of dimension " << currDim++ << "\n";
    if (minIndex >= 0) {
      for (auto index : indicesByDim[minIndex]) {
        print_fil_values(filtValues[index]);
        file << ";\n";
      }
    }
  }
  for (std::size_t i = startIndex; i <= endIndex; ++i) {
    if (!stripComments) {
      file << "# Block of dimension " << currDim << "\n";
      if (reverse)
        ++currDim;
      else
        --currDim;
    }
    for (auto index : indicesByDim[i]) {
      print_fil_values(filtValues[index]);
      file << "; ";
      for (auto b : boundaries[index]) file << shiftedIndices[b] << " ";
      file << "\n";
    }
  }
  if (!reverse) {
    if (!stripComments) file << "# Block of dimension " << currDim << "\n";
    if (maxIndex <= maxDim) {
      for (auto index : indicesByDim[maxIndex]) {
        print_fil_values(filtValues[index]);
        file << ";\n";
      }
    }
  }
}

/**
 * @ingroup multi_persistence
 *
 * @brief Builds a complex from the given bitmap. The bitmap here is a grid where each node contains a 1-critical
 * filtration value, which will be interpreted as a vertex in a cubical complex. The filtration values of the higher
 * dimensional cells are deduced by taking at each parameter the maximal value of its facets at this parameter.
 *
 * Note that for the bitmap to represent a valid multi-parameter filtration, all filtration values have to have the
 * same number of parameters. The behaviour is undefined otherwise.
 *
 * @tparam OneCriticalMultiFiltrationValue Filtration value class respecting the @ref MultiFiltrationValue concept.
 * It will be used as filtration value type of the new complex.
 * @param vertexValues Bitmap with 1-critical filtration values. Represented as a single vector, the next input
 * parameter @p shape indicates the shape of the real bitmap.
 * @param shape Shape of the bitmap. E.g., if @p shape is \f$ {3, 4} \f$, then the bitmap is a \f$ (4 x 3) \f$ grid
 * with four lines and three columns. The vector @p vertexValues should then contain 12 elements: the three first
 * elements will be read as the first line, the three next elements as the second line etc. until having 4 lines.
 */
template <class OneCriticalMultiFiltrationValue>
inline Multi_parameter_filtered_complex<OneCriticalMultiFiltrationValue> build_complex_from_bitmap(
    const std::vector<OneCriticalMultiFiltrationValue>& vertexValues,
    const std::vector<unsigned int>& shape)
{
  using Fil = OneCriticalMultiFiltrationValue;
  using Complex = Multi_parameter_filtered_complex<Fil>;
  using Index = typename Complex::Index;
  using Bitmap_cubical_complex_base = Gudhi::cubical_complex::Bitmap_cubical_complex_base<char>;
  using Bitmap_cubical_complex = Gudhi::cubical_complex::Bitmap_cubical_complex<Bitmap_cubical_complex_base>;

  if (shape.empty() || vertexValues.empty()) return Complex();

  unsigned int numberOfParameters = vertexValues[0].num_parameters();

  Bitmap_cubical_complex cub(shape, std::vector<char>(vertexValues.size()), false);

  const unsigned int numberOfSimplices = cub.num_simplices();

  typename Complex::Dimension_container dimensions(numberOfSimplices);
  typename Complex::Boundary_container boundaries(numberOfSimplices);
  unsigned int i = 0;
  for (unsigned int d = 0; d < shape.size() + 1; ++d) {
    for (auto sh : cub.skeleton_simplex_range(d)) {
      cub.assign_key(sh, i);
      dimensions[i] = d;
      auto& col = boundaries[i];
      for (auto b : cub.boundary_simplex_range(sh)) col.push_back(cub.key(b));
      std::sort(col.begin(), col.end());
      ++i;
    }
  }

  auto get_vertices = [&boundaries](Index i) -> std::set<Index> {
    auto rec_get_vertices = [&boundaries](const auto& self, Index i, std::set<Index>& vertices) -> void {
      if (boundaries[i].empty()) {
        vertices.insert(i);
        return;
      }
      for (auto v : boundaries[i]) self(self, v, vertices);
    };
    std::set<Index> vertices;
    rec_get_vertices(rec_get_vertices, i, vertices);
    return vertices;
  };

  typename Complex::Filtration_value_container filtrationValues(numberOfSimplices, Fil(numberOfParameters));

  for (Index g = 0; g < numberOfSimplices; ++g) {
    if constexpr (Gudhi::multi_filtration::RangeTraits<Fil>::is_dynamic_multi_filtration) {
      // should be faster than doing a proper `push_to_least_common_upper_bound` in the loop after
      filtrationValues[g].force_generator_size_to_number_of_parameters(0);
    }
    for (auto v : get_vertices(g)) {
      for (Index p = 0; p < numberOfParameters; ++p) {
        // 1-critical
        filtrationValues[g](0, p) = std::max(filtrationValues[g](0, p), vertexValues[v](0, p));
      }
    }
  }

  return Complex(std::move(boundaries), std::move(dimensions), std::move(filtrationValues));
}

/**
 * @ingroup multi_persistence
 *
 * @brief Builds a complex from the given simplex tree. The complex will be ordered by dimension.
 *
 * @note The key values in the simplex tree nodes will be overwritten.
 *
 * @tparam MultiFiltrationValue Class following the @ref MultiFiltrationValue concept.
 * @tparam SimplexTreeOptions Class following the @ref SimplexTreeOptions concept. Additionally, if
 * `SimplexTreeOptions::Filtration_value` and `MultiFiltrationValue` are not the same type, there must
 * be a method `as_type` taking `SimplexTreeOptions::Filtration_value` as argument and returning the value as an
 * `MultiFiltrationValue` type. See @ref Gudhi::multi_filtration::as_type for implementations for
 * @ref Gudhi::multi_filtration::Multi_parameter_filtration,
 * @ref Gudhi::multi_filtration::Dynamic_multi_parameter_filtration and
 * @ref Gudhi::multi_filtration::Degree_rips_bifiltration.
 * @param simplexTree Simplex tree to convert. The key values of the simplex tree will be overwritten.
 */
template <class MultiFiltrationValue, class SimplexTreeOptions>
inline Multi_parameter_filtered_complex<MultiFiltrationValue> build_complex_from_simplex_tree(
    Simplex_tree<SimplexTreeOptions>& simplexTree)
{
  // declared here to enable custom `as_type` methods which are not in this namespace.
  using namespace Gudhi::multi_filtration;

  // TODO: is_multi_filtration will discriminate all pre-made multi filtration classes, but not any user made
  // class following the MultiFiltrationValue concept (as it was more thought for inner use). The tests should be
  // re-thought or this one just removed.
  static_assert(RangeTraits<MultiFiltrationValue>::is_multi_filtration,
                "Target filtration value type has to correspond to the MultiFiltrationValue concept.");

  using Complex = Multi_parameter_filtered_complex<MultiFiltrationValue>;

  const unsigned int numberOfSimplices = simplexTree.num_simplices();

  if (numberOfSimplices == 0) return Complex();

  typename Complex::Dimension_container dimensions(numberOfSimplices);
  typename Complex::Boundary_container boundaries(numberOfSimplices);
  typename Complex::Filtration_value_container filtrationValues(numberOfSimplices);

  unsigned int i = 0;
  // keys for boundaries have to be assigned first as we cannot use filtration_simplex_range to ensure that a face
  // appears before its cofaces.
  for (auto sh : simplexTree.complex_simplex_range()) {
    simplexTree.assign_key(sh, i);
    dimensions[i] = simplexTree.dimension(sh);
    ++i;
  }

  // Order simplices by dimension as an ordered Complex is more performant
  std::vector<unsigned int> newToOldIndex(numberOfSimplices);
  std::vector<unsigned int> oldToNewIndex(numberOfSimplices);
  std::iota(newToOldIndex.begin(), newToOldIndex.end(), 0);
  // stable sort to make the new complex more predicable and closer to a lexicographical sort in addition to dimension
  std::stable_sort(newToOldIndex.begin(), newToOldIndex.end(), [&dimensions](unsigned int i, unsigned int j) {
    return dimensions[i] < dimensions[j];
  });
  // Is there a way to directly get oldToNewIndex without constructing newToOldIndex?
  for (unsigned int k = 0; k < numberOfSimplices; ++k) {
    oldToNewIndex[newToOldIndex[k]] = k;
  }

  for (auto sh : simplexTree.complex_simplex_range()) {
    auto index = oldToNewIndex[simplexTree.key(sh)];
    dimensions[index] = simplexTree.dimension(sh);
    if constexpr (std::is_same_v<MultiFiltrationValue, typename SimplexTreeOptions::Filtration_value>) {
      filtrationValues[index] = simplexTree.filtration(sh);
    } else {
      filtrationValues[index] = as_type<MultiFiltrationValue>(simplexTree.filtration(sh));
    }
    typename Complex::Boundary boundary(dimensions[index] == 0 ? 0 : dimensions[index] + 1);
    unsigned int j = 0;
    for (auto b : simplexTree.boundary_simplex_range(sh)) {
      boundary[j] = oldToNewIndex[simplexTree.key(b)];
      ++j;
    }
    std::sort(boundary.begin(), boundary.end());
    boundaries[index] = std::move(boundary);
  }

  return Complex(std::move(boundaries), std::move(dimensions), std::move(filtrationValues));
}

/**
 * @ingroup multi_persistence
 *
 * @brief Builds a slicer for the scc format file given. Assumes that every index appearing in a boundary in the file
 * corresponds to a real line in the file (for example, the lowest dimension has always empty boundaries).
 * See @ref Slicer::write_slicer_to_scc_file "write_slicer_to_scc_file" to write a slicer into a scc format file.
 *
 * @tparam Slicer The @ref Slicer class with any valid template combination.
 * @param inFilePath Path to scc file.
 * @param isRivetCompatible Set to true if the file is written such that Rivet can read it. See TODO ref.
 * Default value: false.
 * @param isReversed Set to true if the cells in the file are written in increasing dimension order instead of
 * the standard decreasing order. Default value: false.
 * @param shiftDimensions Indicates if there is a shift in the dimension written in the file: if the value is 0, it
 * means that the smallest dimension is 0, if the value is positive, the smallest dimension is assumed to be
 * `shiftDimensions` instead of 0, and if the value is negative, the `abs(shiftDimensions)` smallest dimensions in
 * the file are ignored and the smallest remaining dimension is interpreted as 0. Default value: 0.
 */
template <class Slicer>
inline Slicer build_slicer_from_scc_file(const std::string& inFilePath,
                                         bool isRivetCompatible = false,
                                         bool isReversed = false,
                                         int shiftDimensions = 0)
{
  auto cpx = build_complex_from_scc_file<typename Slicer::Filtration_value>(
      inFilePath, isRivetCompatible, isReversed, shiftDimensions);
  return Slicer(std::move(cpx));
}

/**
 * @ingroup multi_persistence
 *
 * @brief Builds a slicer from the given bitmap. The bitmap here is a grid where each node contains a 1-critical
 * filtration value, which will be interpreted as a vertex in a cubical complex. The filtration values of the higher
 * dimensional cells are deduced by taking at each parameter the maximal value of its facets at this parameter.
 *
 * Note that for the bitmap to represent a valid multi-parameter filtration, all filtration values have to have the
 * same number of parameters. The behaviour is undefined otherwise.
 *
 * @tparam Slicer The @ref Slicer class with any valid template combination.
 * @param vertexValues Bitmap with 1-critical filtration values. Represented as a single vector, the next input
 * parameter @p shape indicates the shape of the real bitmap.
 * @param shape Shape of the bitmap. E.g., if @p shape is \f$ {3, 4} \f$, then the bitmap is a \f$ (4 x 3) \f$ grid
 * with four lines and three columns. The vector @p vertexValues should then contain 12 elements: the three first
 * elements will be read as the first line, the three next elements as the second line etc. until having 4 lines.
 */
template <class Slicer>
inline Slicer build_slicer_from_bitmap(const std::vector<typename Slicer::Filtration_value>& vertexValues,
                                       const std::vector<unsigned int>& shape)
{
  auto cpx = build_complex_from_bitmap<typename Slicer::Filtration_value>(vertexValues, shape);
  return Slicer(std::move(cpx));
}

/**
 * @ingroup multi_persistence
 *
 * @brief Builds a slicer from the given simplex tree. The inner complex will be ordered by dimension.
 *
 * @tparam Slicer The @ref Slicer class with any valid template combination.
 * @tparam SimplexTreeOptions Class following the @ref SimplexTreeOptions concept such that
 * @ref SimplexTreeOptions::Filtration_value follows the @ref MultiFiltrationValue concept.
 * @param simplexTree Simplex tree to convert.
 */
template <class Slicer, class SimplexTreeOptions>
inline Slicer build_slicer_from_simplex_tree(Simplex_tree<SimplexTreeOptions>& simplexTree)
{
  auto cpx = build_complex_from_simplex_tree<typename Slicer::Filtration_value, SimplexTreeOptions>(simplexTree);
  return Slicer(std::move(cpx));
}

/**
 * @private
 */
template <bool idx, class U, class Slicer, class F>
std::vector<typename Slicer::template Multi_dimensional_flat_barcode<U>>
persistence_on_slices_(Slicer& slicer, F&& ini_slicer, unsigned int size, [[maybe_unused]] bool ignoreInf = true)
{
  using Barcode = typename Slicer::template Multi_dimensional_flat_barcode<U>;

  if (size == 0) return {};

  std::vector<Barcode> out(size);

  if constexpr (Slicer::Persistence::is_vine) {
    std::forward<F>(ini_slicer)(slicer, 0);
    slicer.initialize_persistence_computation(false);
    out[0] = slicer.template get_flat_barcode<true, U, idx>();
    for (auto i = 1U; i < size; ++i) {
      std::forward<F>(ini_slicer)(slicer, i);
      slicer.vineyard_update();
      out[i] = slicer.template get_flat_barcode<true, U, idx>();
    }
  } else {
#ifdef GUDHI_USE_TBB
    using Index = typename Slicer::Index;
    tbb::enumerable_thread_specific<typename Slicer::Thread_safe> threadLocals(slicer.weak_copy());
    tbb::parallel_for(static_cast<Index>(0), size, [&](const Index& i) {
      typename Slicer::Thread_safe& s = threadLocals.local();
      std::forward<F>(ini_slicer)(s, i);
      s.initialize_persistence_computation(ignoreInf);
      out[i] = s.template get_flat_barcode<true, U, idx>();
    });
#else
    for (auto i = 0U; i < size; ++i) {
      std::forward<F>(ini_slicer)(slicer, i);
      slicer.initialize_persistence_computation(ignoreInf);
      out[i] = slicer.template get_flat_barcode<true, U, idx>();
    }
#endif
  }

  return out;
}

/**
 * @ingroup multi_persistence
 *
 * @brief Returns the barcodes of all the given lines. A line is represented as a pair with the first element being
 * a point on the line and the second element a vector giving the positive direction of the line. The direction
 * container can be empty: then the slope is assumed to be 1.
 *
 * @tparam Slicer Either @ref Slicer or @ref Thread_safe_slicer class with any valid template combination.
 * @tparam T Type of a coordinate element.
 * @tparam U Type of filtration values in the output barcode. Default value: T.
 * @tparam idx If true, the complex indices instead of the actual filtration values are used for the bars. It is
 * recommended to use an integer type for `U` in that case. Default value: false.
 * @param slicer Slicer from which to compute persistence.
 * @param basePoints Vector of base points for the lines. The dimension of a point has to correspond to the number
 * of parameters in the slicer.
 * @param directions Vector of directions for the lines. A direction has to have the same dimension than a point.
 * Can be empty, then the slope is assumed to be 1.
 * @param ignoreInf If true, all cells at infinity filtration values are ignored when computing, resulting
 * potentially in less storage use and better performance. But the parameter will be ignored if
 * PersistenceAlgorithm::is_vine is true.
 */
template <class Slicer, class T, class U = T, bool idx = false>
std::vector<typename Slicer::template Multi_dimensional_flat_barcode<U>> persistence_on_slices(
    Slicer& slicer,
    const std::vector<std::vector<T>>& basePoints,
    const std::vector<std::vector<T>>& directions,
    bool ignoreInf = true)
{
  GUDHI_CHECK(directions.empty() || directions.size() == basePoints.size(),
              "There should be as many directions than base points.");
  GUDHI_CHECK(basePoints.empty() || basePoints[0].size() == slicer.get_number_of_parameters(),
              "There should be as many directions than base points.");

  std::vector<T> dummy;
  auto get_direction = [&](unsigned int i) -> const std::vector<T>& {
    return directions.empty() ? dummy : directions[i];
  };

  return persistence_on_slices_<idx, U>(
      slicer,
      [&](auto& s, unsigned int i) { s.push_to(Line<T>(basePoints[i], get_direction(i))); },
      basePoints.size(),
      ignoreInf);
}

/**
 * @ingroup multi_persistence
 *
 * @brief Returns the barcodes of all the given slices.
 *
 * @tparam Slicer Either @ref Slicer or @ref Thread_safe_slicer class with any valid template combination.
 * @tparam T Type of a slice element.
 * @tparam U Type of filtration values in the output barcode. Default value: T.
 * @tparam idx If true, the complex indices instead of the actual filtration values are used for the bars. It is
 * recommended to use an integer type for `U` in that case. Default value: false.
 * @param slicer Slicer from which to compute persistence.
 * @param slices Vector of slices. A slice has to has as many elements than cells in the slicer.
 * @param ignoreInf If true, all cells at infinity filtration values are ignored when computing, resulting
 * potentially in less storage use and better performance. But the parameter will be ignored if
 * PersistenceAlgorithm::is_vine is true.
 */
template <class Slicer, class T, class U = T, bool idx = false>
std::vector<typename Slicer::template Multi_dimensional_flat_barcode<U>>
persistence_on_slices(Slicer& slicer, const std::vector<std::vector<T>>& slices, bool ignoreInf = true)
{
  GUDHI_CHECK(slices.empty() || slices[0].size() == slicer.get_number_of_cycle_generators(),
              "There should be as many elements in a slice than cells in the slicer.");

  return persistence_on_slices_<idx, U>(
      slicer, [&](auto& s, unsigned int i) { s.set_slice(slices[i]); }, slices.size(), ignoreInf);
}

// Mostly for python
/**
 * @ingroup multi_persistence
 *
 * @brief Returns the barcodes of all the given slices.
 *
 * @tparam Slicer Either @ref Slicer or @ref Thread_safe_slicer class with any valid template combination.
 * @tparam T Type of a slice element.
 * @tparam U Type of filtration values in the output barcode. Default value: T.
 * @tparam idx If true, the complex indices instead of the actual filtration values are used for the bars. It is
 * recommended to use an integer type for `U` in that case. Default value: false.
 * @param slicer Slicer from which to compute persistence.
 * @param slices Pointer to the begining of slices continuously aligned after another in the memory.
 * @param numberOfSlices Number of slices represented by the pointer.
 * @param ignoreInf If true, all cells at infinity filtration values are ignored when computing, resulting
 * potentially in less storage use and better performance. But the parameter will be ignored if
 * PersistenceAlgorithm::is_vine is true.
 */
template <class Slicer, class T, class U = T, bool idx = false, class = std::enable_if_t<std::is_arithmetic_v<T>>>
std::vector<typename Slicer::template Multi_dimensional_flat_barcode<U>>
persistence_on_slices(Slicer& slicer, T* slices, unsigned int numberOfSlices, bool ignoreInf = true)
{
  auto num_gen = slicer.get_number_of_cycle_generators();
  auto view = Gudhi::Simple_mdspan(slices, numberOfSlices, num_gen);

  return persistence_on_slices_<idx, U>(
      slicer,
      [&](auto& s, unsigned int i) {
        T* start = &view(i, 0);
        auto r = boost::iterator_range<T*>(start, start + num_gen);
        s.set_slice(r);
      },
      numberOfSlices,
      ignoreInf);
}

}  // namespace multi_persistence
}  // namespace Gudhi

#endif  // MP_SLICER_HELPERS_H_
