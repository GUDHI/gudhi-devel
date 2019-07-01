/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Vincent Rouvreau
 *
 *    Copyright (C) 2018 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef ALPHA_COMPLEX_OPTIONS_H_
#define ALPHA_COMPLEX_OPTIONS_H_

namespace Gudhi {

namespace alpha_complex {

/**
 * \brief Alpha complex complexity template parameter possible values.
 *
 * \ingroup alpha_complex
 */
enum class complexity : char {
  FAST = 'f',   ///< Fast version.
  SAFE = 's',   ///< Safe version.
  EXACT = 'e',  ///< Exact version.
};

}  // namespace alpha_complex

}  // namespace Gudhi

#endif  // ALPHA_COMPLEX_OPTIONS_H_
