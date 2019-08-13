/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       David Salinas
 *
 *    Copyright (C) 2014 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef UTILS_UI_UTILS_H_
#define UTILS_UI_UTILS_H_

#define UIDBG_VERBOSE

#ifdef UIDBG_VERBOSE
#define UIDBG(a) std::cerr << "UIDBG: " << (a) << std::endl
#define UIDBGMSG(a, b) std::cerr << "UIDBG: " << a << b << std::endl
#define UIDBGVALUE(a) std::cerr << "UIDBG: " <<  #a << ": " << a << std::endl
#define UIDBGCONT(a) std::cerr << "UIDBG: container " << #a << " -> "; for (auto x : a) std::cerr << x << ","; std::cerr << std::endl }
#else
// #define DBG(a) a
// #define DBGMSG(a, b) b
// #define DBGVALUE(a) a
// #define DBGCONT(a) a
#define UIDBG(a)
#define UIDBGMSG(a, b)
#define UIDBGVALUE(a)
#define UIDBGCONT(a)
#endif

#endif  // UTILS_UI_UTILS_H_
