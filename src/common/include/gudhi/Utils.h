 /*    This file is part of the Gudhi Library. The Gudhi library
  *    (Geometric Understanding in Higher Dimensions) is a generic C++
  *    library for computational topology.
  *
  *    Author(s):       David Salinas
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Mediterranee (France)
  *
  *    This program is free software: you can redistribute it and/or modify
  *    it under the terms of the GNU General Public License as published by
  *    the Free Software Foundation, either version 3 of the License, or
  *    (at your option) any later version.
  *
  *    This program is distributed in the hope that it will be useful,
  *    but WITHOUT ANY WARRANTY; without even the implied warranty of
  *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *    GNU General Public License for more details.
  *
  *    You should have received a copy of the GNU General Public License
  *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  */
#ifndef GUDHI_UTILS_H_
#define GUDHI_UTILS_H_


#define PRINT(a) std::cerr << #a << ": " << (a) << " (DISP)"<<std::endl

#define DBG_VERBOSE
#ifdef DBG_VERBOSE
#define DBG(a) std::cerr << "DBG: " << (a)<<std::endl
#define DBGMSG(a,b) std::cerr << "DBG: " << a<<b<<std::endl
#define DBGVALUE(a) std::cerr << "DBG: " <<  #a << ": " << a<<std::endl
#define DBGCONT(a) std::cerr << "DBG: container "<< #a<<" -> "; for(auto x:a) std::cerr<< x << ","; std::cerr<<std::endl
#else
//#define DBG(a) a
//#define DBGMSG(a,b) b
//#define DBGVALUE(a) a
//#define DBGCONT(a) a
#define DBG(a)
#define DBGMSG(a,b)
#define DBGVALUE(a)
#define DBGCONT(a)
#endif




#endif /* UTILS_H_ */
