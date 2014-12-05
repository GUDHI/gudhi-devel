/*
 * Utils.h
 *
 *  Created on: 13 juin 2013
 *      Author: salinasd
 */

#ifndef GUDHI_UTILS_H_
#define GUDHI_UTILS_H_


#define PRINT(a) std::cerr << #a << ": " << (a) << " (DISP)"<<std::endl

//#define DBG_VERBOSE
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
