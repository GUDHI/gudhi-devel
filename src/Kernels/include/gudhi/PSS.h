/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Mathieu Carrière
 *
 *    Copyright (C) 2017  INRIA (France)
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

#ifndef PSS_H_
#define PSS_H_

#define NUMPI 3.14159265359

#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip>

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <limits>
#include <assert.h>

#include <cmath>
#include <math.h>

#include <memory>
#include <stdexcept>
#include <omp.h>

#include <gmp.h>
#include <gmpxx.h>
#include <random>
#include <chrono>
#include <ctime>

#include "figtree.h"
#include "ANN.h"

using PD = std::vector<std::pair<double,double> >;

namespace Gudhi {
namespace persistence_scale_space {

double compute_exact_pss(PD PD1, PD PD2, double sigma = 1){
  double k = 0;
  for(int i = 0; i < PD1.size(); i++){
    for(int j = 0; j < PD2.size(); j++){
      k += exp(  -(    pow(PD1[i].first - PD2[j].first,  2)  +  pow(PD1[i].second - PD2[j].second, 2)     )/(8*sigma))   -\
           exp(  -(    pow(PD1[i].first - PD2[j].second, 2)  +  pow(PD1[i].second - PD2[j].first,  2)     )/(8*sigma));
    }
  }
  return k/(8*NUMPI*sigma);
}

double compute_approximate_pss(PD PD1, PD PD2, double sigma = 1, double error = 1e-2){

  double k = 0;

  int d = 2; int N = PD1.size(); int M = PD2.size(); double h = std::sqrt(8*sigma);
  double* x = new double[2*N]; double* y = new double[2*M]; double* q = new double[N];
  for(int i = 0; i < N; i++){
    q[i] = 1.0/(8*NUMPI*sigma);
    x[2*i] = PD1[i].first; x[2*i+1] = PD1[i].second;
  }
  for(int i = 0; i < M; i++){  y[2*i] = PD2[i].first; y[2*i+1] = PD2[i].second;  }
  double* g_auto = new double[M];
  memset(g_auto, 0, sizeof(double)*M);

  figtree(d, N, M, 1, x, h, q, y, error, g_auto);
  for(int i = 0; i < M; i++)  k += g_auto[i];

  for(int i = 0; i < M; i++){  y[2*i] = PD2[i].second; y[2*i+1] = PD2[i].first;  }

  figtree(d, N, M, 1, x, h, q, y, error, g_auto);
  for(int i = 0; i < M; i++)  k -= g_auto[i];

  delete[] x; delete[] y; delete[] q; delete[] g_auto;
  return k;
}

} // namespace persistence_scale_space

} // namespace Gudhi

#endif // PSS_H_
