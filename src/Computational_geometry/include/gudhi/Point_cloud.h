/*    This file is a prototype for the Gudhi Library.
 *    Author(s):       Clément Maria
 *    Copyright (C) 2021 Inria
 *    This version is under developement, please do not redistribute this software. 
 *    This program is for academic research use only. 
 */

#ifndef POINT_CLOUD_H
#define POINT_CLOUD_H

#include <random>

/** \brief Counts the number of duplicates in a range, using the canonical \f$<\f$ 
 * operator.
 * 
 * ComparableObjectIterator is a LegacyForwardIterator of value_type an object that 
 * is copiable, movable and comparable.
 * 
 * @param[in] beg, end Begin and end pointers on the range.
 */
template< typename ComparableObjectIterator >
size_t duplicates(ComparableObjectIterator beg, ComparableObjectIterator end)  
{
  std::vector rg_copy(beg,end);
  std::sort(rg_copy.begin(), rg_copy.end());
  auto last = std::unique(rg_copy.begin(), rg_copy.end());
  size_t n_dups = rg_copy.end() - last;
  return n_dups;
}
template< typename ComparableObjectRange >
size_t duplicates(ComparableObjectRange &rg) 
{ return duplicates(rg.begin(),rg.end()); }

/** \brief Removes all duplicates in a range, using the canonical \f$<\f$ 
 * operator.
 * 
 * ComparableObjectIterator is a LegacyForwardIterator of value_type an object that 
 * is copiable, movable and comparable.
 * 
 * @param[in] beg, end Begin and end pointers on the range.
 */
template< typename ObjectRange >
size_t remove_duplicates(ObjectRange &rg)  
{
  std::sort(rg.begin(), rg.end());
  auto last = std::unique(rg.begin(), rg.end());
  size_t n_dups = rg.end() - last; 
  rg.erase(last, rg.end());
  return n_dups;
}

/** Computes the maximal and minimal distance between two distinct points in the 
 * range.
 * 
 * Uses a naive quadratic approach.
 * 
 * Scalar is the float type representing distance values.
 * 
 * @param[in] beg, end, iterators repersenting the range of objects.
 * @param[in] dist, distance function taking two objects and returning a Scalar 
 *                  encoding the distance between objects. 
 * 
 * @param[out] a pair of Scalar giving the respectively the maximal and 
 *                      minimal distances between two distinct points.
 */
template< typename Scalar, typename ObjectIterator, typename DistanceFunction >
std::pair<Scalar,Scalar> max_min_distances( ObjectIterator beg, ObjectIterator end
                        , DistanceFunction dist) 
{
  size_t n = end-beg;
  std::vector<Scalar> distances; distances.reserve( n*(n-1)/2 + 1);
  for(auto it1 = beg; it1 != end; ++it1) {
    for(auto it2 = it1+1; it2 != end; ++it2) {
      distances.push_back(dist(*it1, *it2));
    }
  }
  return std::make_pair(*(std::max_element(distances.begin(), distances.end())), *(std::min_element(distances.begin(), distances.end()))); 
}
template< typename Scalar, typename ObjectRange, typename DistanceFunction >
std::pair<Scalar,Scalar> max_min_distances( ObjectRange &rg, DistanceFunction dist) 
{ return max_min_distances<Scalar>(rg.begin(), rg.end(), dist); }

/** \brief Checks whether the pairwise distances between any two distinct points are 
 * distinct.
 * 
 * Specifically, for an input point cloud \f$P\f$ and distance function \f$\mathoperator{d}(\cdot,\cdot)\f$, returns true iff, for any four 
 * points \f$v, w, x, y \in P\f$, such that \f$ v \neq w, x \neq y, (v,w) \neq 
 * (x,y)\f$ (as unordered pairs), one has \f$\mathoperator{d}(v,w) \neq \mathoperator{d}(x,y)\f$.
 * 
 * This is notably useful when garanteeing a unique choice of farthest point ordering of the point cloud, once a seed is chosen.
 * 
 * Scalar is the float type representing distance values.
 * 
 * @param[in] beg, end, iterators repersenting the range of objects.
 * @param[in] dist, distance function taking two objects and returning a Scalar 
 *                   encoding the distance between objects. 
 * 
 * @param[out] a Scalar giving the minimal difference 
 *                      \f$|\mathoperator{d}(v,w) - \mathoperator{d}(x,y)|\f$ 
 *                      between the distance of two pairs of distinct points as 
 *                      above. The point cloud is generic for distance iff this 
 *                      quantity is strictly positive.
 */
template< typename Scalar, typename ObjectIterator, typename DistanceFunction >
Scalar generic_distances( ObjectIterator beg, ObjectIterator end
                        , DistanceFunction dist) 
{
  size_t n = end-beg;
  std::vector<Scalar> distances; distances.reserve( n*(n-1)/2 + 1);
  for(auto it1 = beg; it1 != end; ++it1) {
    for(auto it2 = it1+1; it2 != end; ++it2) {
      distances.push_back(dist(*it1, *it2));
    }
  }
  std::sort(distances.begin(),distances.end());
  Scalar min_delta = std::numeric_limits<Scalar>::infinity();
  for(int i=0; i<distances.size()-1; ++i) {
    if(min_delta > (distances[i+1]-distances[i]) ) { 
      min_delta = (distances[i+1]-distances[i]); 
    }
  }
  return min_delta;
}
template< typename Scalar, typename ObjectRange, typename DistanceFunction >
Scalar generic_distances( ObjectRange &rg, DistanceFunction dist) 
{ return generic_distances<Scalar>(rg.begin(), rg.end(), dist); }

/** \brief Picks uniformly at random a point on the sphere, using Muller's method.
  *  
  * Specifically, the random vector is picked uniformly at random on the 
  * \f$(d-1)\f$-dimensional sphere of radius \f$r\f$, centered on the origin in 
  * \f$\mathbb{R}^d\f$, using Muller's method. 
  * 
  * Scalar is the float type representing coordiantes, RandomGenerator is the random generator as described in the random std.
  * 
  * @param[in] rand_vector, vector in which the coordinates of the output random vector (of the appropriate dimension) are written.
  * @param[in] dim, the dimension of the output rendom vector.
  * @param[in] norm_vector, the Euclidean norm of the output vector.
  * @param[in] distrib, the normal distribution used by Muller's method.
  * @param[in] rand_gen, the random generator.
  * 
  */
template<typename Scalar, typename RandomGenerator>
void uniform_random_on_sphere(  std::vector<Scalar> &rand_vector
                              , int dim
                              , Scalar norm_vector
                              , std::normal_distribution<Scalar> &distrib
                              , RandomGenerator &rand_gen)
{
  rand_vector.clear();
  Scalar norm = 0.;
  for(int i = 0; i < dim; ++i) {
    Scalar rand_num = distrib(rand_gen);
    rand_vector.push_back(rand_num);
    norm += rand_num * rand_num;
  }
  norm = std::sqrt(norm);
  for(int i=0;i < dim; ++i) {
    rand_vector[i] /= norm; rand_vector[i] *= norm_vector;
  }
}

/** \brief Perturbs the input point cloud by adding to each point a vector picked 
  * uniformly at random on a sphere.
  *
  * Scalar is a float type for the coordinates of the points.
  * 
  * @param[in] points, the vector of points to perturb.
  * @param[in] dim, the dimension of the points.
  * @param[in] norm_perturb, the norm of the perturbation vectors.
  */
template<typename Scalar, typename PointRange>
void random_perturbation(PointRange &points, int dim, Scalar norm_perturb) 
{
  std::random_device rand_dev{};
  std::mt19937 rand_gen{rand_dev()};//Mersenne twister 32-bits
  //mean = 0, standard deviation = 2
  std::normal_distribution<Scalar> distrib(0, 2);
  std::vector<Scalar> rand_vector; rand_vector.reserve(dim);
  for(auto &p : points) {
    uniform_random_on_sphere(rand_vector, dim, norm_perturb, distrib, rand_gen);
    int i=0;
    for(auto &coord : p) { coord += rand_vector[i]; ++i; } 
  }
}

#endif //POINT_CLOUD_H
