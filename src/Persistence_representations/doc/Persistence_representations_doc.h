/*    This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
 *    See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2016 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef DOC_GUDHI_STAT_H_
#define DOC_GUDHI_STAT_H_

namespace Gudhi {

namespace Persistence_representations {

/**  \defgroup Persistence_representations Persistence representations
 *
 * \author   Pawel Dlotko
 *
 * @{
 *\section Persistence_representations_idea Idea

 In order to perform most of the statistical tests and machine learning algorithms on a data one need to be able to
 perform only a very limited number of operations on them. Let us fix a representation of
 data of a type A. To perform most of the statistical and machine learning operations one need to be able to compute
 average of objects of type A (so that the averaged object is also of a type A), to
 compute distance between objects of a type A, to vectorize object of a type A and to compute scalar product of a pair
 objects of a type A.

 To put this statement into a context, let us assume we have two collections \f$ c_1,\ldots,c_n\f$ and
 \f$d_1,...,d_n\f$ of objects of a type A. We want to verify if the average of those two collections
 are different by performing a permutation test.
 First of all, we compute averages of those two collections: C average of \f$ c_1,\ldots,c_n \f$ and D average of
 \f$d_1,\ldots,d_n\f$. Note that both C and D are of a type A. Then we compute \f$d(C,D)\f$,
 a distance between C and D.
 Later we put the two collections into one bin:
 \f[B = \{ c_1,...,c_n,d_1,...,d_n \}\f]
 Then we shuffle B, and we divide the shuffled version of B into two classes: \f$B_1\f$ and \f$B_2\f$ (in this case, of
 the same cardinality). Then we compute averages \f$\hat{B_1}\f$ and \f$\hat{B_2}\f$
 of elements in \f$B_1\f$ and \f$B_2\f$. Note that again, \f$\hat{B_1}\f$ and \f$\hat{B_2}\f$ are of a type A.
 Then we compute their distance \f$d(\hat{B_1},\hat{B_2})\f$. The procedure of shuffling and dividing the set \f$B\f$
 is repeated \f$N\f$ times (where \f$N\f$ is reasonably large number).
 Then the p-value of a statement that the averages of  \f$c_1,...,c_n\f$ and \f$d_1,...,d_n\f$ is approximated by the
 number of times \f$d(\hat{B_1},\hat{B_2}) > d(C,D)\f$ divided by \f$N\f$.

 The permutation test reminded above can be performed for any type A which can be averaged, and which allows for
 computations of distances.

 The Persistence\_representations contains a collection of various representations of persistent homology that
 implements various concepts described below:

 \li Concept of a representation of persistence that allows averaging (so that the average object is of the same type).
 \li Concept of representation of persistence that allows computations of distances.
 \li Concept of representation of persistence that allows computations of scalar products.
 \li Concept of representation of persistence that allows vectorization.
 \li Concept of representation of persistence that allows computations of real-valued characteristics of objects.


 At the moment an implementation of the following representations of persistence are available (further details of
 those representations will be discussed later):

 \li Exact persistence landscapes (allow averaging, computation of distances, scalar products, vectorizations and real
 value characteristics).
 \li Persistence landscapes on a grid (allow averaging, computation of distances scalar products, vectorizations and
 real value characteristics).
 \li Persistence heat maps – various representations where one put some weighted or not Gaussian kernel for each point
 of diagram (allow averaging, computation of distances, scalar products,
 vectorizations and real value characteristics).
 \li Persistence vectors (allow averaging, computation of distances, scalar products, vectorizations and real value
 characteristics).
 \li Persistence diagrams / barcodes (allow computation of distances, vectorizations and real value characteristics).


 Note that at the while functionalities like averaging, distances and scalar products are fixed, there is no canonical
 way of vectorizing and computing real valued characteristics of objects. Therefore the
 vectorizations and computation of real value characteristics procedures are quite likely to evolve in the furthering
 versions of the library.

 The main aim of this implementation is to be able to implement various statistical methods, both on the level of C++
 and on the level of python. The methods will operate on the functionalities offered
 by concepts. That means that the statistical and ML methods will be able to operate on any representation that
 implement the required concept (including the ones that are not in the library at the moment).
 That gives provides a framework, that is very easy to extend, for topological statistics.

 Below we are discussing the representations which are currently implemented in Persistence\_representations package:

 \section sec_persistence_landscapes Persistence Landscapes
 <b>Reference manual:</b> \ref Gudhi::Persistence_representations::Persistence_landscape <br>
 Persistence landscapes were originally proposed by Bubenik in \cite bubenik_landscapes_2015. Efficient algorithms to
 compute them rigorously were proposed by Bubenik and Dlotko in \cite bubenik_dlotko_landscapes_2016. The idea of
 persistence landscapes is shortly summarized in below.

 To begin with, suppose we are given a point \f$(b,d) \in \mathbb{R}^2\f$ in a
 persistence diagram. With this point, we associate a piecewise
 linear function \f$f_{(b,d)} : \mathbb{R} \rightarrow [0,\infty)\f$, which is
 defined as

  \f[f_{(b,d)}(x) =
  \left\{ \begin{array}{ccl}
            0     & \mbox{ if } & x \not\in (b, d) \; , \\
            x - b & \mbox{ if } & x \in \left( b, \frac{b+d}{2}
              \right] \; , \\
            d - x & \mbox{ if } & x \in \left(\frac{b+d}{2},
              d \right) \; .
  \end{array} \right.
 \f]

 A persistence landscape of the birth-death
 pairs \f$(b_i , d_i)\f$, where \f$i = 1,\ldots,m\f$, which constitute the given
 persistence diagram is the sequence of functions \f$\lambda_k : \mathbb{R} \rightarrow [0,\infty)\f$ for \f$k \in
 \mathbb{N}\f$, where \f$\lambda_k(x)\f$
 denotes the \f$k^{\rm th}\f$ largest value of the numbers \f$f_{(b_i,d_i)}(x)\f$,
 for \f$i = 1, \ldots, m\f$, and we define \f$\lambda_k(x) = 0\f$ if \f$k > m\f$.
 Equivalently, this sequence of functions can be combined into a single
 function \f$L : \mathbb{N} \times \mathbb{R} \to [0,\infty)\f$ of two
 variables, if we define \f$L(k,t) = \lambda_k(t)\f$.

 The detailed description of algorithms used to compute persistence landscapes can be found in
 \cite bubenik_dlotko_landscapes_2016.
 Note that this implementation provides exact representation of landscapes. That have many advantages, but also a few
 drawbacks. For instance, as discussed
 in \cite bubenik_dlotko_landscapes_2016, the exact representation of landscape may be of quadratic size with respect
 to the input persistence diagram. It may therefore happen
 that, for very large diagrams, using this representation may be memory--prohibitive. In such a case, there are two
 possible ways to proceed:

 \li Use non exact representation on a grid described in the Section \ref sec_landscapes_on_grid.
 \li Compute just a number of initial nonzero landscapes. This option is available from C++ level as a last parameter of
 the constructor of persistence landscape (set by default to std::numeric_limits<size_t>::max()).



 \section sec_landscapes_on_grid Persistence Landscapes on a grid
 <b>Reference manual:</b> \ref Gudhi::Persistence_representations::Persistence_landscape_on_grid <br>
 This is an alternative, not--exact, representation of persistence landscapes defined in the Section \ref
 sec_persistence_landscapes. Unlike in the Section \ref sec_persistence_landscapes we build a
 representation of persistence landscape by sampling its values on a finite, equally distributed grid of points.
 Since, the persistence landscapes that originate from persistence diagrams have slope \f$1\f$ or \f$-1\f$, we have an
 estimate of a region between the grid points where the landscape cab be located.
 That allows to estimate an error make when performing various operations on landscape. Note that for average
 landscapes the slope is in range \f$[-1,1]\f$ and similar estimate can be used.

 Due to a lack of rigorous description of the algorithms to deal with this non--rigorous representation of persistence
 landscapes in the literature, we are providing a short discussion of them in below.

 Let us assume that we want to compute persistence landscape on a interval \f$[x,y]\f$. Let us assume that we want to
 use \f$N\f$ grid points for that purpose.
 Then we will sample the persistence landscape on points \f$x_1 = x , x_2 = x + \frac{y-x}{N}, \ldots , x_{N} = y\f$.
 Persistence landscapes are represented as a vector of
 vectors of real numbers. Assume that i-th vector consist of \f$n_i\f$ numbers sorted from larger to smaller. They
 represent the values of the functions
 \f$\lambda_1,\ldots,\lambda_{n_i}\f$ ,\f$\lambda_{n_i+1}\f$ and the functions with larger indices are then zero
 functions) on the i-th point of a grid, i.e. \f$x + i \frac{y-x}{N}\f$.

 When averaging two persistence landscapes represented by a grid we need to make sure that they are defined in a
 compatible grids. I.e. the intervals \f$[x,y]\f$ on which they are defined are
 the same, and the numbers of grid points \f$N\f$ are the same in both cases. If this is the case, we simply compute
 point-wise averages of the entries of corresponding
 vectors (In this whole section we assume that if one vector of numbers is shorter than another, we extend the shorter
 one with zeros so that they have the same length.)

 Computations of distances between two persistence landscapes on a grid is not much different than in the rigorous
 case. In this case, we sum up the distances between the same levels of
 corresponding landscapes. For fixed level, we approximate the landscapes between the corresponding constitutive
 points of landscapes by linear functions, and compute the \f$L^p\f$ distance between them.

 Similarly as in case of distance, when computing the scalar product of two persistence landscapes on a grid, we sum up
 the scalar products of corresponding levels of landscapes. For each level,
 we assume that the persistence landscape on a grid between two grid points is approximated by linear function.
 Therefore to compute scalar product of two corresponding levels of landscapes,
 we sum up the integrals of products of line segments for every pair of constitutive grid points.

 Note that for this representation we need to specify a few parameters:

 \li Begin and end point of a grid -- the interval \f$[x,y]\f$ (real numbers).
 \li Number of points in a grid (positive integer \f$N\f$).


 Note that the same representation is used in TDA R-package \cite Fasy_Kim_Lecci_Maria_tda.

 \section sec_persistence_heat_maps Persistence heat maps
 <b>Reference manual:</b> \ref Gudhi::Persistence_representations::Persistence_heat_maps <br>
 This is a general class of discrete structures which are based on idea of placing a kernel in the points of
 persistence diagrams.
 This idea appeared in work by many authors over the last 15 years. As far as we know this idea was firstly described
 in the work of Bologna group in \cite Ferri_Frosini_comparision_sheme_1 and \cite Ferri_Frosini_comparision_sheme_2.
 Later it has been described by Colorado State University group in \cite Persistence_Images_2017. The presented paper
 in the first time provide a discussion of stability of the representation.
 Also, the same ideas are used in construction of two recent kernels used for machine learning:
 \cite Kusano_Fukumizu_Hiraoka_PWGK and \cite Reininghaus_Huber_ALL_PSSK. Both the kernel's construction uses
 interesting ideas to ensure stability of the representation with respect to Wasserstein metric. In the kernel
 presented in \cite Kusano_Fukumizu_Hiraoka_PWGK, a scaling function is used to multiply the Gaussian kernel in the
 way that the points close to diagonal got low weight and consequently do not have a big influence on the resulting
 distribution. In \cite Reininghaus_Huber_ALL_PSSK for every point \f$(b,d)\f$ two Gaussian kernels
 are added: first, with a weight 1 in a point \f$(b,d)\f$, and the second, with the weight -1 for a point \f$(b,d)\f$.
 In both cases, the representations are stable with respect to 1-Wasserstein distance.

 In Persistence\_representations package we currently implement a discretization of the distributions described above.
 The base of this implementation is 2-dimensional array of pixels. Each pixel have assigned a real value which
 is a sum of values of distributions induced by each point of the persistence diagram. At the moment we compute the
 sum of values on a center of a pixels. It can be easily extended to any other function
 (like for instance sum of integrals of the intermediate distribution on a pixel).

 The parameters that determine the structure are the following:

 \li A positive integer k determining the size of the kernel we used (we always assume that the kernels are square).
 \li A filter: in practice a square matrix of a size \f$2k+1 \times 2k+1\f$. By default, this is a discretization of
 N(0,1) kernel.
 \li The box \f$[x_0,x_1]\times [y_0,y_1]\f$ bounding the domain of the persistence image.
 \li Scaling function. Each Gaussian kernel at point \f$(p,q)\f$ gets multiplied by the value of this function at the
 point \f$(p,q)\f$.
 \li A boolean value determining if the space below diagonal should be erased or not. To be precise: when points close
 to diagonal are given then sometimes the kernel have support that reaches the region
 below the diagonal. If the value of this parameter is true, then the values below diagonal can be erased.

 In addition to the previous method, we also provide two more methods to perform exact calculations, in the sense that we use functions
 instead of matrices to define the kernel between the points of the diagrams.
 Indeed, in both of these exact methods, the kernel is no longer provided as a square matrix, or a filter (see parameters above), but rather as
 a function assigning a real value to a pair of points in the plane.

 In the first of these exact methods, we aim at obtaining a finite-dimensional representation of the diagram, so we still use a grid of pixels.
 On the other hand, in the second exact method, we represent diagrams implicitly as functions (i.e. infinite-dimensional representations). This way, we no longer require grids,
 but only scalar products and distances are available with these implicit representations. This type of representations is known as
 kernel methods (see \ref sec_persistence_kernels below for more details on kernels).

 Names can be a bit confusing so we recall that, with this second exact method, we implicitly define a kernel representation of diagrams that is built from a kernel between points
 in the plane. Hence, we have two kernels here, which are independent. One is defined between points in the plane (its type in the code is Kernel2D), and is a template parameter,
 whereas the other is defined between persistence diagrams (it is the scalar product of the infinite-dimensional representations of the diagrams).

 \section sec_persistence_vectors Persistence vectors
 <b>Reference manual:</b> \ref Gudhi::Persistence_representations::Vector_distances_in_diagram <br>
 This is a representation of persistent homology in a form of a vector which was designed for an application in 3d
 graphic in \cite Carriere_Oudot_Ovsjanikov_top_signatures_3d. Below we provide a short description of this
 representation.

 Given a persistence diagram \f$D = \{ (b_i,d_i) \}\f$, for every pair of birth--death points \f$(b_1,d_1)\f$ and
 \f$(b_2,d_2)\f$ we compute the following three distances:

 \li \f$d( (b_1,d_1) , (b_2,d_2) )\f$.
 \li \f$d( (b_1,d_1) , (\frac{b_1,d_1}{2},\frac{b_1,d_1}{2}) )\f$.
 \li \f$d( (b_2,d_2) , (\frac{b_2,d_2}{2},\frac{b_2,d_2}{2}) )\f$.

 We pick the smallest of those and add it to a vector. The obtained vector of numbers is then sorted in decreasing
 order. This way we obtain a persistence vector representing the diagram.

 Given two persistence vectors, the computation of distances, averages and scalar products is straightforward. Average
 is simply a coordinate-wise average of a collection of vectors. In this section we
 assume that the vectors are extended by zeros if they are of a different size. To compute distances we compute
 absolute value of differences between coordinates. A scalar product is a sum of products of
 values at the corresponding positions of two vectors.

  \section sec_persistence_kernels Kernels on persistence diagrams
 <b>Reference manual:</b> \ref Gudhi::Persistence_representations::Sliced_Wasserstein <br>
 <b>Reference manual:</b> \ref Gudhi::Persistence_representations::Persistence_heat_maps <br>

 Kernels for persistence diagrams can be regarded as infinite-dimensional vectorizations. More specifically,
 they are similarity functions whose evaluations on pairs of persistence diagrams equals the scalar products
 between images of these pairs under a map \f$\Phi\f$ taking values in a specific (possibly non Euclidean) Hilbert space \f$k(D_i, D_j) = \langle \Phi(D_i),\Phi(D_j)\rangle\f$.
 Reciprocally, classical results of learning theory ensure that such a \f$\Phi\f$ exists for a given similarity function \f$k\f$ if and only if \f$k\f$ is <i>positive semi-definite</i>.
 Kernels are designed for algorithms that can be <i>kernelized</i>, i.e., algorithms that only require to know scalar products between instances in order to run.
 Examples of such algorithms include Support Vector Machines, Principal Component Analysis and Ridge Regression.

 There have been several attempts at defining kernels, i.e., positive semi-definite functions, between persistence diagrams within the last few years. We provide implementation
 for the <i>Sliced Wasserstein kernel</i>---see \cite pmlr-v70-carriere17a, which takes the form of a Gaussian kernel with a specific distance between persistence diagrams
 called the <i>Sliced Wasserstein distance</i>: \f$k(D_1,D_2)={\rm exp}\left(-\frac{SW(D_1,D_2)}{2\sigma^2}\right)\f$. Other kernels such as the Persistence Weighted Gaussian kernel or
 the Persistence Scale Space kernel are implemented in Persistence_heat_maps.

 When launching:

 \code $> ./Sliced_Wasserstein
 \endcode

 the program output is:

 \code $> Approx SW kernel: 0.0693743
 $> Exact  SW kernel: 0.0693218
 $> Distance induced by approx SW kernel: 1.36428
 $> Distance induced by exact  SW kernel: 1.3643
 \endcode

 */
/** @} */  // end defgroup Persistence_representations

}  // namespace Persistence_representations
}  // namespace Gudhi

#endif  // Persistence_representations
