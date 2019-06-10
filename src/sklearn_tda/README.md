# sklearn_tda: a scikit-learn compatible python package for Machine Learning and TDA

**Author**: Mathieu Carrière.

# Description

sklearn_tda is a python package for handling collections of persistence diagrams for machine learning purposes.
Various preprocessing methods, vectorizations methods and kernels for persistence diagrams are implemented in a [scikit-learn](http://scikit-learn.org/) compatible fashion.
Clustering methods from TDA (Mapper and Graph Induced Complex) are also implemented.

### Preprocessing

Currently available classes are: 

  * **BirthPersistenceTransform**: apply the affine transformation (x,y) -> (x,y-x) to the diagrams.

    Parameters: None.

  * **DiagramPreprocessor**: apply a scaler to the diagrams (such as scalers from [scikit-learn](http://scikit-learn.org/)).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    | **use** = False | Whether to use the class or not. |
    | **scalers** = [] | List of scalers to be fit on the diagrams. Each element is a tuple whose first element is a list of coordinates and second element is a scaler (such as sklearn.preprocessing.MinMaxScaler()) for these coordinates.  |

  * **ProminentPoints**: remove points close to the diagonal.

    Parameters:

    | **name** | **description** |
    | --- | --- |
    | **use** = False|     Whether to use the class or not. |
    | **num_pts** = 10|    Cardinality threshold. |
    | **threshold** = -1|  Distance-to-diagonal threshold. |
    | **location** = "upper"|  Whether to keep the points above ("upper") or below ("lower") the previous thresholds. |
    | **point_type** = "finite"|  Specifies the input point type. Either "finite" or "essential". If "finite", the output points are ordered by persistence. |

  * **Padding**: add dummy points to each diagram so that they all have the same cardinality. All points are given an additional coordinate
    indicating if the point was added after padding (0) or already present before applying this class (1).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    | **use** = False|     Whether to use the class or not. |
    

  * **DiagramSelector**: return the finite or essential points of the diagrams.

     Parameters:

    | **name** | **description** |
    | --- | --- |
    | **limit** = np.inf | Diagram points with ordinate equal to **limit** will be considered as essential. |
    | **point_type** = "finite"| Specifies the point type to return. Either "finite" or "essential". |

### Vectorizations


Currently available classes are:

  * **Landscape**: implementation of [landscapes](http://jmlr.org/papers/v16/bubenik15a.html).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**num_landscapes** = 5| Number of landscapes.|
    |**resolution** = 100| Number of sample points of each landscape.|
    |**ls_range** = [np.nan, np.nan]| Range of each landscape. If np.nan, it is set to min and max of x-axis in the diagrams.|

  * **PersistenceImage**: implementation of [persistence images](http://jmlr.org/papers/v18/16-337.html).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**bandwidth** = 1.0 | Bandwidth of Gaussian kernel on the plane.|
    |**weight** = lambda x: 1| Weight on diagram points. It is a python function.|
    |**resolution** = [20,20]| Resolution of image.|
    |**im_range** = [np.nan, np.nan, np.nan, np.nan]| Range of coordinates. If np.nan, it is set to min and max of x- and y-axis in the diagrams.|

  * **BettiCurve**: implementation of [Betti curves](https://www.researchgate.net/publication/316604237_Time_Series_Classification_via_Topological_Data_Analysis).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**resolution** = 100| Number of sample points of Betti curve.|
    |**bc_range** = [np.nan, np.nan]| Range of Betti curve. If np.nan, it is set to min and max of x-axis in the diagrams.|

  * **Silhouette**: implementation of [silhouettes](http://jocg.org/index.php/jocg/article/view/203).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**weight** = lambda x: 1| Weight on diagram points. It is a python function.|
    |**resolution** = 100| Number of sample points of silhouette.|
    |**range** = [np.nan, np.nan]| Range of silhouette. If np.nan, it is set to min and max of x-axis in the diagrams.|

  * **TopologicalVector**: implementation of [distance vectors](https://diglib.eg.org/handle/10.1111/cgf12692).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**threshold** = 10| Number of distances to keep.|

  * **ComplexPolynomial**: implementation of [complex polynomials](https://link.springer.com/chapter/10.1007%2F978-3-319-23231-7_27).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**F** = "R"| Complex transformation to apply on the diagram points. Either "R", "S" or "T". |
    |**threshold** = 10| Number of coefficients to keep. |

### Kernels

Currently available classes are:

  * **PersistenceScaleSpaceKernel**: implementation of [Persistence Scale Space Kernel](https://www.cv-foundation.org/openaccess/content_cvpr_2015/papers/Reininghaus_A_Stable_Multi-Scale_2015_CVPR_paper.pdf).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**bandwidth** = 1.0| Bandwidth of kernel.|
    |**kernel_approx** = None| Kernel approximation method, such as [those in scikit-learn](https://scikit-learn.org/stable/modules/classes.html#module-sklearn.kernel_approximation).| 
  * **PersistenceWeightedGaussianKernel**: implementation of [Persistence Weighted Gaussian Kernel](http://proceedings.mlr.press/v48/kusano16.html).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**bandwidth** = 1.0 | Bandwidth of Gaussian kernel.|
    |**weight** = lambda x: 1| Weight on diagram points. It is a python function.|
    |**kernel_approx** = None| Kernel approximation method, such as [those in scikit-learn](https://scikit-learn.org/stable/modules/classes.html#module-sklearn.kernel_approximation).|
    |**use_pss** = False| Whether to add symmetric of points from the diagonal.|

  * **SlicedWassersteinKernel**: implementation of [Sliced Wasserstein Kernel](http://proceedings.mlr.press/v70/carriere17a.html).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**num_directions** = 10| Number of directions.|
    |**bandwidth** = 1.0| Bandwidth of kernel.|

  * **PersistenceFisherKernel**: implementation of [Persistence Fisher Kernel](papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**bandwidth_fisher** = 1.0| Bandwidth of Gaussian kernel for Fisher distance.|
    |**bandwidth** = 1.0| Bandwidth of kernel.|
    |**kernel_approx** = None| Kernel approximation method, such as [those in scikit-learn](https://scikit-learn.org/stable/modules/classes.html#module-sklearn.kernel_approximation).|

### Metrics

Currently available classes are:

  * **BottleneckDistance**: wrapper for bottleneck distance module of Gudhi. **Requires Gudhi!!**

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**epsilon** = 0.001| Approximation error.|

  * **SlicedWassersteinDistance**: implementation of [Sliced Wasserstein distance](http://proceedings.mlr.press/v70/carriere17a.html).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**num_directions** = 10| Number of directions.|

  * **PersistenceFisherDistance**: implementation of [Fisher Information distance](papers.nips.cc/paper/8205-persistence-fisher-kernel-a-riemannian-manifold-kernel-for-persistence-diagrams).

    Parameters:

    | **name** | **description** |
    | --- | --- |
    |**bandwidth** = 1.0| Bandwidth of Gaussian kernel.|
    |**kernel_approx** = None| Kernel approximation method, such as [those in scikit-learn](https://scikit-learn.org/stable/modules/classes.html#module-sklearn.kernel_approximation).|

### Clustering

Currently available classes are: 

  * **GraphInducedComplex**: wrapper for the [Graph Induced Complex](http://gudhi.gforge.inria.fr/python/latest/nerve_gic_complex_ref.html) module of Gudhi. **Requires Gudhi!!**

    Parameters:

    | **name** | **description** |
    | --- | --- |
    | **graph** = -1 | Radius for the neighborhood graph built on top of the point cloud. If -1, it is computed automatically. |
    | **graph_subsampling** = 100 | Number of bootstrap iterations. Used only if **graph** = -1. |
    | **graph_subsampling_power** = 0.001| Power for the approximation method. Used only if **graph** = -1.|
    | **graph_subsampling_constant** = 10| Constant for the approximation method. Used only if **graph** = -1.|
    | **cover_type** = "functional"| String specifying the cover. Either "functional" or "Voronoi". |
    | **filter** = 0| Filter function. Either an integer, in which case the corresponding coordinate is used, or a numpy array specifying the filter values on each node. Not used if **cover_type** = "Voronoi".|
    | **resolution** = -1| Resolution of intervals. If -1, it is computed automatically. Not used if **cover_type** = "Voronoi".|
    | **gain** = 0.33| Gain of intervals. Not used if **cover_type** = "Voronoi".|
    | **Voronoi_subsampling** = 1000| Number of Voronoi cells. Not used if **cover_type** = "functional".|
    | **mask** = 0| Threshold on the node sizes. |
    | **color** = 0| Color function. Either an integer, in which case the corresponding coordinate is used, or a numpy array specifying the color values on each node. |
    | **verbose** = False| Whether to print info or not. |
    | **input** = "point cloud"| Specifies the input type. Either "point cloud" or "distance matrix". If "distance matrix", some class methods are unavailable.|

  * **MapperComplex**: implementation of the [Mapper](https://research.math.osu.edu/tgda/mapperPBG.pdf). **Requires Gudhi!!**

    Parameters

    | **name** | **description** |
    | --- | --- |
    | **input** = "point cloud" | String specifying input type. Either "point cloud" or "distance matrix". |
    | **filters** = np.array([[0]]) | Numpy array specifying the filter values. Each row is a point and each column is a filter dimension. If only one integer is given per column, the corresponding coordinate is used as filter. |
    | **filter_bnds** = "auto" | Numpy array specifying the lower and upper limits of each filter. If "auto", they are automatically computed.  |
    | **colors** = np.array([[0]]) | Numpy array specifying the color values. Each row is a point and each column is a color dimension. If only one integer is given per column, the corresponding coordinate is used as color. |
    | **resolutions** = -1| List of resolutions for each filter dimension. If -1, they are computed automatically. |
    | **gains** = 0.3| List of gains for each filter dimension. If single number, the same gain is always used. |
    | **clustering** = sklearn.cluster.DBSCAN()| Clustering method. |
    | **mask** = 0| Threshold on the node sizes.|
    | **verbose** = False| Whether to print info or not.|
    | **beta** = 0| Power for the approximation method. Used only if **resolutions** = -1.|
    | **C** = 100| Constant for the approximation method. Used only if **resolutions** = -1.|
    | **N** = 100| Number of bootstrap iterations. Used only if **resolutions** = -1.|

# Installing sklearn_tda

The sklearn_tda library requires:

* python [>=2.7, >=3.5]
* numpy [>= 1.8.2]
* scikit-learn

For now, the package has to be compiled from source. You have to 

* download the code with:
```shell
git clone https://github.com/MathieuCarriere/sklearn_tda
```
* move to the directory:
```shell
cd sklearn_tda
```
* compile with:
```shell
(sudo) pip install .
```

The package can then be imported in a python shell with:
```shell
import sklearn_tda
``` 

Usage
=====

All modules are standard scikit-learn modules: they have fit, transform and fit_transform methods.
Hence, the most common way to use module X is to call X.fit_transform(input).
The input of all modules (except the clustering modules) are lists of persistence diagram, which are represented as lists of 2D numpy arrays.
Various examples can be found [here](example/).

