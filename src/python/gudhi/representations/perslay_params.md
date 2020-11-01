In the following description of PersLay parameters, each parameter, or dictionary key, that contains `_init` in its name is optimized and learned by PersLay during training. If you do not want to optimize the vectorization, set the keys **train_vect** and **train_weight** to False.

  * The following keys are mandatory:

    | **name** | **description** |
    | --- | --- |
    | **layer**              | Either "PermutationEquivariant", "Image", "Landscape", "BettiCurve", "Entropy", "Exponential", "Rational" or "RationalHat". Type of the PersLay layer. "Image" is for [persistence images](https://arxiv.org/abs/1507.06217), "Landscape" is for [persistence landscapes](http://www.jmlr.org/papers/volume16/bubenik15a/bubenik15a.pdf), "Exponential", "Rational" and "RationalHat" are for [structure elements](http://jmlr.org/beta/papers/v20/18-358.html), "PermutationEquivariant" is for the original DeepSet layer, defined in [this article](https://arxiv.org/abs/1703.06114), "BettiCurve" is for [Betti curves](https://www.jstage.jst.go.jp/article/tjsai/32/3/32_D-G72/_pdf) and "Entropy" is for [entropy](https://arxiv.org/abs/1803.08304). |
    | **perm_op**            | Either "sum", "mean", "max", "topk". Permutation invariant operation. |
    | **keep**               | Number of top values to keep. Used only if **perm_op** is "topk". |
    | **pweight**            | Either "power", "grid", "gmix" or None. Weight function to be applied on persistence diagram points. If "power", this function is a (trainable) coefficient times the distances to the diagonal of the points to a certain power. If "grid", this function is piecewise-constant and defined with pixel values of a grid. If "gmix", this function is defined as a mixture of Gaussians. If None, no weighting is applied. |
    | **final_model**        | A Tensorflow / Keras model used to postprocess the persistence diagrams in each channel. Use "identity" if you don't want to postprocess. |

Depending on what **pweight** is, the following additional keys are requested:

  * if **pweight** is "power":

    | **name** | **description** |
    | --- | --- |
    | **pweight_init**         | Initializer of the coefficient of the power weight function. It can be either a single value, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). | 
    | **pweight_power**        | Integer used for exponentiating the distances to the diagonal of the persistence diagram points. |
 
  * if **pweight** is "grid":

    | **name** | **description** |
    | --- | --- |
    | **pweight_size**          | Grid size of the grid weight function. It is a tuple of integer values, such as (10,10). | 
    | **pweight_bnds**          | Grid boundaries of the grid weight function. It is a tuple containing two tuples, each containing the minimum and maximum values of each axis of the plane. Example: ((-0.01, 1.01), (-0.01, 1.01)). | 
    | **pweight_init**          | Initializer for the pixel values of the grid weight function. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.).|

  * if **pweight** is "gmix":

    | **name** | **description** |
    | --- | --- |
    | **pweight_num**           | Number of Gaussian functions of the mixture of Gaussians weight function. |
    | **pweight_init**          | Initializer of the means and variances of the mixture of Gaussians weight function. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |

Depending on what **layer** is, the following additional keys are requested:

  * if **layer** is "PermutationEquivariant":

    | **name** | **description** |
    | --- | --- |
    | **lpeq**                | Sequence of permutation equivariant operations, as defined in [the DeepSet article](). It is a list of tuples of the form (*dim*, *operation*). Each tuple defines a permutation equivariant function of dimension *dim* and second permutation operation *operation* (string, either "max", "min", "sum" or None). Second permutation operation is optional and is not applied if *operation* is set to None. Example: [(150, "max"), (75, None)]. |
    | **lweight_init**        | Initializer for the weight matrices of the permutation equivariant operations. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.).| 
    | **lbias_init**          | Initializer for the biases of the permutation equivariant operations. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **lgamma_init**         | Initializer for the Gamma matrices of the permutation equivariant operations. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.).|

  * if **layer** is "Image":

    | **name** | **description** |
    | --- | --- |
    | **image_size**         | Persistence image size. It is a tuple of integer values, such as (10,10). | 
    | **image_bnds**         | Persistence image boundaries. It is a tuple containing two tuples, each containing the minimum and maximum values of each axis of the plane. Example: ((-0.01, 1.01), (-0.01, 1.01)). |
    | **lvariance_init**      | Initializer for the bandwidths of the Gaussian functions centered on the persistence image pixels. It can be either a single value, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 3.). | 
    
  * if **layer** is "Landscape":

    | **name** | **description** |
    | --- | --- |
    | **lsample_num**        | Number of samples of the diagonal that will be evaluated on the persistence landscapes. |
    | **lsample_init**        | Initializer of the samples of the diagonal. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |

  * if **layer** is "BettiCurve":

    | **name** | **description** |
    | --- | --- |
    | **lsample_num**        | Number of samples of the diagonal that will be evaluated on the Betti curves. |
    | **lsample_init**        | Initializer of the samples of the diagonal. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **theta**              | Sigmoid parameter used for approximating the piecewise constant functions associated to the persistence diagram points. |
    
  * if **layer** is "Entropy":

    | **name** | **description** |
    | --- | --- |
    | **lsample_num**        | Number of samples on the diagonal that will be evaluated on the persistence entropies. |
    | **lsample_init**        | Initializer of the samples of the diagonal. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **theta**              | Sigmoid parameter used for approximating the piecewise constant functions associated to the persistence diagram points. |
    
  * if **layer** is "Exponential":

    | **name** | **description** |
    | --- | --- |
    | **lnum**       | Number of exponential structure elements that will be evaluated on the persistence diagram points. |
    | **lmean_init**          | Initializer of the means of the exponential structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **lvariance_init**      | Initializer of the bandwidths of the exponential structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    
  * if **layer** is "Rational":

    | **name** | **description** |
    | --- | --- |
    | **lnum**       | Number of rational structure elements that will be evaluated on the persistence diagram points. |
    | **lmean_init**          | Initializer of the means of the rational structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **lvariance_init**      | Initializer of the bandwidths of the rational structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    | **lalpha_init**         | Initializer of the exponents of the rational structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    
  * if **layer** is "RationalHat":

    | **name** | **description** |
    | --- | --- |
    | **lnum**      | Number of rational hat structure elements that will be evaluated on the persistence diagram points. |
    | **lmean_init**         | Initializer of the means of the rational hat structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **lr_init**            | Initializer of the threshold of the rational hat structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    | **q**                 | Norm parameter. |
