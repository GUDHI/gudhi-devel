# How to call and use PersLay


The python function `perslay` takes four arguments:

  * the first three ones are: **output** which is the list that will contain the output of PersLay, **name** which is a string defining the name of the PersLay operation for tensorflow, and **diag** which is a numpy array containing the persistence diagrams.

  * the fourth argument is a python dictionary containing the parameters of PersLay. Examples can be found in `tutorialPersLay.ipynb`.
 
In the following description of PersLay parameters, each parameter, or dictionary key, that contains "_init" in its name is optimized and learned by PersLay during training. If you do not want to optimize the vectorization, set the (optional) keys **train_vect** and **train_weight** to False.

  * The following keys are mandatory:

    | **name** | **description** |
    | --- | --- |
    | **layer**              | Either "pm", "im", "ls", "bc", "en", "ex", "rt", "rh". Type of the PersLay layer. "im" is for [persistence images](https://arxiv.org/abs/1507.06217), "ls" is for [persistence landscapes](http://www.jmlr.org/papers/volume16/bubenik15a/bubenik15a.pdf), "ex", "rt" and "rh" are for [structure elements](http://jmlr.org/beta/papers/v20/18-358.html), "pm" is for the original DeepSet layer, defined in [this article](https://arxiv.org/abs/1703.06114), "bc" is for [Betti curves](https://www.jstage.jst.go.jp/article/tjsai/32/3/32_D-G72/_pdf) and "en" is for [entropy](https://arxiv.org/abs/1803.08304). |
    | **perm_op**            | Either "sum", "mean", "max", "topk". Permutation invariant operation. |
    | **keep**               | Number of top values to keep. Used only if **perm_op** is "topk". |
    | **persistence_weight** | Either "linear", "grid", "gmix" or None. Weight function to be applied on persistence diagram points. If "linear", this function is linear with respect to the distance to the diagonal of the point. If "grid", this function is piecewise-constant and defined with pixel values of a grid. If "gmix", this function is defined as a mixture of Gaussians. If None, no weighting is applied. |

Depending on what **persistence_weight** is, the following additional keys are requested:

  * if **persistence_weight** is "linear":

    | **name** | **description** |
    | --- | --- |
    | **coeff_init**         | Initializer of the coefficient of the linear weight function. It can be either a single value, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(1., 1.). | 
    | **coeff_const**        | Boolean specifying if **coeff_init** is initialized with a value (True) or randomly with tensorflow (False). |
 
  * if **persistence_weight** is "grid":

    | **name** | **description** |
    | --- | --- |
    | **grid_size**          | Grid size of the grid weight function. It is a tuple of integer values, such as (10,10). | 
    | **grid_bnds**          | Grid boundaries of the grid weight function. It is a tuple containing two tuples, each containing the minimum and maximum values of each axis of the plane. Example: ((-0.01, 1.01), (-0.01, 1.01)). | 
    | **grid_init**          | Initializer for the pixel values of the grid weight function. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(1., 1.).| 
    | **grid_const**         | Boolean specifying if **grid_init** is initialized with an array (True) or randomly with tensorflow (False). |

  * if **persistence_weight** is "gmix":

    | **name** | **description** |
    | --- | --- |
    | **gmix_num**           | Number of Gaussian functions of the mixture of Gaussians weight function. |
    | **gmix_m_init**        | Initializer of the means of the mixture of Gaussians weight function. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(1., 1.). |
    | **gmix_m_const**       | Boolean specifying if **gmix_m_init** is initialized with an array (True) or randomly with tensorflow (False). |
    | **gmix_v_init**        | Initializer of the bandwidths of the mixture of Gaussians weight function. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(1., 1.). |
    | **gmix_v_const**       | Boolean specifying if **gmix_v_init** is initialized with an array (True) or randomly with tensorflow (False). Used only if **persistence_weight** is "gmix". |

Depending on what **layer** is, the following additional keys are requested:

  * if **layer** is "pm":

    | **name** | **description** |
    | --- | --- |
    | **peq**                | Sequence of permutation equivariant operations, as defined in [the DeepSet article](). It is a list of tuples of the form (*dim*, *operation*). Each tuple defines a permutation equivariant function of dimension *dim* and second permutation operation *operation* (string, either "max", "min", "sum" or None). Second permutation operation is optional and is not applied if *operation* is set to None. Example: [(150, "max"), (75, None)]. |
    | **weight_init**        | Initializer for the matrices of the permutation equivariant operations. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.).| 
    | **weight_const**       | Boolean specifying if **weight_init** is initialized with a value (True) or randomly with tensorflow (False). | 
    | **bias_init**          | Initializer for the biases of the permutation equivariant operations. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **bias_const**         | Boolean specifying if **bias_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **fc_layers**          | Sequence of fully-connected operations to be applied after the permutation invariant operation. It is a list of tuples of the form (*dim*, *pro*, *dropout*). Each tuple defines a fully-connected operation, with dimension *dim* (integer) and processing *pro* (string, e.g. "bdr" ---> batch-norm, dropout, relu). If there is a "d" in string, i.e. dropout, the dropout value can be specified with *dro* (float, default 0.9). Example: [(150,"br"), (75, "bd", 0.85)].|

  * if **layer** is "im":

    | **name** | **description** |
    | --- | --- |
    | **image_size**         | Persistence image size. It is a tuple of integer values, such as (10,10). | 
    | **image_bnds**         | Persistence image boundaries. It is a tuple containing two tuples, each containing the minimum and maximum values of each axis of the plane. Example: ((-0.01, 1.01), (-0.01, 1.01)). |
    | **variance_init**      | Initializer for the bandwidths of the Gaussian functions centered on the persistence image pixels. It can be either a single value, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    | **variance_const**     | Boolean specifying if **variance_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **cv_layers**          | Sequence of convolution operations to be applied after the permutation invariant operation. Used only if **layer** is "im". It is a list of tuples of the form (*num_filters*, *kernel_size*, *pro*, *dropout*). Each tuple defines a convolution operation, with number of filters *num_filters* (integer), kernel size *kernel_size* (integer), and processing *pro* (string, e.g. "bdr" ---> batch-norm, dropout, relu). If there is a "d" in string, i.e. dropout, the dropout value can be specified with *dro* (float, default 0.9). Example: [(10,3,"bd"), (5,3,"dr",0.8)]. | 
    
  * if **layer** is "ls":

    | **name** | **description** |
    | --- | --- |
    | **num_samples**        | Number of samples of the diagonal that will be evaluated on the persistence landscapes. |
    | **sample_init**        | Initializer of the samples of the diagonal. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **sample_const**       | Boolean specifying if **sample_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **fc_layers**          | Sequence of fully-connected operations to be applied after the permutation invariant operation. It is a list of tuples of the form (*dim*, *pro*, *dropout*). Each tuple defines a fully-connected operation, with dimension *dim* (integer) and processing *pro* (string, e.g. "bdr" ---> batch-norm, dropout, relu). If there is a "d" in string, i.e. dropout, the dropout value can be specified with *dro* (float, default 0.9). Example: [(150,"br"), (75, "bd", 0.85)].|

  * if **layer** is "bc":

    | **name** | **description** |
    | --- | --- |
    | **num_samples**        | Number of samples of the diagonal that will be evaluated on the Betti curves. |
    | **sample_init**        | Initializer of the samples of the diagonal. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **sample_const**       | Boolean specifying if **sample_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **theta**              | Sigmoid parameter used for approximating the piecewise constant functions associated to the persistence diagram points. |
    | **fc_layers**          | Sequence of fully-connected operations to be applied after the permutation invariant operation. It is a list of tuples of the form (*dim*, *pro*, *dropout*). Each tuple defines a fully-connected operation, with dimension *dim* (integer) and processing *pro* (string, e.g. "bdr" ---> batch-norm, dropout, relu). If there is a "d" in string, i.e. dropout, the dropout value can be specified with *dro* (float, default 0.9). Example: [(150,"br"), (75, "bd", 0.85)].|

  * if **layer** is "en":

    | **name** | **description** |
    | --- | --- |
    | **num_samples**        | Number of samples on the diagonal that will be evaluated on the persistence entropies. |
    | **sample_init**        | Initializer of the samples of the diagonal. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **sample_const**       | Boolean specifying if **sample_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **theta**              | Sigmoid parameter used for approximating the piecewise constant functions associated to the persistence diagram points. |
    | **fc_layers**          | Sequence of fully-connected operations to be applied after the permutation invariant operation. It is a list of tuples of the form (*dim*, *pro*, *dropout*). Each tuple defines a fully-connected operation, with dimension *dim* (integer) and processing *pro* (string, e.g. "bdr" ---> batch-norm, dropout, relu). If there is a "d" in string, i.e. dropout, the dropout value can be specified with *dro* (float, default 0.9). Example: [(150,"br"), (75, "bd", 0.85)].|

  * if **layer** is "ex":

    | **name** | **description** |
    | --- | --- |
    | **num_elements**       | Number of exponential structure elements that will be evaluated on the persistence diagram points. |
    | **mean_init**          | Initializer of the means of the exponential structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **mean_const**         | Boolean specifying if **mean_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **variance_init**      | Initializer of the bandwidths of the exponential structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    | **variance_const**     | Boolean specifying if **variance_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **fc_layers**          | Sequence of fully-connected operations to be applied after the permutation invariant operation. It is a list of tuples of the form (*dim*, *pro*, *dropout*). Each tuple defines a fully-connected operation, with dimension *dim* (integer) and processing *pro* (string, e.g. "bdr" ---> batch-norm, dropout, relu). If there is a "d" in string, i.e. dropout, the dropout value can be specified with *dro* (float, default 0.9). Example: [(150,"br"), (75, "bd", 0.85)].|

  * if **layer** is "rt":

    | **name** | **description** |
    | --- | --- |
    | **num_elements**       | Number of rational structure elements that will be evaluated on the persistence diagram points. |
    | **mean_init**          | Initializer of the means of the rational structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **mean_const**         | Boolean specifying if **mean_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **variance_init**      | Initializer of the bandwidths of the rational structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    | **variance_const**     | Boolean specifying if **variance_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **alpha_init**         | Initializer of the exponents of the rational structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    | **alpha_const**        | Boolean specifying if **alpha_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **fc_layers**          | Sequence of fully-connected operations to be applied after the permutation invariant operation. It is a list of tuples of the form (*dim*, *pro*, *dropout*). Each tuple defines a fully-connected operation, with dimension *dim* (integer) and processing *pro* (string, e.g. "bdr" ---> batch-norm, dropout, relu). If there is a "d" in string, i.e. dropout, the dropout value can be specified with *dro* (float, default 0.9). Example: [(150,"br"), (75, "bd", 0.85)].|

  * if **layer** is "rh":

    | **name** | **description** |
    | --- | --- |
    | **num_elements**      | Number of rational hat structure elements that will be evaluated on the persistence diagram points. |
    | **mean_init**         | Initializer of the means of the rational hat structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(0., 1.). |
    | **mean_const**        | Boolean specifying if **mean_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **r_init**            | Initializer of the threshold of the rational hat structure elements. It can be either a numpy array of values, or a random initializer from tensorflow, such as tensorflow.random_uniform_initializer(3., 3.). | 
    | **r_const**           | Boolean specifying if **r_init** is initialized with a value (True) or randomly with tensorflow (False). |
    | **q**                 | Norm parameter. |
    | **fc_layers**         | Sequence of fully-connected operations to be applied after the permutation invariant operation. It is a list of tuples of the form (*dim*, *pro*, *dropout*). Each tuple defines a fully-connected operation, with dimension *dim* (integer) and processing *pro* (string, e.g. "bdr" ---> batch-norm, dropout, relu). If there is a "d" in string, i.e. dropout, the dropout value can be specified with *dro* (float, default 0.9). Example: [(150,"br"), (75, "bd", 0.85)].|
