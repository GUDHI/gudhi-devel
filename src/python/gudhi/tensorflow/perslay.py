# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière
#
# Copyright (C) 2021 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import tensorflow as tf
import math

class GridPerslayWeight(tf.keras.layers.Layer):
    """
    This is a class for computing a differentiable weight function for persistence diagram points. This function is defined from an array that contains its values on a 2D grid.
    """
    def __init__(self, grid, grid_bnds, **kwargs):
        """
        Constructor for the GridPerslayWeight class.
  
        Parameters:
            grid (n x n numpy array): grid of values.
            grid_bnds (2 x 2 numpy array): boundaries of the grid, of the form [[min_x, max_x], [min_y, max_y]].
        """
        super().__init__(**kwargs)
        self.grid = tf.Variable(initial_value=grid, trainable=True)
        self.grid_bnds = grid_bnds
    
    def build(self, input_shape):
        return self

    def call(self, diagrams):
        """
        Apply GridPerslayWeight on a ragged tensor containing a list of persistence diagrams.

        Parameters:
            diagrams (n x None x 2): ragged tensor containing n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.

        Returns:
            weight (n x None): ragged tensor containing the weights of the points in the n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.
        """
        grid_shape = self.grid.shape
        indices = []
        for dim in range(2):
            [m,M] = self.grid_bnds[dim]
            coords = tf.expand_dims(diagrams[:,:,dim],-1)
            ids = grid_shape[dim]*(coords-m)/(M-m)
            indices.append(tf.cast(ids, tf.int32))
        weight = tf.gather_nd(params=self.grid, indices=tf.concat(indices, axis=2))
        return weight
    
class GaussianMixturePerslayWeight(tf.keras.layers.Layer):
    """
    This is a class for computing a differentiable weight function for persistence diagram points. This function is defined from a mixture of Gaussian functions.
    """
    def __init__(self, gaussians, **kwargs):
        """
        Constructor for the GridPerslayWeight class.
  
        Parameters:
            gaussians (4 x n numpy array): parameters of the n Gaussian functions, of the form transpose([[mu_x^1, mu_y^1, sigma_x^1, sigma_y^1], ..., [mu_x^n, mu_y^n, sigma_x^n, sigma_y^n]]). 
        """
        super().__init__(**kwargs)
        self.W = tf.Variable(initial_value=gaussians, trainable=True)

    def build(self, input_shape):
        return self
        
    def call(self, diagrams):
        """
        Apply GaussianMixturePerslayWeight on a ragged tensor containing a list of persistence diagrams.

        Parameters:
            diagrams (n x None x 2): ragged tensor containing n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.

        Returns:
            weight (n x None): ragged tensor containing the weights of the points in the n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.
        """
        means     = tf.expand_dims(tf.expand_dims(self.W[:2,:],0),0)
        variances = tf.expand_dims(tf.expand_dims(self.W[2:,:],0),0)
        diagrams  = tf.expand_dims(diagrams, -1)
        dists     = tf.math.multiply(tf.math.square(diagrams-means), 1/tf.math.square(variances))
        weight    = tf.math.reduce_sum(tf.math.exp(tf.math.reduce_sum(-dists, axis=2)), axis=2)
        return weight
    
class PowerPerslayWeight(tf.keras.layers.Layer):
    """
    This is a class for computing a differentiable weight function for persistence diagram points. This function is defined as a constant multiplied by the distance to the diagonal of the persistence diagram point raised to some power.
    """
    def __init__(self, constant, power, **kwargs):
        """
        Constructor for the PowerPerslayWeight class.
  
        Parameters:
            constant (float): constant value.
            power (float): power applied to the distance to the diagonal. 
        """
        super().__init__(**kwargs)
        self.constant = tf.Variable(initial_value=constant, trainable=True)
        self.power = power
        
    def build(self, input_shape):
        return self
    
    def call(self, diagrams):
        """
        Apply PowerPerslayWeight on a ragged tensor containing a list of persistence diagrams.

        Parameters:
            diagrams (n x None x 2): ragged tensor containing n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.

        Returns:
            weight (n x None): ragged tensor containing the weights of the points in the n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.
        """
        weight = self.constant * tf.math.pow(tf.math.abs(diagrams[:,:,1]-diagrams[:,:,0]), self.power)
        return weight
    

class GaussianPerslayPhi(tf.keras.layers.Layer):
    """
    This is a class for computing a transformation function for persistence diagram points. This function turns persistence diagram points into 2D Gaussian functions centered on the points, that are then evaluated on a regular 2D grid.
    """
    def __init__(self, image_size, image_bnds, variance, **kwargs):
        """
        Constructor for the GaussianPerslayPhi class.
  
        Parameters:
            image_size (int numpy array): number of grid elements on each grid axis, of the form [n_x, n_y].
            image_bnds (2 x 2 numpy array): boundaries of the grid, of the form [[min_x, max_x], [min_y, max_y]].
            variance (float): variance of the Gaussian functions. 
        """
        super().__init__(**kwargs)
        self.image_size = image_size
        self.image_bnds = image_bnds
        self.variance   = tf.Variable(initial_value=variance, trainable=True)
        
    def build(self, input_shape):
        return self
        
    def call(self, diagrams):
        """
        Apply GaussianPerslayPhi on a ragged tensor containing a list of persistence diagrams.

        Parameters:
            diagrams (n x None x 2): ragged tensor containing n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.

        Returns:
            output (n x None x image_size x image_size x 1): ragged tensor containing the evaluations on the 2D grid of the 2D Gaussian functions corresponding to the persistence diagram points, in the form of a 2D image with 1 channel that can be processed with, e.g., convolutional layers. The second dimension is ragged since persistence diagrams can have different numbers of points.
            output_shape (int numpy array): shape of the output tensor.
        """
        diagrams_d = tf.concat([diagrams[:,:,0:1], diagrams[:,:,1:2]-diagrams[:,:,0:1]], axis=2)
        step = [(self.image_bnds[i][1]-self.image_bnds[i][0])/self.image_size[i] for i in range(2)]
        coords = [tf.range(self.image_bnds[i][0], self.image_bnds[i][1], step[i]) for i in range(2)]
        M = tf.meshgrid(*coords)
        mu = tf.concat([tf.expand_dims(tens, 0) for tens in M], axis=0)
        for _ in range(2):
            diagrams_d = tf.expand_dims(diagrams_d,-1)
        dists = tf.math.square(diagrams_d-mu) / (2*tf.math.square(self.variance))
        gauss = tf.math.exp(tf.math.reduce_sum(-dists, axis=2)) / (2*math.pi*tf.math.square(self.variance))
        output = tf.expand_dims(gauss,-1)
        output_shape = M[0].shape + tuple([1])
        return output, output_shape
     
class TentPerslayPhi(tf.keras.layers.Layer):
    """
    This is a class for computing a transformation function for persistence diagram points. This function turns persistence diagram points into 1D tent functions (linearly increasing on the first half of the bar corresponding to the point from zero to half of the bar length, linearly decreasing on the second half and zero elsewhere) centered on the points, that are then evaluated on a regular 1D grid.
    """
    def __init__(self, samples, **kwargs):
        """
        Constructor for the GaussianPerslayPhi class.
  
        Parameters:
            samples (float numpy array): grid elements on which to evaluate the tent functions, of the form [x_1, ..., x_n].
        """
        super().__init__(**kwargs)
        self.samples   = tf.Variable(initial_value=samples, trainable=True)
        
    def build(self, input_shape):
        return self
        
    def call(self, diagrams):
        """
        Apply TentPerslayPhi on a ragged tensor containing a list of persistence diagrams.

        Parameters:
            diagrams (n x None x 2): ragged tensor containing n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.

        Returns:
            output (n x None x num_samples): ragged tensor containing the evaluations on the 1D grid of the 1D tent functions corresponding to the persistence diagram points. The second dimension is ragged since persistence diagrams can have different numbers of points.
            output_shape (int numpy array): shape of the output tensor.
        """
        samples_d = tf.expand_dims(tf.expand_dims(self.samples,0),0)
        xs, ys = diagrams[:,:,0:1], diagrams[:,:,1:2]
        output = tf.math.maximum(.5*(ys-xs) - tf.math.abs(samples_d-.5*(ys+xs)), tf.constant([0.]))
        output_shape = self.samples.shape
        return output, output_shape
    
class FlatPerslayPhi(tf.keras.layers.Layer):
    """
    This is a class for computing a transformation function for persistence diagram points. This function turns persistence diagram points into 1D constant functions (that evaluate to half of the bar length on the bar corresponding to the point and zero elsewhere), that are then evaluated on a regular 1D grid.
    """
    def __init__(self, samples, theta, **kwargs):
        """
        Constructor for the FlatPerslayPhi class.
  
        Parameters:
            samples (float numpy array): grid elements on which to evaluate the constant functions, of the form [x_1, ..., x_n].
            theta (float): sigmoid parameter used to approximate the constant function with a differentiable sigmoid function. The bigger the theta, the closer to a constant function the output will be. 
        """
        super().__init__(**kwargs)
        self.samples = tf.Variable(initial_value=samples, trainable=True)
        self.theta   = tf.Variable(initial_value=theta,   trainable=True)
        
    def build(self, input_shape):
        return self
        
    def call(self, diagrams):
        """
        Apply FlatPerslayPhi on a ragged tensor containing a list of persistence diagrams.

        Parameters:
            diagrams (n x None x 2): ragged tensor containing n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.

        Returns:
            output (n x None x num_samples): ragged tensor containing the evaluations on the 1D grid of the 1D constant functions corresponding to the persistence diagram points. The second dimension is ragged since persistence diagrams can have different numbers of points.
            output_shape (int numpy array): shape of the output tensor.
        """
        samples_d = tf.expand_dims(tf.expand_dims(self.samples,0),0)
        xs, ys = diagrams[:,:,0:1], diagrams[:,:,1:2]
        output = 1./(1.+tf.math.exp(-self.theta*(.5*(ys-xs)-tf.math.abs(samples_d-.5*(ys+xs)))))
        output_shape = self.samples.shape
        return output, output_shape

class Perslay(tf.keras.layers.Layer):
    """
    This is a TensorFlow layer for vectorizing persistence diagrams in a differentiable way within a neural network. This function implements the PersLay equation, see `the corresponding article <http://proceedings.mlr.press/v108/carriere20a.html>`_.
    """
    def __init__(self, weight, phi, perm_op, rho, **kwargs):
        """
        Constructor for the Perslay class.

        Parameters:
            weight (function): weight function for the persistence diagram points. Can be either :class:`~gudhi.tensorflow.perslay.GridPerslayWeight`, :class:`~gudhi.tensorflow.perslay.GaussianMixturePerslayWeight`, :class:`~gudhi.tensorflow.perslay.PowerPerslayWeight`, or a custom TensorFlow function that takes persistence diagrams as argument (represented as an (n x None x 2) ragged tensor, where n is the number of diagrams).
            phi (function): transformation function for the persistence diagram points. Can be either :class:`~gudhi.tensorflow.perslay.GaussianPerslayPhi`, :class:`~gudhi.tensorflow.perslay.TentPerslayPhi`, :class:`~gudhi.tensorflow.perslay.FlatPerslayPhi`, or a custom TensorFlow class (that can have trainable parameters) with a method `call` that takes persistence diagrams as argument (represented as an (n x None x 2) ragged tensor, where n is the number of diagrams).
            perm_op (function): permutation invariant function, such as `tf.math.reduce_sum`, `tf.math.reduce_mean`, `tf.math.reduce_max`, `tf.math.reduce_min`, or a custom TensorFlow function that takes two arguments: a tensor and an axis on which to apply the permutation invariant operation. If perm_op is the string "topk" (where k is a number), this function will be computed as `tf.math.top_k` with parameter `int(k)`.
            rho (function): postprocessing function that is applied after the permutation invariant operation. Can be any TensorFlow layer.
        """
        super().__init__(**kwargs)
        self.weight  = weight
        self.phi     = phi
        self.perm_op = perm_op  
        self.rho     = rho

    def build(self, input_shape):
        return self

    def call(self, diagrams):
        """
        Apply Perslay on a ragged tensor containing a list of persistence diagrams.

        Parameters:
            diagrams (n x None x 2): ragged tensor containing n persistence diagrams. The second dimension is ragged since persistence diagrams can have different numbers of points.

        Returns:
            vector (n x output_shape): tensor containing the vectorizations of the persistence diagrams.
        """
        vector, dim = self.phi(diagrams)
        weight = self.weight(diagrams)
        for _ in range(len(dim)):
            weight = tf.expand_dims(weight, -1)
        vector = tf.math.multiply(vector, weight)
          
        permop = self.perm_op
        if type(permop) == str and permop[:3] == 'top':
            k = int(permop[3:])
            vector = vector.to_tensor(default_value=-1e10)
            vector = tf.math.top_k(tf.transpose(vector, perm=[0, 2, 1]), k=k).values
            vector = tf.reshape(vector, [-1,k*dim[0]])
        else:
            vector = permop(vector, axis=1)

        vector = self.rho(vector)
            
        return vector
