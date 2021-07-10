# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carrière
#
# Copyright (C) 2018-2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import tensorflow as tf
import numpy as np

class GridPerslayWeight(tf.keras.layers.Layer):

    def __init__(self, grid, grid_bnds, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.grid = tf.Variable(initial_value=grid, trainable=True)
        self.grid_bnds = grid_bnds
    
    def build(self, input_shape):
        return self

    def call(self, diagrams):
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
    
    def __init__(self, gaussians, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.W = tf.Variable(initial_value=gaussians, trainable=True)

    def build(self, input_shape):
        return self
        
    def call(self, diagrams):
        means     = tf.expand_dims(tf.expand_dims(self.W[:2,:],0),0)
        variances = tf.expand_dims(tf.expand_dims(self.W[2:,:],0),0)
        diagrams  = tf.expand_dims(diagrams, -1)
        dists     = tf.math.multiply(tf.math.square(diagrams-means), 1/tf.math.square(variances))
        weight    = tf.math.reduce_sum(tf.math.exp(tf.math.reduce_sum(-dists, axis=2)), axis=2)
        return weight
    
class PowerPerslayWeight(tf.keras.layers.Layer):
    
    def __init__(self, constant, power, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.constant = tf.Variable(initial_value=constant, trainable=True)
        self.power = power
        
    def build(self, input_shape):
        return self
    
    def call(self, diagrams):
        weight = self.constant * tf.math.pow(tf.math.abs(diagrams[:,:,1]-diagrams[:,:,0]), self.power)
        return weight
    

class GaussianPerslayPhi(tf.keras.layers.Layer):
    
    def __init__(self, image_size, image_bnds, variance, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.image_size = image_size
        self.image_bnds = image_bnds
        self.variance   = tf.Variable(initial_value=variance, trainable=True)
        
    def build(self, input_shape):
        return self
        
    def call(self, diagrams):
        diagrams_d = tf.concat([diagrams[:,:,0:1], diagrams[:,:,1:2]-diagrams[:,:,0:1]], axis=2)
        step = [(self.image_bnds[i][1]-self.image_bnds[i][0])/self.image_size[i] for i in range(2)]
        coords = [tf.range(self.image_bnds[i][0], self.image_bnds[i][1], step[i]) for i in range(2)]
        M = tf.meshgrid(*coords)
        mu = tf.concat([tf.expand_dims(tens, 0) for tens in M], axis=0)
        for _ in range(2):
            diagrams_d = tf.expand_dims(diagrams_d,-1)
        dists = tf.math.square(diagrams_d-mu) / (2*tf.math.square(self.variance))
        gauss = tf.math.exp(tf.math.reduce_sum(-dists, axis=2)) / (2*np.pi*tf.math.square(self.variance))
        return tf.expand_dims(gauss,-1), M[0].shape + tuple([1])
     
class TentPerslayPhi(tf.keras.layers.Layer):
    
    def __init__(self, samples, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.samples   = tf.Variable(initial_value=samples, trainable=True)
        
    def build(self, input_shape):
        return self
        
    def call(self, diagrams):
        samples_d = tf.expand_dims(tf.expand_dims(self.samples,0),0)
        xs, ys = diagrams[:,:,0:1], diagrams[:,:,1:2]
        tent = tf.math.maximum(.5*(ys-xs) - tf.math.abs(samples_d-.5*(ys+xs)), np.array([0.]))
        return tent, self.samples.shape
    
class FlatPerslayPhi(tf.keras.layers.Layer):
    
    def __init__(self, samples, theta, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.samples = tf.Variable(initial_value=samples, trainable=True)
        self.theta   = tf.Variable(initial_value=theta,   trainable=True)
        
    def build(self, input_shape):
        return self
        
    def call(self, diagrams):
        samples_d = tf.expand_dims(tf.expand_dims(self.samples,0),0)
        xs, ys = diagrams[:,:,0:1], diagrams[:,:,1:2]
        flat = 1./(1.+tf.math.exp(-self.theta*(.5*(ys-xs)-tf.math.abs(samples_d-.5*(ys+xs)))))
        return flat, self.samples.shape

class Perslay(tf.keras.layers.Layer):

    def __init__(self, weight, phi, perm_op, rho, **kwargs):
        super().__init__(dynamic=True, **kwargs)
        self.weight  = weight
        self.phi     = phi
        self.pop     = perm_op  
        self.rho     = rho

    def build(self, input_shape):
        return self

    def call(self, diagrams):

        vector, dim = self.phi(diagrams)
        weight = self.weight(diagrams)
        for _ in range(len(dim)):
            weight = tf.expand_dims(weight, -1)
        vector = tf.math.multiply(vector, weight)
          
        permop = self.pop
        if type(permop) == str and permop[:3] == 'top':
            k = int(permop[3:])
            vector = vector.to_tensor(default_value=-1e10)
            vector = tf.math.top_k(tf.transpose(vector, perm=[0, 2, 1]), k=k).values
            vector = tf.reshape(vector, [-1,k*dim[0]])
        else:
            vector = permop(vector, axis=1)

        vector = self.rho(vector)
            
        return vector
