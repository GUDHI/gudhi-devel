"""Module :mod:`perskay.archi` implement the persistence layer."""

# Authors: Mathieu Carriere <mathieu.carriere3@gmail.com>
# License: MIT
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import tensorflow as tf


# DeepSet PersLay
def permutation_equivariant_layer(inp, dimension, perm_op, L_init, G_init, b_init, L_const, G_const, b_const):
    dimension_before, num_pts = inp.shape[2].value, inp.shape[1].value
    lbda = tf.get_variable("L", shape=[dimension_before, dimension], initializer=L_init) if not L_const \
        else tf.get_variable("L", initializer=L_init)
    b = tf.get_variable("b", shape=[1, 1, dimension], initializer=b_init) if not b_const \
        else tf.get_variable("b", initializer=b_init)
    A = tf.reshape(tf.einsum("ijk,kl->ijl", inp, lbda), [-1, num_pts, dimension])
    if perm_op is not None:
        if perm_op == "max":
            beta = tf.tile(tf.expand_dims(tf.reduce_max(inp, axis=1), 1), [1, num_pts, 1])
        elif perm_op == "min":
            beta = tf.tile(tf.expand_dims(tf.reduce_min(inp, axis=1), 1), [1, num_pts, 1])
        elif perm_op == "sum":
            beta = tf.tile(tf.expand_dims(tf.reduce_sum(inp, axis=1), 1), [1, num_pts, 1])
        else:
            raise Exception("perm_op should be min, max or sum")
        gamma = tf.get_variable("G", shape=[dimension_before, dimension], initializer=G_init) if not G_const \
            else tf.get_variable("G", initializer=G_init)
        B = tf.reshape(tf.einsum("ijk,kl->ijl", beta, gamma), [-1, num_pts, dimension])
        return A - B + b
    else:
        return A + b


# Gaussian PersLay
def gaussian_layer(inp, num_gaussians, m_init, s_init, m_const, s_const):
    dimension_before, num_pts = inp.shape[2].value, inp.shape[1].value
    mu = tf.get_variable("m", shape=[1, 1, dimension_before, num_gaussians],
                         initializer=m_init) if not m_const else tf.get_variable("m", initializer=m_init)
    sg = tf.get_variable("s", shape=[1, 1, dimension_before, num_gaussians],
                         initializer=s_init) if not s_const else tf.get_variable("s", initializer=s_init)
    bc_inp = tf.expand_dims(inp, -1)
    return tf.exp(tf.reduce_sum(-tf.multiply(tf.square(bc_inp - mu), tf.square(sg)), axis=2))


# Landscape PersLay
def landscape_layer(inp, num_samples, s_init, s_const):
    # num_pts = inp.shape[1].value
    sp = tf.get_variable("s", shape=[1, 1, num_samples], initializer=s_init) if not s_const \
        else tf.get_variable("s", initializer=s_init)
    return tf.maximum(inp[:, :, 1:2] - tf.abs(sp - inp[:, :, 0:1]), np.array([0]))


# Persistence Image PersLay
def image_layer(inp, im_size, im_bnds, s_init, s_const):
    dimension_before, num_pts = inp.shape[2].value, inp.shape[1].value
    coords = [tf.range(start=im_bnds[i][0], limit=im_bnds[i][1], delta=(im_bnds[i][1] - im_bnds[i][0]) / im_size[i]) for
              i in range(dimension_before)]
    M = tf.meshgrid(*coords)
    mu = tf.concat([tf.expand_dims(tens, 0) for tens in M], axis=0)
    sg = tf.get_variable("s", shape=[1, 1, 1] + [1 for _ in range(dimension_before)],
                         initializer=s_init) if not s_const else tf.get_variable("s", initializer=s_init)
    bc_inp = tf.reshape(inp, [-1, num_pts, dimension_before] + [1 for _ in range(dimension_before)])
    return tf.exp(tf.reduce_sum(-tf.multiply(tf.square(bc_inp - mu), tf.square(sg)), axis=2))
