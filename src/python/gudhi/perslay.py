# This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
# See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
# Author(s):       Mathieu Carriere
#
# Copyright (C) 2019 Inria
#
# Modification(s):
#   - YYYY/MM Author: Description of the modification

import numpy as np
import tensorflow as tf


# Post-processing operation with combination of batch normalization, dropout and relu
def _post_processing(vector, pro, dropout_value=.9):
    for c in pro:
        if c == "b":
            vector = tf.layers.batch_normalization(vector)
        if c == "d":
            vector = tf.nn.dropout(vector, dropout_value)
        if c == "r":
            vector = tf.nn.relu(vector)
    return vector

# Vectorization implementing DeepSet architecture
def permutation_equivariant_layer(inp, dimension, perm_op, L_init, G_init, bias_init, L_const, G_const, bias_const, train_vect):
    """ DeepSet PersLay """
    dimension_before, num_pts = inp.shape[2].value, inp.shape[1].value
    lbda = tf.get_variable("L", shape=[dimension_before, dimension], initializer=L_init, trainable=train_vect)   if not L_const     else tf.get_variable("L", initializer=L_init)
    b    = tf.get_variable("b", shape=[1, 1, dimension], initializer=bias_init, trainable=train_vect)            if not bias_const  else tf.get_variable("b", initializer=bias_init)
    A    = tf.reshape(tf.einsum("ijk,kl->ijl", inp, lbda), [-1, num_pts, dimension])
    if perm_op is not None:
        if perm_op == "max":
            beta = tf.tile(tf.expand_dims(tf.reduce_max(inp, axis=1), 1), [1, num_pts, 1])
        elif perm_op == "min":
            beta = tf.tile(tf.expand_dims(tf.reduce_min(inp, axis=1), 1), [1, num_pts, 1])
        elif perm_op == "sum":
            beta = tf.tile(tf.expand_dims(tf.reduce_sum(inp, axis=1), 1), [1, num_pts, 1])
        else:
            raise Exception("perm_op should be min, max or sum")
        gamma = tf.get_variable("G", shape=[dimension_before, dimension], initializer=G_init, trainable=train_vect) if not G_const else tf.get_variable("G", initializer=G_init)
        B = tf.reshape(tf.einsum("ijk,kl->ijl", beta, gamma), [-1, num_pts, dimension])
        return A - B + b
    else:
        return A + b

# Vectorizations taken from "Learning Representations of Persistence Barcodes"
def rational_hat_layer(inp, num_elements, q, mean_init, r_init, mean_const, r_const, train_vect):
    """ Rational Hat PersLay """
    dimension_before, num_pts = inp.shape[2].value, inp.shape[1].value
    mu = tf.get_variable("m", shape=[1, 1, dimension_before, num_elements], initializer=mean_init, trainable=train_vect)      if not mean_const      else tf.get_variable("m", initializer=mean_init)
    r  = tf.get_variable("r", shape=[1, 1, 1],                              initializer=r_init, trainable=train_vect)         if not r_const         else tf.get_variable("r", initializer=r_init)
    bc_inp = tf.expand_dims(inp, -1)
    norms = tf.norm(bc_inp - mu, ord=q, axis=2)
    return 1/(1 + norms) - 1/(1 + tf.abs(tf.abs(r)-norms)) 

def rational_layer(inp, num_elements, mean_init, variance_init, alpha_init, mean_const, variance_const, alpha_const, train_vect):
    """ Rational PersLay """
    dimension_before, num_pts = inp.shape[2].value, inp.shape[1].value
    mu = tf.get_variable("m", shape=[1, 1, dimension_before, num_elements], initializer=mean_init, trainable=train_vect)      if not mean_const      else tf.get_variable("m", initializer=mean_init)
    sg = tf.get_variable("s", shape=[1, 1, dimension_before, num_elements], initializer=variance_init, trainable=train_vect)  if not variance_const  else tf.get_variable("s", initializer=variance_init)
    al = tf.get_variable("a", shape=[1, 1, num_elements],                   initializer=alpha_init, trainable=train_vect)     if not alpha_const     else tf.get_variable("a", initializer=alpha_init)
    bc_inp = tf.expand_dims(inp, -1)
    return 1/tf.pow(1+tf.reduce_sum(tf.multiply(tf.abs(bc_inp - mu), tf.abs(sg)), axis=2), al)

def exponential_layer(inp, num_elements, mean_init, variance_init, mean_const, variance_const, train_vect):
    """ Exponential PersLay """
    dimension_before, num_pts = inp.shape[2].value, inp.shape[1].value
    mu = tf.get_variable("m", shape=[1, 1, dimension_before, num_elements], initializer=mean_init, trainable=train_vect)      if not mean_const      else tf.get_variable("m", initializer=mean_init)
    sg = tf.get_variable("s", shape=[1, 1, dimension_before, num_elements], initializer=variance_init, trainable=train_vect)  if not variance_const  else tf.get_variable("s", initializer=variance_init)
    bc_inp = tf.expand_dims(inp, -1)
    return tf.exp(tf.reduce_sum(-tf.multiply(tf.square(bc_inp - mu), tf.square(sg)), axis=2))

# Vectorizations implementing persistence landscapes
def landscape_layer(inp, num_samples, sample_init, sample_const, train_vect):
    """ Landscape PersLay """
    sp = tf.get_variable("s", shape=[1, 1, num_samples], initializer=sample_init, trainable=train_vect) if not sample_const else tf.get_variable("s", initializer=sample_init)
    return tf.maximum( .5 * (inp[:, :, 1:2] - inp[:, :, 0:1]) - tf.abs(sp - .5 * (inp[:, :, 1:2] + inp[:, :, 0:1])), np.array([0]))

# Vectorizations implementing Betti curves
def betti_layer(inp, theta, num_samples, sample_init, sample_const, train_vect):
    """ Betti PersLay """
    sp = tf.get_variable("s", shape=[1, 1, num_samples], initializer=sample_init, trainable=train_vect) if not sample_const else tf.get_variable("s", initializer=sample_init)
    X, Y = inp[:, :, 0:1], inp[:, :, 1:2]
    return  1. / ( 1. + tf.exp( -theta * (.5*(Y-X) - tf.abs(sp - .5*(Y+X))) )  )

# Vectorizations implementing persistence entropy
def entropy_layer(inp, theta, num_samples, sample_init, sample_const, train_vect):
    """ Entropy PersLay
    WARNING: this function assumes that padding values are zero
    """
    bp_inp = tf.einsum("ijk,kl->ijl", inp, tf.constant(np.array([[1.,-1.],[0.,1.]], dtype=np.float32)))
    sp = tf.get_variable("s", shape=[1, 1, num_samples], initializer=sample_init, trainable=train_vect) if not sample_const else tf.get_variable("s", initializer=sample_init)
    L, X, Y = bp_inp[:, :, 1:2], bp_inp[:, :, 0:1], bp_inp[:, :, 0:1] + bp_inp[:, :, 1:2]
    LN = tf.multiply(L, 1. / tf.expand_dims(tf.matmul(L[:,:,0], tf.ones([L.shape[1],1])), -1))
    entropy_terms = tf.where(LN > 0., -tf.multiply(LN, tf.log(LN)), LN)
    return  tf.multiply(entropy_terms, 1. / ( 1. + tf.exp( -theta * (.5*(Y-X) - tf.abs(sp - .5*(Y+X))) )  ))

# Vectorizations implementing persistence images
def image_layer(inp, image_size, image_bnds, variance_init, variance_const, train_vect):
    """ Persistence Image PersLay """
    bp_inp = tf.einsum("ijk,kl->ijl", inp, tf.constant(np.array([[1.,-1.],[0.,1.]], dtype=np.float32)))
    dimension_before, num_pts = inp.shape[2].value, inp.shape[1].value
    coords = [tf.range(start=image_bnds[i][0], limit=image_bnds[i][1], delta=(image_bnds[i][1] - image_bnds[i][0]) / image_size[i]) for i in range(dimension_before)]
    M = tf.meshgrid(*coords)
    mu = tf.concat([tf.expand_dims(tens, 0) for tens in M], axis=0)
    sg = tf.get_variable("s", shape=[1], initializer=variance_init, trainable=train_vect) if not variance_const else tf.get_variable("s", initializer=variance_init)
    bc_inp = tf.reshape(bp_inp, [-1, num_pts, dimension_before] + [1 for _ in range(dimension_before)])
    return tf.exp(tf.reduce_sum(  -tf.square(bc_inp-mu) / (2*tf.square(sg[0])),  axis=2)) / (2*np.pi*tf.square(sg[0]))


def perslay_channel(output, name, diag, **kwargs):
    """ PersLay channel for persistence diagrams
        output :   list on which perslay output will be appended
        name :     name of the operation for tensorflow
        diag :     big matrix of shape [N_diag, N_pts_per_diag, dimension_diag (coordinates of points) + 1 (mask--0 or 1)]
    """

    try:
        train_weight = kwargs["train_weight"]
    except KeyError:
        train_weight = True

    try:
        train_vect = kwargs["train_vect"]
    except KeyError:
        train_vect = True

    N, dimension_diag = diag.get_shape()[1], diag.get_shape()[2]
    tensor_mask = diag[:, :, dimension_diag - 1]
    tensor_diag = diag[:, :, :dimension_diag - 1]

    if kwargs["persistence_weight"] == "linear":
        with tf.variable_scope(name + "-linear_pweight"):
            C = tf.get_variable("C", shape=[1], initializer=kwargs["coeff_init"], trainable=train_weight) if not kwargs["coeff_const"] else tf.get_variable("C", initializer=kwargs["coeff_init"])
            weight = C * tf.abs(tensor_diag[:, :, 1:2]-tensor_diag[:, :, 0:1])

    if kwargs["persistence_weight"] == "power":
        with tf.variable_scope(name + "-power_pweight"):
            p = kwargs["power_p"]
            C = tf.get_variable("C", shape=[1], initializer=kwargs["coeff_init"], trainable=train_weight) if not kwargs["coeff_const"] else tf.get_variable("C", initializer=kwargs["coeff_init"])
            weight = C * tf.pow(tf.abs(tensor_diag[:, :, 1:2]-tensor_diag[:, :, 0:1]), p)

    if kwargs["persistence_weight"] == "grid":
        with tf.variable_scope(name + "-grid_pweight"):
            W = tf.get_variable("W", shape=kwargs["grid_size"], initializer=kwargs["grid_init"], trainable=train_weight) if not kwargs["grid_const"] else tf.get_variable("W", initializer=kwargs["grid_init"])
            indices = []
            for dim in range(dimension_diag-1):
                [m, M] = kwargs["grid_bnds"][dim]
                coords = tf.slice(tensor_diag, [0, 0, dim], [-1, -1, 1])
                ids = kwargs["grid_size"][dim] * (coords - m)/(M - m)
                indices.append(tf.cast(ids, tf.int32))
            weight = tf.expand_dims(tf.gather_nd(params=W, indices=tf.concat(indices, axis=2)), -1)

    if kwargs["persistence_weight"] == "gmix":
        with tf.variable_scope(name + "-gmix_pweight"):
            M = tf.get_variable("M", shape=[1,1,2,kwargs["gmix_num"]], initializer=kwargs["gmix_m_init"], trainable=train_weight) if not kwargs["gmix_m_const"] else tf.get_variable("M", initializer=kwargs["gmix_m_init"])
            V = tf.get_variable("V", shape=[1,1,2,kwargs["gmix_num"]], initializer=kwargs["gmix_v_init"], trainable=train_weight) if not kwargs["gmix_v_const"] else tf.get_variable("V", initializer=kwargs["gmix_v_init"])
            bc_inp = tf.expand_dims(tensor_diag, -1)
            weight = tf.expand_dims(tf.reduce_sum(tf.exp(tf.reduce_sum(-tf.multiply(tf.square(bc_inp - M), tf.square(V)), axis=2)), axis=2), -1)

    # First layer of channel: processing of the persistence diagrams by vectorization of diagram points
    if kwargs["layer"] == "pm":  # Channel with permutation equivariant layers
        for idx, (dim, pop) in enumerate(kwargs["peq"]):
            with tf.variable_scope(name + "-perm_eq-" + str(idx)):
                tensor_diag = permutation_equivariant_layer(tensor_diag, dim, pop, kwargs["weight_init"], kwargs["weight_init"], kwargs["bias_init"], kwargs["weight_const"], kwargs["weight_const"], kwargs["bias_const"], train_vect)
    elif kwargs["layer"] == "ls":  # Channel with landscape layer
        with tf.variable_scope(name + "-samples"):
            tensor_diag = landscape_layer(tensor_diag, kwargs["num_samples"], kwargs["sample_init"], kwargs["sample_const"], train_vect)
    elif kwargs["layer"] == "bc":  # Channel with Betti layer
        with tf.variable_scope(name + "-samples"):
            tensor_diag = betti_layer(tensor_diag, kwargs["theta"], kwargs["num_samples"], kwargs["sample_init"], kwargs["sample_const"], train_vect)
    elif kwargs["layer"] == "en":  # Channel with entropy layer
        with tf.variable_scope(name + "-samples"):
            tensor_diag = entropy_layer(tensor_diag, kwargs["theta"], kwargs["num_samples"], kwargs["sample_init"], kwargs["sample_const"], train_vect)
    elif kwargs["layer"] == "im":  # Channel with image layer
        with tf.variable_scope(name + "-bandwidth"):
            tensor_diag = image_layer(tensor_diag, kwargs["image_size"], kwargs["image_bnds"], kwargs["variance_init"], kwargs["variance_const"], train_vect)
    elif kwargs["layer"] == "ex":  # Channel with exponential layer
        with tf.variable_scope(name + "-gaussians"):
            tensor_diag = exponential_layer(tensor_diag, kwargs["num_elements"], kwargs["mean_init"], kwargs["variance_init"], kwargs["mean_const"], kwargs["variance_const"], train_vect)
    elif kwargs["layer"] == "rt":  # Channel with rational layer
        with tf.variable_scope(name + "-bandwidth"):
            tensor_diag = rational_layer(tensor_diag, kwargs["num_elements"], kwargs["mean_init"], kwargs["variance_init"], kwargs["alpha_init"], kwargs["mean_const"], kwargs["variance_const"], kwargs["alpha_const"], train_vect)
    elif kwargs["layer"] == "rh":  # Channel with rational hat layer
        with tf.variable_scope(name + "-bandwidth"):
            tensor_diag = rational_hat_layer(tensor_diag, kwargs["num_elements"], kwargs["q"], kwargs["mean_init"], kwargs["r_init"], kwargs["mean_const"], kwargs["r_const"], train_vect)

    output_dim = len(tensor_diag.shape) - 2

    vector = None  # to avoid warning

    if output_dim == 1:
        # Apply weight and mask
        if kwargs["persistence_weight"] is not None:
            tiled_weight = tf.tile(weight, [1, 1, tensor_diag.shape[2].value])
            tensor_diag = tf.multiply(tensor_diag, tiled_weight)
        tiled_mask = tf.tile(tf.expand_dims(tensor_mask, -1), [1, 1, tensor_diag.shape[2].value])
        masked_layer = tf.multiply(tensor_diag, tiled_mask)

        # Permutation invariant operation
        if kwargs["perm_op"] == "topk":  # k first values
            masked_layer_t = tf.transpose(masked_layer, perm=[0, 2, 1])
            values, indices = tf.nn.top_k(masked_layer_t, k=kwargs["keep"])
            vector = tf.reshape(values, [-1, kwargs["keep"] * tensor_diag.shape[2].value])
        elif kwargs["perm_op"] == "sum":  # sum
            vector = tf.reduce_sum(masked_layer, axis=1)
        elif kwargs["perm_op"] == "max":  # maximum
            vector = tf.reduce_max(masked_layer, axis=1)
        elif kwargs["perm_op"] == "mean":  # minimum
            vector = tf.reduce_mean(masked_layer, axis=1)

        # Second layer of channel: fully-connected (None if fc_layers is set to [], default value)
        for idx, tup in enumerate(kwargs["fc_layers"]):
            # tup is a tuple whose element are
            # 1. dim of fully-connected,
            # 2. string for processing,
            # 3. (optional) dropout value
            with tf.variable_scope(name + "-fc-" + str(idx)):
                vector = tf.layers.dense(vector, tup[0])
            with tf.variable_scope(name + "-bn-" + str(idx)):
                if len(tup) == 2:
                    vector = _post_processing(vector, tup[1])
                else:
                    vector = _post_processing(vector, tup[1], tup[2])

    elif output_dim == 2:

        # Apply weight and mask
        if kwargs["persistence_weight"] is not None:
            weight = tf.expand_dims(weight, -1)
            tiled_weight = tf.tile(weight, [1, 1, tensor_diag.shape[2].value, tensor_diag.shape[3].value])
            tensor_diag = tf.multiply(tensor_diag, tiled_weight)
        tiled_mask = tf.tile(tf.reshape(tensor_mask, [-1, N, 1, 1]), [1, 1, tensor_diag.shape[2].value, tensor_diag.shape[3].value])
        masked_layer = tf.multiply(tensor_diag, tiled_mask)

        # Permutation invariant operation
        if kwargs["perm_op"] == "sum":  # sum
            vector = tf.reduce_sum(masked_layer, axis=1)
        elif kwargs["perm_op"] == "max":  # maximum
            vector = tf.reduce_max(masked_layer, axis=1)
        elif kwargs["perm_op"] == "mean":  # minimum
            vector = tf.reduce_mean(masked_layer, axis=1)

        # Second layer of channel: convolution
        vector = tf.expand_dims(vector, -1)
        for idx, tup in enumerate(kwargs["cv_layers"]):
            # tup is a tuple whose element are
            # 1. num of filters,
            # 2. kernel size,
            # 3. string for postprocessing,
            # 4. (optional) dropout value
            with tf.variable_scope(name + "-cv-" + str(idx)):
                vector = tf.layers.conv2d(vector, filters=tup[0], kernel_size=tup[1])
            with tf.variable_scope(name + "-bn-" + str(idx)):
                if len(tup) == 3:
                    vector = _post_processing(vector, tup[2])
                else:
                    vector = _post_processing(vector, tup[2], tup[3])
        vector = tf.layers.flatten(vector)

    output.append(vector)
    return vector
