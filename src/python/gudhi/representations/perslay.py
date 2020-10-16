import tensorflow as tf
import tensorflow_addons as tfa
from tensorflow import random_uniform_initializer as rui
import numpy as np

def permutation_equivariant_layer(inp, dimension, perm_op, lbda, b, gamma):
    """ DeepSet PersLay """
    dimension_before, num_pts = inp.shape[2], inp.shape[1]
    b = tf.expand_dims(tf.expand_dims(b, 0), 0)
    A = tf.reshape(tf.einsum("ijk,kl->ijl", inp, lbda), [-1, num_pts, dimension])
    if perm_op != None:
        if perm_op == "max":
            beta = tf.tile(tf.expand_dims(tf.math.reduce_max(inp, axis=1), 1), [1, num_pts, 1])
        elif perm_op == "min":
            beta = tf.tile(tf.expand_dims(tf.math.reduce_min(inp, axis=1), 1), [1, num_pts, 1])
        elif perm_op == "sum":
            beta = tf.tile(tf.expand_dims(tf.math.reduce_sum(inp, axis=1), 1), [1, num_pts, 1])
        else:
            raise Exception("perm_op should be min, max or sum")
        B = tf.reshape(tf.einsum("ijk,kl->ijl", beta, gamma), [-1, num_pts, dimension])
        return A - B + b
    else:
        return A + b

def rational_hat_layer(inp, q, mu, r):
    """ Rational Hat PersLay """
    mu, r = tf.expand_dims(tf.expand_dims(mu, 0), 0), tf.expand_dims(tf.expand_dims(r, 0), 0)
    dimension_before, num_pts = inp.shape[2], inp.shape[1]
    bc_inp = tf.expand_dims(inp, -1)
    norms = tf.norm(bc_inp - mu, ord=q, axis=2)
    return 1/(1 + norms) - 1/(1 + tf.math.abs(tf.math.abs(r)-norms)) 

def rational_layer(inp, mu, sg, al):
    """ Rational PersLay """
    mu, sg, al = tf.expand_dims(tf.expand_dims(mu, 0), 0), tf.expand_dims(tf.expand_dims(sg, 0), 0), tf.expand_dims(tf.expand_dims(al, 0), 0)
    dimension_before, num_pts = inp.shape[2], inp.shape[1]
    bc_inp = tf.expand_dims(inp, -1)
    return 1/tf.math.pow(1+tf.math.reduce_sum(tf.math.multiply(tf.math.abs(bc_inp - mu), tf.math.abs(sg)), axis=2), al)

def exponential_layer(inp, mu, sg):
    """ Exponential PersLay """
    mu, sg = tf.expand_dims(tf.expand_dims(mu, 0), 0), tf.expand_dims(tf.expand_dims(sg, 0), 0)
    dimension_before, num_pts = inp.shape[2], inp.shape[1]
    bc_inp = tf.expand_dims(inp, -1)
    return tf.math.exp(tf.math.reduce_sum(-tf.math.multiply(tf.math.square(bc_inp - mu), tf.math.square(sg)), axis=2))

def landscape_layer(inp, sp):
    """ Landscape PersLay """
    sp = tf.expand_dims(tf.expand_dims(sp, 0), 0)
    return tf.math.maximum( .5 * (inp[:, :, 1:2] - inp[:, :, 0:1]) - tf.math.abs(sp - .5 * (inp[:, :, 1:2] + inp[:, :, 0:1])), np.array([0]))

def betti_layer(inp, theta, sp):
    """ Betti PersLay """
    sp = tf.expand_dims(tf.expand_dims(sp, 0), 0)
    X, Y = inp[:, :, 0:1], inp[:, :, 1:2]
    return  1. / ( 1. + tf.math.exp( -theta * (.5*(Y-X) - tf.math.abs(sp - .5*(Y+X))) )  )

def entropy_layer(inp, theta, sp):
    """ Entropy PersLay
    WARNING: this function assumes that padding values are zero
    """
    sp = tf.expand_dims(tf.expand_dims(sp, 0), 0)
    bp_inp = tf.einsum("ijk,kl->ijl", inp, tf.constant(np.array([[1.,-1.],[0.,1.]], dtype=np.float32)))
    L, X, Y = bp_inp[:, :, 1:2], bp_inp[:, :, 0:1], bp_inp[:, :, 0:1] + bp_inp[:, :, 1:2]
    LN = tf.math.multiply(L, 1. / tf.expand_dims(tf.linalg.matmul(L[:,:,0], tf.ones([L.shape[1],1])), -1))
    entropy_terms = tf.where(LN > 0., -tf.math.multiply(LN, tf.math.log(LN)), LN)
    return  tf.math.multiply(entropy_terms, 1. / ( 1. + tf.math.exp( -theta * (.5*(Y-X) - tf.math.abs(sp - .5*(Y+X))) )  ))

def image_layer(inp, image_size, image_bnds, sg):
    """ Persistence Image PersLay """
    bp_inp = tf.einsum("ijk,kl->ijl", inp, tf.constant(np.array([[1.,-1.],[0.,1.]], dtype=np.float32)))
    dimension_before, num_pts = inp.shape[2], inp.shape[1]
    coords = [tf.range(start=image_bnds[i][0], limit=image_bnds[i][1], delta=(image_bnds[i][1] - image_bnds[i][0]) / image_size[i]) for i in range(dimension_before)]
    M = tf.meshgrid(*coords)
    mu = tf.concat([tf.expand_dims(tens, 0) for tens in M], axis=0)
    bc_inp = tf.reshape(bp_inp, [-1, num_pts, dimension_before] + [1 for _ in range(dimension_before)])
    return tf.expand_dims(tf.math.exp(tf.math.reduce_sum(  -tf.math.square(bc_inp-mu) / (2*tf.math.square(sg)),  axis=2)) / (2*np.pi*tf.math.square(sg)), -1)

class PerslayModel(tf.keras.Model):

    def __init__(self, name, diagdim, perslay_parameters, rho):
        super(PerslayModel, self).__init__()
        self.namemodel            = name
        self.diagdim              = diagdim
        self.perslay_parameters   = perslay_parameters
        self.rho                  = rho

        self.vars = [[] for _ in range(len(self.perslay_parameters))]
        for nf, plp in enumerate(self.perslay_parameters):

            weight = plp["pweight"]
            if weight != None:
                Winit, Wtrain, Wname = plp["pweight_init"], plp["pweight_train"], self.namemodel + "-pweight-" + str(nf)
                if not callable(Winit):
                    W = tf.Variable(name=Wname, initial_value=Winit, trainable=Wtrain)
                else:
                    if weight == "power":
                        W = tf.Variable(name=Wname, initial_value=Winit([1]), trainable=Wtrain)
                    elif weight == "grid":
                        Wshape = plp["pweight_size"]
                        W = tf.Variable(name=Wname, initial_value=Winit(Wshape), trainable=Wtrain)
                    elif weight == "gmix":
                        ngs = plp["pweight_num"] 
                        W = tf.Variable(name=Wname, initial_value=Winit([4,ngs]), trainable=Wtrain)
            else:
                W = 0
            self.vars[nf].append(W)

            layer, Ltrain, Lname = plp["layer"], plp["layer_train"], self.namemodel + "-" + str(nf)

            if layer == "PermutationEquivariant":
                Lpeq, LWinit, LBinit, LGinit = plp["lpeq"], plp["lweight_init"], plp["lbias_init"], plp["lgamma_init"]
                LW, LB, LG = [], [], []
                for idx, (dim, pop) in enumerate(Lpeq):
                    dim_before = self.diagdim if idx == 0 else Lpeq[idx-1][0]
                    LWiv = LWinit([dim_before, dim]) if callable(LWinit) else LWinit
                    LBiv = LBinit([dim]) if callable(LBinit) else LBinit
                    LW.append(      tf.Variable(name=Lname+"-W", initial_value=LWiv, trainable=Ltrain))
                    LB.append(      tf.Variable(name=Lname+"-B", initial_value=LBiv, trainable=Ltrain))
                    if pop != None:
                        LGiv = LGinit([dim_before, dim]) if callable(LGinit) else LGinit
                        LG.append(  tf.Variable(name=Lname+"-G", initial_value=LGiv, trainable=Ltrain))
                    else:
                        LG.append([])
                self.vars[nf].append([LW, LB, LG])

            elif layer == "Landscape" or layer == "BettiCurve" or layer == "Entropy":
                LSinit = plp["lsample_init"]
                LSiv = LSinit if not callable(LSinit) else LSinit([plp["lsample_num"]]) 
                LS = tf.Variable(name=Lname+"-S", initial_value=LSiv, trainable=Ltrain)
                self.vars[nf].append(LS)

            elif layer == "Image":
                LVinit = plp["lvariance_init"]
                LViv = LVinit if not callable(LVinit) else LVinit([1])
                LV = tf.Variable(name=Lname+"-V", initial_value=LViv, trainable=Ltrain)
                self.vars[nf].append(LV)

            elif layer == "Exponential":
                LMinit, LVinit = plp["lmean_init"], plp["lvariance_init"]
                LMiv = LMinit if not callable(LMinit) else LMinit([self.diagdim, plp["lnum"]])
                LViv = LVinit if not callable(LVinit) else LVinit([self.diagdim, plp["lnum"]])
                LM = tf.Variable(name=Lname+"-M", initial_value=LMiv, trainable=Ltrain)
                LV = tf.Variable(name=Lname+"-V", initial_value=LViv, trainable=Ltrain)
                self.vars[nf].append([LM, LV])

            elif layer == "Rational":
                LMinit, LVinit, LAinit = plp["lmean_init"], plp["lvariance_init"], plp["lalpha_init"]
                LMiv = LMinit if not callable(LMinit) else LMinit([self.diagdim, plp["lnum"]])
                LViv = LVinit if not callable(LVinit) else LVinit([self.diagdim, plp["lnum"]])
                LAiv = LAinit if not callable(LAinit) else LAinit([plp["lnum"]])
                LM = tf.Variable(name=Lname+"-M", initial_value=LMiv, trainable=Ltrain)
                LV = tf.Variable(name=Lname+"-V", initial_value=LViv, trainable=Ltrain)
                LA = tf.Variable(name=Lname+"-A", initial_value=LAiv, trainable=Ltrain)
                self.vars[nf].append([LM, LV, LA])

            elif layer == "RationalHat":
                LMinit, LRinit = plp["lmean_init"], plp["lr_init"]
                LMiv = LMinit if not callable(LMinit) else LMinit([self.diagdim, plp["lnum"]])
                LRiv = LRinit if not callable(LRinit) else LVinit([1])
                LM = tf.Variable(name=Lname+"-M", initial_value=LMiv, trainable=Ltrain)
                LR = tf.Variable(name=Lname+"-R", initial_value=LRiv, trainable=Ltrain)
                self.vars[nf].append([LM, LR])

    def compute_representations(self, diags, training=False):

        list_v = []

        for nf, plp in enumerate(self.perslay_parameters):

            diag = diags[nf]

            N, dimension_diag = diag.shape[1], diag.shape[2]
            tensor_mask = diag[:, :, dimension_diag - 1]
            tensor_diag = diag[:, :, :dimension_diag - 1]

            W = self.vars[nf][0]

            if plp["pweight"] == "power":
                p = plp["pweight_power"]
                weight = W * tf.math.pow(tf.math.abs(tensor_diag[:, :, 1:2]-tensor_diag[:, :, 0:1]), p)

            elif plp["pweight"] == "grid":
                grid_shape = W.shape  
                indices = []
                for dim in range(dimension_diag-1):
                    [m, M] = plp["pweight_bnds"][dim]
                    coords = tf.slice(tensor_diag, [0, 0, dim], [-1, -1, 1])
                    ids = grid_shape[dim] * (coords - m)/(M - m)
                    indices.append(tf.cast(ids, tf.int32))
                weight = tf.expand_dims(tf.gather_nd(params=W, indices=tf.concat(indices, axis=2)), -1)

            elif plp["pweight"] == "gmix":
                M, V = tf.expand_dims(tf.expand_dims(W[:2,:], 0), 0), tf.expand_dims(tf.expand_dims(W[2:,:], 0), 0) 
                bc_inp = tf.expand_dims(tensor_diag, -1)
                weight = tf.expand_dims(tf.math.reduce_sum(tf.math.exp(tf.math.reduce_sum(-tf.math.multiply(tf.math.square(bc_inp-M), tf.math.square(V)), axis=2)), axis=2), -1)


            lvars = self.vars[nf][1]
            if plp["layer"] == "PermutationEquivariant":
                for idx, (dim, pop) in enumerate(plp["lpeq"]):
                    tensor_diag = permutation_equivariant_layer(tensor_diag, dim, pop, lvars[0][idx], lvars[1][idx], lvars[2][idx])
            elif plp["layer"] == "Landscape":
                tensor_diag = landscape_layer(tensor_diag, lvars)
            elif plp["layer"] == "BettiCurve":
                tensor_diag = betti_layer(tensor_diag, plp["theta"], lvars)
            elif plp["layer"] == "Entropy":
                tensor_diag = entropy_layer(tensor_diag, plp["theta"], lvars)
            elif plp["layer"] == "Image":
                tensor_diag = image_layer(tensor_diag, plp["image_size"], plp["image_bnds"], lvars)
            elif plp["layer"] == "Exponential":
                tensor_diag = exponential_layer(tensor_diag, **lvars)
            elif plp["layer"] == "Rational":
                tensor_diag = rational_layer(tensor_diag, **lvars)
            elif plp["layer"] == "RationalHat":
                tensor_diag = rational_hat_layer(tensor_diag, plp["q"], **lvars)

            # Apply weight
            output_dim = len(tensor_diag.shape) - 2
            if plp["pweight"] != None:
                for _ in range(output_dim-1):
                    weight = tf.expand_dims(weight, -1)
                tiled_weight = tf.tile(weight, [1, 1] + tensor_diag.shape[2:])
                tensor_diag = tf.math.multiply(tensor_diag, tiled_weight)

            # Apply mask
            for _ in range(output_dim):
                tensor_mask = tf.expand_dims(tensor_mask, -1)
            tiled_mask = tf.tile(tensor_mask, [1, 1] + tensor_diag.shape[2:])
            masked_layer = tf.math.multiply(tensor_diag, tiled_mask)

            # Permutation invariant operation
            if plp["perm_op"] == "topk" and output_dim == 1:  # k first values
                masked_layer_t = tf.transpose(masked_layer, perm=[0, 2, 1])
                values, indices = tf.math.top_k(masked_layer_t, k=plp["keep"])
                vector = tf.reshape(values, [-1, plp["keep"] * tensor_diag.shape[2]])
            elif plp["perm_op"] == "sum":  # sum
                vector = tf.math.reduce_sum(masked_layer, axis=1)
            elif plp["perm_op"] == "max":  # maximum
                vector = tf.math.reduce_max(masked_layer, axis=1)
            elif plp["perm_op"] == "mean":  # minimum
                vector = tf.math.reduce_mean(masked_layer, axis=1)

            # Second layer of channel
            vector = plp["final_model"].call(vector, training=training) if plp["final_model"] != "identity" else vector
            list_v.append(vector)

        # Concatenate all channels and add other features
        representations = tf.concat(values=list_v, axis=1)
        return representations

    def call(self, inputs, training=False):

        diags, feats = inputs[0], inputs[1]
        representations = self.compute_representations(diags, training)
        concat_representations = tf.concat(values=[representations, feats], axis=1)
        final_representations = self.rho(concat_representations) if self.rho != "identity" else concat_representations

        return final_representations

