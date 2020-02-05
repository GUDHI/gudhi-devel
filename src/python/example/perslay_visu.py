#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
sys.path.append("expe/")
import numpy             as np
import perslay           as ps
import tensorflow        as tf
import matplotlib.pyplot as plt
from perslay.perslay            import perslay_channel
from perslay.preprocessing      import DiagramScaler, Padding
from sklearn.preprocessing      import MinMaxScaler
from tensorflow                 import random_uniform_initializer as rui


# # Input persistence diagram and tensorflow graph

# ## Plot and preprocess persistence diagram

# In[ ]:


diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]


# In[ ]:


plt.scatter(diag[0][:,0], diag[0][:,1])
plt.plot([0.,6.],[0.,6.])
plt.show()


# In[ ]:


diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
diag = Padding(use=True).fit_transform(diag)


# In[ ]:


plt.scatter(diag[0][:,0], diag[0][:,1])
plt.plot([0.,1.],[0.,1.])
plt.show()


# In[ ]:


D = np.stack(diag, 0)


# ## Initialize tensorflow graph

# In[ ]:


tf.reset_default_graph()


# In[ ]:


diagram = tf.placeholder(tf.float32, shape=D.shape)
feed    = {diagram: D}


# In[ ]:


perslayParameters = {}


# # Persistence weight and permutation invariant operation

# ## Choose persistence weights

# In[ ]:


perslayParameters["persistence_weight"]  = "linear"
perslayParameters["coeff_init"]          = rui(1.0, 1.0)
perslayParameters["coeff_const"]         = False


# In[ ]:


perslayParameters["persistence_weight"]  = "grid"
perslayParameters["grid_size"]           = [100,100]
epsilon = .001
perslayParameters["grid_bnds"]           = ((0. - epsilon, 1. + epsilon), (0. - epsilon, 1. + epsilon))
perslayParameters["grid_init"]           = np.tile(np.arange(0.,100.,1, dtype=np.float32)[np.newaxis,:], [100,1])
perslayParameters["grid_const"]          = True


# In[ ]:


perslayParameters["persistence_weight"]  = "gmix"
perslayParameters["gmix_num"]            = 3
perslayParameters["gmix_m_init"]         = rui(0., 1.)
perslayParameters["gmix_m_const"]        = False
perslayParameters["gmix_v_init"]         = rui(5, 5)
perslayParameters["gmix_v_const"]        = False


# In[ ]:


perslayParameters["persistence_weight"]  = None


# ## Choose permutation invariant operation

# In[ ]:


perslayParameters["perm_op"] = "sum"


# In[ ]:


perslayParameters["perm_op"] = "topk"
perslayParameters["keep"]    = 3


# In[ ]:


perslayParameters["perm_op"] = "max"


# In[ ]:


perslayParameters["perm_op"] = "mean"


# # Persistence representation

# ## Persistence image

# In[ ]:


perslayParameters["layer"]          = "im"
perslayParameters["image_size"]     = (100, 100)
epsilon = .001
perslayParameters["image_bnds"]     = ((-.5-epsilon, 1.5+epsilon), (-.5-epsilon, 1.5+epsilon))
perslayParameters["variance_init"]  = rui(10., 10.) 
perslayParameters["variance_const"] = False
perslayParameters["cv_layers"]      = []


# In[ ]:


list_v = []
ps.perslay.perslay_channel(output=list_v, name="perslay", diag=diagram, **perslayParameters)
vector = tf.concat(list_v, 1)


# In[ ]:


init = tf.global_variables_initializer()
with tf.Session() as sess:
    sess.run(init)
    
    # Plot representation
    V = vector.eval(feed_dict=feed)[0,:]
    V = np.flip(np.reshape(V, [int(np.sqrt(V.shape[0])), int(np.sqrt(V.shape[0]))]), 0)
    plt.figure()
    plt.imshow(V, cmap="Purples")
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    plt.show()
    
    # Plot weight
    if perslayParameters["persistence_weight"] == "grid":
        W = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,"perslay-grid_pweight/W")[0]).T
        weights = np.flip(W, 0)
        plt.figure()
        plt.imshow(weights, cmap="Purples", zorder=1)
        ((xm,xM),(ym,yM)) = perslayParameters["grid_bnds"]
        [xs, ys] = perslayParameters["grid_size"]
        plt.scatter([int(xs*(x-xm)/(xM-xm)) for x in diag[0][:,0]], 
                    [ys-int(ys*(y-ym)/(yM-ym)) for y in diag[0][:,1]], 
                    s=10, color="red", zorder=2)
        plt.show()
        
    if perslayParameters["persistence_weight"] == "gmix":
        means = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "perslay-gmix_pweight/M")[0])
        varis = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "perslay-gmix_pweight/V")[0])
        x = np.arange(-.5, 1.5, .001)
        y = np.arange(-.5, 1.5, .001)
        xx, yy = np.meshgrid(x, y)
        z = np.zeros(xx.shape)
        for idx_g in range(means.shape[3]):
            z += np.exp(-((xx-means[0,0,0,idx_g])**2 * (varis[0,0,0,idx_g])**2 
                        + (yy-means[0,0,1,idx_g])**2 * (varis[0,0,1,idx_g])**2 ))
        plt.contourf(xx, yy, z)
        plt.scatter(diag[0][:,0], diag[0][:,1], s=50, color="red")
        plt.show()


# ## Persistence landscape / entropy / Betti curve

# In[ ]:


perslayParameters["layer"]          = "ls"
#perslayParameters["layer"]          = "bc"
#perslayParameters["layer"]          = "en"


# In[ ]:


perslayParameters["num_samples"]    = 3000
perslayParameters["sample_init"]    = np.array([[ np.arange(-1.,2.,.001) ]], dtype=np.float32)
perslayParameters["sample_const"]   = True
perslayParameters["theta"]          = 100 # used only if layer is "bc" or "en"
perslayParameters["fc_layers"]      = []


# In[ ]:


list_v = []
perslay_channel(output=list_v, name="perslay", diag=diagram, **perslayParameters)
vector = tf.concat(list_v, 1)


# In[ ]:


init = tf.global_variables_initializer()
with tf.Session() as sess:
    sess.run(init)
    
    #Plot representation
    V = vector.eval(feed_dict=feed)[0,:]
    plt.figure()
    if perslayParameters["perm_op"] == "topk":
        V = np.reshape(V, [-1, perslayParameters["keep"]])
        for k in range(perslayParameters["keep"]):
            plt.plot(V[:,k], linewidth=5.0)
    else:
        plt.plot(V, linewidth=5.0)
    plt.show()
    
    # Plot weight
    if perslayParameters["persistence_weight"] == "grid":
        W = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES,"perslay-grid_pweight/W")[0]).T
        weights = np.flip(W, 0)
        plt.figure()
        plt.imshow(weights, cmap="Purples", zorder=1)
        ((xm,xM),(ym,yM)) = perslayParameters["grid_bnds"]
        [xs, ys] = perslayParameters["grid_size"]
        plt.scatter([int(xs*(x-xm)/(xM-xm)) for x in diag[0][:,0]], 
                    [ys-int(ys*(y-ym)/(yM-ym)) for y in diag[0][:,1]], 
                    s=10, color="red", zorder=2)
        plt.show()
        
    if perslayParameters["persistence_weight"] == "gmix":
        means = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "perslay-gmix_pweight/M")[0])
        varis = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "perslay-gmix_pweight/V")[0])
        x = np.arange(-.5, 1.5, .001)
        y = np.arange(-.5, 1.5, .001)
        xx, yy = np.meshgrid(x, y)
        z = np.zeros(xx.shape)
        for idx_g in range(means.shape[3]):
            z += np.exp(-((xx-means[0,0,0,idx_g])**2 * (varis[0,0,0,idx_g])**2 
                        + (yy-means[0,0,1,idx_g])**2 * (varis[0,0,1,idx_g])**2 ))
        plt.contourf(xx, yy, z)
        plt.scatter(diag[0][:,0], diag[0][:,1], s=100)
        plt.show()


# In[ ]:




