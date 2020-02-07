import numpy             as np
import tensorflow        as tf
import matplotlib.pyplot as plt

from gudhi.perslay              import perslay_channel
from gudhi.representations      import DiagramScaler, Padding
from sklearn.preprocessing      import MinMaxScaler
from tensorflow                 import random_uniform_initializer as rui

diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
diag = Padding(use=True).fit_transform(diag)
plt.scatter(diag[0][:,0], diag[0][:,1])
plt.plot([0.,1.],[0.,1.])
plt.show()

D = np.stack(diag, 0)


# Initialize tensorflow graph
tf.reset_default_graph()

diagram = tf.placeholder(tf.float32, shape=D.shape)
feed    = {diagram: D}

perslayParameters = {}

#####################
# Persistence image #
#####################

epsilon = .1
perslayParameters["persistence_weight"]  = "grid"
perslayParameters["grid_size"]           = (100,100)
perslayParameters["grid_bnds"]           = ((0.-epsilon, 1.+epsilon), (0.-epsilon, 1.+epsilon))
perslayParameters["grid_init"]           = np.array(np.random.uniform(size=[100,100]), dtype=np.float32)
perslayParameters["grid_const"]          = True
perslayParameters["perm_op"]             = "sum"
perslayParameters["layer"]               = "im"
perslayParameters["image_size"]          = (100, 100)
perslayParameters["image_bnds"]          = ((0.-epsilon, 1.+epsilon), (0.-epsilon, 1.+epsilon))
perslayParameters["variance_init"]       = rui(.1, .1) 
perslayParameters["variance_const"]      = False
perslayParameters["cv_layers"]           = []

list_v = []
perslay_channel(output=list_v, name="perslay", diag=diagram, **perslayParameters)
vector = tf.concat(list_v, 1)

init = tf.global_variables_initializer()
with tf.Session() as sess:
    sess.run(init)
    
    # Plot representation
    V = vector.eval(feed_dict=feed)[0,:]
    V = np.flip(np.reshape(V, [int(np.sqrt(V.shape[0])), int(np.sqrt(V.shape[0]))]), 0)
    plt.figure()
    plt.imshow(V, cmap="Purples")
    cb = plt.colorbar()
    plt.show()
    
    # Plot weight
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

#################################################
# Persistence landscape / entropy / Betti curve #
#################################################

perslayParameters["persistence_weight"]  = "gmix"
perslayParameters["gmix_num"]            = 3
perslayParameters["gmix_m_init"]         = np.array(np.random.uniform(size=[2,3]), dtype=np.float32)
perslayParameters["gmix_m_const"]        = True
perslayParameters["gmix_v_init"]         = np.array(10 * np.ones([2,3]), dtype=np.float32)
perslayParameters["gmix_v_const"]        = True
perslayParameters["perm_op"]             = "topk"
perslayParameters["keep"]                = 3
perslayParameters["layer"]               = "ls"
perslayParameters["num_samples"]         = 3000
perslayParameters["sample_init"]         = np.array([[ np.arange(-1.,2.,.001) ]], dtype=np.float32)
perslayParameters["sample_const"]        = True
perslayParameters["fc_layers"]           = []

list_v = []
perslay_channel(output=list_v, name="perslay", diag=diagram, **perslayParameters)
vector = tf.concat(list_v, 1)

init = tf.global_variables_initializer()
with tf.Session() as sess:
    sess.run(init)
    
    #Plot representation
    V = vector.eval(feed_dict=feed)[0,:]
    plt.figure()
    V = np.reshape(V, [-1, perslayParameters["keep"]])
    for k in range(perslayParameters["keep"]):
        plt.plot(V[:,k], linewidth=5.0)
    plt.show()
    
    # Plot weight
    means = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "perslay-gmix_pweight/M")[0])
    varis = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "perslay-gmix_pweight/V")[0])
    x = np.arange(-.1, 1.1, .001)
    y = np.arange(-.1, 1.1, .001)
    xx, yy = np.meshgrid(x, y)
    z = np.zeros(xx.shape)
    for idx_g in range(means.shape[3]):
        z += np.exp(-((xx-means[0,0,0,idx_g])**2 * (varis[0,0,0,idx_g])**2 
                    + (yy-means[0,0,1,idx_g])**2 * (varis[0,0,1,idx_g])**2 ))
    plt.contourf(xx, yy, z)
    plt.scatter(diag[0][:,0], diag[0][:,1], s=100)
    plt.show()
