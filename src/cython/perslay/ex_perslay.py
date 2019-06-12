from tda_deep_set import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

diag = {}
diag["test"] = [ np.array([[0.,1.],[.5,1.5],[2.,4.]]), np.array([[0.,.5],[1.5,2.5],[3.,5.],[0.,3.]]) ]
num_pts = len(diag["test"])

filtration   = "test"
idx_diagram  = 1

diag_example = tda.DiagramSelector(limit=np.inf).fit_transform(diag[filtration])
pre = tda.DiagramPreprocessor(use=True, scalers=[([0,1], MinMaxScaler())]).fit(diag_example)
[mx,my],[Mx,My] = pre.scalers[0][1].data_min_, pre.scalers[0][1].data_max_
#print("Minimum x = " + str(mx) + ", Maximum x = " + str(Mx) + ", Minimum y = " + str(my) + ", Maximum y = " + str(My))

xs, ys = diag_example[idx_diagram][:,0], diag_example[idx_diagram][:,1]
#plt.scatter(xs,ys)
#plt.plot([min(xs)-.5,max(ys)+.5],[min(xs)-.5,max(ys)+.5])
#plt.axis([min(xs),max(xs),  min(ys),max(ys)])
#plt.gca().set_aspect('equal', adjustable='box')
#plt.axis("equal")
#plt.axis("off")
#plt.show()
#plt.savefig("examplePD", format="pdf")

# Remark: no need to consider finite points in homological dimension 1 for graphs
prm = {"test":  {"ProminentPts__num_pts":10}}
num_filt = len(prm.keys())

scalers3d = [([0,1],  Pipeline([("1",tda.BirthPersistenceTransform()),("2",MinMaxScaler())])), ([2], Pipeline([("1",log_transform()), ("2",MinMaxScaler())]))]
scalers2d = [([0,1],  Pipeline([("1",tda.BirthPersistenceTransform()),("2",MinMaxScaler())]))]

preprocess = Pipeline([
    ("ProminentPts",  tda.ProminentPoints(use=True, point_type="finite")),
    ("Scaler",        tda.DiagramPreprocessor(use=True,  scalers=scalers2d)),
    ("NuSeparator",   tda.DiagramPreprocessor(use=False, scalers=[([0,1],nu_separator(nu=.1))])),
    ("Padding",       tda.Padding(use=True)),
                      ])

D = []
for dt in prm.keys():
    param = prm[dt]
    preprocess.set_params(**param)
    D.append(preprocess.fit_transform(diag[dt]))

filtration   = 0
idx_diagram  = 1

diag_example = D[filtration]
pre = tda.DiagramPreprocessor(use=True, scalers=[([0,1,2], MinMaxScaler())]).fit(diag_example)
[mx,my,mz],[Mx,My,Mz] = pre.scalers[0][1].data_min_, pre.scalers[0][1].data_max_
#print("Minimum x = " + str(mx) + ", Maximum x = " + str(Mx) + ", Minimum y = " + str(my) + ", Maximum y = " + str(My))

xs, ys = diag_example[idx_diagram][:,0], diag_example[idx_diagram][:,1]
#plt.scatter(xs,ys)
#plt.axis([min(xs)-.5,max(xs)+.5,min(ys)-.5,max(ys)+.5])
#plt.gca().set_aspect('equal', adjustable='box')
#plt.axis("off")
#plt.show()
#plt.savefig("examplePDpro", format="pdf")

D_pad = []
for dt in range(num_filt):
    D_pad.append(np.concatenate([D[dt][i][np.newaxis,:] for i in range(len(D[dt]))], axis=0))
    #print(D_pad[dt].shape)

layer               = "pm"
idx_diag            = 1

# Specific of "pm"
peq                 = [(2,None)]
weight_init         = rui(0.,1.)
weight_const        = False
bias_init           = rui(0.,1.)
bias_const          = False

#Specific of "la"
laguerre_weight     = [np.array([[0,0,0,0],[0,0,0,0]], dtype=np.float32)]
laguerre_const      = 10
laguerre_coords     = True

# Specific of "gs"
num_gaussians       = 500
mean_init           = rui(0.,1.)
mean_const          = True

# Specific of "im"
image_size          = [1000,1000]
image_bnds          = [[-1.,2.],[-1.,2.]]

# Used for "gs" and "im"
variance_init       = rui(5.,5.)
variance_const      = False

# Used for "ls" and "la"
num_samples         = 3000
sample_init     = np.array([[ np.arange(-1.,2.,.001) ]], dtype=np.float32) # specific example for "ls"
if layer == "la":
    X,Y             = np.meshgrid(np.linspace(-1,2,1000), np.linspace(-1,2,1000))
    pos             = np.reshape(np.transpose(np.stack([X,Y], axis=0), (1,2,0)), [num_samples,2]).T
    sample_init     = np.array([[ pos ]], dtype=np.float32) # specific example for "la"
sample_const        = True

perm_op             = "sum"
keep                = 3

persistence_weight  = None

# Specific of "grid"
grid_size           = [100,100]
grid_bnds           = [[-1.,2.],[-1.,2.]]
grid_init           = np.tile(np.arange(0.,100.,1, dtype=np.float32)[np.newaxis,:], [100,1])
#grid_init           = np.array(np.random.rand(100,100), dtype=np.float32)
grid_const          = True

#Specific of "linear"
coeff_init          = rui(1.,1.)
coeff_const         = False

# Reset graph
tf.reset_default_graph()

# Neural network input
indxs =   tf.placeholder(shape=[None,1], dtype=tf.int32)
diags = [ tf.placeholder(shape=[None, D_pad[dt].shape[1], D_pad[dt].shape[2]], dtype=tf.float32) for dt in range(num_filt) ]

# Compute vectorization

#ls_diag = landscape_layer(inp=tf.convert_to_tensor(diags[0], dtype=tf.float32), num_samples=num_samples, s_init=sample_init, s_const=sample_const)
im_diag = image_layer(inp=tf.convert_to_tensor(diags[0], dtype=tf.float32)[:,:,:-1], im_size=image_size, im_bnds=image_bnds, s_init=variance_init, s_const=variance_const)

list_v = []
vect = deep_set_channel(list_v, "fin-0", diags[0], layer=layer, peq=peq, num_samples=num_samples, num_gaussians=num_gaussians, 
        mean_init=mean_init, variance_init=variance_init, mean_const=mean_const, variance_const=variance_const,
        weight_init=weight_init, bias_init=bias_init, weight_const=weight_const, bias_const=bias_const, coeff_init=coeff_init, coeff_const=coeff_const,
        sample_init=sample_init, sample_const=sample_const, image_size=image_size, image_bnds=image_bnds, 
        laguerre_weight=tf.gather_nd(laguerre_weight[0], indxs), laguerre_const=laguerre_const, laguerre_coords=laguerre_coords,
        perm_op=perm_op, keep=keep, persistence_weight=persistence_weight, grid_size=grid_size, grid_bnds=grid_bnds, grid_init=grid_init, grid_const=grid_const)

init = tf.global_variables_initializer()
with tf.Session() as sess:
        
    sess.run(init)
    V = vect.eval(feed_dict={diags[0]: D_pad[0], indxs: np.array([[idx_diag]])})

    plt.rcParams["axes.grid"] = False

    if persistence_weight == "grid":
        weight_func = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "fin-0-grid_pweight/W")[0])
        #plt.imshow(np.flip(weight_func.T,0), cmap="Purples")
        #plt.colorbar()
        #plt.axis("off")
        #plt.savefig("exampleWE", format="pdf")
        #plt.show()
    
    if layer == "im":
        IM = im_diag.eval(feed_dict={diags[0]: D_pad[0], indxs: np.array([[idx_diag]])})
        
        #fig = plt.figure()
        #ax = fig.gca(projection='3d')

        #X = np.arange(-1., 2., .01)
        #Y = np.arange(-1., 2., .01)
        #idxx = np.squeeze(np.argwhere(abs(X-D_pad[0][idx_diag,0,0]) <= 0.5))
        #idxy = np.squeeze(np.argwhere(abs(Y-D_pad[0][idx_diag,0,1]) <= 0.5))
        #X, Y = np.meshgrid(X, Y)
        #IMd = IM[idx_diag,0,:,:]
        #ax.plot_surface(X[idxy, :][:,idxx], Y[idxy,:][:,idxx], IMd[idxy,:][:,idxx], alpha=.5)

        #X = np.arange(-1., 2., .01)
        #Y = np.arange(-1., 2., .01)
        #idxx = np.squeeze(np.argwhere(abs(X-D_pad[0][idx_diag,1,0]) <= 0.5))
        #idxy = np.squeeze(np.argwhere(abs(Y-D_pad[0][idx_diag,1,1]) <= 0.5))
        #X, Y = np.meshgrid(X, Y)
        #IMd = IM[idx_diag,1,:,:]
        #ax.plot_surface(X[idxy, :][:,idxx], Y[idxy,:][:,idxx], IMd[idxy,:][:,idxx], alpha=.5)

        #X = np.arange(-1., 2., .01)
        #Y = np.arange(-1., 2., .01)
        #idxx = np.squeeze(np.argwhere(abs(X-D_pad[0][idx_diag,2,0]) <= 0.5))
        #idxy = np.squeeze(np.argwhere(abs(Y-D_pad[0][idx_diag,2,1]) <= 0.5))
        #X, Y = np.meshgrid(X, Y)
        #IMd = IM[idx_diag,2,:,:]
        #ax.plot_surface(X[idxy, :][:,idxx], Y[idxy,:][:,idxx], IMd[idxy,:][:,idxx], alpha=.5)

        #X = np.arange(-1., 2., .01)
        #Y = np.arange(-1., 2., .01)
        #idxx = np.squeeze(np.argwhere(abs(X-D_pad[0][idx_diag,3,0]) <= 0.5))
        #idxy = np.squeeze(np.argwhere(abs(Y-D_pad[0][idx_diag,3,1]) <= 0.5))
        #X, Y = np.meshgrid(X, Y)
        #IMd = IM[idx_diag,3,:,:]
        #ax.plot_surface(X[idxy, :][:,idxx], Y[idxy,:][:,idxx], IMd[idxy,:][:,idxx], alpha=.5)

        #plt.show()
      
        V = np.flip(np.reshape(V[idx_diag,:], image_size), 0)
        plt.figure()
        plt.imshow(V, cmap="Purples")
        #plt.scatter(xs,ys,s=100.)
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=14)
        plt.axis("off")
        plt.show()
        plt.savefig("examplePI", format="jpeg")
        
    if layer == "ls":
        plt.figure()
        if perm_op == "topk":
            #V = np.reshape(V[idx_diag,:], [num_samples,keep])
            #for k in range(keep):
            #    plt.plot(V[:,k], linewidth=5.0)
            LS = ls_diag.eval(feed_dict={diags[0]: D_pad[0], indxs: np.array([[idx_diag]])})
            for k in range(3):
                plt.plot(LS[idx_diag,k,:], linewidth=3.)
        else:
            plt.plot(V[idx_diag,:], linewidth=5.0)

        plt.axis([0.,3000.,-0.2,1.2])
        plt.xticks([0,500,1000,1500,2000,2500,3000], [-1.,-.5,0.,.5,1.,1.5,2.])
        plt.xticks(size=20)
        plt.yticks(size=20)
        #plt.scatter(1000*(xs+1.),ys,s=100.,zorder=10)

        #plt.axis([0.,1000.,0.,1.])
        #plt.xticks([0,200,400,600,800,1000], [0.0,0.2,0.4,0.6,0.8,1.0])
        #plt.legend(["1st landscape","2nd landscape","3rd landscape"], loc="upper right")
        #plt.savefig("exampleSI", format="pdf")
        plt.show()

    if layer == "la":
        V = np.reshape(V[idx_diag,:],[1000,1000,2]) if laguerre_coords else np.reshape(V[idx_diag,:],[1000,1000])
        plt.figure()
        plt.imshow(np.flip(V[:,:,1],0)) if laguerre_coords else plt.imshow(np.flip(laguerre_const-V,0))
        plt.colorbar()
        
    if layer == "pm":
        weight = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "fin-0-perm_eq-0/L")[0])
        biases = sess.run(tf.get_collection(tf.GraphKeys.TRAINABLE_VARIABLES, "fin-0-perm_eq-0/b")[0])
        for idx in range(weight.shape[1]):
            [l1, l2], b = weight[:,idx], biases[:,:,idx]
            plt.plot([-.5,2.], [-.5*l2/l1,2.*l2/l1], linewidth=2., zorder=5, c="black")
        plt.scatter(xs,ys,s=10.,zorder=10)
        plt.axis([-.5,2.,-.5,2.])
        plt.show()
        #plt.savefig("examplePM", format="pdf")

