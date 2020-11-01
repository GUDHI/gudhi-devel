import sys
import numpy             as np
import tensorflow        as tf
import matplotlib.pyplot as plt

from sklearn.preprocessing import MinMaxScaler
from tensorflow            import random_uniform_initializer as rui

from gudhi.representations import DiagramScaler, Padding, PerslayModel


def test_perslay_image():

    diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
    diag = Padding(use=True).fit_transform(diag)
    D = np.stack(np.array(diag, dtype=np.float32), 0)
    diagrams, empty_feats = [D], np.empty([1,0], dtype=np.float32)
    perslayParameters = {}

    perslayParameters["pweight"]         = None
    perslayParameters["perm_op"]         = "sum"
    perslayParameters["layer"]           = "Image"
    perslayParameters["layer_train"]     = False
    perslayParameters["image_size"]      = (2,2)
    perslayParameters["image_bnds"]      = ((-.501, 1.501), (-.501, 1.501))
    perslayParameters["lvariance_init"]  = .1
    perslayParameters["final_model"]     = tf.keras.Sequential([tf.keras.layers.Flatten()])
    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()

    assert vector.shape == (1,4)
    assert np.abs(vector-np.array([[0,0,5.6e-5,3.3668644]])).sum() <= 1e-6
 
def test_perslay_landscape():

    diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
    diag = Padding(use=True).fit_transform(diag)
    D = np.stack(np.array(diag, dtype=np.float32), 0)
    diagrams, empty_feats = [D], np.empty([1,0], dtype=np.float32)
    perslayParameters = {}

    perslayParameters["pweight"]        = None
    perslayParameters["perm_op"]        = "topk"
    perslayParameters["keep"]           = 3
    perslayParameters["layer"]          = "Landscape"
    perslayParameters["layer_train"]    = False
    perslayParameters["lsample_num"]    = 3
    perslayParameters["lsample_init"]   = np.array(np.arange(-.1,1.1,.5), dtype=np.float32)
    perslayParameters["final_model"]    = "identity"
    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()

    assert vector.shape == (1,9)
    assert np.abs(vector-np.array([[0.,0.,0.,0.1,0.025,0.,0.1,0.1,0.]])).sum() <= 1e-6

