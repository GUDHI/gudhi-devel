import sys
import numpy             as np
import tensorflow        as tf
import matplotlib.pyplot as plt

from sklearn.preprocessing import MinMaxScaler
from tensorflow            import random_uniform_initializer as rui

my_devices = tf.config.experimental.list_physical_devices(device_type='CPU')
tf.config.experimental.set_visible_devices(devices=my_devices, device_type='CPU')
tf.config.experimental.set_visible_devices([], 'GPU')

from gudhi.representations import DiagramScaler, Padding, PerslayModel

np.random.seed(0)
gauss_init = np.array(np.vstack([np.random.uniform(0.,10.,[2,3]), 1e-5*np.ones([2,3])]), dtype=np.float32)

def test_perslay_image():

    diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
    diag = Padding(use=True).fit_transform(diag)
    D = np.stack(np.array(diag, dtype=np.float32), 0)
    diagrams, empty_feats = [D], np.empty([1,0], dtype=np.float32)
    perslayParameters = {}

    perslayParameters["layer"]           = "Image"
    perslayParameters["layer_train"]     = False
    perslayParameters["image_size"]      = (2,2)
    perslayParameters["image_bnds"]      = ((-.501, 1.501), (-.501, 1.501))
    perslayParameters["lvariance_init"]  = .1

    perslayParameters["pweight"]         = None
    perslayParameters["perm_op"]         = "sum"

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

    perslayParameters["layer"]          = "Landscape"
    perslayParameters["layer_train"]    = False
    perslayParameters["lsample_num"]    = 3
    perslayParameters["lsample_init"]   = np.array(np.arange(-.1,1.1,.5), dtype=np.float32)
    perslayParameters["final_model"]    = "identity"

    perslayParameters["pweight"]        = None
    perslayParameters["perm_op"]        = "topk"
    perslayParameters["keep"]           = 3

    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()
    assert vector.shape == (1,9)
    assert np.abs(vector-np.array([[0.,0.,0.,0.1,0.025,0.,0.1,0.1,0.]])).sum() <= 1e-6

    perslayParameters["pweight"]        = "power"
    perslayParameters["pweight_power"]  = 2
    perslayParameters["pweight_init"]   = 1.
    perslayParameters["pweight_train"]   = False
    perslayParameters["perm_op"]        = "sum"

    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()
    assert vector.shape == (1,3)
    assert np.abs(vector-np.array([[0., 0.03476562, 0.04531251]])).sum() <= 1e-6

def test_perslay_betti():

    diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
    diag = Padding(use=True).fit_transform(diag)
    D = np.stack(np.array(diag, dtype=np.float32), 0)
    diagrams, empty_feats = [D], np.empty([1,0], dtype=np.float32)
    perslayParameters = {}

    perslayParameters["layer"]          = "BettiCurve"
    perslayParameters["layer_train"]    = False
    perslayParameters["lsample_num"]    = 3
    perslayParameters["lsample_init"]   = np.array(np.arange(-.1,1.1,.5), dtype=np.float32)
    perslayParameters["theta"]          = 1.
    perslayParameters["final_model"]    = "identity"

    perslayParameters["pweight"]        = "grid"
    perslayParameters["pweight_size"]   = [100,100]
    perslayParameters["pweight_bnds"]   = ((-.001, 10.001), (-.001, 10.001))
    perslayParameters["pweight_init"]   = np.tile(np.arange(0.,100.,1, dtype=np.float32)[np.newaxis,:], [100,1])
    perslayParameters["pweight_train"]   = False
    perslayParameters["perm_op"]        = "sum"
    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()
    assert vector.shape == (1,3)
    assert np.abs(vector-np.array([[10.091741, 12.746357, 13.192123]])).sum() <= 1e-6

def test_perslay_entropy():

    diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
    diag = Padding(use=True).fit_transform(diag)
    D = np.stack(np.array(diag, dtype=np.float32), 0)
    diagrams, empty_feats = [D], np.empty([1,0], dtype=np.float32)
    perslayParameters = {}

    perslayParameters["layer"]          = "Entropy"
    perslayParameters["layer_train"]    = False
    perslayParameters["lsample_num"]    = 3
    perslayParameters["lsample_init"]   = np.array(np.arange(-.1,1.1,.5), dtype=np.float32)
    perslayParameters["theta"]          = 1.
    perslayParameters["final_model"]    = "identity"

    perslayParameters["pweight"]        = "gmix"
    perslayParameters["pweight_num"]    = 3
    perslayParameters["pweight_init"]   = gauss_init
    perslayParameters["pweight_train"]  = False
    perslayParameters["perm_op"]        = "sum"

    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()
    assert vector.shape == (1,3)
    assert np.abs(vector-np.array([[1.4855406, 1.7884576, 1.6987829]])).sum() <= 1e-6

def test_perslay_rational():

    diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
    diag = Padding(use=True).fit_transform(diag)
    D = np.stack(np.array(diag, dtype=np.float32), 0)
    diagrams, empty_feats = [D], np.empty([1,0], dtype=np.float32)
    perslayParameters = {}

    perslayParameters["layer"]          = "Rational"
    perslayParameters["layer_train"]    = False
    perslayParameters["lnum"]           = 3
    perslayParameters["lmean_init"]     = gauss_init[:2,:]
    perslayParameters["lvariance_init"] = gauss_init[2:,:]
    perslayParameters["lalpha_init"]    = rui(1., 1.)

    perslayParameters["pweight"]         = "power"
    perslayParameters["pweight_power"]   = 2
    perslayParameters["pweight_init"]    = 1.
    perslayParameters["pweight_train"]   = False
    perslayParameters["perm_op"]         = "sum"

    perslayParameters["final_model"]     = tf.keras.Sequential([tf.keras.layers.Flatten()])
    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()
    assert vector.shape == (1,3)
    assert np.abs(vector-np.array([[0.7186792, 0.7186759, 0.718668]])).sum() <= 1e-6

    perslayParameters["layer"]          = "RationalHat"
    perslayParameters["layer_train"]    = False
    perslayParameters["lnum"]           = 3
    perslayParameters["q"]              = 1.
    perslayParameters["lmean_init"]     = gauss_init[:2,:]
    perslayParameters["lr_init"]        = rui(1., 1.)

    perslayParameters["pweight"]         = "power"
    perslayParameters["pweight_power"]   = 2
    perslayParameters["pweight_init"]    = 1.
    perslayParameters["pweight_train"]   = False
    perslayParameters["perm_op"]         = "sum"

    perslayParameters["final_model"]     = tf.keras.Sequential([tf.keras.layers.Flatten()])
    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()
    assert vector.shape == (1,3)
    assert np.abs(vector-np.array([[-0.00675799, -0.00620097, -0.00510298]])).sum() <= 1e-6

def test_perslay_exponential():

    diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
    diag = Padding(use=True).fit_transform(diag)
    D = np.stack(np.array(diag, dtype=np.float32), 0)
    diagrams, empty_feats = [D], np.empty([1,0], dtype=np.float32)
    perslayParameters = {}

    perslayParameters["layer"]          = "Exponential"
    perslayParameters["layer_train"]    = False
    perslayParameters["lnum"]           = 3
    perslayParameters["lmean_init"]     = 1e3 * gauss_init[:2,:]
    perslayParameters["lvariance_init"] = gauss_init[2:,:]

    perslayParameters["pweight"]        = None
    perslayParameters["perm_op"]        = "max"

    perslayParameters["final_model"]     = tf.keras.Sequential([tf.keras.layers.Flatten()])
    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()
    assert vector.shape == (1,3)
    assert np.abs(vector-np.array([[0.9940388, 0.99311596, 0.99222755]])).sum() <= 1e-6

def test_perslay_peq():

    diag = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diag = DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diag)
    diag = Padding(use=True).fit_transform(diag)
    D = np.stack(np.array(diag, dtype=np.float32), 0)
    diagrams, empty_feats = [D], np.empty([1,0], dtype=np.float32)
    perslayParameters = {}

    perslayParameters["layer"]          = "PermutationEquivariant"
    perslayParameters["layer_train"]    = False
    perslayParameters["lpeq"]           = [(5, "sum"), (5, "sum")]
    perslayParameters["lweight_init"]   = rui(1e-1, 1e-1)
    perslayParameters["lbias_init"]     = rui(0.1,   0.1)
    perslayParameters["lgamma_init"]    = rui(1e-1, 1e-1)

    perslayParameters["pweight"]        = None
    perslayParameters["perm_op"]        = "topk"
    perslayParameters["keep"]           = 3

    perslayParameters["final_model"]     = tf.keras.Sequential([tf.keras.layers.Flatten()])
    model = PerslayModel(name="perslay", diagdim=2, perslay_parameters=[perslayParameters], rho="identity")
    vector = model([diagrams, empty_feats]).numpy()
    assert vector.shape == (1,15)
    assert np.abs(vector-np.array([[0.4375, 0.41875, 0.375, 0.4375, 0.41875, 0.375, 0.4375, 0.41875, 0.375, 0.4375, 0.41875, 0.375, 0.4375, 0.41875, 0.375]])).sum() <= 1e-6
