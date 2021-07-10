import sys
import numpy             as np
import tensorflow        as tf
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import gudhi.representations as gdr 

def test_gaussian_perslay():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = gdr.GaussianPerslayPhi((100, 100), ((-.5, 1.5), (-.5, 1.5)), .1)
    weight = gdr.PowerPerslayWeight(1.,0.)
    perm_op = tf.math.reduce_sum
    
    perslay = gdr.Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

def test_tent_perslay():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = gdr.TentPerslayPhi(np.array(np.arange(-1.,2.,.001), dtype=np.float32))
    weight = gdr.PowerPerslayWeight(1.,0.)
    perm_op = 'top3'

    perslay = gdr.Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

def test_flat_perslay():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity
    phi = gdr.FlatPerslayPhi(np.array(np.arange(-1.,2.,.001), dtype=np.float32), 100.)
    weight = gdr.PowerPerslayWeight(1.,0.)
    perm_op = tf.math.reduce_sum
    
    perslay = gdr.Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

def test_gmix_weight():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = gdr.FlatPerslayPhi(np.array(np.arange(-1.,2.,.001), dtype=np.float32), 100.)
    weight = gdr.GaussianMixturePerslayWeight(np.array([[.5],[.5],[5],[5]], dtype=np.float32))
    perm_op = tf.math.reduce_sum

    perslay = gdr.Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

def test_grid_weight():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = gdr.FlatPerslayPhi(np.array(np.arange(-1.,2.,.001), dtype=np.float32), 100.)
    weight = gdr.GridPerslayWeight(np.array(np.random.uniform(size=[100,100]),dtype=np.float32),((-0.01, 1.01),(-0.01, 1.01)))
    perm_op = tf.math.reduce_sum
    
    perslay = gdr.Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

