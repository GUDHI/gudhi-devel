import numpy             as np
import tensorflow        as tf
from sklearn.preprocessing import MinMaxScaler
from gudhi.tensorflow.perslay import *
import gudhi.representations as gdr

def test_gaussian_perslay():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = GaussianPerslayPhi((5, 5), ((-.5, 1.5), (-.5, 1.5)), .1)
    weight = PowerPerslayWeight(1.,0.)
    perm_op = tf.math.reduce_sum
    
    perslay = Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

    print(vectors.shape)

    assert np.linalg.norm(vectors.numpy() - np.array(
[[[[1.7266072e-16],
   [4.1706043e-09],
   [1.1336876e-08],
   [8.5738821e-12],
   [2.1243891e-14]],

  [[4.1715076e-09],
   [1.0074080e-01],
   [2.7384272e-01],
   [3.0724244e-02],
   [7.6157507e-05]],

  [[8.0382870e-06],
   [1.5802664e+00],
   [8.2997030e-01],
   [1.2395413e+01],
   [3.0724116e-02]],

  [[8.0269419e-06],
   [1.3065740e+00],
   [9.0923014e+00],
   [6.1664842e-02],
   [1.3949171e-06]],

  [[9.0331329e-13],
   [1.4954816e-07],
   [1.5145997e-04],
   [1.0205092e-06],
   [7.8093526e-16]]]]) <= 1e-7)

test_gaussian_perslay()

def test_tent_perslay():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = TentPerslayPhi(np.array(np.arange(-1.,2.,.1), dtype=np.float32))
    weight = PowerPerslayWeight(1.,0.)
    perm_op = 'top3'

    perslay = Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

    assert np.linalg.norm(vectors-np.array([[0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.09999999, 0.,         0.,
                                             0.2,        0.05,       0.,         0.19999999, 0.,         0.,
                                             0.09999999, 0.02500001, 0.,         0.125,      0.,         0.,
                                             0.22500002, 0.,         0.,         0.3,        0.,         0.,
                                             0.19999999, 0.05000001, 0.,         0.10000002, 0.10000002, 0.,
                                             0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.,         0.,         0.,
                                             0.,         0.,         0.,         0.,         0.,         0.        ]])) <= 1e-7

def test_flat_perslay():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity
    phi = FlatPerslayPhi(np.array(np.arange(-1.,2.,.1), dtype=np.float32), 100.)
    weight = PowerPerslayWeight(1.,0.)
    perm_op = tf.math.reduce_sum
    
    perslay = Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

    assert np.linalg.norm(vectors-np.array([[0.0000000e+00, 0.0000000e+00, 1.8048651e-35, 3.9754645e-31, 8.7565101e-27,
                                             1.9287571e-22, 4.2483860e-18, 9.3576392e-14, 2.0611652e-09, 4.5398087e-05,
                                             5.0000376e-01, 1.0758128e+00, 1.9933071e+00, 1.0072457e+00, 1.9240967e+00,
                                             1.4999963e+00, 1.0000458e+00, 1.0066929e+00, 1.9933071e+00, 1.9999092e+00,
                                             1.0000000e+00, 9.0795562e-05, 4.1222914e-09, 1.8715316e-13, 8.4967405e-18,
                                             3.8574998e-22, 1.7512956e-26, 7.9508388e-31, 3.6097302e-35, 0.0000000e+00]]) <= 1e-7)

def test_gmix_weight():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = FlatPerslayPhi(np.array(np.arange(-1.,2.,.1), dtype=np.float32), 100.)
    weight = GaussianMixturePerslayWeight(np.array([[.5],[.5],[5],[5]], dtype=np.float32))
    perm_op = tf.math.reduce_sum

    perslay = Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

    assert np.linalg.norm(vectors-np.array([[0.0000000e+00, 0.0000000e+00, 1.7869064e-35, 3.9359080e-31, 8.6693818e-27,
                                             1.9095656e-22, 4.2061142e-18, 9.2645292e-14, 2.0406561e-09, 4.4946366e-05,
                                             4.9502861e-01, 1.0652492e+00, 1.9753191e+00, 9.9723548e-01, 1.9043801e+00,
                                             1.4844525e+00, 9.8947650e-01, 9.9604094e-01, 1.9703994e+00, 1.9769192e+00,
                                             9.8850453e-01, 8.9751818e-05, 4.0749040e-09, 1.8500175e-13, 8.3990662e-18,
                                             3.8131562e-22, 1.7311636e-26, 7.8594399e-31, 3.5682349e-35, 0.0000000e+00]]) <= 1e-7)

def test_grid_weight():

    diagrams = [np.array([[0.,4.],[1.,2.],[3.,8.],[6.,8.]])]
    diagrams = gdr.DiagramScaler(use=True, scalers=[([0,1], MinMaxScaler())]).fit_transform(diagrams)
    diagrams = tf.RaggedTensor.from_tensor(tf.constant(diagrams, dtype=tf.float32))

    rho = tf.identity 
    phi = FlatPerslayPhi(np.array(np.arange(-1.,2.,.1), dtype=np.float32), 100.)
    weight = GridPerslayWeight(np.array(np.random.uniform(size=[100,100]),dtype=np.float32),((-0.01, 1.01),(-0.01, 1.01)))
    perm_op = tf.math.reduce_sum
    
    perslay = Perslay(phi=phi, weight=weight, perm_op=perm_op, rho=rho)
    vectors = perslay(diagrams)

    assert np.linalg.norm(vectors-np.array([[0.0000000e+00, 0.0000000e+00, 1.5124093e-37, 3.3314498e-33, 7.3379791e-29,
                                             1.6163036e-24, 3.5601592e-20, 7.8417273e-16, 1.7272621e-11, 3.8043717e-07,
                                             4.1902456e-03, 1.7198652e-02, 1.2386327e-01, 9.2694648e-03, 1.9515079e-01,
                                             2.0629172e-01, 2.0210314e-01, 2.0442720e-01, 5.4709727e-01, 5.4939687e-01,
                                             2.7471092e-01, 2.4942532e-05, 1.1324385e-09, 5.1413016e-14, 2.3341474e-18,
                                             1.0596973e-22, 4.8110000e-27, 2.1841823e-31, 9.9163230e-36, 0.0000000e+00]]) <= 1e-7)
