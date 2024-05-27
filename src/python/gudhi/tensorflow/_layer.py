import tensorflow
# The sklearn parse_version function is taken from 'packaging' and this dependency is not mandatory
from sklearn.utils.fixes import parse_version

class TensorflowKerasLayer(tensorflow.keras.layers.Layer):
    """
    TensorFlow Keras layer
    """
    def __init__(self, **kwargs):
        """
        Constructor for the TensorflowKerasLayer class
        """
        # On tensorflow < 2.16.0 set dynamic argument to True
        if parse_version(tensorflow.__version__) < parse_version('2.16.0'):
            kwargs['dynamic'] = True
        super().__init__(**kwargs)
