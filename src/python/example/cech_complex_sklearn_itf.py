# Standard data science imports
import numpy as np
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

# Import TDA pipeline requirements
from gudhi.sklearn.cech_persistence import CechPersistence
from gudhi.representations.vector_methods import PersistenceLengths
# To build the dataset
from gudhi.datasets.generators import points

no_plot = False

# Mainly for tests purpose - removed from documentation
import argparse
parser = argparse.ArgumentParser(description='Plot average landscapes')
parser.add_argument("--no-plot", default=False, action="store_true")
args = parser.parse_args()
no_plot = args.no_plot

# Build the dataset
dataset_size = 1000
# Noise is expressed in percentage of radius - set it to 0. for no noise
noise = 0.1
# target is a list of 1000 random circle radiuses (between 1. and 10.)
target = 1. + 9. * np.random.rand(dataset_size,1)
# use also a random number of points (between 100 and 300)
nb_points = np.random.randint(100,high=300, size=dataset_size)
X = []

for idx in range(dataset_size):
    pts = points.sphere(nb_points[idx], 2, radius=target[idx])
    ns = noise * target[idx] * np.random.rand(nb_points[idx], 2) - (noise / 2.)
    X.append(pts + ns)

# Cech filtration are squared radius, so transform targets with squared radius values for train/predict purposes
target = target * target

# Split the dataset for train/predict
Xtrain, Xtest, ytrain, ytest = train_test_split(X, target, test_size = 0.25)

pipe = Pipeline(
    [
        ("cech_pers", CechPersistence(homology_dimensions=1, n_jobs=-2)),
        ("max_pers", PersistenceLengths(num_lengths=1)),
        ("regression", LinearRegression()),
    ]
)

model = pipe.fit(Xtrain, ytrain)

# Let's see how our model is predicting squared radiuses from points on random circles
predictions = model.predict(Xtest)
_, ax = plt.subplots()
ax.set_xlabel('target')
ax.set_ylabel('prediction')
ax.scatter(ytest, predictions)
ax.set_aspect("equal")
if no_plot == False:
    plt.show()
