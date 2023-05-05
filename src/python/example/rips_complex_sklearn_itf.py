# Standard data science imports
import numpy as np
from scipy.stats import bootstrap
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt

# Import TDA pipeline requirements
from gudhi.sklearn.rips_persistence import RipsPersistence
from gudhi.representations import DiagramSelector, Landscape
# To fetch the dataset
from gudhi.datasets import remote

import argparse
parser = argparse.ArgumentParser(description='Plot average landscapes')
parser.add_argument("--no-diagram", default=False, action="store_true")
args = parser.parse_args()

# constants for subsample
nb_times = 80
nb_points = 200

def subsample(array):
    sub = []
    # construct a list of nb_times x nb_points
    for sub_idx in range(nb_times):
        sub.append(array[np.random.choice(array.shape[0], nb_points, replace=False)])
    return sub

# constant for plot_average_landscape
landscape_resolution = 1000
# Nothing interesting after 0.4, you can set this filter to 1 (number of landscapes) if you want to check
filter = int(0.4 * landscape_resolution)

def plot_average_landscape(landscapes, color, label):
    landscapes = landscapes[:,:filter]
    rng = np.random.default_rng()
    res = bootstrap((np.transpose(landscapes),), np.std, axis=-1, confidence_level=0.95, random_state=rng)
    ci_l, ci_u = res.confidence_interval
    plt.fill_between(np.arange(0,filter,1), ci_l, ci_u, alpha=.3, color=color, label=label)

walking_ds = remote.fetch_daily_activities(subset="walking")
walking = subsample(walking_ds)

stepper_ds = remote.fetch_daily_activities(subset="stepper")
stepper = subsample(stepper_ds)

cross_ds = remote.fetch_daily_activities(subset="cross_training")
cross = subsample(cross_ds)

jumping_ds = remote.fetch_daily_activities(subset="jumping")
jumping = subsample(jumping_ds)

pipe = Pipeline(
    [
        ("rips_pers", RipsPersistence(homology_dimensions=1, n_jobs=-2)),
        ("finite_diags", DiagramSelector(use=True, point_type="finite")),
        ("landscape", Landscape(num_landscapes=1,resolution=landscape_resolution)),
    ]
)

walking_landscapes = pipe.fit_transform(walking)
plot_average_landscape(walking_landscapes, 'black', 'walking')

stepper_landscapes = pipe.fit_transform(stepper)
plot_average_landscape(stepper_landscapes, 'green', 'stepper')

cross_landscapes = pipe.fit_transform(cross)
plot_average_landscape(cross_landscapes, 'red', 'cross tr.')

jumping_landscapes = pipe.fit_transform(jumping)
plot_average_landscape(jumping_landscapes, 'blue', 'jumping')

plt.title('Average landscapes')
plt.legend()
if args.no_diagram == False:
    plt.show()