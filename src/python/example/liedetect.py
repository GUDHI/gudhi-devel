"""This file is part of the Gudhi Library - https://gudhi.inria.fr/ - which is released under MIT.
See file LICENSE or go to https://gudhi.inria.fr/licensing/ for full license details.
Author(s):       Henrique Ennes & Raphaël Tinarrage

Copyright (C) 2016 Inria & Institute of Science and Technology Austria
"""

import matplotlib.pyplot as plt
import sklearn
from gudhi.datasets import linear_orbits
from gudhi.liedetect import OrbitFitter

rep_type = ((1, 2),)
group = "torus"
seed = 42
nb_points = 500
print("Target representation type:", rep_type, "\n")

pts = linear_orbits.sample_from_lie_group_rep(
    group=group,
    rep_type=rep_type,
    nb_points=nb_points,
    seed=seed,
)

orbit_fitter = OrbitFitter(pts)
nb_neighbors = 10
orbit_dim = 1

# Print eigenvalues, should be only one small
_ = orbit_fitter.lie_pca(
    nb_neighbors=nb_neighbors, orbit_dim=orbit_dim, method="PCA", verbose=False
)
orbit_fitter.print_lie_pca_eigenvalues(return_vals=False)

# Computes closest Lie algebra
method = "abelian"
frequency_max = 4
la, _ = orbit_fitter.closest_algebra(
    group="torus",
    group_dim=1,
    frequency_max=frequency_max,
    method=method,
    verbose=False,
)
print("\n\nClosest representation type:", la)

# Sample from estimated orbit
pred = orbit_fitter.sample_orbit(nb_points=1000)
print("\n\nMinimal Hausdorff distance:", min(orbit_fitter.hausdorff_distances_))

# Visualize the point cloud and the estimated orbit in 3D using PCA
pca = sklearn.decomposition.PCA(n_components=3).fit(pts)
pts_pca = pca.transform(pts)
orbit_pca = pca.transform(pred)
_, ax = plt.subplots(figsize=(5, 5), subplot_kw={"projection": "3d"})
ax.scatter(
    pts_pca[:, 0],
    pts_pca[:, 1],
    pts_pca[:, 2],
    c="black",
    s=10,
    label="Original point cloud",
)
ax.plot(
    orbit_pca[:, 0],
    orbit_pca[:, 1],
    orbit_pca[:, 2],
    c="magenta",
    lw=10,
    alpha=0.5,
    label="Estimated orbit",
)
plt.legend()
plt.show()
