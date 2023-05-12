# Standard scientific Python imports
import numpy as np

# Standard scikit-learn imports
from sklearn.datasets import fetch_openml
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn import metrics

# Import TDA pipeline requirements
from gudhi.sklearn.cubical_persistence import CubicalPersistence
from gudhi.representations import PersistenceImage, DiagramSelector

X, y = fetch_openml("mnist_784", version=1, return_X_y=True, as_frame=False)

# Take only the 500 images from this 70.000 images dataset
# Remove these 2 lines if you want to test the full dataset
X = X[:500]
y = y[:500]
# Reshapes the 1d inputs as images in 28 x 28.
X = X.reshape((-1, 28, 28))

# Target is: "is an eight ?"
y = (y == "8") * 1
print("There are", np.sum(y), "eights out of", len(y), "numbers.")

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.4, random_state=0)
pipe = Pipeline(
    [
        ("cub_pers", CubicalPersistence(homology_dimensions=0, n_jobs=-2)),
        # Or for multiple persistence dimension computation
        # ("cub_pers", CubicalPersistence(homology_dimensions=[0, 1])),
        # ("H0_diags", DimensionSelector(index=0), # where index is the index in homology_dimensions array
        ("finite_diags", DiagramSelector(use=True, point_type="finite")),
        (
            "pers_img",
            PersistenceImage(bandwidth=50, weight=lambda x: x[1] ** 2, im_range=[0, 256, 0, 256], resolution=[20, 20]),
        ),
        ("svc", SVC()),
    ]
)

# Learn from the train subset
pipe.fit(X_train, y_train)
# Predict from the test subset
predicted = pipe.predict(X_test)

print(f"Classification report for TDA pipeline {pipe}:\n" f"{metrics.classification_report(y_test, predicted)}\n")

report = metrics.classification_report(y_test, predicted, output_dict=True)
assert report["0"]["f1-score"] > 0.65
assert report["1"]["f1-score"] > 0.65
