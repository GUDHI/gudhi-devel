from .knn import KNN


class DTM:
    def __init__(self, k, q=2, **kwargs):
        """
        Args:
            q (float): order used to compute the distance to measure. Defaults to the dimension, or 2 if input_type is 'distance_matrix'.
            kwargs: Same parameters as KNN, except that metric="neighbors" means that transform() expects an array with the distances to the k nearest neighbors.
        """
        self.k = k
        self.q = q
        self.params = kwargs

    def fit_transform(self, X, y=None):
        return self.fit(X).transform(X)

    def fit(self, X, y=None):
        """
        Args:
            X (numpy.array): coordinates for mass points
        """
        if self.params.setdefault("metric", "euclidean") != "neighbors":
            self.knn = KNN(self.k, return_index=False, return_distance=True, **self.params)
            self.knn.fit(X)
        return self

    def transform(self, X):
        """
        Args:
            X (numpy.array): coordinates for query points, or distance matrix if metric is "precomputed", or distances to the k nearest neighbors if metric is "neighbors" (if the array has more than k columns, the remaining ones are ignored).
        """
        if self.params["metric"] == "neighbors":
            distances = X[:, : self.k]
        else:
            distances = self.knn.transform(X)
        distances = distances ** self.q
        dtm = distances.sum(-1) / self.k
        dtm = dtm ** (1.0 / self.q)
        return dtm
