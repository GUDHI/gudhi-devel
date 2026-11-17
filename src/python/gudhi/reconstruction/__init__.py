try:
    # if no CGAL or no Eigen
    from .tangential_complex import TangentialComplex

    __all__ = [
        'TangentialComplex',
    ]
except ImportError:
    pass
