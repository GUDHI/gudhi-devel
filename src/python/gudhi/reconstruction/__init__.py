from .witness_complex import WitnessComplex
from .strong_witness_complex import StrongWitnessComplex
from .nerve_gic import CoverComplex

__all__ = [
    'WitnessComplex',
    'StrongWitnessComplex',
    'CoverComplex',
]

try:
    # if no CGAL or no Eigen
    from .euclidean_strong_witness_complex import EuclideanStrongWitnessComplex
    from .euclidean_witness_complex import EuclideanWitnessComplex

    from .tangential_complex import TangentialComplex

    __all__ += [
        'EuclideanStrongWitnessComplex',
        'EuclideanWitnessComplex',
        'TangentialComplex',
    ]
except ImportError:
    pass
