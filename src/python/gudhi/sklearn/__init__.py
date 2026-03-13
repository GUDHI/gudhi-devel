from .cubical_persistence import CubicalPersistence
from .rips_persistence import RipsPersistence

__all__ = [
    'CubicalPersistence',
    'RipsPersistence',
]

try:
    # if no CGAL
    from .cech_persistence import CechPersistence, WeightedCechPersistence
    
    __all__ += [
        'CechPersistence',
        'WeightedCechPersistence',
    ]
except ImportError:
    pass
