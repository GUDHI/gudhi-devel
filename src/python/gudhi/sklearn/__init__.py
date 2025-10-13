from .cubical_persistence import CubicalPersistence
from .rips_persistence import RipsPersistence

try:
    # if no CGAL
    from .cech_persistence import CechPersistence, WeightedCechPersistence
except ImportError:
    pass
