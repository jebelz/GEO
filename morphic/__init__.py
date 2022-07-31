"""Tools for geomorphology:

monge.MongePatch takes the components of an array-like vector (with optional
weights). From that, it fits the 2nd-order surface at [0, 0, 0] with a
taylor.Quadric instance. That then knows how to call the functions
in gauss.py, which compute geomorphological curvatures; it can also
compute the differential geometry curvatures too.

numpy and numpy.linalg are required.
"""
## \namespace geo::morphic Geomorphology: the shapes of topographic
# surfaces.
from . import gauss
from . import monge
from . import taylor
