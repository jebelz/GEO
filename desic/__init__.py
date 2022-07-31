"""The goal of geo.desic is to represent the algebra of the Grand
Unified Theory of the Standard Model in a Poincare invariant manner--
possibly with supersymmetry.


Just kidding. It does dabble in Minkowski-space (minkowski/) and
spinors (hilbert/).

Full quaternions (hamilton.py) and their extensions (cayley_dickson.py)
have their own modules.

The plucker/ package is only an idea: Projective geometry, with a little
mobius.py transform. """
## \namespace geo.desic  A deeper look into algebra of FLAT space(time)--for
# fun.
# @image html gpb.gif
# @image pdf gpb.gif
from . import cayley_dickson
from . import clifford
from . import minkowski
from . import plucker

