"""Geometric Objects.

Core:

euclid/   Rank-N Tensors and Affine Transformations
euler/    Rotations of all kinds.

Extras:

frenet_serret/  Space Curves

schur_weyl/ Schur-Weyl Duality explained (and implemented)

Bonus:

wigner/   Irreducible Representations of Tensors

SO3.py    Rotation Generators, from scratch.
riemann.py An empty module for non-euclidean geometry.
levi_civita.py An alternating symbol factory."""
## \namespace geo.metric Geometry of Plain Euclidean Space
# (\f$ \mathbf{R}^3\f$)
from .euclid import *
from .euler import *
from .wigner import *
from .schur_weyl import *
