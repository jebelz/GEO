"""Landau and Lifshitz Style (QED) style interface to Stokes
Parameters, Jones & Mueller Calculus. We don't need their formulation,
because we have complex 3 (and 4) vectors and tensors. That's all you
need:

With propagation in the Zed directions, pure polarization can be
described by a complex vector with e, with e.z = 0. Mixed states
are described by a density matrix (that is, outer products of vectors
with 0 zed component, unless you hvae some kind of longitudinal polarization,
which could appear, for instance, If your photon turns into a...wait for it...
vector meson).
"""

## \namespace geo.magnetic.faraday Polarization Vectors,
# <a href="http://en.wikipedia.org/wiki/Stokes_parameters">Stokes Parameters
# </a> and the
# <a href="http://en.wikipedia.org/wiki/Jones_calculus">Jones</a> and
# <a href="http://en.wikipedia.org/wiki/Mueller_calculus">Mueller</a>
# calculi.
# \bug complex behavior is in flux.

ROOT2 = 2 ** 0.5

from numpy import exp, pi
from ..metric.wigner import eckart
from .. import Tensor

# Create Spin Eigenstates, convert to Cartesian
Rket = -eckart.Vector(1, 0, 0).tocart()
Lket = -eckart.Vector(0, 0, 1).tocart()

# Create Linear Polarized Combos
Hket = (Rket - Lket)/ROOT2

# V
Vket = (1j)*(Rket + Lket)/ROOT2

## Conjugate states
Rbra = Rket.C
Lbra = Lket.C
Hbra = Hket.C
Vbra = Vket.C

## Density Matricies (Are also operators)
R = Rbra & Rket  # Right polarizer
L = Lbra & Lket  # Left polarizer
H = Hbra & Hket  # H polarizer
V = Vbra & Vket  # V polarizer

## Unpolarized Density Matrix (sum H,V and/or R, L)
U = ((H + V)/2 + (R + L)/2)/2

## 1/4 wave plate (fast horizontal)
QWPH = Tensor(1, 0, 0, 0, 1j, 0, 0, 0, 0) * exp(1j*pi/4)
QH = ((1+1j) * H - (1-1j) * V)/ROOT2
## 1/4 wave plate (fast veritcall)
QWPV = Tensor(1, 0, 0, 0, -1j, 0, 0, 0, 0) * exp(1j*pi/4)

QV = ((H+V) + 1j * (H-V))/ROOT2
