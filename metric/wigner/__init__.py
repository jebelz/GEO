"""Representations of the rotation group and their tensor products and sums.

eckart/is a subpackage devoted to implementing the nitty gritty of
irreducible representations of Vectors and Tensors-- that is, the hard
part of actually converting them to and from Cartesian form.

wigner.py is in progress. It does Wigner matrices, and other necessities.

clebsch_gordan.py computed Clebsch-Gordan (CG) coefficients (which, while are
design for combining quantum angular momentum, they are in fact indispensable
in combining tensor products of the fundamental SO(3) rep and computing
tensor sums of irreps.

racah.py As in Racah Algebra, puts (CG) to work in Kets(j, m) that can be
multiplied and converted to sums of irreps Kets(J, M). For example, the
product of 2 vectors is a tensor (dyad) that can have Scalar, Vector, or
pure rank-2 tensors parts. That is 100% analogous to combining the spin of
2 identical spin-1 particles (vector Bosons) into states of total angular
momentum 0, 1, and 2... I mean *exactly*.

casimir.py is a utilities module.

Any tensor may be expressed as a sum of irreducible tensors transforming
as irreps of SO(3). For a general rank n Cartesian tensor (with 3**n
components), the irreducible tensors are characterized by a weight
0 <= j <= n with 2j+1 components each. When n > 2, each weight may be
represented more than once, requiring another label, q, called the
seniority index.

As an example, consider the 27-component rank-3 tensor (see eckart/three.py)--
which of course get all mixed up under rotations. However, there are linear
combinations with the following dimensions that do not mix with other linear
combinations (they're called Weyl-modules):

10
8
8
1


The Weyl modules can be further decomposed into projections that transform
as follows:

10 = 7 + 3    (natural form rank 3 plus a vector)
8  = 5 + 3    (natural form rank 2 plus a vector)
8  = 5 + 3    (natural form rank 2 plus a vector)
1  = 1        (scalar)

That means a general rank-3 Cartesian tensor has a 7-dimensions (2j+1 for
j=3) symmetric, trace-free pure rank-3 part, 2 pure rank 2 parts, 3 vector
parts and one scalar part. Since various ranks are repeated, the
idea of a senority index (q) is introduced. It counts from
the maximal j on down--and in the interest of python, it's indexed from 0.
So in summary:


Weyl Module      Weights &
Dimensions       Seniority
--------------------------------------------------------------
10               (j=3)       +   (j=1, q=0)
8                (j=2, q=0)  +   (j=1, q=1)
8                (j=2, q=1)  +   (j=1, q=2)
1                (j=0)

Note the ordering of the j=2 seniority is derived from the ordering
of young tableaux.

What this means is that the j=0 part transforms like a scalar, the j=1
parts transform like vectors, the j=2 pieces transform like a symmetric trace
free rank 2 tensors. Finally the j=3 piece is called the Natural Form:
it is symmetric in all indices, and trace-free on all pairs of indices.

For rank N, the Natural Form is the j=N part of the tensors, with 2j+1
components labeled by the order m.

The |j=N, m=+/-N> parts can be computed by finding the polyadic product:

|+>|+>...N times |+>  or
|->|->...N times |->, where

|+/-> are non-trivial vector eigenvalues of an azimuthal rotation (e.g.,
|+/-> = |j=1, m=+/-1>), which can be found in eckart.py.

The other m-states can be computed by applying the raising (lower) operators.

In theory, all the weights, orders, and senorities can be computed from
the Clebsch Gordan coefficients (see clebsch_gordan.py) found in the 3**N
N-adic products of |->, |0>, |+>. Needless-to-say, that is involved, and
hence there is racah.py: this module help do the algebra (which is called
Racah Algebra)

For a simpler case, consider the Rank 2 tensor. It's Weyl modules are
familiar:

6  totally symmetric S = (T_ij + T_ji) / 2
3  ....antisymmetric A = (T_ij - T_ji) / 2

Those are linear combinations. The symmetric module can be further reduced
into:

5 + 1

where the "1" is scalar like--that is, an isotropic tensor. There is only
one rank-2 isotopic tensor, the Kronecker Delta. Hence, the projection
of T onto the Kronecker Delta is "1":

S_2^0 = <S|delta> = delta * Tr(S).

Subtracting that from T_S given the pure j=2 tensor, which has the following
properties:

5 = 2j+1 components
Symmetric
Trace-Free.

How are those 5 components distributed among:

|2, +2>
|2, +1>
|2,  0>
|2, -1>
|2, -2>.....?


Take for instance

|+> = -(x + iy) / sqrt(2)

Its dyadic product, |+>|+>, is:

(xx - yy + i(xy+yx)) / 2

and that solves the problem:

|2, 2> = S.xx - S.yy + 1j*(S.xy + S.yx) = S.xx - S.yy + 2*1j*S.xy

Meanwhile, the racah.py tells is:

>>>print racah.Ket(1, 1)*racah.Ket(1, 0)
[sqrt(1/2) | 1.0, 1> +
sqrt(1/2) | 2.0, 1>]

= (|1, 1> + |2, 1>)/sqrt(2)

So that |+>|0> is split between the j=2 and j=1 parts. With |0> = z,
they dyadic becomes:

|+>|0> = -(x + iy)z / sqrt(2) = -(xz - iyz) / sqrt(2), so that:

|2, 1> = -(S.xz - 1j*S.yz) / 2.

and so on. For N=3 (three.py) and N=4 (four.py), this gets out-of-hand, fast.
"""
## \namespace geo.metric.wigner Irreducible Tensor Representations:
#\f${\bf T}_{(n)}=\sum_{j=0}^n{\sum_{q=1}^{N_n^{(j)}}{{\bf T}_{(n)}^{(j;q)}}}\f$

from . import casimir
from . import clebsch_gordan
from . import racah
from . import eckart
from . import wigner
