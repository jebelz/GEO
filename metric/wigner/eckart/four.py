"""This Module has tools for converting rank-4 from Cartesian to
the Spherical basis and back. The method is the same as wigner/three.py,
so go see that for the outline.

3 x 3 x 3 x 3 breaks up as follow:

Weyl Module  dimW  Irreps  j, q
Tableau
=====================================================================
[0][1][2][3]  15     1    (0; 0)         w = F.iijj
                     5    (2; 0)         T = F.iijk - (1/3) w * DELTA
                     9    (4; 0)         N = F - 6/7T - (w&D&D.S)/5


[0][1][2]     15     3    (1; 1)         V = F.ijjk.dual()
[3]                  5    (2; 1)         T = F.iijk - V.dual()/2
                     7    (3; 0)  4--->3 S4= F - 4(D&T)/3 - 2*(DV.d)/5
                                         S3= (EPSILON & S3).ijrstij

[0][1][3]     15     3    (1; 2)
[2]                  5    (2; 2)
                     7    (3; 1)

[0][2][3]     15     3    (1; 3)
[1]                  5    (3; 2)
                     7    (2; 3)

[0][1]         6     1    (0; 1)         w = F.iijj
[2][3]               5    (2; 4)         T = F.iijk - (DELTA * w) /3

[0][2]         6     1    (0; 2)
[1][3]               5    (2; 5)

[0][1]         3     3    (1; 4)     T.iijk.dual()
[2]
[3]

[0][2]         3     3    (1; 5)     T.ijik.dual()
[1]
[3]

[0][3]         3     3    (1; 6)     T.jkii.dual()
[1]
[2]

[0]            0     None
[1]
[2]
[3]
"""
import functools
import operator

## \namespace geo.metric.wigner.eckart.four \f$ \otimes^4{\bf 3} \f$:
# NOT Implemented.
from ....utils.trig import sqrt
from ...schur_weyl.young import Diagram
from ...euclid.syzygy import DELTA, EPSILON
from . import three

## __9__ + __5__ + __1__
Y_0123 = Diagram(4).normalize()

## 3 * (__7__ + __5__ + __3__)
Y_012_3, Y_013_2, Y_023_1 = Diagram(3, 1).standard_tableaux()

## 2 * (__5__ + __1__)
Y_01_23, Y_02_13 = Diagram(2, 2).standard_tableaux()

## 2 * __3__
Y_01_2_3, Y_02_1_3, T_03_1_2 = Diagram(2, 1, 1).standard_tableaux()


## Extract the vector associated with a Young Tableau.
# \param T A Cartesian tensor, \f$ T_{ijk} \f$
# \param tab A Young Tableau, \f$ \lambda \f$
# \returns \f$ \vec{v_i} = [T^{(\lambda)}]_{ijj} \f$
def _vq(T, tab):
    return T.symmetrize(tab).ijj


## Put all the module function together into a class that knows who to call
class JQM(three.JQM):
    """jqm = JQM(T)

    T is a rank 3 Cartesian tensor.

    jqm is a container and a function:

    T[j, q, m] = jqm[j, q, m] = jqm(q, j, m)

    Note: call reverses the container ordering of q and m (so that they
    can have default values of 0).
    """

    ## Construct from a rank 4 tensor
    def __init__(self, T):
        assert T.rank == 4, "Just checking."
        ## Rank 4 cartesian tensor
        self.T = T

    ## Nominal format: T[j, q, m]
    def __getitem__(self, index):
        j, q, m = index
        return self(j, m, q)

    ## Default Function format:  T(j, m=0, q=0)
    # \throws ValueError for unruly combinations of j, m, and q.
    def __call__(self, j, m=0, q=0):
        """Welcome to branchville.
        if j == 0:
            if q == 0:
                s = s0(self.T)
            elif q == 1:
                s = s1(self.T)
            elif q == 2:
                s = s2(self.T)
            else:
                self.throw(j, q, m)
            return s[m]

        elif j == 1:
            try:
                # using a hash-table since partial function aren't tracked in
                # Doxygen's call-graph anyway.
                v = {0: V10, 1: V11, 2: V12,
                     3: V13, 4: V14, 5: V15}[q](self.T)
            except KeyError:
                raise ValueError("j={}, m={}, q={} is Invalid".format(j, q, m))
            return v.spherical()[m]
        elif j == 2:
            if q == 0:
                t = NotImplemented
            elif q == 1:
                t = NotImplemented
            elif q == 2:
                t = NotImplemented
            elif q == 3:
                t = NotImplemented
            elif q == 4:
                t = NotImplemented
            elif q == 5:
                t = NotImplemented
            else:
                self.throw(j, q, m)
            return t.spherical()[j, m]
        elif j == 3:
            if q == 0:
                t = NotImplemented
            elif q == 1:
                t = NotImplemented
            elif q == 2:
                t = NotImplemented
            else:
                self.throw(j, q, m)
            # let three to the work
            return three.Three.fromcart(t).irrep(j=3, q=0, m=m)
        elif j == 4:
            if q == 0:
                if m == -4:
                    return t4_m4(self.T)
                elif m == -3:
                    return t4_m3(self.T)
                if m == -3:
                    return t4_m3(self.T)
                elif m == -2:
                    return t4_m2(self.T)
                elif m == -1:
                    return t4_m1(self.T)
                elif m == 0:
                    return t4_0(self.T)
                elif m == 1:
                    return t4_1(self.T)
                elif m == 2:
                    return t4_2(self.T)
                elif m == 3:
                    return t4_3(self.T)
                elif m == 4:
                    return t4_3(self.T)"""
        self.throw(j, q, m)


def T44(m, S):
    """<4, +/-4 | S>, for S a rank-4 natural form Cartesian tensor."""

def T43(m, S):
    """<4, +/-3 | S>, for S a rank-4 natural form Cartesian tensor."""

def T42(m, S):
    """<4, +/-2 | S>, for S a rank-4 natural form Cartesian tensor."""

def T41(m, S):
    """<4, +/-1 | S>, for S a rank-4 natural form Cartesian tensor."""

def T40(m, S):
    """<4, 0 | S>, for S a rank-4 natural form Cartesian tensor."""
    return None


PROJECTIONS = {4: functools.partial(T44, 1),
               3: functools.partial(T43, 1),
               2: functools.partial(T42, 1),
               1: functools.partial(T41, 1),
               0: functools.partial(T40, 1),
               -1: functools.partial(T41, -1),
               -2: functools.partial(T42, -1),
               -3: functools.partial(T43, -1),
               -4: functools.partial(T44, -1)}
