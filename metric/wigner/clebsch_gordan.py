"""Clebsch Gordan Coefficients, algorithm derived from:

# angmom.py  written by Ragnar Stroberg
#                     Feb 5 2013

with additional support added for python's fractions library. (Obviously, when
dealing with 1/2-integer quanta, the fractions.Fraction(1, 2) is far
superior to the float 0.5).


This module has a numpy dependency- that needs to be R E M O V E D.
"""

#pylint: disable=C0103,R0913,E1103

## \namespace geo.metric.wigner.clebsch_gordan
# \f$ \langle j_1m_1,j_2m_2|JM \rangle \f$
from functools import partial
import operator

#from numpy import sqrt, array, arange
from ...utils.trig import sqrt

from ..schur_weyl.pascal import factorial

from . import casimir

prod = partial(reduce, operator.mul)


## Range function for integer and 1/2 integer inputs
def hrange(start, stop, step):
    """[j1, ..., j2] = hrange(j1, j2, step)"""
    from fractions import Fraction
    f = lambda dum: int(round(2*dum))
    return [Fraction(item, 2) for item in range(f(start), f(stop), f(step))]

## Generate valid J values
# \param j1 \f$ j_1\f$
# \param j2 \f$ j_2\f$
# \yields \f$ j \in (|j_1-j_2|, |j_1-j_2|+1, \ldots, j_1+j_2-1, j_1+j_2)\f$
def jrange(j1, j2):
    """jrange(j1, j2) will yield valid j combinations between the two."""
    two_jmax = int(2*(j1 + j2))
    j = int(2*abs(j1 - j2))
    while j <= two_jmax:
        yield float(j) / 2
        j += 2


## Raising Operator scale factor
# \param j degree
# \param m order
# \return \f$ \sqrt{j(j+1) - m(m+1)} \f$
def Jplus(j, m):
    """Jplus(j, m) is really c+_jm's scale factor."""
    return sqrt(casimir.invariant(j) - m*(m+1))


## Lowering Operator scale factor
# \param j degree
# \param m order
# \return \f$ \sqrt{j(j+1) - m(m-1)} \f$
def Jminus(j, m):
    """Jminus(j, m) is really c-_jm's scale factor."""
    return sqrt(casimir.invariant(j) - m*(m-1))


## <a href="http://en.wikipedia.org/wiki/Clebsch-Gordan_coefficients">
#  Clebsch-Gordan coefficient</a>.
# \param j1 1st degree \f$ j \f$ 1/2 or whole integer
# \param m1 1st order  \f$ m \f$  j-like
# \param j2 2nd degree  \f$ j' \f$  1/2 or whole integer
# \param m2 2nd order  \f$ m' \f$  j'-like
# \param J Irrep degree  \f$ J \f$  1/2 or whole integer
# \param M Irrep Order  \f$ M \f$  J-like
# \returns \f$ c_{j_1m_1j_2m_2,JM} = \frac{
#  2i^{j+j'-J}
# \sqrt{\pi} (\frac{J+j-j'}{2})!  (\frac{J-j+j'}{2})!
#  (\frac{j+j'-J'}{2})! \sqrt{j+j'+J+1}  } {
# \sqrt{(\frac{J+j-j'}{2})!}
# \sqrt{(\frac{J-j+j'}{2})!}
# \sqrt{(\frac{j+j'-J}{2})!}
# \sqrt{(\frac{j+j'+J}{2})!}
# \sqrt{2j+1}\sqrt{2j+1'}
# }
# \oint_{d\Omega}{Y_j^m(\Omega)Y_{j'}^{m'}(\Omega)\bar{Y}_J^M(\Omega) d\Omega}
#    \f$ for whole integer arguments.
def clebsch_gordan(j1, m1, j2, m2, J, M):
    """<j1, m2, j2, m2|J, M> = clebsch_gordan(j1, m1, j2, m2, J, M)

    for 1/2 or whole integer arguments."""

    if m1 + m2 != M or abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:  # guard
        return 0

    if j1 > j2 + J or j1 < abs(j2 - J):  # guard
        return 0

    # Step 1:  Find all overlaps < (j1, k1)(j2, J-k1) | JJ >
    # They are related by applying J+ to the j1 and j2 components
    # and setting the sum equal to J+|JJ> = 0.
    # Normalization is from completeness.
    # Start with < (j1, j1)(j2, J-j1) | JJ >, with weight 1
    n = 1
    Norm = [1]
    K = hrange(j1, J-j2-1, -1)
    
    for k in K[1:]:
        n *= -Jplus(j2, J-k-1) / Jplus(j1, k)
        Norm.append(n)
    Norm /= sqrt(sum(item**2 for item in Norm))
    
    # Step 2: Apply J- to get from |JJ> to |JM>, and do the same
    # with j1 and j2. Do this for all the overlaps found in Step 1
    Jm1 = partial(map, partial(Jminus, j1))
    Jm2 = partial(map, partial(Jminus, j2))
    JmM = partial(map, partial(Jminus, J))
        
    cg = 0
    for i, k1 in enumerate(K):
        k2 = J - k1
        if k1 < m1 or k2 < m2:  # filter
            continue
        # multiply by a factor F to account for all the ways
        # you can get there by applying the lowering operator
        F = factorial(k1-m1+k2-m2) / (factorial(k1-m1)*factorial(k2-m2))
 
        c1 = prod(Jm1(hrange(k1, m1, -1)), 1)
        c2 = prod(Jm2(hrange(k2, m2, -1)), 1)
        C  = prod(JmM(hrange(J, M, -1)), 1)
            
        cg += c1*c2/C * Norm[i] * F

    return cg


## Wigner 3-J symbol
# \param J1
# \param M1
# \param J2
# \param M2
# \param J
# \param M
# \return \f$ \left[\begin{array}{ccc} J_1 & J_2 & J \\
#  M_1 & M_2 & M  \end{array} \right] =
# \frac{(-1)^{J_1-J_2+M}}{\sqrt{2J+1}}c_{J_1M_1J_2M_2,JM} \f$
def ThreeJ(J1, M1, J2, M2, J, M):
    """ThreeJ(J1, M1, J2, M2, J, M)"""
    return (pow(-1, J1-J2+M) / sqrt(2*J+1) *
            clebsch_gordan(J1, M1, J2, M2, J, -M))


## Triangle coefficient
# \param j1
# \param j2
# \param j3
# \returns \f$ \Delta(j_1,j_2,j_3) = \frac{(j_1+j_2-j_3)!(j_1-j_2+j_3)!
#  (-j_1+j_2+j_3)!}{(j_1+j_2+j_3+1)!} \f$
def Tri(j1, j2, j3):
    """Tri(j1, j2, j3)"""
    return (factorial(j1+j2-j3) * factorial(j1-j2+j3) * factorial(-j1+j2+j3) /
            float(factorial(j1+j2+j3+1)))


##  Wigner 6-J symbol
# \param j1
# \param j2
# \param j3
# \param J1
# \param J2
# \param J3
# \returns \f$ \left[
# \begin{array}{ccc}j_1 & j_2 & j_3 \\ J_1 & J_2 & J_3 \end{array}\right] \f$
def SixJ(j1, j2, j3, J1, J2, J3):
    """SixJ(j1, j2, j3, J1, J2, J3)"""
    from numpy import array
    triads = (j1, j2, j3), (j1, J2, J3), (J1, j2, J3), (J1, J2, j3)
    for (a, b, c) in triads:
        if a+b < c or abs(a-b) > c or (a+b+c)%1 > 0:  # guard
            return 0

    f = lambda t: factorial(array(
        [t-j1-j2-j3, t-j1-J2-J3, t-J1-j2-J3, t-J1-J2-j3,
         j1+j2+J1+J2-t, j2+j3+J2+J3-t, j3+j1+J3+J1-t])).prod()

    sixj = 0
    for t in hrange(j1+j2+j3, j1+j2+j3+J1+J2+j3+1):
        ff = f(t)
        if ff > 0:  # filter
            sixj += pow(-1, t) * factorial(t+1)/ff

    for (a, b, c) in triads:
        sixj *= sqrt(Tri(a, b, c))
    return sixj


##  Racah W
# \param j1
# \param j2
# \param j3
# \param J1
# \param J2
# \param J3
# \returns \f$ (-1)^{(j_1+j_2+j_3+J_1)}\left[
# \begin{array}{ccc}j_1 & j_2 & j_3 \\ J_1 & J_2 & J_3 \end{array}\right] \f$
def RacahW(j1, j2, j3, J1, J2, J3):
    """RacahW(j1, j2, j3, J1, J2, J3)"""
    return pow(-1, j1+j2+j3+J1) * SixJ(j1, j2, j3, J1, J2, J3)


##  Wigner 9-J symbol TODO: factor out Jrange.
# \param j1
# \param j2
# \param J12
# \param j3
# \param j4
# \param J34
# \param J13
# \param J24
# \param J
# \returns \f$ \left[
# \begin{array}{ccc}j_1 & j_2 & J_{12} \\ j_3 & j_4 & J_{34} \\
# J_{13} & J_{24} & J     \end{array} \right]\f$
def NineJ(j1, j2, J12, j3, j4, J34, J13, J24, J):
    """NineJ(j1, j2, J12, j3, j4, J34, J13, J24, J)"""
    ninej = 0
    for g in range(abs(J-j1), J+j1+1):
        ninej += (pow(-1, 2*g) * (2*g+1)
                    * SixJ(j1, j2, J12, J34, J, g)
                    * SixJ(j3, j4, J34, j2, g, J24)
                    * SixJ(J13, J24, J, g, j1, j3))
    return ninej


##  Normalized Wigner 9-J
# \param j1
# \param j2
# \param J12
# \param j3
# \param j4
# \param J34
# \param J13
# \param J24
# \param J
# \returns \f$
# \langle (J_1, J_3, J_{13})(J_2, J_4, J_{24})J |
#    (J_1, J_2, J_{12})(J_3, J_4, J_{34})J \rangle =
#  \sqrt{(2J_{13}+1)(2J_{24}+1)(2J_{12}+1)(2J_{34}+1)}
#  \left[
# \begin{array}{ccc}j_1 & j_2 & J_{12} \\ j_3 & j_4 & J_{34} \\
# J_{13} & J_{24} & J     \end{array} \right]\f$
def NormNineJ(j1, j2, J12, j3, j4, J34, J13, J24, J):
    """NormNineJ(j1, j2, J12, j3, j4, J34, J13, J24, J)"""
    return (sqrt((2*J13 +1) * (2*J24+1) * (2*J12+1) * (2*J34+1)) *
            NineJ(j1, j2, J12, j3, j4, J34, J13, J24, J))
