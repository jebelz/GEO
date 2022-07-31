"""This Module has tools for converting rank-3 from Cartesian to
the Spherical basis. See:

Note: it does NOT follow:
Proceedings of the World Congress on Engineering 2013 Vol I, WCE 2013,
July 3 - 5, 2013, London, U.K.

too many typos in there, instead, it follows my own algorithm where one
1st looks at Weyl modules, and then operates on those tensors to get
a detraced max weight, and a lower weight thingy that can have its "m" value
selected using the rank-1 and 2 machinery (that is, do work with code that
has already been written elsewhere).

This module does implement the pure rank-3 breakdown from scratch, but with
the fully symmetric Young tableau and racah.py, much of that work has already
been done too.


3 x 3 x 3 breaks up as follow:

Weyl Module  dimW Irreps (j; q)         T_lambda (all tensors are
Tableau                                            symmetrized)
=====================================================================
[0][1][2]     10     3   (1; 0)              v = T.ijj
                     7   (3; 0)           T - (3/5)(DELTA & v)

[0][1]         8     3   (1; 1)              v = T.ijj
[2]                  5   (2; 0)           T + (3/4)(DELTA & v)

[0][2]         8     3   (1; 2)              v = T.ijj
[1]                  5   (2; 1)           T - (3/2)(DELTA & v)

[0]            1     1   (0; 0)            (EPSILON & T).ijkijk
[1]
[2]


Note to the maintainer: This implementation challagens the Zen:

Simple is better than complex.
Complex is better than complicated

However, it evolved from "simple" implementations that turned out to
be less so.

Start with:

>>>PROJECTIONS

which is a dictionary mappong "m" to projecting functions from the
natural form rank-3 symmetric tensor, N_ijk, to |3, m>:

T_(3, m) = PROJECTIONS[m](N_ijk)
"""
## \namespace geo.metric.wigner.eckart.three \f$ \otimes^3{\bf 3} \f$:
# Implemented.
import functools
import operator

from ....utils.trig import sqrt
from ...schur_weyl.young import Tableau
from ...euclid.syzygy import DELTA, EPSILON

## Hash table implements \f$ \pm x\f$ 
SGN = {1: operator.pos, -1: operator.neg}

## Hash table implements \f$ x \pm y\f$ 
ADD = {1: operator.add, -1: operator.sub}



## The symmetrizing tableau.
Y_012 = Tableau([0, 1, 2])

## The 1st mixed symmetry (standard) tableau.
Y_01_2 = Tableau([0, 1], [2])

## The 2nd mixed symmetry (semi-standard) tableau.
Y_02_1 = Tableau([0, 2], [1])

## The antisymmetrizing tableau.
Y_0_1_2 = Tableau([0], [1], [2])


## Extract the vector associated with a Young Tableau.
# \param T A Cartesian tensor, \f$ T_{ijk} \f$
# \param tab A Young Tableau, \f$ \lambda \f$
# \returns \f$ \vec{v_i} = [T^{(\lambda)}]_{ijj} \f$
def _vq(T, tab):
    return T.symmetrize(tab).ijj


## _vq() with tab = ::Y_012 \n
V10 = functools.partial(_vq, tab=Y_012)


## _vq() with tab = ::Y_01_2
V11 = functools.partial(_vq, tab=Y_01_2)


## _vq() with tab = ::Y_02_1
V12 = functools.partial(_vq, tab=Y_02_1)


## Extract tensor associated with vector associated with a Young Tableau.
# \param T A Cartesian tensor, \f$ T_{ijk} \f$
# \param tab A Young Tableau, \f$ \lambda \f$
# \param factor A float
# \returns \f$ A_{ijk} = [\delta_{ij}(T^{(\lambda)})_{kmm}]^{(\lambda)} \f$
def _t1q(T, tab, factor):
    return factor * (DELTA & _vq(T, tab=tab)).symmetrize(tab)


## _t1q() with tab = ::Y_012 and factor = 3/5
T10 = functools.partial(_t1q, tab=Y_012, factor=operator.truediv(3, 5))


## _t1q() with tab = ::Y_01_2 and factor = -3/4
T11 = functools.partial(_t1q, tab=Y_01_2, factor=-operator.truediv(3, 4))


## _t1q() with tab = ::Y_02_1 and factor = 3/2
T12 = functools.partial(_t1q, tab=Y_02_1, factor=operator.truediv(3, 2))


## Extract heighest weight natural form tensor assiciated with a Tableau.
# \param T A Cartesian tensor, \f$ T_{ijk} \f$
# \param tab A Young Tableau, \f$ \lambda \f$
# \param detracer A partial invocation of _t1q()
# \return A natural form tensor: \f$ T^{(\lambda)} - {\mathrm{detracer}}(T) \f$
def _tJq(T, tab, detracer):
    return T.symmetrize(tab) - detracer(T)


## _tJq() with tab = ::Y_012 and detracer = ::T10 \n
# \f$ T_{ijk}^{(3; 0)} =
# T^{(\lambda)} - \frac{3}{5}[\delta_{ij}v^{(0)}_k]^{(\lambda)} \f$
T30 = functools.partial(_tJq, tab=Y_012, detracer=T10)


## _tJq() with tab = ::Y_01_2 and detracer = ::T11 \n
# \f$ T_{ijk}^{(2; 0)} =
# T^{(\lambda)} + \frac{3}{4}[\delta_{ij}v^{(1)}_k]^{(\lambda)} \f$
T20 = functools.partial(_tJq, tab=Y_01_2, detracer=T11)


## _tJq() with tab = ::Y_02_1 and detracer = ::T12 \n
# \f$ T_{ijk}^{(2; 1)} =
# T^{(\lambda)} - \frac{3}{2}[\delta_{ij}v^{(2)}_k]^{(\lambda)} \f$
T21 = functools.partial(_tJq, tab=Y_02_1, detracer=T12)


## Convert rank 3's ::T20 part to a rank-2 spherical tensor (by fiat)
# \param t A general rank-3 Cartesian tensor
# \returns \f$ M_{ij} = \epsilon_{ijk}T^{(2; 0)}_{jkl} \f$
def T20_as2(t):
    """Returns |2,m>_0 part"""
    return EPSILON.inner(T20(t), 2)


## Convert rank 3's ::T21 part to a rank-2 spherical tensor (by fiat)
# \param t A general rank-3 Cartesian tensor
# \returns \f$ M_{ij} = \frac{1}{2}\epsilon_{ijk}T^{(2; 1)}_{jkl} \f$
def T21_as2(t):
    """Returns |2,m>_1 part"""
    return EPSILON.inner(T20(t), 2) / 2.


## Get Scalar part of a rank-3 Cartesian tensor
# \param T A Cartesian tensor, \f$ T_{ijk} \f$
# \returns \f$ \frac{1}{6}\epsilon_{ijk} T^{(j=0)}_{ijk} \f$
def s0(T):
    """scalar = |0,0>_0 weight"""
    return EPSILON.inner(T00(T), 3).w / 6.


## Get Pure Scalar rank-3 Cartesian tensor from T.
# \param T A Cartesian tensor, \f$ T_{ijk} \f$
# \returns \f$ T^{(j=0)}_{ijk} = {\ bf{T}}^{(\lambda_{[0], [1], [2]})} \f$
def T00(T):
    """symmetrize with Y_0_1_2"""
    return T.symmetrize(Y_0_1_2)


## Put all the module functions together into a class that knows whom to call
# when finding |j, m>_q.
class JQM(object):
    """jqm = JQM(T)

    T is a rank 3 Cartesian tensor.

    jqm is a container and a function:

    T[j, q, m] = jqm[j, q, m] = jqm(q, j, m)

    Note: call reverses the container ordering of q and m (so that they
    can have default values of 0). One may also build a tensor from:

    spherical = eckart.Three.fromlist(list(jqm(euclid.three)))
    """

    ## Construct from a rank 3 tensor
    def __init__(self, T):
        assert T.rank == 3, "Just checking."
        ## Rank 3 cartesian tensor
        self.T = T
        print "There are bugs in this object"

    ## iter(jqm) applies an order to jqm[j, q, m], but it must get that
    # order from another class (inappropiate intamacy)--on the other hand,
    # not doing that requires Three to know how JQM works, and that is
    # inverted.
    # \yields \f$ ^qT_j^m \f$ In the order that builds the spherical tensor.
    def __iter__(self):
        from . import Three
        return (self[j, q, m] for j, q, m in Three.ijqm())
        
    ## Nominal format: T[j, q, m]
    def __getitem__(self, index):
        try:
            j, q, m = index
        except TypeError:
            raise ValueError("Expected 3 indices: jqm[j, q, m]")
        return self(j, m, q)

    ## Default Function format:  T(j, m=0, q=0)
    # \param j (manadatory) weight
    # \param m =0, azimuth order
    # \param q =0, seniority index
    # \returns \f$ \langle j, q|_q {\bf T} |j, m\rangle_q \f$
    def __call__(self, j, m=0, q=0):
        """qT_jm = jqm(j, m, q) = jqm[j, q, m]
        .                 ^
        NB: call orders q *last*."""
        if j == 0:
            if q == 0:
                if m == 0:
                    return s0(self.T)
        elif j == 1:
            try:
                v = {0: V10, 1: V11, 2: V12}[q](self.T)
            except KeyError:
                self.throw(j, q, m)
            return v.spherical()[m]
        elif j == 2:
            if q == 0:
                t = T20_as2(self.T)
            elif q == 1:
                t = T21_as2(self.T)
            else:
                self.throw(j, q, m)
            return t.spherical()[j, m]
        elif j == 3:
            if q == 0:
                try:
                    # A straight highest-weight projection..
                    func = PROJECTIONS[m]
                except KeyError:
                    pass
                else:
                    # ...is applied to the natural form rank-3 tensor.
                    return func(self.T.natural_form())
        self.throw(j, q, m)

    ## Generic Error Message.
    # \throws ValueError for unruly combinations of j, m, and q.
    @staticmethod
    def throw(j, q, m):
        """throw(j, q, m) throws a ValueError"""
        raise ValueError("j={}, m={}, q={} is Invalid".format(j, q, m))

    ## This decorator preconditions rank-3 tensors for the extraction of
    # of the highest weight: it projects out the normal form before calling
    # the decorated function.
    @staticmethod
    def func_of_natural_form(func):
        """Convert func(T) to func(T.natural_form())"""

        @functools.wraps(func)
        def nf_func(t):
            """apply func to arg"""
            return func(t.natural_form())

        return nf_func


## Project out 3, -3
# \param m \f$ \pm 1 \f$
# \param S Natural Form Rank-3 Cartesian Tensor
# \returns \f$\pm \sqrt{\frac{1}{8}}
# [(-S_{xxx} + 3S_{xyy}) \mp i(+S_{yyy} - 3S_{xxy})] \f$
def T33(m, S):
    """<3, +/-3 | S> (S is a natural form rank-3 tensor)"""
    # V&V'd
    return SGN[m](
        ADD[m](-S.xxx + 3*S.xyy, -1j*(S.yyy - 3*S.xxy)) / 8**0.5)
    

## Project out j=3 ||m||=2
# \param m \f$ \pm 1 \f$
# \param S Natural Form Rank-3 Cartesian Tensor
# \returns \f$ \frac{1}{2}
# [(S_{xxz} - S_{yyz}) \mp 2iS_{xyz}] \f$
def T32(m, S):
    """<3, +/-2 | S> (S is a natural form rank-3 tensor)"""
    # V&V'd
    return ADD[m](S.xxz - S.yyz, -2*1j*S.xyz) / 2
    

## Project out J=3, ||M||=2.
# \param m \f$ \pm 1 \f$
# \param S Natural Form Rank-3 Cartesian Tensor
# \returns \f$ 3\{
# \sqrt{\frac{8}{15}}[\pm \frac{1}{2}(S_{xzz}-iS_{yzz})] +
# \sqrt{\frac{2}{15}}[\mp(S_{xxx}+S_{xyy})-\frac{i}{4}(S_{yyy}+S_{xxy})]\}\f$
def T31(m, S):
    """<3, +/-1 | S> (S is a natural form rank-3 tensor)"""
    return 3 * (
        (sqrt(2*4/15.) * (SGN[m](S.xzz) - 1j*S.yzz)/2. +
         sqrt(2./15.) * (SGN[-m](S.xxx + S.xyy) - 1j*(S.yyy + S.xxy))/4))


## Project out J=3, M=0.
# \param S Natural Form Rank-3 Cartesian Tensor
# \returns \f$ 3\{
# \sqrt{\frac{1}{10}}\frac{1}{2}[-S_{xxz} - S_{yyz}]
# - \sqrt{\frac{2}{5}}S_{zzz}\} \f$
def T30(S):    
    """<3, 0|S> (S is a rank-3 natural form)."""
    return 3 * (sqrt(1./10.) * (-S.xxz - S.yyz) / 2 +
                sqrt(2./5.) * (-S.zzz))


## Natural Form Projectors keyed by tuple: (j, m)
PROJECTIONS = {3: functools.partial(T33, 1),
               2: functools.partial(T32, 1),
               1: functools.partial(T31, 1),
               0: T30,
               -1: functools.partial(T31, -1),
               -2: functools.partial(T32, -1),
               -3: functools.partial(T33, -1)}
