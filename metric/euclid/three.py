"""Rank-3 and Rank 4 Tensors.

Module Constants:

EPSILON


In dev, so some functions may be deprecated.
"""
#pylint: disable=C0301, R0902, R0913, R0904
## \namespace geo::metric::euclid::three Rank Three cartesian tensors.
import operator
import itertools

from ...utils import keydefaultdict
from .. import levi_civita
from . import euclid
from .tensor import DELTA

__all__ = ('HigherOrder', 'Three', 'EPSILON')


## Experimental generalized transformer
# \param R A Rotation (or other transform) \f$ R_{ij} \f$
# \param T A tensor of any rank \f$ T_{ij\ldots mn} \f$
# \returns T' \f$ T'_{i'j' \ldots m'n'} = R_{i'i}R_{j'j}\cdots R_{m'm}
# R_{n'n} T_{ij \ldots mn} \f$
def transformer(R, T):
    """R(T_ijklm) = T_i'j'k'l'm'
    R_m'm *
    R_l'l *
    R_k'k *
    R_j'j *
    R_i'i *
    T_ijklm

    (R * (R * T)._left())._left()
    """
    #pylint: disable=W0212
    Tprime = T
    for dummy in xrange(T.rank):
        Tprime = (R * Tprime)._left()
    return Tprime


## Not sure wether this is used.
def _dilation(t3, alpha):
    """dilate"""
    return type(t3)(*[alpha*item for item in t3.iter()])


## see _dilation().
def _dilation0(t3, s):
    """dilate"""
    return _dilation(t3, s.w)


## Defaults for rank > 2
class HigherOrder(euclid.Tensor_):
    """Do computed __init__ from components"""
    # lazy: direct import of a method


    ## Dispatch things with no index, and ranked things.
    _dispatch_mul = keydefaultdict(lambda rank_: euclid.Tensor_.inner,
                                   {None: _dilation, 0: _dilation0})

    ## Arbitrary Init.
    def __init__(self, *args):
        for attr, arg in itertools.izip(self.components, args):
            setattr(self, attr, arg)

    ## Recursive str
    def __str__(self):
        return "\n".join(itertools.imap(str, itertools.imap(
            self.e, xrange(euclid.DIMENSIONS))))

    ## Slot fillinf multilinear function
    # \param *args \f$ {\bf v_i} \f$
    # \returns \f$ T({\bf v_1}, {\bf v_2}, \ldots, {\bf v_n}) \f$
    def __call__(self, *args):
        from collections import deque
        inputs = deque(args)
        result = self
        while inputs:
            vector_ = inputs.pop()
            if vector_ is None:  # raise
                raise NotImplementedError("Skipped slots is not implemented")
            assert vector_.rank == 1
            result = result * vector_
        return result
        
    ## Dot contracts on ALL indices.
    @euclid.wrank    
    def dot(self, other):
        """T_ijk U_ijk..."""
        return self.inner(other, min(self.rank, other.rank))

    ## Reverse indices.
    # \return \f$ T_{ijklm} \rightarrow T_{mlkji} \f$
    def _reverse_indices(self):
        return self.transpose(*reversed(range(self.rank)))

    ## transpose indices to the right
    # \return \f$ T_{ijklm} \rightarrow T_{mijkl} \f$
    def _right(self):
        args = [self.rank-1] + range(self.rank-1)
        return self.transpose(*args)

    ## tranpose indices to the left
    # \return \f$ T_{ijklm} \rightarrow T_{jklmi} \f$
    def _left(self):
        args = range(1, self.rank) + [0]
        return self.transpose(*args)

    ## Maximum right contraction (unused, in practice).
    @euclid.wrank    
    def _dot_right(self, other):
        return self.inner(other.transpose(*reversed(range(other.rank))),
                          count=other.rank)

    ## Maximum left contraction (unused, in theory).
    @euclid.wrank
    def _dot_left(self, other):
        #pylint: disable=W0212
        return other._dot_right(self).transpose(
            *reversed(range(abs(self.rank - other.rank))))

    ## Multilinear Algebra's
    # <a href="http://en.wikipedia.org/wiki/Higher-order_singular_value_decomposition">
    # SVD</a>.
    # \throws NotImplementedError
    def svd(self):
        """multilinear algebra SVD is not implemented."""
        raise NotImplementedError("for {}".format(type(self)))

    

##  <a href="http://en.wikipedia.org/wiki/Tensor">Rank-3 Tensor</a>
#  class transforms as \f$ T_{ijk}' = M_{ii'}M_{jj'}M_{kk'}T_{i'j'k'} \f$
class Three(HigherOrder):
    """Three(*27)"""

    __metaclass__ = euclid.ranked

    ## The rank is 3
    rank = 3

    ## Transpose and descend rank to make a nice string.
    def __str__(self):
        """something like this:
        [0, 1, 2] [9, 10, 11]  [18, 19, 20]
        [3, 4, 5] [12, 13, 14] [21, 22, 23]
        [6, 7, 8] [15, 16, 17] [24, 25, 26]
        """
        return "\n".join([" ".join(item) for item in
                          map(str.split,
                              map(str,
                                  map(self.transpose(1, 0, 2).e,
                                      xrange(self.rank))))])

    ## Right Index Symmetry
    # \returns \f$ T_{i\{jk\}} \f$
    def Sright(self):
        """T.i{jk}"""
        return getattr(self, 'i{jk}')

    ## Left Index Symmetry
    # \returns \f$ T_{\{ij\}k} \f$
    def Sleft(self):
        """T.{ij}k"""
        return getattr(self, '{ij}k')

    ## The Natural Form (Symmetric and Trace Free in All indices)
    # \return \f$ {\bf T}^S- \frac{3}{5}[T_{iij}\delta_{kl}]^S \f$
    def natural_form(self):
        """Symmetric and Trace Free in all indices."""
        return (self -
                (1/5.)*(
                    (self.iij + self.iji + self.jii) & DELTA
                )).S()

    ## Dual is the pseudoscalar (which is a scalar in geo)
    # \returns \f$ \frac{1}{6}\epsilon_{ijk}T_{ijk} \f$
    def dual(self):
        """Pseudoscalar full contraction with Alternator."""
        return self.dot(EPSILON) / 6

    def spherical(self):
        """As a spherical tensor."""
        from ..wigner import eckart
        return eckart.Three.fromcart(self)

    ## <a href="http://am.ippt.pan.pl/am/article/viewFile/v63p383/pdf">
    # Some Invariants</a> of the  rank-3 tensor.
    # \param order
    # \yields L^(order) Invariants of specified order. 
    def invariants(self, order):
        if order == 1:  # arg
            yield self.dot(EPSILON)
        elif order == 2:
            yield self.iij ** 2
            yield self.iji ** 2
            yield self.jii ** 2
            yield (self & self).ijkijk
            yield (self & self).ijkikj
            yield (self & self).ijkjik
            yield (self & self).ijkjki
            yield (self & self).iijjkk
            yield ((self & self).iijkjk + (self & self).ijikkj) / 2
            yield ((self & self).ijijkk + (self & self).ijjkik) / 2
            yield ((self & self).ijkkij + (self & self).ijkjki) / 2
        elif order == 3:
            yield abs(self.iij ^ self.iji ^ self.jii)  # i amde this up
        else:
            raise ValueError
            
## \f$ \epsilon_{ijk} \f$ \n
# <a href="http://en.wikipedia.org/wiki/Levi-Civita_symbol">
# The fuly antisymmetric tensor</a> in 3-dimensions
EPSILON = Three(*levi_civita.epsilon(euclid.DIMENSIONS))

## The Alternating Tensor
ALTERNATOR = EPSILON

## Also, the Levi-Civta Symbol
LEVI_CIVITA = EPSILON

from collections import namedtuple
## Basis of Cartesian Triads
BASIS = namedtuple("Basis", sorted(EPSILON.__dict__.keys()))(*Three.basis())


