"""compute tuple levi-civita elements in n-dimensions

epsilon(n)"""
## \namespace geo.metric.levi_civita The
# <a href="http://en.wikipedia.org/wiki/Levi-Civita_symbol">Levi-Civita
# Symbol</a>, \f$\epsilon\f$, in N dimensions.
import itertools
import operator
import collections
import math

__all__ = ('parity', 'epsilon', 'sgn')


## Traditional Sign Functions
# \param x a number
# \returns sgn(x) with 0-->0.
def sgn(x):
    """traditional sgn function."""
    return int(math.copysign(1, x)) if x else 0


## Parity of a sequence relative to (1, ..., len(args))
#  \param *args a permutation of (1, .., n)
#  \retval +/-1 sign of permutation
def parity(*args):
    """opc parity finder"""
    lst = list(args)
    parity_ = 1
    for i in range(0, len(lst)-1):
        if lst[i] != i:
            parity_ *= -1
            mn = min(range(i, len(lst)), key=lst.__getitem__)
            lst[i], lst[mn] = lst[mn], lst[i]
    return parity_


## True If no elements are repeated
#  \param *args A list or arugments
#  \retval bool IFF args elements are all unique
def _isunique(*args):
    """True If *args is unique"""
    return max(collections.Counter(args).viewvalues()) == 1


## Levi-Civita Component (in _N_-dimensions)
# \param *args The indices in question: \f$ (a_1, \cdots, a_n) \f$
# \returns \f$ \epsilon_{a_1, \cdots, a_n} = \Pi_{1 \le i < j \le n} {\rm sgn}{(a_j-a_i)} \f$
def eps_ijk(*args):
    """eps_ijk(*args):

    Take pair-wise combinations of *args (indices) and difference them,
    negate that and get the sign (-1, 0, 1): then take the product."""
    return reduce(
        operator.mul,
        itertools.imap(
            sgn,
            itertools.imap(
                operator.neg,
                itertools.starmap(
                    operator.sub,
                    itertools.combinations(
                        args,
                        2)))))


## Antisymetric Permutataions
# \param n Dimension
# \returns
# \f$ \epsilon_{a_1 \ldots a_i \ldots a_n} \forall i \in \sigma \cdot
# (1, \ldots, n) \forall \sigma \in S_n \f$
def epsilon(n):
    """n**n element tuple of alternating symbols."""
    return tuple(itertools.starmap(eps_ijk,
                                   itertools.product(xrange(n), repeat=n)))
