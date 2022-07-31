## \namespace geo.politics.units._entropy
# <a href="http://en.wikipedia.org/wiki/Entropy">Entropy</a>

from ._bases import *

__all__ = ('entropy', 's2omega', 'nat', 'shannon', 'hartley', 'shannon2nat',
           'nat2shannon', 'hartley2nat', 'nat2hartley')


## Theromodynamic (Statistical) Entropy
# \param rho Density Matrix
def entropy(rho):
    from scipy.linalg import logm
    return constants.k * (rho * logm(rho)).trace()


## Number of microstates.
def s2omega(s):
    return scipy.exp(s/constants.k)


## <a href="http://en.wikipedia.org/wiki/Nat_(unit)">Nat</a>
# \param p_i An array
# \returns \f$ H = -\sum_i{p_i\ln{p_i}} \f$
def nat(p_i):
    """nat(p_i)... sum if performed on the last axis:

    H = -(p_i * scipy.log(p_i)).sum(axis=-1)"""
    return -(p_i * scipy.log(p_i)).sum(axis=-1)


## <a href="http://en.wikipedia.org/wiki/Shannon_(unit)">Shannon</a>
def shannon(p_i):
    return nat2shannon(nat(p_i))


## <a href="http://en.wikipedia.org/wiki/Hartley_(unit)">Hartley</a>.
def hartley(p_i):
    return nat2hatley(nat(p_i))

## shannon() --> nat()
shannon2nat = SI(scipy.log(2), 'Sh', 'nat')

## nat() --> shannon()
nat2shannon = ~shannon2nat

## hartley() --> nat()
hartley2nat = SI(scipy.log(10), 'Hart', 'nat')

## nat() --> hartley().
nat2hartley = ~hartley2nat
