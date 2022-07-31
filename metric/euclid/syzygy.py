"""Isotropic Tensor Catalog:

This module does not address the so called Syzygy Problem, rather, it gets its
name from the fact that syzygies restrict the number of isotropic tensors
from rank 5, 7, 8, ...and higher. That is: any combination of EPSILON
tensors and Kronecker DELTA tensors with a total of N indices *is* a rank-N
isotropic tensor, but: not all combinations are linearly independent b/c of
syzygies.

For example:
g
>>>print abs(
(EPSILON & DELTA) -
(EPSILON & DELTA).jklim +
(EPSILON & DELTA).klijm -
(EPSILON & DELTA).lijkm
)

Scalar(0)


The isotropic tensors of rank-n are n-tuple values in the module constant
dictionary:

ISO[n]
"""
## \namespace geo.metric.euclid.syzygy Catalog of
# <a href="http://nvlpubs.nist.gov/nistpubs/jres/79b/jresv79bn1-2p49_a1b.pdf">
# Isotropic Tensors</a>.

from .scalar import ONE
from .tensor import DELTA
from .three import EPSILON


## Count FICTS:
# \param n Rank
# \returns For n even: \f$ \frac{n!}{(\frac{n}{2})!2^{n/2}} \f$ \n
# For n odd: \f$ \frac{n!}{3!(\frac{n-3}{2})!2^{(n-3)/2}} \f$ \n
def ficts(n):
    """Number of FICTs of rank n = ficts(n)"""
    from math import factorial
    result = factorial(n)
    if n % 2:
        result //= factorial(3) * factorial((n-3)//2) * pow(2, (n-3)//2)
    else:
        result //= factorial(n//2) * pow(2, n//2)
    return result


## <a href="http://mathworld.wolfram.com/IsotropicTensor.html">Count</a>.
# <a href="http://mathworld.wolfram.com/Syzygy.html">independent</a>
# isotropic tensors of rank-n.
# \param n Rank
# \returns int \f$ -\frac{1}{2}(-3)^{n+1}\ _2F_1(
# -n-1,\frac{1}{2},\frac{1}{2};-\frac{1}{3}) \f$
def motzkin(n):
    """number of independent tensors of rank n = motzkin(n)"""
    if n in ISO:  # seed the recursion
        return len(ISO[n])
    return int(float(n-1)/float(n+1)*(2*motzkin(n-1) + 3*motzkin(n-2)))


## \f$ \delta_{ij}\delta_{kl}\f$
DD = DELTA & DELTA


## \f$ \epsilon_{ijk}\delta_{lm}\f$
ED = EPSILON & DELTA


## \f$ \delta_{ij}\delta_{kl}\delta_{mn}\f$
D3 = DD & DELTA


## \f$ \epsilon_{ijk}\delta_{lm}\delta_{no}\f$
EDD = ED & DELTA


## A starting point for Isotropic Tensors (up to rank-6).
ISO = {0: tuple([ONE]),
       1: tuple(),
       2: tuple([DELTA]),
       3: tuple([EPSILON]),
       4: tuple([DD, DD.ikjl, DD.iljk]),
       5: tuple([ED, ED.ijlkm, ED.ijmkl, ED.ikljm, ED.ikmlj, ED.ilmjk]),
       6: tuple([D3, D3.ikjnlm, D3.imjlkn,
                 D3.ijkmln, D3.iljkmn, D3.imjnkl,
                 D3.ijknlm, D3.iljmkn, D3.injklm,
                 D3.ikjlmn, D3.iljnkm, D3.injlkm,
                 D3.ikjmln, D3.imjkln, D3.injmkl])}

## String of 7 attritubes copied from
# http://nvlpubs.nist.gov/nistpubs/jres/79b/jresv79bn1-2p49_a1b.pdf
SEVEN = """ijkmpqr ijkmqpr ijkmrpq \
ijmkpqr ijmkqpr ijmkrpq \
ijpkmqr ijpkqmr ijpkrmq \
ijqkmpr ijqkpmr ijqkrmp
ijrkmpq ijrkpmq ijrkqmp \
ikmjpqr ikmjqpr ikmjrpq \
ikpjmqr ikpjqmr ikpjrmq \
ikqjmpr ikqjpmr ikqjrmp \
ikrjmpq ikrjpmq ikrjqmp \
impjkqr impjqkr impjrkq
imqjkpr imqjpkr imqjrkp \
imrjkpq imrjpkq imrjqkp \
ipqjkmr ipqjmkr ipqjrkm \
iprjkmq iprjmkq iprjqkm \
iqrjkmp iqrjmkp iqrjpkm"""
#.replace(
#    'r', 'o').replace(
#        'q', 'n').replace(
#            'm', 'l').replace(
#                'p', 'm').split()


