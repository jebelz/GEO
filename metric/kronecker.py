"""The Generalized Kronecker Symbol:

>>>delta_2n = kronecker.delta(n)

The rank 2N Generalized Kroneker Delta:

     u1 u2 u3 ... uN
delta
     v1 v2 v3 ... vN

epsilon(n). Since geo does not support mixed indices (it's Euclidean,
after all), the top indices appear on the left, c.f:

delta_(v1, v2, u1, u2)

As such, the rank-2N Kronecker symbol can be unqiuely described as the
tensor with entries of 0, 1, or -1 whose symmetry is defined by the
Young Tableaux:

[0][N+1]
[1][N+2]
 .
 .
 .
[N][2N]
"""

## \namespace geo.metric.kronecker The
# <a href="https://en.wikipedia.org/wiki/Kronecker_delta#generalized_Kronecker_delta">Gernaeralized Kronecker Delta</a>.
import itertools
import operator

__all__ = ('generalized_kronecker_delta')

from .euclid.euclid import ZOO, AXES


## Generalized Kronecker Delta factory function
# \param n \f$ n \in {\bf N} \f$
# \returns \f$ \delta_{\nu_1, \ldots, \nu_n}^{\mu_1, \ldots, \mu_n} =
# \sum_{\sigma \in S_n}{{\mathrm{sgn}}{\sigma}\Pi_{i=1}^n
# \delta_{\nu_{\sigma(i)}}^{\mu_i}} \f$
def delta(n):
    """delta = delta(n) is the rank 2n generalied Kronecker delta
    tensor."""
    from .schur_weyl import monte    
    # Get the group of all permutations on n-letters
    sym_n = monte.Sym(n)
    # This next line does many things:
    # 1) Get/create rank-2n tensor class
    # 2) Loop over all tensor attributes
    # 3) for each attribute, loop over all permutations in Sym(n)
    # 4) split the into top and bottom indices, transpose bottom and
    # 5) compute product pairwise kronecker deltas
    # 6) multiply by parity of the permutation
    # 7) unpack the whole thing into the class.
    args = []
    for indices in itertools.product(AXES, repeat=2*n):
        value = 0
        for sig in sym_n:
            dval = sig.sign()*reduce(
                operator.mul,
                itertools.imap(operator.eq,sig(indices[:n]), indices[n:]))
            value += dval
        args.append(value)
    return ZOO[2*n](*args)

        

## Unqiue synonym
generalized_kronecker_delta = delta
