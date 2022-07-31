"""Voigt and Love Multiindex wrapper:

Maps symmetric tensors to 6D vectors, and rank-4 tensors to
6 x 6 symmetric (or not) matrices.

The point is to deal with the fully anisotropic rank-4 elaticisty
tensor in a reduced space of 21 DoF.

Likewise for electric susecptability, in 36 DoFs.
"""
## \namespace geo.magnetic.voigt
# <a href="http://en.wikipedia.org/wiki/Voigt_notation">Voigt Notation</a>.
from geo.metric.euclid.euclid import AXES


## Reversibale Mapping from cartesian pairs to (0-5) 6 DoF indices
# @image html Voigt_notation_Mnemonic_rule.png
# @image latex Voigt_notation_Mnemonic_rule.png
MULTIINDEX = tuple("{}{}".format(AXES[i], AXES[j]) for i, j in
                   zip((0, 1, 2, 1, 2, 0), (0, 1, 2, 2, 0, 1)))


## Voigt to Cartesian index converter
v2c = MULTIINDEX.__getitem__


## Heavily indented conversion from an attr look-up on an even rank tensor
def c2v(attr):
    """(i_1, ..., i_n) = multiindex('xyz ..etc.. yzx')

    Convert rank 2n cartesian component names to numpy-style array
    indices tuple, per Voigt / Mandel / Nye convention.
    """
    if len(attr) % 2:
        raise ValueError("TODO")
    result = []
    for chunk in map("".join, map(None, *([iter(attr)] * 2))):
        try:
            index = MULTIINDEX.index(chunk)
        except ValueError:
            raise ValueError("TODO: fancy pants exception")
        else:
            result.append(index)
    return tuple(result)


## Convert Rank-2 Tensor to 6-element Column
# \param S \f$ \sigma_{ij} \f$
# \returns \f$ \left[\begin{array}{c}
# \sigma_{xx} \\ \sigma_{yy} \\ \sigma_{zz} \\
# \sigma_{yz} \\ \sigma_{zx} \\ \sigma_{xy} \end{array} \right ] \f$
# \todo add weight
def voigt2(S, w=2):
    """convert rank 2 tensor to 6 element numpy column vector."""
    from numpy import matrix
    return matrix([getattr(S, i) for i in map(v2c, range(6))]).T


## Convert Rank-4 Tensor to 6-by-6 matrix
# \param D \f$ D_{ijkl} \f$
# \returns \f$ \left[\begin{array}{cccccc}
# D_{xxxx} & D_{xxyy} & D_{xxzz} & D_{xxyz} & D_{xxxz} & D_{xxxy} \\
# D_{yyxx} & D_{yyyy} & D_{yyzz} & D_{yyyz} & D_{yyxz} & D_{yyxy} \\
# D_{zzxx} & D_{zzyy} & D_{zzzz} & D_{zzyz} & D_{zzxz} & D_{zzxy} \\
# D_{yzxx} & D_{yxyy} & D_{yzzz} & D_{yzyz} & D_{yzxz} & D_{yzxy} \\
# D_{xzxx} & D_{xzyy} & D_{xzzz} & D_{xzyx} & D_{xzxz} & D_{xzxy} \\
# D_{xyxx} & D_{xyyy} & D_{xyzz} & D_{xyyz} & D_{xyxz} & D_{xyxy} \\
#\end{array} \right ] \f$
# \todo add weight
def voigt4(D, w=2**0.5):
    """Conver rank 4 tensor to 6x6 numpy matrix."""
    from numpy import matrix
    return matrix([[getattr(D, i+j) for j in map(v2c, range(6))]
                   for i in map(v2c, range(6))])


## Convert rank 2 or 4 tesnor to Voigt notation.
# \param T a cartesian tensor (rank 2 or rank 4)
# \returns matrix from voigt2() or voigt4().
# \throws TypeError
def voigt(T):
    """rank 2 or 4 tensor to Voight notation."""
    if T.rank == 2:  # raise
        f = voigt2
    elif T.rank == 4:
        f = voigt4
    else:
        raise TypeError("Rank Error")
    return f(T)
