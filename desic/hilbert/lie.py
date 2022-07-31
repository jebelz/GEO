"""Lie algebra generator...is probably intractable, although rumor has it
that it is not:


M. A. A. van Leeuwen, A. M. Cohen, and B. Lisser, Information about LiE,"
http: //www-math.univ-poitiers.fr/~maavl/LiE/ """
## \namespace geo.desic.hilbert.lie
# <a href="http://en.wikipedia.org/wiki/Lie_group">Lie
# Algebra</a> wish list.
# @image html lie.png
# @image pdf lie.png

import collections
import itertools

from ...metric.levi_civita import parity


## <a href="http://en.wikipedia.org/wiki/Lie_algebra">Lie Algebra</a> base.
class Lie(type):
    """The base-metaclass for Lie Algebra"""


## <a href="http://en.wikipedia.org/wiki/Special_unitary_group">Special
# Unitary</a>  Metaclass.
class SU(Lie):
    """SU(N) Metaclass"""


## <a href="http://en.wikipedia.org/w/index.php?title=Orthogonal_group&redirect=no">Special Orthogonal</a> metaclass.
class SO(Lie):
    """SO(N) metaclass."""


## <a href="http://en.wikipedia.org/wiki/Spin_group">Spin</a> metaclass.
class Spin(Lie):
    """Spin(N) metaclass"""


f = [[[0]*8]*8]*8


## Permute arguments and assign parity: should be superceeded by shur_weyl/.
def permute(*args):
    """d = permute(*args)
    {d[pi(arg)]: +/-1}
    """
    return {tuple(args[i] for i in i_perm): parity(*i_perm) for i_perm
            in itertools.permutations(range(len(args)))}


_SCHelper = collections.namedtuple("_SCHelper", ("indices", "weight"))


## Curried Function
class SCHelper(_SCHelper):
    """Help Permute."""

    def __call__(self):
        perms = {}
        for perm, par in permute(*self.indices).iteritems():
            perms[perm] = par * self.weight
            continue
        return perms


## Structure Constant Maker (curried function).
class _StructureConstantMaker(object):
    """*ssc"""

    def __init__(self, *ssc):
        self.ssc = ssc

    def __call__(self):
        struct_constants = StructureConstant()
        for ssc in self.ssc:
            struct_constants.update(ssc())
            continue
        return struct_constants


## <a href="http://en.wikipedia.org/wiki/Structure_constants">Structure
# Constant</a> Dictionary
class StructureConstant(dict):
    """Structure Constants"""

    ## Dimension of Adjoint Rep, computed ad hoc.
    def dimension(self):
        """Adjoit Rep Dimension?"""
        return max(map(max, self.iterkeys()))

    ## See dimension().
    def __int__(self):
        return self.dimension()

    # Note offset of index in array lingo versur group lingo.
    def __array__(self):
        from numpy import zeros
        result = zeros(tuple(itertools.repeat(int(self), 3)), dtype=complex)
        for index, value in self.iteritems():
            i, j, k = [item-1 for item in index]
            result[i, j, k] = value
        return result


## SU(2) Structure Constants
f_su2 = _StructureConstantMaker(SCHelper((1, 2, 3), 1j))()


## SU(3) Structure Constants
f_su3 = _StructureConstantMaker(
    SCHelper([1, 2, 3], 1.),
    SCHelper([1, 4, 7], 0.5),
    SCHelper([2, 4, 6], 0.5),
    SCHelper([2, 5, 7], 0.5),
    SCHelper([3, 4, 5], 0.5),
    SCHelper([1, 5, 6], -0.5),
    SCHelper([3, 6, 7], -0.5),
    SCHelper([4, 5, 8], 0.5 * 3**0.5),
    SCHelper([6, 7, 8], 0.5 * 3**0.5)
    )()


## SU(3) anti-commutator structure constants.
d_su3 = _StructureConstantMaker(
    SCHelper([1, 1, 8], 1./3**0.5),
    SCHelper([2, 2, 8], 1./3**0.5),
    SCHelper([3, 3, 8], 1./3**0.5),
    SCHelper([8, 8, 8], -1./3**0.5),
    SCHelper([4, 4, 8], -1./2/3**0.5),
    SCHelper([5, 5, 8], -1./2/3**0.5),
    SCHelper([6, 6, 8], -1./2/3**0.5),
    SCHelper([7, 7, 8], -1./2/3**0.5),
    SCHelper([1, 4, 5], 0.5),
    SCHelper([1, 5, 7], 0.5),
    SCHelper([2, 4, 7], -0.5),
    SCHelper([2, 5, 6], 0.5),
    SCHelper([3, 4, 4], 0.5),
    SCHelper([3, 5, 5], 0.5),
    SCHelper([3, 6, 6], -0.5),
    SCHelper([3, 7, 7], -0.5)
    )()
