"""A (hashable) multiset class.

The need for a hashable collections.Counter arises b/c each permutation
in Sym(N) has a cycle_structure(), e.g:

>>>pi = monte.Perm(1, 0, 4, 2, 3)
>>>print repr(pi)
(01)(243)

So that is 1 length 2 and 1 length 3 cycle:

>>>print pi.cycle_structure()
Multiset(2=1, 3=1)

The question then is, how many permutations in Sym(5) have this cycle
structure (and likewise for other cycle structures)?

>>>S5 = monte.Sym(5)
>>>print S5.conjugacy_table()
Multiset(
...Multiset(1=5)=1,
...Multiset(5=1)=24,
...Multiset(1=1, 2=2)=15,
...Multiset(1=1, 4=1)=30,
...Multiset(1=2, 3=1)=20,
...Multiset(1=3, 2=1)=10,
...Multiset(2=1, 3=1)=20
)

which is a multiset of multisets.
"""
## \namespace geo.metric.schur_weyl.de_bruijn Multisets
import collections


#pylint: disable=W0613

## A hashable version of the collections.Counter builtin.
class Multiset(collections.Counter):
    """This is a collections.Counter that is hashable."""

    def __key(self):
        return tuple(sorted(self.items()))

    def __repr__(self):
        return "{0}({1})".format(
            type(self).__name__, ", ".join("{0}={1}".format(
                str(i[0]), repr(i[1])) for i in self.__key()))

    def __hash__(self):
        return hash(self.__key())

    def __delitem__(self, key):
        raise TypeError("{0} does not support item assignment"
                        .format(type(self).__name__))

    ## \throws TypeError
    def clear(self):
        """Not Allowed"""
        raise TypeError("{0} does not support item assignment"
                        .format(type(self).__name__))

    ## \throws TypeError
    def pop(self, *args, **kwargs):
        """Not Allowed"""
        raise TypeError("{0} does not support item assignment"
                        .format(type(self).__name__))

    ## \throws TypeError
    def popitem(self, *args, **kwargs):
        """Not Allowed"""
        raise TypeError("{0} does not support item assignment"
                        .format(type(self).__name__))

    ## \throws TypeError
    def setdefault(self, *args, **kwargs):
        """Not Allowed"""
        raise TypeError("{0} does not support item assignment"
                        .format(type(self).__name__))
