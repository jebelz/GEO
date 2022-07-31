"""monte.py: Permutations classes:

Cycle : a cycle representation of permutations
Perm  : A permutation
Sym    : A permutation group (and subgroups Alt and Zn).

An additional helper class:
StandardRep computes the standard representation of Sym(n) on GL(n-1)."""
## \namespace geo.metric.schur_weyl.monte
# <a href="http://en.wikipedia.org/wiki/Three-card_Monte">Permutations</a>,
# their Group, and their Cycles.
import collections
import functools
import itertools
import math
import operator
import types

from .de_bruijn import Multiset

#pylint: disable=R0201,W0613

__all__ = ('Cycle', 'Sym', 'Perm', 'randor', 'Alt')


## \f$ \sum_i{\lambda_i^(n)} \ne n \f$
class PartitionError(UserWarning, ValueError):
    """Invalid partition of an integer, where a valid one is required."""


## Random ordered list of length n
def randor(n):
    """Random Order of Length n"""
    a = Sym(n).neutral()
    a.shuffle()
    return list(a)


## Map integer to a character in\n
# 0, 1, ---, 9, a, b, ---, z, A, B, ---, X, Y, Z, [62], ---
def char(i):
    """char(i) will map 0-87 to a unique single character
    (for cycle notation)."""
    if i < 10:  # sort
        o = 48
    elif i < 36:
        o = 87
    elif i < 62:
        o = 29
    else:
        return "[{}]".format(i)
    return chr(i+o)


## Greatest Common Divisor of 2 arguments
# \param a
# \param b
# \returns GCD(a, b)
def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:
        a, b = b, a % b
    return a


## Least Commom Mutiple of any number of arguments
# \param args some number of integers (forced)
# \return LCM
def lcm(*args):
    """Return least common multiple."""
    return reduce(lambda a, b: a*b // gcd(a, b), map(int, args))


## A decorator to ensure a "binary_op"erator is applied to elements of the
# SAME permutation group.
def ingroup(binary_op):
    """checked_method(binary_op) ensures binary_op is called with
    elements of the same group. It uses functools.wrap to preserved to
    docstring/signature of binary_op, which is otherwise overwritten
    by the process of decoration."""

    ## Add a guard clause before calling the method.
    @functools.wraps(binary_op)  # for users: it preserves the signature.
    def checked_method(self, other):
        """This binary of gets type checked"""
        if not isinstance(other, type(self)):  # poly
            return NotImplemented  # this allows other.__rop__ take over.
        if len(self) != len(other):  # raise for elements of different groups
            raise GroupError(len(self), len(other), binary_op)
        return binary_op(self, other)

    return checked_method


## Apply a binary operation on the tuple of the arguments.
def ontuple(binary_op):
    """binary_op must be function. The decorator applies it to
    tuple(self) and tuple(other)."""
    assert isinstance(binary_op, types.FunctionType), 'got: {}'.format(
        type(binary_op)
        )

    @functools.wraps(binary_op)
    def method(self, other):
        """Tuple of self"""
        return binary_op(self, other)(tuple(self), tuple(other))

    return method


## Raised when combining permutations from different groups
class GroupError(UserWarning, ValueError):
    """GroupError(n, m, op) means you mixed S_n and S_m in op"""

    ## Use call-super to make error message, so clients don't have to.
    def __init__(self, n, m, op):
        # assert that it was used properly.
        assert n != m, "Go debug your code."
        super(GroupError, self).__init__(
            "Can't mix S{} and S{} in operation {}".format(n, m, op)
            )


## A cycle stored in a tuple.
#@functools.total_ordering
class Cycle(tuple):
    u"""c = Cycle(n, \u03c3(n), \u03c3^2(n-1), ...).

    Acts as an impure function on a list (see __call__.__doc__)"""

    ## Tuple constructor
    # \param iterable just like a tuple
    # \returns self The instance
    # \throws ValueError on data that doesn't qualify as a cycle.
    def __new__(mcs, iterable):
        for item in iterable:
            if iterable.count(item) != 1:  # raise
                raise ValueError("Repeat Entry: {}".format(item))
        try:
            self = super(Cycle, mcs).__new__(mcs, map(int, iterable))
        except (ValueError, TypeError):
            raise ValueError("integers only, pls")
        return self

    def __str__(self):
        return "({})".format("".join(map(char, self)))

    ## Combine 2 cycles into a \f$S_n\f$ permutation.
    # \throws ValueError
    def __mul__(self, other):
        """sigma = c1 * c2"""
        if self & other:  # raise
            raise ValueError("Cycles are not DisJoint!")
        return Perm.fromcycles(self, other)

    ## OO-wise: inappropriate intimacy, but Cycles and Perms are INTIMATE:
    # \param other must be an Sym() of cycles disjoint to self.
    # \returns Sym() that combines all cycles.
    def __rmul__(self, other):
        """sigma = c1 * c2 * c3 * ... cN"""
        return other.fromcycles(*(list(other.cycles()) + [self]))

    ## Unpack cycle's portion into a list
    # \param list_ A list to be modified
    # \returns list_ modifiefd
    # \sideeffect list_ is modified! But it is idempotent.
    # \throws ValueError If list is wrong length.
    def __call__(self, list_):
        """cycle(list_) modifies list_ to match cycle."""
        for start, end in itertools.izip(self[:-1], self[1:]):
            try:
                list_[start] = end
            except IndexError as err:
                if len(list_) < max(self):  # raise
                    msg = "List length = {} can't map from cycle's {}"
                    raise ValueError(msg.format(len(list_), max(self)))
                raise err
        list_[self[-1]] = self[0]
        return list_

    ## Set-like AND
    # \param other a Cycle
    # \return set
    def __and__(self, other):
        return set(self) & set(other)

    ## Set-like OR
    # \param other a Cycle
    # \return set
    def __or__(self, other):
        return set(self) | set(other)

    ## Set-like XOR
    # \param other a Cycle
    # \return set
    def __xor__(self, other):
        return set(self) | set(other)

    ## Test Disjointness
    # \param other A cycle
    # \returns bool
    def disjoint(self, other):
        """True Iff argument is a disjoint cycle."""
        return not self & other

    ## Parity from cycle length
    # \param _hash is a private map from even (odd) --> 1 (-1).
    # \return \f$ \pm 1 \f$ Even (Odd).
    def parity(self):
        """+/- = cycle.parity()"""
        return -pow(-1, len(self))

    def parity_true(self):
        """true parity is the number of cycles in the decomposition/"""
        return len(list(self.decompose()))

    def sign_true(self):
        """true sign in -((-1)**len(decomp))"""
        return -pow(-1, self.parity_true())

    sgn = sign_true
    
    ## Is it an (adjacent) transposition
    # \kwd adjacent boolean keyword
    # \returns bool If it's an (adjacent) transposition
    def is_transposition(self, adjacent=True):
        """True iff cycle us an (adjacent=True) tranposition"""
        return ((len(self) == 2) and
                (not adjacent or abs(operator.sub(*self)) == 1))

    ## Convert to the smallest permutation containing the cycle.
    # \returns Sym instance
    def asperm(self):
        """Convert to a minimal element of Sym(n)"""
        return Perm.fromcycles(self)

    ## Decompose cycle into transpositions
    # \yields transpositions (cycles of length 2)
    def decompose(self):
        """Generate a tranposition decomposition If cycle."""
        for target in self[1:]:
            yield type(self)([self[0], target])

    ## cmp(c, c') overrides tuple behavior so that 2 things happen:\n
    # (1) rows (cycles) are sorted by length
    # (2) equal length rows are sorted in "standard" tableau order.
    def __cmp__(self, other):
        """Use for sorted(Sym.cycles(), cmp) to for a tableau ordering."""
        # cmp len or reversed max (which NOT min).
        return cmp(len(self), len(other)) or  cmp(max(other), max(self))


## Symmetric Group in n letters: \f$ S_n \f$
@functools.total_ordering
class Sym(object):
    u"""a = Sym(\u03c3(0), \u2026, \u03c3(n-1)).

    Operations:

    int(sn) = n
    len(sn) = n!

    list(sn) = [list of all p in Sym]


    a * b   What?


    a == b  test equality
    a != b  in equality"""
    ## n-letter Identity
    # \returns \f$ e_n \f$, the neutral elements on \f$ S_n \f$.
    def neutral(self):
        """n-letter neutral element."""
        return self(*xrange(int(self)))

    ## The Order Reversing Permutation
    # \returns \f$ \sigma \f$ The maximal element wrt to the Bruhat order().
    def order_reversing(self):
        """n-letter Order Reversing Permutation."""
        return self(*reversed(range(int(self))))

    ## Yield Perm() 's in the npr().
    def __iter__(self):
        return self.npr()

    ## Order by rank-- this implementation needs to invert Perm.rank()
    def __getitem__(self, index):
        return list(self)[index]

    ## Not allowed
    # \throws TypeError
    def __setitem__(self, index, value):
        raise TypeError("Group does not support item assignment.")

    ## Not allowed
    # \throws TypeError
    def __delitem__(self, index):
        raise TypeError("Group does not support item deletion.")

    ## The natural permutation representation
    # \param n \f$ S_n \f$
    # \yields \f$ \sigma \forall \sigma \in S_n \f$
    def npr(self):
        """Yield all  n-element permutations.

        min(Sym.perms(n)) => Sym.neutral(n)
        max(Sym.perms) == -min(Sym.perms) = Sym.order_reversing(n)
        """
        return itertools.starmap(self,
                                 itertools.permutations(xrange(int(self))))

    ## The Trvial Rep
    def trivial_rep(self):
        """Trivial Irreducible Representation"""
        result = {1: set()}
        for k, v in itertools.imap(operator.pos, self):
            result[k].add(v)
        return result

    ## Returns Sign Rep (partiy:{perms})
    def sign_rep(self):
        """Sign Rep"""
        result = {-1: set(), 1: set()}
        for k, v in itertools.imap(operator.neg, self):
            result[k].add(v)
        return result

    def trivial_rep_chi(self):
        from numpy import matrix
        for perm in self:
            yield 1

    def sign_rep_chi(self):
        from numpy import matrix
        for perm in self:
            yield perm.parity()

    ## Standard Rep in (n-1) square matrix.
    # \returns list of (n-1)-square numoy matrices.
    def standard_rep(self):
        """Standard Rep (in matrix form) (TBD)"""
        return list(StandardRep(int(self)))

    ## Cycle rep as a dictionary
    def cycle_rep(self):
        """Cycle Rep"""
        result = collections.defaultdict(lambda: [])
        for p in self:
            result[p.transposition_number()].append(p)
        return dict(result)

    ## Matrix Rep on \f$ {\bf R}^N \f$
    # \returns list of \f$ {\bf M_{\pi}} \in {\bf R}^{N=n} \forall \pi \in
    # S_n \f$
    def Rn_rep(self):
        """Matrix Rep of NPR. Also called the "defining representation"."""
        from numpy import matrix
        return map(matrix, self)

    ## Defining representation.
    defining_rep = Rn_rep
    
    ## See Rn_rep()
    def __array__(self):
        from numpy import array
        return array(self.Rn_rep())

    ## Generate Transpositions
    # \param n \f$ S_n \f$
    # \kwd adjacent =False
    # \yields \f$ \tau \in S_n \f$
    def tranpositions(self, adjacent=True):
        """Yield all n-element (adjacent=True) transpositions."""
        for tp in (item for item in itertools.imap(Cycle,
                                                   itertools.combinations(
                                                       range(self.n), 2))
                   if item.is_transposition(adjacent=adjacent)):  # filter
            root = range(self.n)
            root[tp[0]], root[tp[1]] = root[tp[1]], root[tp[0]]
            yield self(*root)

    ## Construct elements from an explicit parition of n
    # \param args who sum to \f$ n \f$
    # \yields \f$ \sigma \in S_n \f$ with cycle structure given by partition
    # \throws PartitionError
    def frompartition(self, *args):
        """Generate elements of Sym whos cycle structure is *args:

        >>>for sigma in Sym(lam1, lam2, lamN):   # n = sum(lam_i)
        ........<suite>
        """
        if (list(args) != sorted(args,  # raise
                                 key=operator.neg) or sum(args) != self.n):
            raise PartitionError(*args)
        counter = collections.Counter(args)
        return (sigma for sigma in self
                if counter == sigma.cycle_structure()) # filter

    ## refactor
    def factorial(self):
        """deprecated"""
        print "Use order, not factorial"
        return self.order()

    ## The degree of a group
    # \return \f$ n \f$
    def degree(self):
        """The degree is n"""
        return int(self)

    ## Total number of elements of \f$Sym(n)\f$
    # \returns \f$ n! \f$
    def order(self):
        """n! is the number of elements."""
        return math.factorial(int(self))

    ## Construct from -all- arugments (no cycle shorthand, yet).
    # \param *args
    # \throws ValueError for invalid permutations.
    def __init__(self, n):
        ## Private elements-- no one needs to know this name.
        self.n = int(n)

    ## Degree of group Sym(n)
    # \return n
    def __int__(self):
        """int(Sym) == n"""
        return self.n

    ## Order of group Sym(n)
    # \returns n!
    def __len__(self):
        """len(Sym) == n!"""
        return self.order()

    ## Compares int(sym)
    # \returns \f$ n == n' \f$
    def __eq__(self, other):
        """Test equality."""
        return int(self) == int(other)

    ## Compares int(sym)
    # \returns \f$ n > n' \f$
    def __gt__(self, other):
        return int(self) > int(other)

    ## See __and__()
    # \returns NotImplemented
    # \sideeffect print deprecation warning
    def __mul__(self, other):
        if isinstance(other, type(self)):
            print 'Deprecate: use __and__'
            return self & other
        return NotImplemented

    ## S_n
    def __str__(self):
        return "S_{}".format(self.n)

    ## Apply to sequence, with offset: A branchy mess - not to mention
    # tensor app.
    def __call__(self, *args):
        if len(args) != self.n:  # raise
            msg = "{} got {} items to permute."
            raise ValueError(msg.format(str(self), len(args)))
        return Perm(*args)

    ## The Perfect Shuffle
    # \returns Perm that is a perfect shuffle for even order.
    # \throws ValueError for odd order.
    def perfect_shuffle(self):
        """Perfect shuffle permutation for even n."""
        h, r = divmod(self.n, 2)
        if r:  # raise
            raise ValueError("order must be even")
        z = range(self.n)
        return self(*itertools.chain(*zip(z[:h], z[h:])))

    ## Brute for algorithm- so disappointing.
    # \param list_ A length n inversion vector
    # \returns Perm in Sym(n) with inversion_vector() == list_
    # \throws ValueError for an invalid inversion vector.
    def from_inversion_vector(self, list_):
        """sn = Sym([i0, ..., i_(n-1)])

        builds the permutations from its inversion vector."""
        for item in self:
            if item.inversion_vector() == list(list_):  # filter
                return item
        raise ValueError("Invalid Inversion Vector")

    ## Generate Conjugacy Classes
    # \yields list of Perm() object in the same conjugacy classes for all
    # classes.
    def conjugacy_classes(self):
        """conjugacy_classes()

        yields lists of elements of Sym that are in the same conjugacy class.
        """
        perms = list(self)
        while perms:
            result = [perms.pop()]
            func = operator.methodcaller("isconjugateto", result[0])
            map(result.append, itertools.imap(
                perms.pop, itertools.imap(perms.index,
                                          itertools.ifilter(func, perms[:]))))
            yield result
        
    ## A Multiset of cycle structure Multiset.
    # \returns Multiset() of Multisets counting multiplicity of
    # conjgacy classes.
    def conjugacy_table(self):
        """A Mutliset of conjugacy classes."""
        return Multiset(p.cycle_structure() for p in self)

    ## Yield Rows in the Cayley Table
    # \yields Rows in the multiplication tables.
    def cayley_table(self):
        """Multiplication table, as a..."""
        Gn = iter(itertools.tee(iter(self), 1+len(self)))
        return ([a*b for a in next(Gn)] for b in next(Gn))

    ## Array version of cayley_table()
    # \return array (n+1)-square representing the cayley table.
    def cayley_array(self):
        """Cayley table as a numy array."""
        from numpy import array
        return array(map(functools.partial(map, int), self.cayley_table()))

    ## The class number is
    # <a href="https://oeis.org/A000041">the parition function.</a>
    # \returns int The number of conjugacy classes.
    def class_number(self):
        """The number of conjugacy classes is the parition function."""
        return len(map(len, self.conjugacy_classes()))
        
    ## Get Alternating Subgroup
    # \returns Alt()
    def alternating(self):
        """Alt(n) = Sym(n).alternating()"""
        return Alt(self)

    ## Get Cyclic subgroup
    # \returns Zn()
    def cyclic(self):
        """Zn(n) = Sym(n).cyclic()

        is the cyclic subgroup of order n"""
        return Zn(self)


## The Alternating Group
class Alt(Sym):
    """Alt(n) is the Alternaing group of degree n."""

    ## Generate Perm() -utations in Alt(n)
    # \yields Perm() for even permutations in Sym(n).
    def __iter__(self):
        return itertools.ifilterfalse(bool, super(Alt, self).__iter__())

    ## Order
    # \returns \f$ n!/2 \f$
    def order(self):
        """Total number of elements"""
        return super(Alt, self).order() / 2

    ## The class number is
    # <a href="https://oeis.org/A046682">this</a>.
    def class_number(self):
        """The number of conjugacy classes is the divisor function."""
        return super(Alt, self).class_number()


## \f$ {\bf Z}_n = {\bf Z}/n{\bf Z} \f$
class Zn(Sym):
    """Z(n) is the Cyclic group of degree n."""

    ## Gererate the Cyclic subgroup
    # \yields Perm for cyclic permutations.
    def __iter__(self):
        root = collections.deque(range(self.n))
        for _ in xrange(self.n):
            yield self(*root)
            root.rotate(-1)


    ## Number of elements:
    # \returns \f$ n \f$
    def order(self):
        """Total number of elements"""
        return self.n

    ## The class number is
    # <a href="https://oeis.org/A000005">the divisor function</a>.
    # \returns \f$ \sigma_0(n)\f$
    def class_number(self):
        """The number of conjugacy classes is the divisor function."""
        return super(Zn, self).class_number()


## A Permutation
class Perm(object):
    u"""a = Sym(\u03c3(0), \u2026, \u03c3(n-1)).

    Operations:

    n = len(a) for Sym.

    a[i]    \u03c3(i)

    bool(a) is True (False) for Even (Odd) Parity

    ~a      is the inverse
    a * b   composition: (a*b)(x) = a(b(x))
    a ** n  for n in (..., -2, -1, 0, 1, 2, ...)

    a == b  test equality
    a != b  in equality

    int(a)  see rank

             cmp() is based on __int__
    a < b   int(a) < int(b)
    a <= b   int(a) <= int(b)
    a > b   int(a) > int(b)
    a >= b   int(a) >= int(b)

    +a      copy
    -a      reversed

    a >> n   pad right  [a0,a1,..aM] --> [a0,a1,..aM, M, M+1, ...M+n]
    a << n   pad left   --> [0, 1, .... n-1, a0+n, ... aM+n]

    a & b    combine into larger group

    array(a) numpy array of \u03c3_i

    list = a(sequence) applies a to a sequence, returning a list.

    cycle in a

    str(a)  Cauchy Two Line Format
    repr(a) Cycle Format String."""
    ## Construct from cycles.
    # \param args Some number of cycles (or tuples)
    # \returns Sym instance.
    # \throws ValueError for invalid args.
    @classmethod
    def fromcycles(cls, *args):
        """pi_n = fromcycles(c1, ...)"""
        vals = range(max(map(max, args)) + 1)
        for item in itertools.imap(Cycle, args):
            item(vals)
        if None in vals:  # raise
            raise ValueError("Invalid cycles: {}".format(vals))
        return cls(*vals)


    ## Construct from -all- arugments (no cycle shorthand, yet).
    # \param *args
    # \throws ValueError for invalid permutations.
    def __init__(self, *args):
        ## Private elements-- no one needs to know this name.
        self._args = collections.deque(args)

    ## Pair permutation with "1"
    # \return tuple \f$(1, \sigma)\f$
    def __pos__(self):
        """+s --> (1, s)"""
        return (1, self)

    ## Pair permutation with parity
    # \return tuple \f$({\mathrm {sgn}}{\sigma}, \sigma)\f$
    def __neg__(self):
        """-s --> (sgn(s), s)"""
        return (self.sign(), self)

    ## This is python's len- not group theory's.
    # \returns \f$ n \f$
    def __len__(self):
        return len(self._args)

    ## treat permuted values as a list.
    def __getitem__(self, index):
        if not isinstance(index, slice):  # external polymorphism.
            return self._args[index]
        return list(self)[index]

    ## You can't set/del individual items - a well defined sequence always
    # defines setitem and delitem, even If you can't use them.
    # \throws TypeError
    def __setitem__(self, index, value):
        raise TypeError("Can't set {}[{}]={}".format(self, index, value))

    ## You can't set/del individual items - a well defined sequence always
    # defines setitem and delitem, even If you can't use them.
    # \throws TypeError
    def __delitem__(self, index):
        raise TypeError("Can't delete {}[{}]".format(self, index))

    ## See rank()
    def __int__(self):
        return self.rank()

    ## Order: positive exponent required to reach neutral element.
    def order(self):
        """m = sigma.order() <==> simga**m == sigma ** 0 == 1"""
        return lcm(*self.cycle_lengths())

    ## Length of cycles
    # \yields int Lengths of each cycle.
    def cycle_lengths(self):
        return itertools.imap(len, self.cycles())

    ## True (False) is Odd (Even) Permutation - a bit cryptic, but:
    # -1 ** bool(self) is the parity.
    def __nonzero__(self):
        print "Deprecate this lame overload"
        return self.sign() == -1

    ## Test equality for elements of the same group; 1st decorator is a macro
    # that ensure that condition, 2nd decorator applies operator overload
    # to tuple of arguments.
    @ingroup
    @ontuple
    def __eq__(self, other):
        """Test equality."""
        return operator.eq

    ## Just like __eq__, but not.
    @ingroup
    @ontuple
    def __ne__(self, other):
        """not eq."""
        return operator.ne

    ## Order by inversion number.
    def __gt__(self, other):
        return int(self) > int(other)

    ## Order by inversion number.
    def __lt__(self, other):
        return int(self) > int(other)

    ## Order by inversion number.
    def __ge__(self, other):
        return int(self) >= int(other)

    ## Order by inversion number.
    def __le__(self, other):
        return int(self) <= int(other)

    ## Inverse Permutation
    # Use self ** 0 to get range(len(self)) DRY
    # Use tuple(self) to get iter(self) and the index method--it is the inverse
    def __invert__(self):
        return type(self)(*map(tuple(self).index, self ** 0))

    ## Composition of Permutations (from the same group): it is left-to-right
    # (somewhat arbitrary).
    @ingroup
    def __mul__(self, other):
        return type(self)(*map(other.__getitem__, self))

    ## Solve c = a/b, so that: b * (a/b) == b.
    # \param other \f$ b \f$
    # \returns \f$ ab^{-1} \f$
    @ingroup
    def __div__(self, other):
        return self * ~(type(self)(*other))

    ## Concatenation / outer product
    # \param other \f$ \sigma' \in S_m \f$
    # \returns \f$ \sigma'' \in S_{n+m} \f$
    def __and__(self, other):
        return (self >> len(other)) * (other << len(self))

    ## n-fold composition with self
    def __pow__(self, n):
        """Pardon the multiple returns, but pow is pretty generic like that."""
        if n < 0:
            # Inverse to a positive power.
            return (~self) ** -n
        elif n == 0:  # can't use 'is' because "np.unit32(0) is not 0" e.g.
            # Neutral Element.
            return type(self)(*xrange(len(self)))
        else:
            # repeat(x, n) is better than [x] * n, esp. for large n.
            return reduce(operator.mul, itertools.repeat(self, n))

    ## Pad Left
    def __lshift__(self, n):
        return type(self)(*(range(n) + [item + n for item in self]))

    ## Pad Left, e.g:  3 << sigma >> 5
    def __rlshift(self, n):
        return self << n

    ## Pad Right
    def __rshift__(self, n):
        return type(self)(*(list(self) + range(len(self), len(self) + n)))

    ## <a href="http://en.wikipedia.org/wiki/Permutation_matrix">
    # Matrix Notation </a>
    # \returns numpy.matrix An orthogonal matrix in __SO__(n) with
    # \f$ \frac{\pi}{2} \f$ rotations (i.e., axis swapping)- so really,
    # __GL__(n, Z_2).
    def __array__(self):
        from numpy import zeros
        result = zeros((len(self),)*2, dtype=int)
        result[range(len(self)), list(self)] = 1
        return result

    ## Get the symmetry group to which the permutation belongs:
    # \return Sym in with Perm should reside.
    def group(self):
        """Get the permutations Sym(N) group."""
        return Sym(len(self))

    ## Charteristic Polynomial
    # \returns \f$ p(x) = \sum_{k=0}^n{a_k x^k} \f$ with a_k in (-1, 0, 1)
    def poly(self):
        """Characteristic Polynomial: a numpy.poly1d"""
        from numpy import poly, poly1d, matrix
        return poly1d(poly(matrix(self)))

    ## Image of permutation in the  Standard Rep
    # \returns matrix (n-1) square matrix.
    def standard_rep(self):
        """This may be in the wrong spot."""
        return StandardRep(len(self))(self)

    ## Only numpy can handle it as an index, but it calls __int__
    # \throws TypeError
    def __index__(self):
        raise TypeError("use: x[list(s)], not x[s]")

    ## Permuations are orthogonal matrices (think 90 degree rotations in
    # n-dimensions).
    # \returns Inverse permutation.
    @property
    def T(self):
        """sigma.T is sigma.I"""
        return ~self

    ## Inverse is Transpose
    I = T

    ## Cauchy 2 line format.
    def __str__(self):
        return self.cauchy_two_line_format()

    ## Cauchy 2-line format.
    # \kwd _prop8 private code determines details.
    # \return str
    def cauchy_two_line_format(self, _prop8=" | {} | \n | {} | ".format):
        """Returns:
        |(0 ,  1, ...  n-1) |
        | (s0, s1, ..., sn_1) | """
        return _prop8(tuple(range(len(self))), tuple(self))

    ## Apply to sequence, with offset: A branchy mess - not to mention
    # tensor app.
    def __call__(self, seq, offset=0):
        if hasattr(seq, 'permute'):  # external poly, it's a Tableau
            return seq.permute(self)

        if hasattr(seq, 'transpose'): # then this is a tensor, short circuit?
            return seq.transpose(*self)

        if offset > 0:
            result = (self << offset)(seq)
        elif offset < 0:
            raise ValueError("Offset is negative: {}".format(offset))
        else:
            delta = len(seq) - len(self)
            if delta < 0:
                raise TypeError("Sequence is {} item(s) short.".format(-delta))
            elif delta > 0:
                result = (self >> delta)(seq)
            else:
                result = [seq[n] for n in self]
        if isinstance(seq, basestring):  # external poly, post facto
            result = "".join(result)
        return result

    ## Use a generator to yield cycles.
    def cycles(self):
        """for cycle in self.cycles():
                <suite>
        cycle is a tuple linking permutations."""
        # list repping the permutation
        deal = list(self)
        # function pops items off deal based on their index
        pop = lambda i: deal.pop(deal.index(i))
        while deal:
            # start with the smallest value
            index = min(deal)
            # save it
            result = [pop(index)]
            while True:
                # now follow it
                index = self[index]
                try:
                    # pop it off
                    index = pop(index)
                except ValueError:
                    # pop failed? cycle is complete
                    if len(result) > False: # kwd skip fixed points on request.
                        yield Cycle(result)
                    # go start with the rest of deal
                    break
                else:
                    # append target to list
                    result.append(index)

    ## The discriminant (or the number of disjoint cycles).
    def discriminant(self):
        """det(Pi) is the number of disjoint cycles."""
        return len(list(self.cycles()))

    __abs__ = discriminant

    ## Rotate permuated elements (in place).
    # \kwd m =1, steps to rotate
    # \returns None
    # \sideffect See collections.dequeue.rotate
    def rotate(self, m=1):
        """Rotate perm by m=1 steps."""
        self._args.rotate(m)

    ## Decompose cycle into transpositions
    # \yields transpositions (cycles of length 2)
    def decompose(self):
        """Generate transpositions in to which permuation can
        be decomposed:

        In theory:

        sigma == reduce(mul, sigma.decompose()).

        In practice: see compose.
        """
        return itertools.chain(*itertools.imap(
            operator.methodcaller("decompose"), self.cycles()))

    ## Transposition Number (True Parity)
    def transposition_number(self):
        """NUmber of transpositions..."""
        return len(list(self.decompose()))

    ## Some say the parity is the number of transpositions
    parity_true = transposition_number

    ## The sign:
    # \returns \f$ (-1)^{P(\sigma)} \f$
    def sign_true(self):
        return pow(-1, self.parity_true())
        
    ## The Cycle Structure
    # \returns Multiset() of length of cycles.
    def cycle_structure(self):
        """A Multiset of the cycle structure:

        hence:

        {n:m, p:q} = perm.cycle_structure()

        means there are "m" length-n cycles and "q" length "p" cycles."""
        return Multiset(list(self.cycle_lengths()))

    ## Cycle Format.
    def __repr__(self):
        return self.cycle_str(reduced=True)

    ## A string in cycle-notation
    # \kwd reduced bool ignores trivial cycles.
    # \returns str
    def cycle_str(self, reduced=False):
        """cycle_str([reduced=False])

        is a string in cycle structure: reduced flag eliminates
        trivial cycles."""
        return "".join([str(item) for item in self.cycles()
                        if len(item) > reduced])  # leaky abstraction

    ## Parity (from cycles() decomposition).
    def parity(self):
        """Parity, as a product of the cycle's parity."""
        return reduce(operator.mul,
                      (cycle.parity() for cycle in self.cycles()))

    ## \f$ {\mathrm{sgn}}(\sigma) \f$, computed from cycles,
    # not transpositions
    # \returns \f$ (-1)^{P(\sigma)} \f$
    def sign(self):
        """Parity, as a product of the cycle's parity."""
        return reduce(operator.mul,
                      (cycle.parity() for cycle in self.cycles()))

    ## The decrement of P
    # \returns len(p)-abs(p)
    def decrement(self):
        """The decrement is the length minus number of cycles"""
        return len(self) - abs(self)

    ## Conjuagte Product:
    # \param self \f$ a \f$
    # \param other \f$ b \f$
    # \returns \f$ a b a^{-1} \f$
    def conjugate(self, other):
        """a * b * (~a) == a.conjugate(b)  """
        return self * other * (~self)

    ## TBD/
    def __mod__(self, other):
        return self.cycle_structure() & other.cycle_structure()

    ## Test if argument is conjugate
    # \param other A permutation of same order
    # \return bool If other is conjugate to self...
    # \throws GroupError If order differs.
    @ingroup
    def isconjugateto(self, other):
        """True iff other is conjugate to self."""
        return self.cycle_structure() == other.cycle_structure()

    ## Generate Conjugacy class
    # \yields \f$ p \forall p \in Cl(p) \f$
    def conjugates(self):
        """p.conjugates() yields all elements of Sym conjugate to p."""
        return itertools.ifilter(self.isconjugateto, Sym(len(self)))

    ## Fixed points are not permuted
    # \returns \f$ \{i  | \pi(i) = i \}\f$
    def fixed_points(self):
        """Return a list of fixed points., e.g:

        >>>[2,1,0,3].fixed_points()
        [1, 3]"""
        return set(count for count, item in enumerate(self) if count == item)

    ## is <a href="http://en.wikipedia.org/wiki/Derangement"> Degrangement.
    # \returns bool Iff it a degrangement.
    def isderangement(self):
        """bool = sn.isderangement()"""
        return not tuple(self.fixed_points())

    ## Randomly shuffle in-place.
    # \returns None
    def shuffle(self):
        """sn.shuffle() randomly shuffles inplace."""
        from random import shuffle
        shuffle(self._args)

    ## <a href="http://mathworld.wolfram.com/InversionVector.html">Inversion
    # Vector</a>
    # \returns list
    def inversion_vector(self):
        """n-length list = sn.inversion_vector()"""
        return map(int,
                   (sum(t > s for t in self[:i]) for i, s in enumerate(self)))

    # Inversion number.
    # \returns int Total number of inversions in inversion_vector().
    def inversions(self):
        """Number of inversions."""
        return sum(self.inversion_vector())

    ## Is an inversion.
    # \returns bool Iff inversions() is 1.
    def is_inversion(self):
        """Predicate: is 1 inversion."""
        return self.inversions() == 1

    ## Brute for algorithm- so disappointing.
    @classmethod
    def from_inversion_vector(cls, list_):
        """sn = Sym([i0, ..., i_(n-1)])

        builds the permutations from its inversion vector."""
        for item in Sym(len(list_)):
            if item.inversion_vector() == list_:  # filter
                return item

    ## Find Records as: (index, value)
    # \return list of records (pairs of (when, what)).
    def records(self, result=(0,)):
        """list (of records) = sn.records()"""
        #pylint: disable=E1103
        for count, i in enumerate(self):
            if count:  # init
                if i > max(self[:count]):  # filter
                    result.append((count, i))
            else:
                result += (i,)
                result = [result]
        return result

    ## Count Number of Records
    # \returns int (number of records)
    def record(self):
        """n = sn.records() is the number of records"""
        return len(self.records())

    ## Lexigraphical Rank, computed left-to-right with MSB having a value
    # of (n-1)!, (n-2)!, ... 1, 1; with the value determined by the sorted
    # order.
    # \returns int \f$ 0 \le r < n! \f$
    def rank(self):
        """Lexigraphical rank, un related to cmp()-- which compares
        structure."""
        n = len(self)
        ref = range(n)  # reference list of "bit" values.
        rank = 0
        for position, value in enumerate(self):
            order = ref.index(value)  # count wrt remaining numbers
            rank += order * math.factorial(n - position - 1)
            ref.pop(order)  # remove "used" number
        return rank

    ## <a href="http://en.wikipedia.org/wiki/Lehmer_code">Lehmer Code<a/>
    lehmer_code = rank

    ## Convert to a Young's Tableaux, ensuring cycles of equal length are
    # as close to standard as possible
    # \returns young.Tableaux
    def tableau(self):
        """young.Tableau = sn.tableau()"""
        from .young import Tableau
        return Tableau(*sorted(self.cycles(), cmp=cmp, reverse=True))

    ## Get corresponding Young Diagram
    # \returns young.Diagram
    def diagram(self):
        """young.diagram = sn.diagram()
        see Sym.tableau."""
        return self.tableau().diagram

    def transpose_index(self):
        from geo.metric.euclid.einstein import RUN
        return "".join()

    ## This is the character with respect to a representation.
    def character(self, rep='defining_rep'):
        n = len(self)
        sym_n = Sym(n)
        matrep = operator.methodcaller(rep)(sym_n)
        for perm, mat in itertools.izip(sym_n, matrep):
            if perm == self:
                return mat.trace()[0, 0]

## Compose cycles
# \param args
# \returns Sym (after right padding so all elements have the same 'n')
# \bug args should be reversed to match convention?
def compose(*args):
    """Sym = compose(*args)

    Each arg is a Cycle (or tuple rep), e.g:

    >>>print compose((1,3),(2,3,4),)
    |(0, 1, 2, 3, 4)|
    |(0, 4, 3, 1, 2)|

    """
    size = max(map(max, args))
    return reduce(operator.mul, ((item.asperm() >> size - max(item))
                                 for item in itertools.imap(Cycle,
                                                            args) if
                                 len(item) > 1))


## This class computes the Standard Representation of Sym(N)
class StandardRep(object):
    """srep = StandardRep(n)

    is a tool for computing the standard representation of Sym(n), which
    is a homomorphism to GL(n-1, R); that is:

    >>>for mat in StandardRep(n):
    ....print mat

    is the (n-1, n-1) sqaure numpy matrix associated with each permutation
    in Sym(n).

    As an (n, n-1, n-1) array:

    >>>array(StandardRep(n))

    Computation time rises as 10 ** ( N**2/11 + N/3 - 4),...  N**N.

    NOTE: Matrix multiplcation proceeds in the OPPOSITE order w.r.t
    permutation multiplicaiton.
    """

    def __init__(self, n):
        self.n = int(n)

    def __int__(self):
        return self.n

    ## Get the m-th basis array
    # \param m  for \f$ 0 \le m < n-1 \f$
    # \returns \f$ e_m - e_{m+1} \f$ as an n-dimensional array.
    # \throws ValueError for m out of bounds.
    def make_basis(self, m):
        """make_basis(m) --> [0, 0, ..., 1, -1, 0, ..., 0] (with m
        leading zero pads."""
        from numpy import array
        n = int(self)
        if not 0 <= m < n-1:  # raise
            raise ValueError("Failed: 0 <= {} < {}".format(m, n-1))
        return array([0]*m + [1, -1] + [0]*(n-m-2))

    ## This is an n-1 length list of all valid bases from make_basis().
    @property
    def bases(self):
        """A list of all valid bases (memoized)"""
        while True:
            try:
                return self._bases
            except AttributeError:
                self._bases = map(self.make_basis, range(int(self)-1))

    ## Project permuted basis on original basis (by brute force)
    # \param b A permuted basis
    # \returns tuple of cooefficients multiplying bases()
    # \throws StopIteration for a bum argument.
    def _project(self, b):
        return next(
            itertools.ifilter(
                lambda x: all(
                    reduce(
                        operator.add,
                        (itertools.starmap(
                            operator.mul,
                            itertools.izip(
                                x,
                                self.bases)))) == b),
                itertools.product((-1, 0, 1), repeat=int(self)-1)))

    ## Solve a permutation's projection onto bases().
    # \param perm \f$ \pi \in Sym(n) \f$
    # \returns matrix (n-1) square
    def __call__(self, perm):
        """SR(perm)--> matrix representation of a perm."""
        from numpy import matrix
        return matrix(map(self._project,
                          itertools.imap(
                              perm,
                              map(list,
                                  self.bases))))

    ## Yield standard reps of each permutations in Sym(n)
    # \yields \f$ \{R_{ij}^{(\pi)}\in GL(n-1,R)\}\forall \pi \in Sym(n) \f$
    def __iter__(self):
        return itertools.imap(self, Sym(int(self)))

    ## Array of square arrays
    # \returns array (n!, n-1, n-1).
    def __array__(self):
        from numpy import array
        return array(map(array, self))
