"""Young Tablaux and the Permutation Group (for acting on tensor indices)
and other needs. Classes are:

Diagram
Tableau (Skew-, SemiStandard-, Standard-, Normal-)


This module raison detre is to implement tensor symmetries through the
Schur-Weyl duality formalism, and it works (!), via:

young.Diagram(n1, n2, ..., nm).symmetrize(T)

where rank(T) == n1 + n2 +... nm.

Here are the steps to creating a tensor with the symmetry defined by
a diagram, for example, we'll do 3 = 2 + 1, the mixed symmetry part
of a rank 3 tensor.

(1) Draw the diagram:

>>>d = Diagram(2, 1)

[ ][ ]
[ ]

(2) Now fill it to make the normal tableau

>>>t = d.fill(0, 1, 2) = d.normalize()

[0][1]
[2]

That makes a Tableau. Now get the Row and Colum symmetrizers:

(3) The Row symmetrizer is the set of permutations that leave the tableau
row-equivalent, e.g.,

   >>>t' = t.permute(sigma)

which changes the 'i' in t into sigma(i)--- it's by value, not position,
and ensure that

   >>>t'.is_row_equivalent(t) is True


That's done via:
>>>row = list(t.Row())
[, (01)]

So the 1st permutation is the identity (which in cycle notation is empty);
The second permutation swaps 0 and 1:
[1][0]    <==> [0][1]
[2]       <==> [2]

which is row equivalent to the original tableau (on the left).

(4) Now do the Column Symmetrizer, likewise (this includes the parity)

>>>col = t.Col()
[(1,), (-1, 02)]

>>>print t.permute(col[1])
[2][1]
[0]


Now compute the Symmetrizer, it is the product of the column * row, with
the sign on the column included:

S ~ (sign)col * row
  = [1 - (02)] * [1 + (01)]
  =  1 + (01) - (02) - (021)

We must include a weight, given by the dimension of irrep divided by
the dimension of the permutation group:

d.dim() / int(d)! = 2 / 6 = 1/3

So we now apply those permutations to the tensor indices, so that:

T_(2,1) = 1/3 * [ T_ijk +
                  T_jik -
                  T_kji -
                  T_kij ]

Tensors with this symmetry form the Weyl Module for the diagram, with
dimension (in N=3 dimensions):

>>>d.dimW(3)
8

The remaining 8 dimensions with (2,1) symmetry are the projection of the
other standard tableau. Add the 1 antisymmetric dimensions and the 10
symmetric, and you get 27 = 3 x 3 x 3.

Wow.


Furthermore:

Diagrams are extended to allow skew diagrams (though, nothing is implemented
with them). The general idea begin a skew diagram is that is a normal diagram
combined with another diagram that 'fits' inside it. It is represented
by division, e.g:

>>>print x
[][][][]
[][][]
[][]
[]

>>>print y
[][]
[]

>>>print x/y
    [][]
  [][]
[][]
[]



Tableaux are extended to allow representation theory beyond the
Robinson-Schensted bijection.
"""
## \namespace geo.metric.schur_weyl.young Youngs Diagram and Tableau
import collections
import functools
import itertools
import math
import operator

from ...utils.exceptions import RankError
from ..euclid.einstein.albert import FREE

from .monte import Cycle, Sym, Perm


__todo__ = """
(1) Skew Tableau don't exisit (Diagrams do).
(2) How do you order Tableau? DO you use the lehmer code of the permutation,
or the Yamamachi symbol
(3) Construction from a Yamamachi symbol doesn't exits.
(4) Do we need a class Yamamachi(tuple)?

"""


## Is a list non-decreasing?
# \ param *args
# \return True IFF list_ is non-decreasing.
def is_nondecreasing(*args):
    """True IFF arguments are non-decreasing."""
    return list(args) == sorted(args)

## Is a list increasing?
# \ param list_ A list
# \return True IFF list_ is increasing.
def is_increasing(*args):
    return (is_nondecreasing(*args) and
            max(map(args.count, args)) == 1)


## Convert index number to an alphabetic index as seen in Einstein summation
int2alpha = lambda n: chr(n + ord(FREE[0]))

#pylint: disable=W0613,C0103

## Partition Error of any kind
class PartitionError(UserWarning, ValueError):
    """Invalid Parition (e.g. Young Tableau) Format"""

    def __init__(self, *args):
        super(PartitionError, self).__init__(
            "Valid parition format is not {}".format(sorted(args,
                                                            reverse=True)))


## Skew formation failure
class SkewValueError(UserWarning, ValueError):
    """Skew diagram fail: shape-wise, or numerator or demoniator is
    already skew."""


class SemiStandardError(UserWarning, TypeError):
    """Raised when an operation needs at least a SemiStandardTableau"""

class StandardError(SemiStandardError):
    """Raised when an operation needs at least a StandardTableau"""

class NoramlError(StandardError):
    """Raised when an operation needs a Noraml Tableau"""


## Check if arguments are a valid, ordered partion
# \param seq A sequence
# \returns bool Iff seq is reversed sorted.
def ispartition(seq):
    """ispartitions(list_) Iff list_ *is* reverse sorted."""
    return list(seq) == sorted(seq, reverse=True)
        

## Tensor symmetrizer
# \f$ \frac{1}{N!}\sum_{\pi \in \Pi(N)}{T_{\pi(i, j, \ldots)}}\f$
# \param t a Tensor \f$T\f$ of rank \f$ k \f$
# \returns function \f$ P_{\lambda}(T) = \frac{1}{k!}
# \sum_{\sigma \in G_k}{\tau_{\sigma}T} \f$
def P_lambda(n):
    """The Symmetrizer for rank n (factory function):

    P = P_lambda(n)"""
    return Tableau(range(n)).symmetrize


## Antisymmetrizer
# \f$ \frac{1}{N!}\sum_{\pi \in \Pi(N)}{sgn(\pi)T_{\pi(i, j, \ldots)}}\f$
# \param t a Tensor \f$T\f$ of rank \f$ k \f$
# \returns function
# \f$ Q_{\lambda}(T) = \frac{1}{N!}
# \sum_{\pi \in \Pi(N)}{(-1)^{P(\pi)}{T_{\pi(i, j,\ldots)}}}\f$
def Q_lambda(n):
    """The Antisymmetrizer for rank n (factory function):

    Q = Q_lambda(n)
    """
    return (~Tableau(range(n))).symmetrize


## The information that makes a cell---this is just a description for
# informational purposes.
class Cell(collections.namedtuple("Cell", "i j arm leg hook value")):
    """Cell(i, j, arm, leg, hook, value) is all the information you
    need to know about a cell in a Young Tableau."""
    
    def yamanouchi(self):
        return self.value, self.i
    

## <a href="http://en.wikipedia.org/wiki/Young_tableau#Diagrams">
# Young Diagram</a>.
class Diagram_(object):
    """d = Diagram(*row_lengths), e.g.:

    d = Diagram(5,5,2,1)

    unicode(d) draws the diagram.

    d is a generator function:
    for t_tab in d(tensor):
    ...dumps the weyl module
    """
    
    ## Count the number of standard tableau with n-boxes
    # \returns long from pascal.Pascal.telphone()
    @staticmethod
    def involution_number(n):
        """involution_number(n) is the number of standard tableau on n-boxes.
        It is the telephone number from pascal.Pascal."""
        from .pascal import Pascal
        return Pascal(n).telephone()

    ## Construct from a multiset (Counter)
    # \param b A collection.Counter
    # \returns Diagram
    @classmethod
    def frommultiset(cls, b):
        return cls(*sum(([[k]*b[k] for k in reversed(b.keys())])))
    
    ## Construct from Row Lengths
    # \param args decreasing row lengths
    # \throws PartitionError If args do not make an ORDERED partition of N.
    def __init__(self, *args, **kwargs):
        if not ispartition(args):  # raise
            raise PartitionError(*args)
        ## rows are row length.
        self.rows = args
        #self.skew = kwargs.get('skew') or ()

    ## iterate plus skew part    
    def _view(self):
        for diag, skew in itertools.izip_longest(self, self.skew, fillvalue=0):
            yield diag, skew
        
    def __div__(self, other):
        if self.skew:  # raise
            raise SkewValueError("Diagram is already skew")
        return type(self)(*self, **dict(skew=other))
        
    ## There are many ways to order diagrams, here it goes:\n
    # sum --> len --> min
    def __cmp__(self, other):
        return (self._cmp_size(other) or
                self._cmp_len(other) or
                self._cmp_min(other))

    ## The weight of a diagram is the number of boxes, which is the sum of
    # row lengths.
    # \returns int The number of cells.
    def weight(self):
        """Count the number of boxes in the diagram."""    
        return sum(self)

    ## The length (need reference) is the number of rows.
    # \returns int The number of rows.
    def length(self):
        """Count the number of rows in the diagram."""    
        return len(self)

    ## The shortest row
    # \returns int The length of the shortest row.
    def shortest_row(self):
        """Length of the shortest row."""
        return min(self)

    ## The longest row is the length of the complementary ::Diagram.
    # \returns int The length of the longest row.
    def length_complement(self):
        """Length of the longest row."""
        return max(self)    
    
    ## Compare number of boxes
    # \param other Diagram of list
    def _cmp_size(self, other):
        return cmp(self.weight(), sum(other))

    ## Compare number of rows
    # \param other Diagram of list
    def _cmp_len(self, other):
        return cmp(self.length(), len(other))

    ## Compare shortest row
    # \param other Diagram of list
    def _cmp_min(self, other):
        return cmp(self.shortest_row(), min(other))

        
    ## Number of rows.
    def __len__(self):
        return len(self.rows)

    ## Get length of row[n]
    def __getitem__(self, index):
        return self.rows[index]

    ## Number of boxes in main
    def __int__(self):
        return self.weight()

    ## Number of configurations when filled.
    def __long__(self):
        return long(self.num_fillings())

    def num_fillings(self):
        return math.factorial(int(self))
        
    ## Yield rows
    def __iter__(self):
        return iter(self.rows)

    ## Check If a box is in the diagram, or a diagram fits in another.
    # \param pair (i, j) OR a diagram
    # \returns bool ith-row has a box at jth column.
    def __contains__(self, pair):
        """(row, col) in diagram <==> row, col has a box."""
        if isinstance(pair, Diagram_):  # guard for external poly
            if self.skew or pair.skew:  # raise
                raise TypeError("Skew doesn't compare")
            return all([s <= p for p, s in itertools.izip_longest(self, pair,
                                                                  fillvalue=0)])
        # 
        try:
            # unpack row, col (counting from 1) argument.
            i, j = pair
        except (TypeError, ValueError):
            raise ValueError(
                "Expected a container of length 2 (row, col), got: {}".format(
                    type(pair)))
        try:
            # get size of row (If it exists)
            row = self[i-1]
        except IndexError:
            # the row did not exist
            result = False
        else:
            # the column existence is a straight compare.
            result = 0 < j <= row
        return result
            
    ## Complementary Diagram
    # \returns Diagram that is complement (row <--> col).
    def __invert__(self):
        """~diagram is the complementary diagram."""
        # need to do skew part, without kicking off recursion
        skew = ~Diagram(*self.skew) if self.skew else ()
        return type(self)(
            *[1 + self.leg(1, j) for j in range(1, 1 + list(self)[0])],
            **dict(skew=skew))

    ## ASCII Boxes
    # \returns ASCII Boxes
    def __str__(self, _box = "[{}]"):
        return "\n".join([("   "*y + "[{}]"*(x-y)).format(*(" " * (x-y)))
            for x, y in self._view()])
        
    ## Unicode Boxes
    # \param _box -> u2610 (private)
    # \returns unicode diagram with boxes
    def __unicode__(self):
        box = u"\u2610"
        return "\n".join([(" "*y + box*(x-y)) for x, y in self._view()])

    ## Faithful Rep. (not like that)
    # \returns str that can be eval'd
    def __repr__(self):
        return type(self).__name__ + repr(self.rows)

    ## The length is the same as the python len().
    length = __len__

    ## Conjugate (TODO: have invert call conjuagte).
    conjugate = __invert__

    ## Diagram Addition
    # \param other Diagram
    # \returns Diagram with sum of row lengths
    def __add__(self, other):
        return type(self)(*itertools.starmap(
            operator.add,
            itertools.izip_longest(self, other, fillvalue=0)))

    ## TODO: invert addition
    def __sub__(self, other):
        raise TypeError("subtraction not implemented")

    ## The <a href="http://en.wikipedia.org/wiki/Multiset">multiset</a>
    # version of a Young Diagram is a set with multiplicties of elements.
    # \returns de_bruijn.Multiset as {row length: multiplicity}
    def multiset(self):
        """Convert to a multiset (or bag, or collections.Counter)"""
        from .de_bruijn import Multiset
        return Multiset(self)

    ## Multiplication is like multiset union
    # \returns Diagram that is the multiset union of row lengths.
    def __mul__(self, other):
        """Union: see collection.Counter.__or__"""
        return self.frommultiset(self.multiset() | other.multiset())

    ## Shape *is* the diagram.
    def shape(self):
        return tuple(self)
    
    ## Loop over all combos and yield possible diagrams
    # \yields list Candidate tensor products
    def _gen_tensor_prod0(self, other):
        """The idea is to add the "other" to the right side of self,
        and then just yield the valid diagrams. It seems to work, however,
        testing is not exhaustive due to lack of expertize."""
        # now step the right diagram down 1 row at a time
        for offset in range(len(self)+1):
            # start with a copy of 1 diagram
            base_diagram = list(self)
            # add the right's boxs to left's, rowwise.
            for count, item in enumerate(other):
                index = offset + count
                # check If the result is longer
                if index > len(base_diagram) - 1:
                    # result *is* longer
                    base_diagram.append(item)
                else:
                    # result is NOT longer, just add the rows.
                    base_diagram[index] += item
            yield base_diagram

    ## Generate tensor product graphically.
    # \yields Diagram in the tensor product            
    def _gen_tensor_prod(self, other):
        """Filter valid partitions and construct Diagram from
        _gen_tensor_prod0"""
        return itertools.starmap(
            Diagram,
            itertools.ifilter(ispartition,
                              self._gen_tensor_prod0(other)))
        
    ## Choose ordering for diagram generator
    def _tensor_product_with_diagram(self, other):
        # The alogrithm needs the "smaller" diagram on the left
        # This does NOT agree with cmp.
        if other > self:  # order arguments
            g = other._gen_tensor_prod(self)
        else:
            g = self._gen_tensor_prod(other)
        return list(g)

    ## Solve each element in a list and concatenate
    def _tensor_product_with_list(self, other):
        return reduce(list.__add__,
                      map(self._tensor_product_with_diagram, other))

    ## Hashtable selects method for computing tensor product base on type.
    _tp_choser = {False: _tensor_product_with_list,
                  True: _tensor_product_with_diagram}
                  
    def _tensor_product(self, other):
        """tensor_product(other)

        other is a Diagram, or list thereof. It returns a list of
        Diagrams."""
        method = self._tp_choser[isinstance(other, Diagram)]
        return sorted(method(self, other))

    ## Tensor product of n diagrams is a set of new diagrams
    # This is a __MAJOR__ part of representation theory- please
    # compare to/with racah.py.
    # \param other A diagram of sequence of diagrams
    # \returns a list of diagrams
    def tensor_product(self, *args):
        """diagram.tensor_products(d1, ..., d_n)"""
        return reduce(operator.xor, (self,) + args)

    ## This shows the power of the Young Diagram, lexigraphically
    @staticmethod
    def _lexical_tensor_product(r, *args):

        # get dimensions of diagrams
        dim_r = operator.methodcaller("dimW", r=r)

        left = u" \u2297 ".join(["{}"]*len(args)).format(*map(dim_r,
                                                            args))
        
        # compute tensor sum (filter out zero dimension results)
        prod = filter(dim_r, reduce(operator.xor, args))

        # variable length right side, with filling
        right = u" \u2295 ".join(["{}"]*len(prod)).format(*map(dim_r, prod))
                           
        return u"{} = {}".format(left, right)

    @classmethod
    def ltp(cls, r):
        return functools.partial(cls._lexical_tensor_product, r)
        
    
    ## Direct assigment
    __xor__ = _tensor_product

    ## This does: list ^ Diagram
    __rxor__ = __xor__
        
    ## <a href="http://en.wikipedia.org/wiki/Young_tableau#Arm_and_leg_length">
    # Arm length</a> of cell.
    # \param i row
    # \param j column
    # \returns \f$ a_{\lambda}(i, j) \f$
    def arm(self, i, j):
        """arm length = diagram.arm(i, j)"""
        return list(self)[i-1] - j if (i, j) in self else 0

    ## <a href="http://en.wikipedia.org/wiki/Young_tableau#Arm_and_leg_length">
    # Leg length</a> of cell.
    # \param i row
    # \param j column
    # \returns \f$ l_{\lambda}(i, j) \f$
    def leg(self, i, j):
        """leg length = diagram.leg(i, j)"""
        return sum(row >= j for row in list(self)[i:])

    ## <a href="http://en.wikipedia.org/wiki/Hook_length_formula">
    # Dimension of irreducible rep of \f$ \pi_{\lambda} \f$
    # \returns \f$ d_{\lambda} = \frac{n!}{\prod_{(i, j)}{ h(i, j)_{\lambda}}}
    # \f$
    def dim(self):
        """n = diagram.dim()

        Is the number of Standard Tableau with the same shape as the diagram,
        and that is the number of generators need to generate the
        irrep."""
        return math.factorial(int(self)) / self.hookproduct()
    
    ## Hook Number
    # \param i Row index (from 1)
    # \param j Col index (from 1)
    # \returns \f$ h_{\lambda}(i, j) \f$
    def hook(self, i, j):
        """hook number = diagram.hook(i, j)"""
        return 1 + self.arm(i, j) + self.leg(i, j) if (i, j) in self else 0

    ## The Hooktableau
    # \returns Tableau filled with hook lengths of each cell.
    def hooktableau(self):
        """tableau = diagram.hook_tableau() fills each cell with its
        hook length."""
        return self.fill(*[cell.hook for cell in self.cell_info()])

    ## The Hookproduct
    # \returns int The product of all the hook lengths.
    def hookproduct(self):
        return reduce(operator.mul, self.hooktableau().cells())
        
    ## <a href="http://en.wikipedia.org/wiki/Young_tableau#Dimension_of_a_representation">
    # Dimension of W
    # \param r =3 Dimension of euclidean space
    # \returns int Dimension of tensor subspace irrep.
    def dimW(self, r=3):
        """D = dimW(r=3)
        What is it:

        The dimension, D, of the tableaus representaton in r dimensions
        (spatial, or internal symmetry).

        Apply the tableau's symmetrizer to the indices --> t'
        D is the dimension of the t's subspace.

        For exmaple:
        >>>d = Diagram(2,)  # A totally symmetric rank-2 object
        >>>print d
        [ ][ ]

        >>>print d.dimW(r=3)  # Symmetric tensor in R3 has 6 dimensions
        6

        c.f: the inertia tensor
        
        >>>d = Diagram(1, 1)  # Antisymmetric rank-2 object
        >>>print d
        [ ]
        [ ]

        >>>print d.dimW(r=4)  # Antusymm. tensor in 4D has 6 dimensions
        6

        c.f.: Electromagentic Field Stregnth tensor in 3+1 spacetime.
        """
        # get 2 copies of hook lengths in diagram
        g = iter(
            itertools.tee(((i, j, h) for (i, j, h) in
                           ((i, j, self.hook(i, j)) for i, j in
                            itertools.product(
                                range(1, 1 + int(self)),
                                repeat=2)) if h), 2))
        # take ratio of products in numerator and denominator
        return (reduce(operator.mul,
                       ((r+j-i) for i, j, h in next(g)))/
                reduce(operator.mul, (h for i, j, h in next(g))))

    ## <a href="http://en.wikipedia.org/wiki/Young_tableau#Restricted_representations">
    # Restricted Representations</a>.
    # \yields \f$S^i_{n-1} \f$ such that \f$S_n = \oplus_i S^i_{n-1} \f$
    def restricted_reps(self):
        """generate restricted representations."""
        list_ = list(self)
        for count, row in enumerate(self):
            partition = list_[:]
            if count != len(self) - 1:
                if list_[count+1] == row:
                    continue
            partition[count] -= 1
            yield type(self)(*[i for i in partition if i > 0])

    ## Permutations in irreducible rep (push to class method call on
    # Sym)
    # \yields \f$ \sigma \in \pi_n \f$
    def specht(self):
        """Yield permutations in diagram's irrep
        (aka it's Specht module)."""
        list_ = list(self)
        for sig in Sym(self):
            if list_ == list(sig.diagram()):
                if sig.tableau().isstandard():
                    yield sig

    ## \f$ {\mathrm Tab}(\lambda) \f$
    # \yields Tableau() in \f$ \pi_n \f$
    def standard_tableaux(self, _func=operator.methodcaller("tableau")):
        """Generate all tableaux in the Specht module."""
        return itertools.imap(_func, self.specht())
    
    ## Fill Diagram with arguments (row-wise)
    # \param args to sequetially fill cells
    # \returns Tableau()
    # \throws ValueError Iff wrong number of arguments.
    def fill(self, *args):
        """tableau = diagram(n1, ..., n-N)."""
        if self.skew:  # guard
            raise NotImplementedError("Skew diagram fill is not implemented")
        nargs = collections.deque(args)
        try:
            return Tableau(*[[nargs.popleft() for dummy in range(count)]
                             for count in self])
        except IndexError:
            raise ValueError("Not enough arguments to fill Diagram")
        finally:
            if nargs:  # raise (there are leftover arguments).
                raise ValueError("Too many arguments for Diagram.")

    ## Construct normal Tableax.
    # \returns \f$ A(\lambda) \f$
    def normalize(self):
        """fill diagram with increasing # running down columns from
        L -> R."""
        return self.fill(*xrange(int(self)))  # self.weight()

    ## The diagram's <a href="https://en.wikipedia.org/wiki/Young_symmetrizer">
    # symmetrizers</a>.
    # \yields tupple (tab, tab's symmetrizer) for
    def S(self):
        """yield pairs of (sign, sigma) that compose the diagram's
        symmetrizer."""
        return ((tab, list(tab.S())) for tab in self.standard_tableaux())

    ## Symmetrize a tensor."""
    # \param T
    # \yields \f$ T_{(...)} \forall T(\lambda)\in D(\lambda) \f$
    # that are standard
    def symmetrize(self, tensor_):
        """T_(n1, n2, ...) for all n = list(Diagram(n1, n2, ..)(T))

        for a rank n1+n2+... tensor, T.

        That is: call each standard tableaux.
        """
        # note: non Demeter friendly implementation?
        return (tab.symmetrize(tensor_) for tab in self.standard_tableaux())

    __call__ = symmetrize

    ## Show the content and context of each cell
    # \yields Cell() objects for each cell, in order.
    def cell_info(self):
        """Yields a Cell object describing each cell."""
        for i, j in itertools.chain(
                *itertools.starmap(
                    itertools.product,
                    [([j+1], r) for j, r in
                     enumerate([range(1, 1+x) for x in self.rows])])):
            yield Cell(i, j, self.arm(i, j), self.leg(i, j), self.hook(i, j),
                       None)

    ## Induced Representation
    def Ind(self, other):
        """Induced Representation"""
        return self ^ other

    ## Restricted Representation
    def Res(self, other):
        """Restricted Representation"""
        return NotImplemented
        
## The Skew Diagram is not even thought out-- but this class can
# create them (by adding None(s) to the left side of a rows).
class Diagram(Diagram_):
    """d = SkewDiagram(*row_lengths, skew=tuple), e.g.:

    d = Diagram(5,5,2,1, skew=(2,1))

    unicode(d) draws the diagram.
    """

    skew = ()
    
    ## Construct from row lengths, and possible skew factor
    # \param *args Row lengths
    # \param **kwargs accept 1 value: skew = Diagram
    # \throws SkewValueError
    def __init__(self, *args, **kwargs):
        super(Diagram, self).__init__(*args)
        ## Tuple of skew offsets, per rows.
        skew = kwargs.get('skew', ())
        # check validity (sorry- I hate work in init)
        if skew:  # refactor this mess.
            if skew not in self:
                self.skew = skew
                assert any([d < s for d, s in self._view()])
                raise SkewValueError("Skew diagram doesn't fit in main diagram")
            else:
                self.skew = skew
                assert not any([d < s for d, s in self._view()])
            
    def __repr__(self):
        result = super(Diagram, self).__repr__()
        if self.skew:
            result = "{}, skew={})".format(result.rstrip(')'),
                                           repr(self.skew))
        return result


## This decorator helps Tableau use Diagram's method--the inheritance
# is not entirely justifiable, as the structure are different enough-
# so this patches that problem-- an antipattern, perhaps--this
# relation is not Liskov safe.
def _usediagram(method):

    @functools.wraps(method)
    def _dmethod(self, *args):
        return method(self.diagram, *args)

    return _dmethod


## <a href="http://en.wikipedia.org/wiki/Young_tableau">
# Young Tableau</a>
class Tableau_(object):
    """tab = Tableau(*args)

    where each arg is a list of entries in a row, e.g:

    >>>tab = Tableau([0,2,3], [1])

    >>>print tab
    [0][2][3]
    [1]


    >>>print ~tab  # complement
    [0][1]
    [2]
    [3]


    The symmetrizers are computed as follows:

    >>>print list(tab.Row())  #  row-equivalent permutations
    [, (23), (02), (023), (032), (03)]

    >>>print list(tab.Col())  #  col-equivalent permutations (and their signs)
    [(1, ), (-1, (01))]

    >>>print list(tab.S())  # COl() * Row() is the total symmetrizer
    [(1, ), (1, (23)), (1, (02)), (1, (023)), (1, (032)), (1, (03)),
    (-1, (01)), (-1, (01)(23)), (-1, (012)), (-1, (0123)), (-1, (0132)),
    (-1, (013))]

    >>>print tab.lexigraphicS() # how it is applied to a tensor
    (+T.ijkl +T.ijlk +T.kjil +T.kjli +T.ljik +T.ljki -T.jikl -T.jilk -T.jkil
    -T.jkli -T.jlik -T.jlki) / 8

    # AS A FUNCTION, the Tableuax symmetrizes a tensor:
    >>>print tab(DELTA & X & Y*8)
    =================================================
    [0.0, 0.0, 0.0] [-2.0, 0.0, 0.0] [0.0, 0.0, 0.0]
    [0.0, 0.0, 0.0] [0.0, 2.0, 0.0] [0.0, 0.0, 1.0]
    [0.0, 0.0, 0.0] [0.0, 0.0, 0.0] [0.0, 1.0, 0.0]
    -------------------------------------------------
    [2.0, 0.0, 0.0] [0.0, 0.0, 0.0] [0.0, 0.0, 1.0]
    [0.0, -2.0, 0.0] [0.0, 0.0, 0.0] [0.0, 0.0, 0.0]
    [0.0, 0.0, 0.0] [0.0, 0.0, 0.0] [1.0, 0.0, 0.0]
    -------------------------------------------------
    [0.0, 0.0, 0.0] [0.0, 0.0, -1.0] [0.0, 0.0, 0.0]
    [0.0, 0.0, -1.0] [0.0, 0.0, 0.0] [0.0, 0.0, 0.0]
    [0.0, -1.0, 0.0] [-1.0, 0.0, 0.0] [0.0, 0.0, 0.0]
    =================================================


    Finally, the R/S Correspondance matches the tableaux with a permutation:

    >>>print tab.pi()
    | (0, 1, 2, 3) | 
    | (2, 1, 3, 0) | 
    """

    ## Permutations start here-- this may be deprecated
    start = 0

    involution_number = staticmethod(Diagram_.involution_number)
    
    ## Find a standard Tableau's Yamanouchi symbol (see
    # https://docs.python.org/2/howto/sorting.html#key-functions for
    # algorithm).
    # \returns tuple matching tuple index to value (to row).
    # \throws TypeError
    def yamanouchi(self, _func=operator.itemgetter(1)):
        """Convert standard tableau to a Yamanouchi symbol. Returns
        a tuple maping Yamanouchi index to Tableau row (which uniquely
        defines a standard Tableau)."""
        if not self.isstandard():  # raise
            raise ValueError("Only standard tableaux can be converted")
        return tuple(
            map(_func,
                sorted(item.yamanouchi() for item in self.cell_info())))
    
    ## Construct from rows
    # \param args iterates rows
    def __init__(self, *args):
        ## labled rows that are self.
        self.rows = map(list, args)  # assure it's made of lists.

    ## \yields rows
    def __iter__(self):
        return iter(self.rows)

    ## \returns Number of rows.
    def __len__(self):
        return len(self.rows)

    ## Should we sort by yamamochi number or lethem rank?
    def __cmp__(self, other, func=lambda dum: dum.pi().rank()):
        return (cmp(self.diagram, other.diagram) or
                cmp(func(self), func(other)))
    
    def __repr__(self):
        return type(self).__name__ + repr(self.rows)

    ## Functional skewness-enabled string function.
    def __str__(self):
        return "\n".join(
            map("".join,
                map(functools.partial(
                    map,
                    lambda n: "   " if n is None else "[{}]".format(
                        int2alpha(n))),
                    self)))

    ## Get a row, or a cell.
    def __getitem__(self, index):
        try:
            i, j = index
        except TypeError:
            if index == 0:  # raise
                raise IndexError("count from 1")
            return self.rows[index-1]
        if i == 0:  # raise
            raise IndexError('count from 1')
        return self.rows[i-1][j-1]

    ## Complementary Table (may not perserve standard-ness)
    # \param _func f(x) = filter(y != None, x)
    # \returns Tableau (the complementary one).
    def __invert__(self,
                   _func=functools.partial(filter,
                                           functools.partial(operator.is_not,
                                                             None))):
        """This flips the Tableau as follows:

        >>>g = partial(is_not, None)  # Is boolean, True For not None
        >>>f = partial(filter, f1) # keeps not None values in a list

        So start with a Tableau, e.g.:
        >>>t = Diagram(5,3,1).normalize()

        >>>print t
        [i][j][k][l][m]
        [n][o][p]
        [q]

        >>>list(t)
        [[0, 1, 2, 3, 4], [5, 6, 7], [8]]

        Now flip it, while filling with Nones:
        >>>list(izip_longest(*t))
        [(0, 5, 8),
         (1, 6, None),
         (2, 7, None),
         (3, None, None),
         (4, None, None)]

        And map that list to "f" (which is _func in the private keyword):

        >>>map(f2, izip_longest(*t))
        [(0, 5, 8), (1, 6), (2, 7), (3,), (4,)]

        then unpack that into Tableau's constructor:

        >>>print Tableau((0, 5, 8), (1, 6), (2, 7), (3,), (4,))
        [i][n][q]
        [j][o]
        [k][p]
        [l]
        [m]

        Boom: that is the complementary Tableau."""
        return Tableau(*map(_func, itertools.izip_longest(*self)))

    ## The weight of a Tableaux counts the multiplicity of cell values.
    # \returns \f$ \mu \f$ is a tuple of multiplicities
    def weight(self):
        """The weight is a tuple of cell multiplicities"""
        return tuple(
            reversed(sorted(collections.Counter(self.cells()).keys())))

    def weight(self):
        """The weight is a tuple of cell multiplicities"""
        return tuple(
            reversed(sorted(collections.Counter(self.cells()).keys())))
    
    ## Test for skewness
    # \returns bool IFF Tableau is skew.
    def isskew(self):
        """Bool tests skewness"""
        return None in list(self.cells())

    ## Semistandard: \n Weakly Inceasing Rows \n Strictly Incearsing Columns
    # \returns bool Iff Tableau is semistandard
    def issemistandard(self):
        """True Iff tableau is semistandard:

        Weakly increasing rows,
        strictly increasing columns."""
        return self.is_row_weak() and self.is_col_strict()

    ## Are all rows non-decreasing?
    # \return True IFF all rows are non-decreasing.
    def is_row_weak(self):
        """True IFF rows are non-descreasing."""
        return all(list(itertools.starmap(is_nondecreasing, self)))
    
    ## Are all rows increasing?
    # \return True IFF all rows are increasing.
    def is_row_strict(self):
        """True IFF rows are increasing."""
        return all(list(itertools.starmap(is_increasing, self)))

    ## Are all columns increasing?
    # \return True IFF all columns are increasing.
    def is_col_strict(self):
        """True IFF columns are increasing."""
        return (~self).is_row_strict()
        
    ## Standard: \n Stricty Inceasing Rows/Cols \n Filling is of (0,...,n-1)
    # \returns bool Iff Tableau is standard
    def isstandard(self):
        """True Iff tableau is standard."""
        # 1: verify T and ~T are semistandard.
        # 2: verify filling numbers are (0,...,n-1)
        return (self.issemistandard() and (~self).issemistandard() and
                sorted(list(self.cells())) == range(int(self)))

    ## Normal: Ordered Filling.
    # \returns bool Iff Tableau is normal (by comparing cells to range(N)).
    def isnormal(self):
        """True Iff tableau is normal"""
        # all is better, as direct comparison to range is type (list/tuple)
        # sensitive, and may lead to errors: we care about the VALUES in the
        # cells-- or one can check for standard and then verify that
        # yamanouchi() is sorted-- but it is 50 times slower.
        return all(map(operator.eq, self.cells(), range(int(self))))

    ## Robinson Schested Correspondence:
    # \yields sn.Cycles
    def cycles(self):
        """Yields cycles from rows."""
        return itertools.imap(Cycle, self)

    ## Robinson Schested Correspondence:
    # \returns \f$ \sigma \f$ corresponding monte.Sym permutations
    def pi(self):
        """Returns Sym permutation corresonding to tableau."""
        return Perm.fromcycles(*self.cycles())

    ## Iterate over cell values
    # \yields cells values.
    def cells(self):
        """yields cell values (raster-scan)."""
        return itertools.chain(*self)

    ## Add Tableua value to Diagram.cell_info()
    # \yields Cell
    def cell_info(self):
        """detailed infor about each cell"""
        for value, cell in itertools.izip(self.cells(),
                                          self.diagram.cell_info()):
            yield Cell(*cell[:-1], value=value)  # can't set value in tuple.

    ## Permute with a permutation
    # \param sigma \f$ \sigma_n\f$ or a callable permutation operator
    # \returns Tableau
    def permute(self, sigma):
        """Permute VALUES in Tableau, NOT positions."""
        c = list(self.cells())
        s = Perm(*c)
        # self.diagram.fill(*s(sig((~s)(c)))) ... is the same as:
        return self.diagram.fill(*s(list(sigma)))

    ## Row Symmetrizer
    # \yields \f$ \sigma \f$ for \f$ \sigma \in R \f$
    def Row(self):
        """Row Symmetrizer:

        [s_i, ...] = tableau.Row()

        are the permutations that leave tableau row equivalent."""
        return itertools.ifilter(self.isinR, Sym(int(self)))

    ## Column Symmetrizer (signed)
    # \yields \f$ ({\mathrm {sgn}}(\sigma), \sigma) \f$
    # for \f$ \sigma \in C \f$
    def Col(self):
        """Column Symmetrizer:

        [(p, s_i), ...] = tableau.Col()

        are the permutations that leave tableau col equivalent.
        Unlike the row symmetrizer, which yields a perm 'sig',
        the yields a tuple:
        -sig --> (sgn(sig), sig)"""
        return itertools.imap(operator.neg,
                              itertools.ifilter(self.isinC,
                                                Sym(int(self))))

    ## <a href="http://en.wikipedia.org/wiki/Young_symmetrizer">The
    # Symmetrizer</a>
    # \yields tuples of (sgn(pi), pi) for \f$ \pi \in S_n \f$.
    def S(self):
        """The Young Symmetrizer: S_A => (signed-Col) * Row

        So this yields (sign, sigma) pairs.
        """
        return ((p, c*r) for (p, c), r in
                itertools.product(self.Col(), self.Row()))

    ## Print Tensor equation
    # \kwd variable =T is the variable in the string.
    # \return str
    def lexigraphicS(self, variable='T', _sign=' +-', latex=False):
        """T_ijk + T... = tableau.lexigraphicS(
        variable='T',
        _sign=' +-',
        latex=False)."""
        indices = "".join(chr(105+n) for n in range(self))
        if latex:  # kwd
            dot = "_"
            front = "{"
            back = "}"
        else:
            dot = "."
            front = back = ''
        T = "{}".format(variable)
        T += dot
        
        return ("(" +
                " ".join(_sign[s] + T + front + ''.join(p(indices)) + back
                         for s, p in self.S()) +
                ") / {}".format(math.factorial(self) / self.diagram.dim()))

    ## Is a permutation in the Row symmetrizer?
    # \param sig
    # \returns bool
    def isinR(self, sig):
        """True IFF sig is in tableaux's Row symmetrizer"""
        #return self | self.permute(sig)
        return self.is_row_equiv_to(self.permute(sig))

    ## Is a permutation in the Col symmetrizer?
    # \param sig
    # \returns bool
    def isinC(self, sig):
        """True IFF sig is in tableaux's Col symmetrizer"""
        return self.is_col_equiv_to(self.permute(sig))

    ## Row Equivalence
    # \param self (implicit) \f$ t_1 \f$
    # \param other (explicit) \f$ t_2 \f$
    # \return \f$ t_1 ~ t_2 \f$
    def is_row_equiv_to(self, other):
        """bool = t.is_row_equivalent_to(t')"""
        return (len(self) == len(other) and
                map(sorted, self) == map(sorted, other))

    ## Col Equivalence
    # \param self (implicit) \f$ t_1 \f$
    # \param other (explicit) \f$ t_2 \f$
    # \return \f$ \bar{t}_1 ~ \bar{t}_2 \f$
    def is_col_equiv_to(self, other):
        """bool = t.is_col_equivalent_to(t')"""
        return (~self).is_row_equiv_to(~other)

    ## Row equivalence: q | p
    # \return bool Iff q is row equivlent to p.
    def __or__(self, other):
        """is_row_equiv_to"""
        return self.is_row_equiv_to(other)
    
    arm = _usediagram(Diagram.arm)
    leg = _usediagram(Diagram.leg)
    hook = _usediagram(Diagram.hook)
    shape = _usediagram(Diagram.shape)
    __unicode__ = _usediagram(Diagram.__unicode__)
    __long__ = _usediagram(Diagram.__long__)
    __int__ = _usediagram(Diagram.__int__)
    __len__ = _usediagram(Diagram.__len__)

    ## Tensor symmetrizer w/o weight factor.
    # \param self \f$ \lambda \f$
    # \param tensor_ \f$ {\bf T} \f$
    # \returns \f$ P_{\lambda}({\bf T}) \f$.
    # \throws RankError
    def _symmetrize(self, tensor_):
        """It appears we take the sum over all the symmetrizers
        in the diagrams tableaux?

        but we may need to count dim() again
        """
        if int(self) != tensor_.rank:  # raise
            raise RankError(int(self), tensor_.rank)
        return reduce(operator.add, (sgn * tensor_.ntranspose(sig)
                                     for sgn, sig in self.S()))

    ## Symmetrize a euclid.tensor.Tensor
    # \param tensor_ \f$ {\bf T} \f$
    # \returns \f$ {\bf T}_{(\lambda_1, ..., \lambda_n)} \f$
    # \throws RankError
    def symmetrize(self, tensor_):
        """T_(tab) = tab.symmetrize(tensor_)

        symmetrized the tensor.
        """
        if tensor_.rank != int(self):  # raise
            raise RankError(int(self), tensor._rank)
        return self.diagram.dim() * self._symmetrize(tensor_) / long(self)

    ## As a function, Tableau are Tensor symmetrizing operators.
    __call__ = symmetrize


## The Most General tableau
class SkewTableau(type):
    """SkewTableau is TBD"""

    ## Construtor caller.
    # \returns young.Tableau_ instance.
    def __new__(mcs, *args):
        # Instantiate
        self = Tableau_(*args)
        ## young.Diagram matching Tableuax
        self.diagram = Diagram(*map(len, args))
        # verify shape (may require diagram to exist).
        mcs._verify(self)
        return self

    ## Method to ensure arguments to constructor are valid.
    # \param tab The Tableau under construction.
    # \throws ValueError Iff tableau tries to be skew, but fails.
    @classmethod
    def _verify(mcs, tab):
        """Verify that a skew tableau has the proper shape."""
        if tab.isskew():  # input validate (pre-raise).
            nones = [row.count(None) for row in tab]
            if list(reversed(nones)) != sorted(nones):  # raise
                raise ValueError("Skew shape is not a 'parallelogram'")


## General Young Tableaux, without Skewness.
class Tableau(SkewTableau):
    """Tableau([0,1,2],[3,4]) (for example).

    Enter list of integers to rep rows of cells."""

    ## VerIfy that it is NOT skew.
    # \throws TypeError on skew tableaux.
    @classmethod
    def _verify(mcs, tab):
        if tab.isskew():  # raise
            raise TypeError("Tableau is Skew")

        
## The Tableau metaclass's sole purpose is to enforce rules on the
# filling of Diagrams.
class SemiStandardTableau(Tableau):
    """SemiStandardTableau([0,1,2],[3,4]) (for example).

    Enter list of integers to rep rows of cells.

    Entries *must* form a semi-standard tableaux."""

    ## VerIfy that it is semi-standard
    # \throws ValueError non-semi-standard tableaux.
    @classmethod
    def _verify(mcs, tab):
        # call super: verify that it is not skew
        super(SemiStandardTableau, mcs)._verify(tab)
        if not tab.issemistandard():  # raise
            raise TypeError("Tableau is not semistandard")


## Tableau metaclass ensures it is Standard
class StandardTableau(SemiStandardTableau):
    """StandardTableau([0,1,2],[3,4]) (for example).

    Enter list of integers to rep rows of cells.

    Entries *must* form a standard tableaux."""

    ## VerIfy that it is standard
    # \throws TypeError non-standard tableaux.
    @classmethod
    def _verify(mcs, tab):
        # call super: verify semi-standard
        super(StandardTableau, mcs)._verify(tab)
        if not tab.isstandard():  # raise
            raise TypeError("Tableau is not Standard")
        
    @staticmethod
    def fromyamanouchi(args):
        # create empty (and DIFFERENT) rows
        rows = [list() for dum in range(max(args))]
        for x, y in enumerate(args):
            rows[y-1].append(x)
        return Tableau(*rows)


## Normal Tableau factory
class NormalTableau(StandardTableau):
    """NormalTableau([0,1,2],[3,4]) (for example).

    Enter list of integers to rep rows of cells.

    Entries *must* form a normal tableaux."""

    ## VerIfy that it is normal
    # \throws TypeError non-normal tableaux.
    @classmethod
    def _verify(mcs, tab):
        # call super: verify standard
        super(NormalTableau, mcs)._verify(tab)
        if not tab.isnormal():  # raise
            raise TypeError("Tableau is not Normal")



## Generate rows of
# a <a href="https://en.wikipedia.org/wiki/Young%27s_lattice">
# Young's Lattice</a>.
def lattice():
    """yield rows of Young's Lattice."""
    y, basis = itertools.repeat(Diagram(1), 2)
    yield [basis]
    while True:
        y ^= basis
        while any(map((1).__sub__, map(y.count, y))):
            print len(y)
            for item in y:
                if y.count(item) > 1:
                    y.pop(y.index(item))
        yield y


## The minimnal Young Diagram (using lattice() generator)
BOX = next(lattice())[0]
        
def otest():
    yield Diagram(1, 1, 1)
    yield Diagram(1)
    yield Diagram(3, 1)
    yield Diagram(2, 2)

    yield Diagram(2)
    yield Diagram(1, 1)
    yield Diagram(2, 1, 1)
    yield Diagram(3)
    yield Diagram(2, 1)
    yield Diagram(1, 1, 1, 1)

    yield Diagram(4)

        
