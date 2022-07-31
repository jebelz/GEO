"""This module has the python structure of spherical tensors.

They are containers of the following

j     L: weight (as in 2j+1)-- maybe degree, as in Y_lm (with l=j)
q     seniority index: after rank 2, there are degeneracies.
m     the order (as per usual).

and this module is about shuffling values to the right spot.

In general, you get items via:

T[j, q, m] or T[j, q] for all m or T[j] for all m and q.

but at rank-2 (or less), there is only 1-q, so that:

T[j, m] or T[j] for all m

while for a vector, there's also only 1 j:

V[m]

and for the lowly, scalar, j=q=m=0:

S[0]

The construtor is a bit more annoying for 2-reasons:

(1) All the different j, q, m configs are a mess

(2) -j <= m <= j so that negative m is rampant. To make this match python's
indexing convention, representations (that's what you call a j and all its
m's) are filled like this:

[0, 1, 2, ..., j, -j, ..., -1]

which is not in numerical order. It is ordered a lot like the frequency in
the standard FFT package, so it should not be totally unfamilar.

Let's work it in reverse:

s = Scalar(w)

v = Vector(zero, plus, minus)

At rank-2 we have to add depth to out containers:

t = Tensor(w, [zero, plus, minus], [0, 1, 2, -2, -1])

At higher rank, another layer is added to deal with degeneracy (the seniority
index).

There are 2 methods to help you along:

jqm() which returns a nested list in the shape that the contructor expects.

Meanwhile:

ijqm()

is a generator, and it yields triplets:

(j, q, m)

in the ORDER they are expected, flattened. All tensor classes have a
pack classmethod that will construct, as needed.
"""

## \namespace geo.metric.wigner.eckart.ec Python Structure of Spherical Tensors
import abc
import itertools

from ....utils.exceptions import RepresentationError, NilpotencyError
from ..casimir import jiter, miter, is_in_N
from ..racah import qn, str2unicode, Nice

def c2r(method):
    def c2r_method(self, *args):
        result = method(self, *args)
        try:
            im = result.imag
        except AttributeError:
            pass
        else:
            if im == 0:
                result = result.real
        return result
    return c2r_method



## Out-of bounds Seniority Index (counting from ZERO).
class SeniorityIndexError(IndexError, UserWarning):
    """Raised when seniority index request is garb."""


## Error for \f$ l < 0 \f$
class CasimirOperatorError(UserWarning, ValueError):
    """Garbage Total Angular Momentum"""


## The base class for container.
class SphericalTensor(object):
    """The spherical tensors as a container of (a container of (
    a container of)) Z-rotation eigenstates."""

    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def L(self):
        """All concretions need an weight (ala total angular momentum)."""

    ## For now
    def __repr__(self):
        return "{}=\n{}".format(type(self), str(self))

    ## Yields valid weights
    # \yields j for all weights
    @classmethod
    def jiter(cls):
        """yield weights, j."""
        return jiter(cls.L)

    ## Generate flattened (j, q, m) triplets: this is a VERY IMPORTANT
    # method, as it tells us the order in which spherical tensors are
    # built: compare with jqm().
    # \yields (j, q, m) in order.
    @classmethod
    def ijqm(cls):
        """Yields ordered (j, q, m) triplets."""
        for j in map(int, cls.jiter()):
            for q in xrange(cls.degeneracy(j)):
                for m in miter(j):
                    yield j, q, m

    ## Iterate azimuth order parameter for heightest weight
    # \yields m according to casimir.miter()
    @classmethod
    def miter(cls):
        """Yields m for maximum j (=cls.L) according to correct order.""" 
        return miter(cls.L)
        
    ## Get [j[q[m]]] nested struture for any rank.
    # \kwd container=tuple Chooses the inner-most container.
    # \returns list[j][q][m] (nested)
    @classmethod
    def jqm(cls, container=tuple):
        """Returns a j-wise list of q-wise list of tuples of m values"""
        J = []
        for j in map(int, cls.jiter()):
            Q = []
            for dummy in xrange(cls.degeneracy(j)):
                M = container(miter(j))
                Q.append(M)
            J.append(Q)
        return J

    ## Linear Iteration in the order of ijqm().
    def __iter__(self):
        """imap(self.__getitem__,  self.ijqm())"""
        return (self[j, q, m] for j, q, m in self.ijqm())

    ## Construct from j-wise lists of q-wise lists of miter-wise tuples.
    # \param args is a nested j->q->m container of (j, q, m)
    def __init__(self, *args):
        for j, jargs in enumerate(args):
            for q, qargs in enumerate(jargs):
                for m, marg in zip(miter(j), qargs):
                    self.__dict__[(j, q, m)] = marg

    def __getslice__(self, *args):
        print args
        raise Exception("slice is not impelemted")

    ## If you ask for (j, q, m), it's easy peasy; you want special help:
    # you get cyclomatic cmoplexity and non-local flow.
    @c2r
    def __getitem__(self, index):
        try:
            result = self.__dict__[index]
        except KeyError:
            try:
                j, q, m = index
            except TypeError:  # dump j=index
                # this is __getslice__ request: self[j, :qmax, miter(j)]
                result = [[self[index, q, m] for m in miter(index)]
                          for q in range(self.degeneracy(index))]
                if not result:  # raise
                    raise RepresentationError("j={}?".format(index))
            except ValueError: # dump j, q = index
                try:
                    j, q = index
                except TypeError:
                    raise IndexError(
                        "T[j, q, m] is not T[{}]".format(", ".join(
                            map(str, index))))
                # this is a __getslice__ request running over m
                result = [self[j, q, m] for m in miter(j)]
            else:
                if not all(map(is_in_N, index)):  # raise
                    raise ValueError("Common Man! Integer only.")
                if j < 0 or j > self.L:  # raise
                    raise RepresentationError(
                        "N={} j={} is nonsense".format(self.L, j)
                    )
                if not 0 <= q < self.degeneracy(j):  # raise
                    raise SeniorityIndexError(
                        "Seniority index violates: 0 <= {} < {}".format(
                            q, self.degeneracy(j)))
                if not -j <= m <= j:  # raise
                    raise NilpotencyError(" | j={}, m={}>... really?".format(
                        j, m))
                raise ValueError(
                    "Failed to Dx rank {} with  [j, q, m]={}".format(self.L,
                                                                     index))
        try:  # TODO FIX THIS MESS
            return result.w
        except AttributeError:
            return result

    _memo_str = ""

    def __str__(self):
        return self._memo_str
    
    def super_str(self):
        if not self._memo_str:
            self._memo_str = self._str(R=lambda _s: str(Nice(_s)))
        return self._memo_str
            
    def __unicode__(self):
        return str2unicode(self.super_str())

    ## Pack a stream into correct format:
    # \param args
    # \returns instance of cls, with args reformated
    @classmethod
    def _pack(cls, *args):
        result = cls.jqm(container=list)
        for value, (j, q, m) in zip(args, cls.ijqm()):
            result[j][q][m] = value
        return cls(*result)

    ## See _pack().
    # \param args Each one is a |j, m>_q value, in jqm() order.
    # \returns instance of cls
    # \throws ValueError If the number of args is not \f$3^L\f$.
    @classmethod
    def pack(cls, *args):
        """t == T._pack(*t)

        That is, pack reforms a 1 deep list  in ijqm() order to look
        like jqm(), which is what the constructor expects."""
        if len(args) != pow(3, cls.L):  # raise
            msg = "Expected {} arguments, got: {}"
            raise ValueError(msg.format(pow(3, cls.L), len(args)))
        return cls._pack(*args)

    ## pack() a list.
    # \param list_ of |j, m>_q values.
    # \returns instance of class.
    @classmethod
    def fromlist(cls, list_):
        """instance = cls.fromlist(list_)

        list is unraveled."""
        return cls._pack(*list_)

    ## Construct a basis (eigen-) state for rank L.
    # \param j weight or degree (total angular momentum quantum number)
    # \param m order (Z-projection angular momentum quantum numner)
    # \param q seniority index
    # \returns \f$ |j, m\rangle_q \f$
    @classmethod
    def eigenstate(cls, j, q, m):
        """|j, m>_q = basis(j, m, q)"""
        return cls.fromlist(map(int, map((j, q, m).__eq__, cls.ijqm())))

    ## Construct ALL Basis states.
    # \returns JQM List of List of List of eigenstate()
    @classmethod
    def basis(cls):
        """basis()[j][q][m] == eigenstate(j, q, m)"""
        result = cls.jqm(list)
        for j, q, m in cls.ijqm():
            result[j][q][m] = cls.eigenstate(j, q, m)
        return result

    ## Get each polyadic (summing to self)
    # \returns lists \f$ T_j^{m; q} \forall (j, q, m) \f$
    def polyadics(self):
        """polyadics()[j][q][m] == self[j, q, m]*eigenstate(j, q, m)"""
        result = self.jqm(list)
        for j, q, m in self.ijqm():
            result[j][q][m] = self[j, q, m] * self.eigenstate(j, q, m)
        return result
        

## Handle j = q = m = 0
class Scalar(SphericalTensor):
    """Scalar contrainer"""

    ## Scalars are still not numbers.
    L = 0

    ## Scalar(w) --> Super([[w]])
    def __init__(self, w):
        # pack w into m=0 list and then a q=0 list.
        super(Scalar, self).__init__([[w]])

    ## scalar[0] == scalar[0, 0] = scalar[0, 0, 0]
    def __getitem__(self, index):
        try:
            length = len(index)
        except TypeError:
            # index is [m]
            index = (index, 0, 0)
        else:
            # index is [j, m]
            if length == 2:  # usage
                index = (index[0], 0, index[1])
        # call super with full j, q, m.
        return super(Scalar, self).__getitem__(index)

    ## String is a scalar in brackets.
    def __str__(self):
        return "{}".format(str(Nice(self[0])))

    
    ## Pack 1 value into a Scalar.
    @classmethod
    def _pack(cls, *args):
        return cls(args[0])

    ## Complex value
    def __complex__(self):
        return complex(self[0])

    ## float value
    def __float__(self):
        return float(self[0])

    ## Integer value.
    def __int__(self):
        return int(float(self))

        
## Handle j=1, q=0.
class Vector(SphericalTensor):
    """Vector Container"""

    ## Vectors have 1-unit of angular momentum.
    L = 1

    ## Vector(0, 1, -1) --> Super([[]], [[0, 1, -1]])
    def __init__(self, e0, e_plus, e_minus):
        # pack empty scalar, and then put vector into q=0 list.
        super(Vector, self).__init__([[]], [[e0, e_plus, e_minus]])

    ## Get "m", since j=1 and q=0, always.
    def __getitem__(self, index):
        try:
            length = len(index)
        except TypeError:
            # index is m, so insert j=1 and q=0.
            index = (1, 0, index)
        else:
            if length == 2:  # usage
                # index is [j, m], so insert q=0.
                index = (index[0], 0, index[1])
        # call super with j=1, q=0, m=m.
        return super(Vector, self).__getitem__(index)

    ## String is a column of kets with decreasing m
    # \param _template="{}|+{}>\n{}|{}>\n{}|{}>" is a private template
    # \returns \n
    # <v|1>|+1> \n
    # <v|0>|0>\n
    # <v|1>|-1>
    def __str__(self):
        return super(Vector, self).__str__() or self._str()

    def _str(self, _template="{}|+{}>\n{}|{}>\n{}|{}>", R=lambda _s: _s):
        return _template.format(
            *itertools.chain(
                *[(R(self[m]), m)
                  for m in reversed(sorted(list(self.miter())))]))

        
    ## Pack arguments into a Vector.
    @classmethod
    def _pack(cls, *args):
        return cls(args[0], args[1], args[-1])


## Handle j = (0, 1, 2) with q = 0.
class Tensor(SphericalTensor):
    """Tensor Container: 3 x 3 = 5 + 3 + 1"""

    ## 3 x 3
    L = 2

    ## Tensor(0, list_3, list_5) --> Super([[0]], [list_3], [list_5]
    def __init__(self, rank0, rank1, rank2):
        # pack scalar in m=0 list and then pack them all in q=0 list.
        super(Tensor, self).__init__([[rank0]], [rank1], [rank2])

    ## get [j, m], with q=0, always, or get all m-values for fixed j.
    def __getitem__(self, index):
        try:
            length = len(index)
        except TypeError:  # process a skipped q
            # dump all m at constant j=index.
            return [self[index, 0, m] for m in miter(index)]
        ## add q=0 layer
        if length == 2:
            index = (index[0], 0, index[1])
        return super(Tensor, self).__getitem__(index)

    ## complex formatting function shows J=0, J=1, J=2 Irreps.
    def _str(self,
             template=" | {}> ",
             col=('m \\', '+2: ', '+1: ', ' 0: ', '-1: ', '-2: '),
             R = lambda _s: _s):
        k = [""] * 5
        k[2] = template.format(R(self[0, 0]))
        def just(n, sign):
            """local justify function"""
            from operator import methodcaller
            name = ('r' if sign > 1 else 'l') + 'just'
            func = methodcaller(name, n)
            for count, item in enumerate(k):
                k[count] = func(item)

        just(len(k[2]), 1)
        for index, m in enumerate((1, 0, -1)):
             k[index + 1] += template.format(R(self[1, m]))
        just(max(map(len, k)), -1)
        just(max(map(len, k)), 1)
        for index, m in enumerate((2, 1, 0, -1, -2)):
            k[index] += template.format(R(self[2, m]))
        just(max(map(len, k)), -1)
        just(max(map(len, k)), 1)
        top = " J=0"
        top = top.ljust(k[1].index(template[0]))
        top += "     J=1"
        top = top.ljust(k[0].index(template[0]))
        top += "     J=2"
        k.insert(0, top)
        return "\n".join(
            itertools.starmap(
                "{}{}".format,
                zip(col,
                    [type(item).rstrip(item) for item in k])))

    ## Ylm format
    def __str__(self):
        """print Tensor([1], [3,4,2], [7,8,9,5,6])
        m \ J=0  J=1  J=2
        +2:        |9>
        +1:    |4> |8>
         0:|1> |3 >|7>
        -1:    |2> |6>
        -2:        |5>
        """
        return self._str()

    @classmethod
    def _pack(cls, *args):
        return cls(args[0], args[1:4], args[4:])


## Rank 3 spherical tensor
class Three(SphericalTensor):
    """Rank 3 container:
    3 x 3 x 3 = 7 + 5 + 5 + 3 + 3 + 3 + 1"""

    ## Rank 3
    L = 3

    ## complex formatting function (deprecated)
    def _str(self,
             template=" | {}> ",
             col=('m \\', '+3 ', '+2: ', '+1: ', ' 0: ', '-1: ', '-2: ', '-3')):
        k = [""] * 7
        k[2] = template.format(self[0, 0])
        def just(n, sign):
            """local justify function"""
            from operator import methodcaller
            name = ('r' if sign > 1 else 'l') + 'just'
            func = methodcaller(name, n)
            for count, item in enumerate(k):
                k[count] = func(item)


                
        just(len(k[2]), 1)
        for index, m in enumerate((1, 0, -1)):
            k[index + 1] += template.format(self[1, 0, m])
        just(max(map(len, k)), -1)
        just(max(map(len, k)), 1)

        for index, m in enumerate((2, 1, 0, -1, -2)):
            k[index] += template.format(self[2, 0, m])
        just(max(map(len, k)), -1)
        just(max(map(len, k)), 1)

        for index, m in enumerate((3, 2, 1, 0, -1, -2, -3)):
            k[index] += template.format(self[3, 0, m])
        just(max(map(len, k)), -1)
        just(max(map(len, k)), 1)

        top = " L=0"
        top = top.ljust(k[1].index(template[0]))
        top += " L=1"
        top = top.ljust(k[0].index(template[0]))
        top += " L=2"
        k.insert(0, top)
        return "\n".join(
            itertools.starmap(
                "{}{}".format,
                zip(col,
                    [type(item).rstrip(item) for item in k])))
        top += " L=3"
        k.insert(0, top)

    ## Ylm format
    def __str__(self):
        s = "[][][]\n"
        s = ''
        # [][][] module
        for j in (3, 1):
            s += "\nL={}, q=0\n".format(j)
            s += ("{} "*(2*j+1)).format(*self[j, 0]) + '\n'
        # [] []
        # []
        #        s += "[][]\n[]:\n"
        for q0 in (0, 1):
            for count, j in enumerate([2, 1]):
                q = q0 + count
                s += "\nL={}, q={}\n".format(j, q)
                s += ("{} "*(2*j+1)).format(*self[j, q]) + '\n'

        s += "\nL=0: {}".format(self[0, 0, 0])
        return s


## The rank-4 spherical tensor.
class Four(SphericalTensor):
    """Rank 4 container:
    3 x 3 x 3 x 3 = 9 + .... + 1"""

    ## Rank 4
    L = 4
