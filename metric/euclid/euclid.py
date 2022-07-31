"""This module is full of abstract base classes for Tensors and maps on
E3.

The ZOO and its VIEW are here. The ZOO is full of geometric animals living
in E3.
"""
## \namespace geo::metric::euclid::euclid Base classes for Geometric Animals
#  living in
#  <a href="http://en.wikipedia.org/wiki/Euclidean_space">\f$R^3\f$</a>
import abc
import collections
import functools
import itertools
import operator

from ...utils import cauchy
from ...utils import arraylike
from ...utils import keydefaultdict
from ...utils import exceptions

from . import einstein

__all__ = ('AXES', 'DIMENSIONS', 'rank', 'ranked', 'ZOO', 'VIEW', 'Tensor_')

AXES = einstein.AXES

DIMENSIONS = len(AXES)


## Catch bad calls for chart type conversion.
class ChartError(UserWarning, TypeError):
    """Raised if you try to convert a chart to a non-existent type"""


## get the rank from the class of the argument, or None.
# \param tensor_
# \returns rank or None
def rank(tensor_):
    """get rank attribute or None"""
    try:
        return tensor_.rank
    except AttributeError:
        pass


## Wrapper for binary ops between ranked objects
# \param method A binary operation
# \returns method
# \throws GeometricException
def wrank(method):
    """run a tensor method that takes tensor arg with
    exception handling."""
    
    @functools.wraps(method)
    def method_rank(self, other):
        try:
            return method(self, other)
        except AttributeError:
            raise exceptions.GeometricException(
                "Expected Ranked Object, got: {}".format(type(other)))
    return method_rank
    

## <a href="http://en.wikipedia.org/wiki/Norm_(mathematics)#p-norm">
# p-norm</a>
# \param p \f$ p \f$
# \param x \f$ x_{i\ldots k} \f$ (any rank)
# \returns scalar.Scalar \f$ ||{\bf x}||_p= ( \sum_{i\ldots k}
# { |x_{i\ldots k}|^p } )^{\frac{1}{p}} \f$
def pnorm(p, x):
    """s = pnorm(p, x) = ||x||_p"""
    return ZOO[0](reduce(operator.add,
                         (abs(item)**p for item in x.iter())) ** (1./p))


## Left Dilation of an arbitrary Tensor
# \param s a Scalar
# \param t a Tensor with a dilation methods
# \returns result dilated tensor
# \throws exceptions.GeometricException
def _scalar_times_tensor(s, t):
    """dilation via arguments dilation method"""
    try:
        result = t.left_dilation(s.w)
    except AttributeError:
        raise exceptions.GeometricException(
            "arguments must be a scalar and a dilatable")
    return result


## Detracer Operator \n
# \f$ D_n:T^{(n)}_{i_0\cdots i_{n-1}} = \frac{1}{(2n-1)!!}
# \sum_{m=0}^{\lfloor{\frac{n}{2}}\rfloor}{(-1)^m(2n-2m-1)!!
# \sum_{Sym_n(i)}{\delta_{i_0i_1}\cdots\delta_{i_{2m-2}i_{2m-1}}
# T^{(n)}_{j_0j_0\ldots j_{m-1}j_{m-1}i_{2m\ldots i_{n-1}}}}} \f$
class DeTracer(object):
    """dt = DeTracer(T)
    will compute the detraced T, where T is a symmetric tensor

    of rank n.

    >>>dt[m]  # for 0 <= m <= n//2 is the m-th term in the sum, so that:

    >>>list(dt)  # is a list of all the terms (unnormalized)

    >>>T_m = dt()   # is sum(dt) / (2n-1)!!
    """
    
    def __repr__(self):
        return repr(self.tensor)
        
    def __init__(self, tensor_):
        self.tensor = tensor_

    ## Rank
    def __int__(self):
        return self.tensor.rank

    ## Outer products of Kronecker ::DELTA
    # \param m
    # \return \f$ \delta_{i_1j_1} \ldots \delta_{i_m j_m} \f$
    @staticmethod
    def delta(m):
        """D&D...(m-times)..&D = delta(m)"""
        from operator import and_
        from itertools import repeat
        from .tensor import DELTA
        return reduce(and_, repeat(DELTA, m)) if m else 1

    ## Double factorial short-hand
    # \param k
    # \returns \f$ (2k-1)!! \f$
    @staticmethod
    def func(k):    
        from geo.metric.schur_weyl.pascal import factorial2
        return float(factorial2(2*k - 1))
        
    ## Compute free indices
    def ifree(self, m):
        from geo.metric.euclid.einstein import RUN
        result = RUN[m:][:int(self)]
        return result[:2*m], result[2*m:]
        
    ## Compute trace indices
    @staticmethod
    def itrace(m):
        from geo.metric.euclid.einstein import RUN
        return "".join([RUN[i]*2 for i in range(m)])

    ## Permutations
    # \returns \f$Sym(n) \f$
    def sym(self):
        from geo.metric.schur_weyl.monte import Sym
        return Sym(int(self))

    ## Get m-th term in detrace operations
    def __getitem__(self, m):
        if m > int(self)//2:
            return None
        factor = self.func(int(self)-m)
        sign = pow(-1, m)
        if m == 0:
            return self.tensor, factor, sign
        deltas = self.delta(m)

        trace = self.itrace(m)
        d_free, t_free = self.ifree(m)

        n_left = len(d_free)
        n_right = len(t_free)

        result = 0

        attr0 = d_free + t_free

        done = []
        result = 0 * self.tensor

        print 'attr0=', attr0
        
        for pi in self.sym():
            attr = pi(attr0)
            d_samp = attr[:n_left]
            if d_samp in done:
                continue
            else:
                done.append(d_samp)
                
            t_samp = attr[n_left:]
            full = d_samp + trace + t_samp

            print 'D_{}'.format(d_samp), trace, t_samp
            
            result += getattr(deltas & self.tensor, full)
        
        return result, factor, sign

    def __iter__(self, result=1):
        for m in itertools.count(0):
            result = self[m]
            if result:  # filter
                yield reduce(operator.mul, result)
            else:
                break
            
    def __call__(self):
        return sum(list(self), 0*self.tensor) / self.func(int(self))
        


## metaclass computes indices from rank and assigns them to components,
# This is not for users -- most likely, this need to be scrapped for
# proper components implementation.... (TODO)- so a metaclass
# for derivative ranked objects and a sub that includes
# submission to the zoo- plus a local definition of AXES, so that
# object defined on spaces other than Euclidean can have
# natural access- e.g. rgb, ud, etc...
class ranked(abc.ABCMeta):
    """A metaclass-- used for classes with a rank static attribute that
    need components derived from it.

    See the rank2coordinates static method"""

    ## Keyed Default Dict Maps rank to class name for rank > 3.
    NAMER = keydefaultdict("Tensor{}".format,
                           {count + 5: name for count, name in
                            enumerate(
                                ('Five', 'Six', 'Seven', 'Eight', 'Nine'))})

    ## Extend type.__new__ so that it add components using rank2coordinates()
    # \sideeffect Updates ZOO.
    def __new__(mcs, *args, **kwargs):
        cls = super(ranked, mcs).__new__(mcs, *args, **kwargs)
        cls.components = mcs.rank2coordinates(cls.rank)
        ## The global ZOO of geometric animals
        if cls.rank in cls.ZOO: # guard
            ## This is a development path that need to be refactored away
            print "ZOO CONTAMINATION", cls
        else:
            # add cls to the ZOO
            cls.ZOO[cls.rank] = cls
            # and add scalar dilation to Scalar's "*" overload (this is
            # totally stovepipish).
            #pylint: disable=W0613
            cls.ZOO[0]._dispatch_mul[cls.rank] = _scalar_times_tensor
        return cls

    ## A function that computes a tensor's attributes from it's rank
    # Starting with "xyz" or "w".
    @staticmethod
    def rank2coordinates(rank_):
        """Convert non-negative integer rank into a tuple of
        tensor index attributes:
        0 --> ('w',)
        1 --> ('x', 'y', 'z')
        2 --> ('xx', 'xy', ..., 'zz')
        3 --> ('xxx', ..., 'zzz'
        e
        t
        c
        .
        .
        N --> ('x'*N, ...<3**N-2>, 'z'*N)
        .
        .
        """
        return tuple(
            map("".join, itertools.product(AXES, repeat=rank_))
        ) if rank_ else ('w',)  # guard on degenerate indices

    ## This method need to be refactored into __new__.
    # \param rank_ Rank
    # \returns class A new class.
    @classmethod
    def tensor_factory(mcs, rank_):
        """tensor_factory(rank_) makes a new tensor class."""
        from .three import HigherOrder
        dict_ = {'rank': rank_, '__metaclass__': mcs}
        return mcs(mcs.NAMER[rank_], (HigherOrder,), dict_)


## A temporary class until I decide If ops are composed or convolved
#  (brain freeze).
class LinearMap(object):
    """This class ensure linear map's compose method calls the convolve
    method"""
    __metaclass__ = abc.ABCMeta

    ## Transpose
    # \pure
    # \bug Is the linear map really where this should be?
    @abc.abstractproperty
    def T(self):
        """Transpose"""

    @abc.abstractmethod
    def alibi_transform(self, other):
        """The money method"""

    ## Chain together a serious of instances -- it's a static method
    # \param args like minded arguments so compose
    # \bugs reversed of arg might be a mistake, which works, possible b/c of
    # cancelation with calling order?
    # \returns type(arg)
    @staticmethod
    def chain(*args):
        """chain(*args)-->reduce(compose, args)"""
        compose = operator.mul
        return reduce(compose, reversed(args))

    ## Linear Map of Basis Vectors, as a Tensor
    # \return tensor of basis vectors, transformed.
    def frame(self):
        """'Tensor' of zipped ordered frame."""
        from .vector import BASIS
        from .tensor import zipcols ## Circular import
        return zipcols(*map(self, BASIS))

    ## This is the mother ship of external polymorphism.
    # \param name A type, and instance of a type, or a
    # \_\_name\_\_ of a type
    # \returns instance of type type.
    def aschart(self, name):
        """aschart(name) calls self.name to convert to type <name> as follows:
        self.aschart(instance) -->
        self.aschart(type(instance)) --> self.aschart(type)
        self.aschart(type) --> self.aschart(type.__name__.lower())
        --> self.aschart(name)
        self.(name) --> self.name

        SO: you can pass this method an Instance, a Type, or a name

        it it will find and call the method that converts self to and
        instance of that type (which by converntion, is called 'name' or
        the type).

        For exmaple:
        >>>tensor_ = versor_.astype('tensor') --> versor_.tensor
        """
        if isinstance(name, basestring):  # guard clause on type
            method_caller = operator.methodcaller(name)
            try:
                return method_caller(self)
            except AttributeError:
                msg = "No chart is associated with type {}".format(name)
                raise ChartError(msg)
        from types import TypeType
        if isinstance(name, TypeType):  # guard clause on type
            return self.aschart(name.__name__.lower())
        return self.aschart(type(name))

    ## alibi_transform() another Animal-- this is not implement nicely, b/c
    # linear map should follow the slotted tensor construction
    # \param other Some animal to transform.
    # \returns type(other) alibi transformed other.
    def __call__(self, other):
        """call *is* alibi_transform, and operators on a single
        argument."""
        return self.alibi_transform(other)

    ## Inverse's alibi_transform()
    ## alibi_transform() another Animal.
    # \param other Some animal to transform.
    # \returns type(other) alias transformed other.
    def alias_transform(self, other):
        """Alibi Transform."""
        return (~self)(other)

    ##  \f$ (f\circ g)(v) \rightarrow g(f(v))  \f$
    # \bug this also is revered (should be Alibi).
    def compose(self, other):
        """Composition is __mul__"""
        return self * other


## Base class for animals living in \f$ R^3 \f$, and python.
class Geometric(arraylike.ABCNumpy):
    """Abstract Base class Geometric objects.

    They are iterable over their components (but not pythonically),
    you have to use the iter() method.

    Some basic operations are also defined here:

    >>>-g  # negation
    >>>g.C  # complex conjugation
    >>>g.T  # generic transposition
    >>>g.H  # Hermitian conjugate

    There is also the "clean" method--which is for debugging.
    """
    __metaclass__ = abc.ABCMeta

    ## neg is the same class with all components negated, use map and
    #  operator.neg with iter().
    def __neg__(self):
        """-tensor maps components to neg and builds a new tensor"""
        return type(self)(*map(operator.neg, self.iter()))

    ## Repr is as repr does
    def __repr__(self):
        guts = "({})".format(", ".join(map(str, self.iter())))
        return repr(type(self)).lstrip("<class '").rstrip("'>") + guts

    ## just for clean up numeric precision woes when testing.
    # \kwd{1e-14} tol How clean?
    # \kwd digits how much to round
    # \returns type(self) rounded.
    def clean(self, tol=1e-14, digits=2):
        """round results for repping and debugging: that is,
        this should clean up numerical errors when comparing, say,
        DCMs that are computed from different sources."""
        #pylint: disable=E1101
        # Zero out real part's numerical error
        R = type(self)(*[(item * (abs(item) > tol))
                         for item in self.real.iter()])
        # Zero out imaginary... IBID....
        I = type(self)(*[(item * (abs(item) > tol))
                         for item in self.imag.iter()])
        # If the imaginary part exists, then include it.
        if abs(I) > tol:
            R += I
        # now round-- in a non traditional manner.
        scale = 10 ** digits
        return (scale*R).broadcast(round) / scale

    ## \f$ A^* \f$
    @property
    def C(self):
        """conjugate"""
        return self.conjugate()

    ## \f$ A^T \f$
    @property
    def T(self):
        """transpose"""
        #pylint: disable=E1101
        return self.transpose(*reversed(range(self.rank)))

    ## \f$ (A^*)^T \f$
    @property
    def H(self):
        """Hermitian conjugate-transpose"""
        #pylint: disable=E1101
        return self.T.C


## This decorator decorates element-wise operations on tensors
# \param op A binary operation
# \returns wrapped_op A method decorator
# \throws geo.utils.exceptions.NonCovariantOperation
# \throws TypeError
# \throws AttributeError
def _elementwise(op):
    """func = elementwise(op):

    op is a binary arithmetic operator.

    So is func, except it works on elements of the Tensor, returning a new
    tensor of the same type.
    """
    @functools.wraps(op)
    def wrapped_op(self, other):
        """wrapped op deals with possible array behavior."""
        try:
            # apply op to each pair of item, and build a new instance
            result = type(self)(*[
                op(*items) for items in itertools.izip_longest(
                    self.iter(), other.iter())])
        except TypeError as err:
            if None in items:  # raise Dx
                raise exceptions.NonCovariantOperation(
                    exceptions.error_message(op, self, other)
                    )
            print "This is a debug message: ", items
            raise err
        except AttributeError as err:
            # Forensic check on exception....
            x = (exceptions.NonCovariantOperation if
                 isinstance(other, arraylike.ABCNumpy) else
                 err)
            raise x(exceptions.error_message(op, self, other))
        return result
    return wrapped_op


## Base class for ranked tensors: \f$ {\bf T} \in V^{\otimes n} \f$ for
# \f$ V \in (R, C)\f$ and \f$ n \in (0,1,2,\ldots) \f$.
class Tensor_(Geometric, einstein.GymClass):
    """Base class For Any Rank Tensor:

    It *is* a Geometric object and it *has* a mixin for index
    gymnastics."""

    ## Invoke np.fromile to take construct.
    # \param args and kwargs are passed straight to np.fromfile
    # \kwd shape will override 1-D shape default
    # \return instance of cls.
    @classmethod
    def fromfile(cls, *args, **kwargs):
        """fromfile(*args, **kwargs) --> np.fromfile (with reshaping)

        If shape is in kwargs, than that overrides 1-D behavior."""
        from numpy import fromfile
        shape = kwargs.get(shape) or (-1, len(cls.components))
        return cls(*fromfile(*args, **kwargs).reshape(shape))

    ## Generalized Symmetrizer, projects onto \f$S^n(V)\f$m which
    # has dimension \f$ {n+3-1 \choose 3 - 1} \f$ (for 3=3).
    # \returns \f$ {\bf T}_{\{\}} = P_{\lambda}({\bf T}) \f$
    def S(self):
        """Symmetrized via Young Tableaux."""
        from ..schur_weyl.young import P_lambda
        return P_lambda(self.rank)(self)

    ## Generalized Antisymmetrizer:
    # \returns \f$ {\bf T}_{[]} = Q_{\lambda}({\bf T}) \f$
    def A(self):
        """Antisymmetrized via Young Tableaux."""
        from ..schur_weyl.young import Q_lambda
        return Q_lambda(self.rank)(self)

    ## Abstract Properties - is there a better way
    rank = 0
    components = None
    _dispatch_mul = {}

    ## The Zoo of geometric animals living in R3 is of extreme importance:
    # it is where they all live, keyed by their rank--and it is where
    # they are found by other animals-- both static and dynamic, but wait:
    # it is also a metaclass that makes new tensors at runtime.
    ZOO = keydefaultdict(ranked.tensor_factory)

    ## \f$ \sum_{\sigma(i,j,\ldots) \in S_n}{T_{\sigma}T^*_{\sigma}} \f$
    def __float__(self):
        """float(T) or map(float, T) computes real <T|T*>"""
        return sum(a*a.conjugate() for a in self.iter()).real

    ## \f$ \sum_{\sigma(i,j,\ldots) \in S_n}{T_{\sigma}T_{\sigma}} \f$
    def __complex__(self):
        """complex(T) or map(complex, T) computes real <T|T>"""
        return sum(a*a for a in self.iter())

    @staticmethod
    ## \param rank_
    # \yields tuples of indices (ordered)
    def iter_index(rank_):
        """iter_index(rank) yields indices for any rank."""
        return itertools.product(xrange(DIMENSIONS), repeat=rank_)

    ## \returns [0, 1, 2]
    @classmethod
    def _list_index(cls):
        """return [index for index, in cls.iter_index(1)]"""
        return reduce(operator.add, map(list, cls.iter_index(1)))

    ## \param cls
    # \yields Tensor indices.
    @classmethod
    def iter_indices(cls):
        """yield tensor indices."""
        return cls.iter_index(cls.rank)

    ## Antisymmetric Product:
    # \param self \f$ A_{efghi} \f$
    # \param other \f$ B_{klmn} \f$
    # \returns \f$  A_{efghi} \epsilon_{ijk} B_{klmn} \f$
    def wedge(self, other):
        """not recommended"""
        from .three import EPSILON
        # \todo refactor using einstein
        return self.inner(EPSILON.transpose(2, 1, 0)).inner(other)

    ## Wedge Product
    # \return wedge()
    def __xor__(self, other):
        return self.wedge(other)

    ## Outer product
    def __and__(self, other):
        if rank(other):
            return self.outer(other)
        # Trivial outer product: one or both are Scalar/number.
        return self * other

    ## Transpose any slots (currently reversed order...)
    # \param *args None or a permutation of indices from (0, ..., rank-1)
    # \returns Tensor with indices transposed.
    # \throws geo.utils.exceptions.IndicesError
    # \throws RuntimeError
    def transpose(self, *args):
        """transpose(*args):

        () --> (0,1,2,rank-1) to (rank-1, ..., 2, 1, 0)

        Otherwise:

        args are some permutations of (0, 1, ..., rank-1).

        The permutation sends the implicit index (by order) to
        the explicit index (by value), i.e:

        T.kijlmn == T.transpose(2, 0, 1, 3, 4, 5)
        """
        # set default
        if not args:  # kwd formatting
            args = list(reversed(xrange(self.rank)))
        # Ensure indices make sense
        if len(args) != self.rank:
            raise exceptions.IndicesError(
                "Rank {} Tensor got {} indices.".format(self.rank,
                                                        len(args)))
        # Convert negative slots
        args = [(self.rank+item if item < 0 else item) for item in args]
        # check for duplicate
        if len(set(args)) != self.rank:  # raise
            bag = collections.Counter(
                [item for item in args if args.count(item) != 1])
            raise exceptions.IndicesError(
                "Duplicate Index: {}".format(bag.keys()))
        # Now do it:
        try:
            index_combos = [[item[index] for index in args] for item in
                            map(list, self.iter_indices())]
        except IndexError:
            for n in args:
                if 0 <= n < self.rank-1:  # raise
                    continue
                raise exceptions.IndicesError(
                    "No index slot = {}".format(n))
            else:
                msg = "Undiagnosed Error for rank={}, indices={}"
                raise RuntimeError(msg.format(self.rank, args))
        return type(self)(*list(itertools.starmap(self.e, index_combos)))

    ## numpy-style transpose (inverts permutation)
    # \param args
    def ntranspose(self, args):
        """Numpy style transpose: takes an iterable,
        and inverts definition, e.g:

        T.jkilmn == T.ntranspose([2, 0, 1, 3, 4, 5])"""
        from ..schur_weyl.monte import Perm
        return self.transpose(*(~Perm(*args)))

    ## All Transpositions
    # \yields \f$ T_{\pi} \forall \pi \in S_n \f$
    def itranspositions(self):
        """Yield all transpositions of indices."""
        from ..schur_weyl.monte import Sym
        return itertools.imap(self.ntranspose,
                              itertools.imap(list, Sym(self.rank)))

    ## Tensor Symmetry of a Weyl module (standard tableaux)
    # \param tableaux A standard tableau \f$
    # {\bf \lambda} =(\lambda_1,\ldots,\lambda_n)\f$
    # \yields \f$ {\bf T_{\lambda}} \forall\f$ Young tableaux in the
    # Young diagram's Weyl module.
    def weyl_module(self, tableau):
        """for T_[n1, ..., nk] in tensor.weyl_module(n1, ..., nk):
        .....<yields all tensors for diagram>...
        where the n_i are the (ordered) length of rows in
        the young.Diagram for the symmetrizing partition of the rank."""
        return tableau(self)
    
    ## Symmetrize according to a standard Young tableaux, takes
    # 0.00064, 0.0012, 0.0043, 0.022, 0.25, 4.37, 96, 2532, 75441 seconds.
    # \param rows a list of rows, which are lists of cells.
    # \returns tensor
    def symmetrize(self, rows):
        """T_{n1, nN} = tensor.symmetrize(n1, ..., nN}

        where n1 + ... + nN is an integer partition of the tensor rank."""
        from ..schur_weyl.young import StandardTableau
        return StandardTableau(*rows)(self)

    ## Integer Partitions of the Tensor's Rank.
    # \yields \f$ \lambda^{(n_i)} \f$
    # with \f$ \sum_i{n_i} = {\mathrm{rank}}(T)\f$
    def partitions(self):
        """Yield partitions of the tensor's rank, e.g.:
        [3],
        [2, 1]
        [1, 1, 1,]

        for rank 3."""
        from ..schur_weyl.pascal import P
        return P(self.rank).partitions()

    ## Diagrams for each partition in partitions().
    # \yields young.Diagram() for each partition in pascal.P.partitions().
    def diagrams(self):
        """Yield young.Diagram for each partition in self.partitions(),
        e.g., for rank 3:

        [ ][ ][ ]

        [ ][ ]
        [ ]

        [ ]
        [ ]
        [ ]"""
        from ..schur_weyl.young import Diagram
        return itertools.starmap(Diagram, self.partitions())

    ## All Standard Tableaux in the Weyl Modules
    # \kwd _func Private method caller (young.Diagram_.standard_tableaux()).
    def standard_tableaux(self,
                          _func=operator.methodcaller('standard_tableaux')):
        """Yield all standard tableaux for tensor's rank, e.g., for rank 3:

        [i][j][k]

        [i][j]
        [k]

        [i][k]
        [j]

        [i]
        [j]
        [k]"""
        return itertools.chain(*itertools.imap(_func, self.diagrams()))
    
    ## Generate all possible weyl_module() results-- this is the pinacle
    # of tensor mathematics.
    # \yields tuple of (geo.metric.schur_weyl.young.Tableau,
    # symmetrized tensor)
    def weyl_modules(self):
        """Generates tuples of:

        (young.Tableau, Tensor symmetrized wrt to the tableau).
        """
        # get 2 copies of the standard tableaux iterator
        tab1, tab2 = itertools.tee(self.standard_tableaux(), 2)
        # yield the tableaux paired with its weyl module symetry
        return itertools.izip(tab1, itertools.imap(self.weyl_module, tab2))
                
    ## Young symmetrized tensors (without the tableuax).
    # \yields \f$ T_{\lambda^{\alpha}} \f$
    def symmetries(self):
        """generate all tensors in the weyl_modules()"""
        return itertools.imap(self.weyl_module, self.standard_tableaux())

    ## This is just a helper for those dealing with high rank tensors,
    # \kwd tol=1e-7
    # \returns None
    # \sideeffect Print data to stdout.
    # \throws excpetions.PolymorphismError for array_like tensors.
    def show_modules(self, tol=1e-7):
        """show_modules([stream=sys.stdout, tol=1.e-7])
        prints numerically significant weyl modules' Young tableaux and
        L2-norm (abs)."""
        try:
            for tab, ten in self.weyl_modules():
                if abs(ten).w > tol:  # filter
                    print tab, '\n', abs(ten), '\n'*2
        except ValueError:
            raise exceptions.PolymorphismError("method is for singletons.")

    ## Decompose into scalars and tuples of basis vectors
    # \yields \f$(T_{ij\ldots k},\ ({\bf \hat{e}_i, \hat{e}_j, \cdots,
    # \hat{e}_k})) \forall (i, j, \cdots, k) \f$
    def idecompose(self):
        """Iterative decomposition of rank N T into pairs of
        a scalar and decomposed polyadics (i.e, N unit vectors):

        T_00...0, (X, X, ..., X)
        T_00...1, (X, X, ..., Y)
        .          .  .       .
        .          .  .       .
        .          .  .       .
        T_22...2, (Z, Z, ..., Z)"""
        #pylint: disable=E1102
        from .vector import BASIS
        for units in itertools.imap(list,
                                    itertools.product(BASIS,
                                                      repeat=self.rank)):
            # note order: Tensor.__call__ conforms to "slotted" notation.
            yield self(*units), units


    ## Polyadic decomposition of rank N tensors:
    # \yields \f$T_{ij\ldots k}{\bf \hat{e}_i\hat{e}_j\cdots \hat{e}_k} \f$
    def polyadics(self):
        """Generator of Cartesian components.

        For T[i, j, ...k], yield:

        T_ij...k * (e_i & e_j & ... & e_k)

        for all i, j, ..., k (in order)."""
        for t_ijk, etuple in self.idecompose():
            yield reduce(operator.and_, etuple) * t_ijk.w

    ## Unit Polyadic Basis Iterator:
    # \yields \f$ {\bf \hat{e}_i\hat{e}_j\cdots \hat{e}_k} \f$
    @classmethod
    def ibasis(cls):
        u"""Generator of basis polyadics.

        For Rank N, yields:

        X\u2026(N times)\u2026X
        X\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026XY
        X\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026XZ
        X\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026YX
        \u22EE           \u22EE
        Z\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026\u2026ZY
        Z\u2026(N times)\u2026Z"""
        return cls(*itertools.repeat(1, pow(DIMENSIONS, cls.rank))).polyadics()

    ## Instantiate a basis vector
    # \param cls (implicit) A rank-N tensor class
    # \param args N axis labels: \f$ (i_0, j_1, \cdots, m_n )\f$
    # \returns \f$ {\bf \hat{e}}_{i_0j_1\cdots m_n} \f$
    # \throws exceptions.IndicesError Wrong number of indices
    # \throws execeptions.TensorIndexValueError Invalid axis
    @classmethod
    def baseunit(cls, *args):
        """e_ij..k = T.baseunit(i, j, ..., k)

        is the unit polyadic for indices i, j, ..., k."""
        if len(args) != cls.rank:  # raise
            msg = "Expect {} indices, got: {}"
            raise exceptions.IndicesError(msg.format(cls.rank, len(args)))
        values = [0.] * pow(DIMENSIONS, cls.rank)
        index = 0
        for count, arg in enumerate(reversed(args)):
            if not 0 <= arg < DIMENSIONS:  # raise
                raise exceptions.TensorIndexValueError(arg)
            index += pow(DIMENSIONS, count) * arg
        values[index] += 1.
        return cls(*values)

    ## Cartesian Basis elements
    # \returns tuple of \f$ e_i \otimes e_j \otimes \cdots e_m \forall
    # (i, j, \ldots, m) \f$
    @classmethod
    def basis(cls):
        """T.basis() --> (xx...x, xx...y, ..., zzz..z)"""
        return tuple(itertools.starmap(cls.baseunit, cls.iter_indices()))

    ## The Null Element
    # \returns \f$ {\bf 0} \f$
    @classmethod
    def null(cls):
        """The Null Element"""
        return cls(*([0.]*len(cls.components)))

    ## Generalize Proj/Ref/Ref-tion from any rank tensor to a vector
    # \param func operator to use: ">>", "<<", or "|".
    # \param index slot on which to work: (0, .., RANK-1)
    # \param other A vector that will be in the polyadic (at index).
    # \returns Tensor like self, but func'd on index and other.
    # @image html cauchy.png
    # @image latex cauchy.png
    def _jector(self, func, other, index):
        """This is complicated:

        func    result
        rshift projection
        lshift rejection
        or_    reflection

        other: a vector

        index: integer in (0, ..., RANK-1)

        For instance:
        ((T >> X)(0) >> T)(1) --> T.xy only non-zero element

        ((T << X)(0) << T)(1) --> zeros T.x* and T.*y.

        It implements a generalized tensor-->tensor' Proj, Rej, Ref ection
        for any rank combo (If you can even figure the geometric significance
        of that-- it transforms, so it must mean something.)
        """
        # initialize result as a null rank(self) tensor
        m = self * 0
        # loop over idecompose scalar and polyadic components.
        for t_ijk, basis_vecs in self.idecompose():
            for count, e in enumerate(basis_vecs):
                # convert 'e' (basis vector) to a Xjection onto other.
                if count == index % self.rank:  # filter
                    e = func(e, other)
                # initialize or continue tensor product construction
                #pylint: disable=E0601
                ee = ee & e if count else e  # init
                continue
            # add basic_vec's image to result
            m += t_ijk * ee
        return m

    # \param index \f$ i_k \f$ with \f$ 0\le k < rank \f$
    # \returns function \f$ f(\vec{v}) \f$ projecting slot k onto v.
    def __rshift__(self, other):
        """Tensor projection:
        T" = T >> T'  projects T onto T' to make a new tensor T".
        f  = T >> v   returns a function:

        f(0) projects first index onto v
        f(-1) projects last index onto v."""
        if rank(other) == 1:  # process (todo decorate)
            return functools.partial(self._jector, operator.rshift, other)
        return self.projection(other)

    def __lshift__(self, other):
        """Tensor rejection:
        T" = T >> T'  rejects T onto T' to make a new tensor T".
        f  = T >> v   returns a function:

        f(0) rejects first index onto v
        f(-1) rejects last index onto v."""
        if rank(other) == 1:  # process (todo decorate)
            return functools.partial(self._jector, operator.lshift, other)
        return self.rejection(other)

    def __or__(self, other):
        """Tensor reflection:
        T" = T >> T'  reflects T onto T' to make a new tensor T".
        f  = T >> v   returns a function:

        f(0) reflects first index onto v
        f(-1) reflects last index onto v."""
        if rank(other) == 1:  # process (todo decorate)
            return functools.partial(self._jector, operator.or_, other)
        return self.reflection(other)

    ## Rank independent Projection method via direct injection.
    projection = cauchy.projection

    ## Rank independent Rejection method.
    rejection = cauchy.rejection

    ## Rank independent Reflection method.
    reflection = cauchy.reflection

    ## Get basis componets based TENSOR indices (NOT array indices).
    # \param indices or Ellipsis, e.g: \f$ (0, i, 1, j) \f$
    # \returns \f$ T_{0ij1} \f$
    # \throw geo.utils.exceptions.TensorIndexValueError
    def e(self, *args):
        """e(*args) runs on indices."""
        # Do we need deeper index processing?
        if len(args) != self.rank or any(map(self._isrunning, args)):  # kwd
            result = self._eparse(*args)
        else:
            index = sum((arg * DIMENSIONS**count) for count, arg in
                        enumerate(reversed(args)))
            try:
                result = getattr(self, self.components[index])
            except IndexError:
                raise exceptions.TensorIndexValueError(
                    "Invalid dimension in {}".format(args))
        return result

    ## helper decides whether or not an index request is running or fixed.
    _isrunning = staticmethod((Ellipsis, None).__contains__)

    ## Parse index-notation with running indices.
    def _eparse(self, *args):
        # 'zero' pad, truncate junk
        args = (list(args) + [Ellipsis] * (self.rank - len(args)))[:self.rank]
        # process request type
        index_type = map(self._isrunning, args)
        return self.ZOO[sum(index_type)](
            *itertools.starmap(
                self.e,
                itertools.product(*[
                    (self._list_index() if i_type else [i_value])  # kwd
                    for i_type, i_value in zip(index_type, args)])))

    ## Contract on inner indices:
    # \param self \f$ T_{i \ldots k} \f$
    # \param other \f$ T'_{i \ldots k} \f$
    # \kwd count=1 Number of indices to contract
    # \returns \f$ T_{i \ldots k} T'_{k \ldots m} \f$
    # \throws exceptions.NonGeometricTypeError ill-fated argument
    # \throws exceptions.IndicesError No index on geometric object
    # \throws TypeError
    def inner(self, other, count=1):
        """Generalized inner product, contracts on 'count' pairs of
        of indices:

        result = self.inner(other, count=1)

        result.rank = self.rank + other.rank - 2 * count
        """
        ## Guard on non-geometric arguments.
        try:
            orank = other.rank
        except AttributeError:
            raise exceptions.NonGeometricTypeError(
                "argument is not an appropriate geometric object.")

        ## Running over Scalar 'w' make  garbage, so don't do it.
        if not (self.rank and orank):  # raise
            raise exceptions.IndicesError(
                "No indicies, no contraction!\n{}".format(other))

        ## Right stub is in the call signature
        left_stub = [None] * (self.rank-count)

        ## Need to do this to avoid crash after pulling floats out of self.
        func = operator.mul if self.rank == count else operator.and_

        # If rank is count, then each arg to func is a pure component (
        # number or ndarray), and must use their "*"- otherwise,
        # at least 1 is a tensor, so Tensor.__and__ must be invoked)
        # The, then summation is over indices, and you must use reduce(add,
        # b/c sum doesn't treat ndarrays right-- though it does work for
        # singletons--a messy method, but it is all for the dear user:
        # who just write code as equations.
        try:
            result = reduce(
                operator.add,
                [func(self.e(*(left_stub + list(indexes))), other.e(*indexes))
                for indexes in self.iter_index(count)])
        except AttributeError:
            raise TypeError("inner expected a tensor, not: {}".format(
                            type(other)))
            
        # If contraction was complete (V*V, or T*T count=2, etc),  then
        # result is now a number or ndarray-- so make it a Scalar
        if rank(result) is None:  # review this- should it be a scalar?
            result = self.ZOO[0](result)
        return result

    ## Outer product for all orders
    # \param self \f$ T_{i \ldots k} \f$
    # \param other \f$ T'_{i' \ldots k'} \f$
    # \returns \f$ T_{i \ldots k} T'_{i' \ldots k'} \f$
    def outer(self, other):
        """generalized outer product:

        result.rank = self.rank + other.rank"""
        # Get the class (i.e. rank) of the result (note: ZOO will make it at
        # runtime, If needed- ZOO is a dictionary and a class factory)
        cls = self.ZOO[self.rank + other.rank]

        # mix and match indices:
        return cls(*[(self.e(*item[:self.rank])*other.e(*item[-other.rank:]))
                     for item in cls.iter_indices()])

    ## Contract a pair of indices
    # \param i1 An index, counted from left-to-right starting a 0
    # \param i2 _ditto_
    # \returns Tensor contraction in \f$i_1\f$ and \f$i_2\f$.
    # \throws exceptions.TensorIndexValueError Axis to large or repeated
    def contract(self, i1, i2):
        """contract(i1, i2)

        T_0123..(N-1)

        will contract on indices (i1,i2), or If i < 0, i --> N-i
        (which is not a tensor algebra thing, but it is a python thing)"""
        if i1 < 0:  # kwd (negative indices)
            i1 += self.rank
        if i2 < 0:  # kwd (negative indices)
            i2 += self.rank
        # sort puts left index on the left, and right on the right.
        i1, i2 = sorted([i1, i2])
        if i2 >= self.rank or i1 == i2:  # raise
            # chose message
            raise exceptions.TensorIndexValueError(
                i2,
                comment="it is " +
                {True: 'repeated', False: 'too big'}[i1 == i2])
        # Guard ends: method starts.
        items = []
        # behold, the itertools magic
        # (1) loop over the results indices
        for fixed in itertools.imap(list, self.iter_index(self.rank-2)):
            item = 0.  # initialize sum
            # now do some Einstein summation over x, y, z:
            for i in self._list_index():
                index = fixed[:]  # copy final indices to work on
                index.insert(i1, i) # now insert axis values into index
                index.insert(i2, i)
                item += self.e(*index)  # go get the component, increment sum
            items.append(item)  # now add this component to results list.
        # unpack results list into a tensor with the correct rank
        return self.ZOO[self.rank-2](*items)

    ## Unpack a numpy array into constructor
    # \param cls (implicit)
    # \param array_ A numpy array with len correct
    # \returns cls instance with array_'s items..
    # \throws TypeError If len(array_) doesn't match cls.__init__
    @classmethod
    def fromarray(cls, array_):
        """From a numpy array"""
        return cls(*array_) if cls.rank else cls(array_)  # guard: degenerate

    ## Create a real (or complex) tensor with random components.
    @classmethod
    def random(cls, n, imag=False):
        """Create an instance with random normal components:
        imag=True: adds a complex component."""
        from numpy import random
        result = cls(*random.randn(len(cls.components), n or None))
        if imag:  # kwd
            result += (1j)*cls.random(n=n, imag=(not imag))
        return result

    ## Multimethod-style implementation.
    # \param self (implicit) a tensor
    # \param other (explicit) another tensor
    # \returns tensor Adjacent index contraction
    def __mul__(self, other):
        """Note to user: __mul__ is inherited for Tensor. self's _dispatch_mul
        dictionary is keyed by other's rank in order to get the correct
        function to multiply the 2 objects.

        This default "*" does adjacent index contraction via inner.
        """
        return self._dispatch_mul.get(rank(other),
                                      type(self).inner)(self, other)

    ## reflected mul is always a non-Tensor, ans is left_dilation().
    # \param self \f$ {\bf T} \f$
    # \param other \f$ \alpha \f$ A scalar or matrix 
    # \returns \f$ {\bf T} = \alpha{\bf T}\f$ Components are left multiplied.
    def __rmul__(self, other):
        """Left Dilation:

        T' = alpha * T.
        """
        return self.left_dilation(other)

    ## elementwise() decorated addition
    @_elementwise
    def __add__(self, other):
        """t3_i...k = t1_i...k + t2_i...k  (elementwise add only)"""
        return operator.add(self, other)

    ## elementwise() decorated subtraction
    @_elementwise
    def __sub__(self, other):
        """t3_i...k = t1_i...k - t2_i...k  (elementwise sub only)"""
        return operator.sub(self, other)

    ## Division should be pretty straight forward-- but its not always for
    # arrays-- that is, supporting singleton and array components can get
    # tricky.
    def __div__(self, other):
        if rank(other) == 0:  # guard on scalar division
            return self.__div__(other.w)
        try:
            elements = [item/other for item in self.iter()]
        except TypeError as err:
            try:
                # does other know why it can't be a divisor?
                raise other.DivisionError
            except AttributeError:
                pass
            raise err
        return type(self)(*elements)

    ## What?
    def __str__(self):
        return ", ".join("{}={}".format(attr, getattr(self, attr))
                         for attr in self.components)

    ## The L2 norm is a Scalar, from normsq()**0.5 -d e p r e c a t e d
    def L2norm(self):
        """normsq()**0.5"""
        return self.norm(p=2)

    ## \f$L^p\f$-norm.
    # \kwd p=2 Exponent in L-norm
    # \returns \f$ L_p(T) = \big[
    # \sum_{i=\ldots} |T_i|^p\big]^{\frac{1}{p}} \f$
    def norm(self, p=2):
        """Norm-- should work for R3 and C3"""
        return pnorm(p, self)

    ## General "unit" tesnor.
    hat = cauchy.hat

    ## Abs value is always the L2-norm, inspite of "det".
    __abs__ = L2norm

    ## Left dilation of components via broadcast().
    # \param self \f$ {\bf T} \f$
    # \param other \f$ \alpha \f$ A scalar or matrix 
    # \returns \f$ {\bf T}=\alpha{\bf T}\f$ Components are left multiplied.
    def left_dilation(self, other):
        """left_dilation does left multiplication of components."""
        return self.broadcast(functools.partial(operator.mul, other))

    ## Right dilation of components
    # \returns broadcast with lambda x: x * other
    def right_dilation(self, other):
        """right_dilation does right multiplication of components."""
        try:
            result = self.broadcast(lambda x: x * other)
        except TypeError:
            # unknown external polymorphism (e.g. polar * tensor).
            print "EXTERNAL POLY 1020"
            result = other.__rmul__(self)
        return result

    ## Default right dilation is called by "__mul__"
    dilation = right_dilation

    ## Return True/False for singleton-like tensors, and bool arrays for
    # array_like input-- this is a troublesome overload.
    def __nonzero__(self):
        """
        ## Note: isSequenceType doesn't work here
        try:
            for dummy in self:
                break
        except TypeError:
            # deal with singelton tensor/coordinates
            result = any(map(bool, self.iter()))
        else:
            # Now deal with numpy array)like attributes
            from numpy import array
            result = array(map(bool, self))
        return result"""
        return any(map(bool, self.iter()))

    ## Generalized transformer
    # \param R A Rotation (or other transform) \f$ R_{ij} \f$
    # \param T A tensor of any rank \f$ T_{ij\ldots mn} \f$
    # \returns T' \f$ T'_{i'j' \ldots m'n'} = R_{i'i}R_{j'j}\cdots R_{m'm}
    # R_{n'n} T_{ij \ldots mn} \f$
    def alibi_transform_any(self, other):
        """R(T_ijklm) = T_i'j'k'l'm'
        R_m'm *
        R_l'l *
        R_k'k *
        R_j'j *
        R_i'i *
        T_ijklm
        (R * (R * T)._left())._left()"""
        #pylint: disable=W0613
        T = other
        for dummy in xrange(other.rank):
            T = (self * T)._left()
        return T

    ## Angle between another vectors, or whatever.
    # \returns
    # \f$ \cos^{-1}{\hat{a}\cdot\hat{b}} \f$ as a Scalar instance \n this
    # farms out the trig call so the developer doesn't have to worry about it.
    # \bug Is round-off error weak.
    def angle(self, other):
        """a.theta(b) --> Scalar(a*b / |a||b|)
        -- for real, it's a rank-0 obj"""
        from ...utils.trig import arccos
        return arccos(cauchy.scalar_projection(self.hat(), other).w)

    ## Convert to a numpy array
    # \returns numpy.ndarray
    def asarray(self):
        """Convert to a numpy.ndarray ([len(self),] len(self.components))"""
        from numpy import array
        return array(self.tolist()).T

    ## Convert to array and write to a file
    # \sideeffect Writes a file.
    def tofile(self, *args, **kwargs):
        """tofile(*args, **kwargs) goes straight to numpy.tofile"""
        self.asarray().tofile(*args, **kwargs)


    ## Get Cartesian form of eigenstate
    # \param j Weight
    # \param q Seniority Index
    # \parma m azimuth order.
    # \returns Cartesian Tensor part this is |j, m>_q eigenstate.
    def irrep(self, j, q, m):
        """tensor.irrep(j, q, m) projects the weight j, senority index q,
        and azimuth order m eigenstate out, in a Cartesian basis."""
        stensor = self.spherical()
        value = stensor[j, q, m]
        return stensor.fromlist(
            [value*(j == j0 and q == q0 and m == m0)
             for j0, q0, m0 in stensor.ijqm()]).tocart()
    
    def hyperdet(self):
        """sum([reduce(mul, [t.e(*item) for item in [[p[n] for p in perms] for n in range(3)]]) * reduce(mul, [p.parity() for p in perms]) for perms in combinations(S(3),2)])/6
        """
        from ..schur_weyl.monte import Sym
        sn = list(Sym(3))
        k = self.rank
        result = 0
        for perms in itertools.combinations(sn, k):
            sign = reduce(operator.mul, [p.sign() for p in perms])

            term = 1.

            for i in range(DIMENSIONS):
                indices = [p[i] for p in perms]
                term *= self.e(*indices)

            result += sign * term
        return result / 6.
    
    def natural_form2(self):
        print "This is in developement"
        return DeTracer(self)
        
    
                
## The ZOO of geometric animal: it's a dictionary and a metaclass factory
ZOO = Tensor_.ZOO

## A view into the ZOO (no touchy!).
VIEW = Tensor_.ZOO.viewitems()
