"""Generalized base class for non-Cartesian vectors."""
## \namespace geo.metric.euclid.kant._bases ABC for Non-Cartesian Coordinates
import abc
import itertools
import operator
import functools

from ..euclid import rank
#from .. import vector


## Remove Before Flight
def debug(method):
    """There is a decorator that will fix all overloads. I can't find it."""

    @functools.wraps(method)
    def dmethod(self, *args, **kwargs):
        """TBD decorator."""
        return method(self, *args, **kwargs)

    return dmethod


## These geo.metric.euclid.vector.Vector methods are added to non-Cartesian
# classes at import-time.
TRAITS = ('dual', 'versor', 'hat', 'right_versor', 'quaternion',
          'projector', 'reflector', 'rejector', 'householder',
          'polar', 'cylinder', 'parabola', 'spherical', 'var')


## Add some instancemethod that Cartesian Vector do.
def kant(cls):
    """kant(cls) decorator adds TRAITS as bound methods operating on
    'vector' property."""
    from functools import partial
    from new import instancemethod
    for func in TRAITS:
        setattr(cls, func, instancemethod(partial(_closure, func),
                                          None,
                                          cls))
        continue
    return cls


## Closure calls method in self.vector, and converts vector results to
# correct Non Cartesian class.
def _closure(func, self):
    # Call self.vector.func()
    result = operator.methodcaller(func)(self.vector)
    # If it's a vector, then convert to type(self)
    if isinstance(result, vector.Vector):  # type check
        result = type(self).fromvector(result)
    return result


## Non Cartesian Object
@kant
class NonCartesianBase(object):
    """This is base class, duh. It has a loft duty:

    For any operation on any non-Cartesian vector (self) with any other
    vector of any type:

    do that op polymorphicaly and return the result in self's coordinates.

    Ideally: a decorator will do all that, and the overloaded methods
    will just return from the operator module. We're here yet.
    """
    __metaclass__ = abc.ABCMeta

    rank = 1

    ## Construct from 3 parameters.
    # \param p0 1st coordinate
    # \param p1 2nd coordinate
    # \param p2 3rd coordinate
    def __init__(self, p0, p1, p2):
        for attr, value in zip(self.components, (p0, p1, p2)):
            setattr(self, attr, value)

    ## \yields components
    def iter(self):
        """Iterate components."""
        return itertools.imap(self.__getattribute__, self.components)

    ## \returns list of components
    def tolist(self):
        """List of components"""
        return list(self.iter())

    ## \returns vector.Vector
    @property
    def vector(self):
        """As a cartesian vector."""
        from ..vector import Vector
        return Vector(*self._tocart(*self.iter()))

    @abc.abstractmethod
    def _tocart(self):
        """This is the conversion to a Cartesian triplet"""

    @classmethod
    def fromvector(cls, vector_):
        """Construct from Cartesian vector."""
        #pylint: disable=E1101
        data = list(vector_.iter())
        data2 = cls._fromcart(*data)
        return cls(*data2)

    @classmethod
    def fromeuclid(cls, self):
        return locals()
    
    
    @classmethod
    def basis(cls):
        """Cartesian Basis in local basis (at origin)."""
        from ..vector import BASIS
        return map(cls.fromvector, BASIS)

    ## \return \f$ +{\bf \vec{v}} \f$
    @debug
    def __pos__(self):
        return self

    ## \returns \f$ -{\bf \vec{v}} \f$
    @debug
    def __neg__(self):
        return type(self).fromvector(self.vector.__neg__())

    ## \param other a vector
    # \returns \f$ \vec{a} + \vec{b} \f$ in a's coordinates.
    @debug
    def __add__(self, other):
        return type(self).fromvector(operator.add(self.vector, other.vector))

    @debug
    def __radd__(self, other):
        return type(other).fromvector(operator.add(other.vector, self.vector))

    ## \param other a vector
    # \returns \f$ \vec{a} - \vec{b} \f$ in a's coordinates.
    @debug
    def __sub__(self, other):
        return type(self).fromvector(operator.sub(self.vector, other.vector))

    @debug
    def __rsub__(self, other):
        return type(other).fromvector(operator.sub(other.vector, self.vector))

    ## External Polymorphism (-1 point).
    @debug
    def __mul__(self, other):
        rank_ = rank(other)
        if not rank_:  # external polymorphism
            result = self.dilation(other)
        elif rank_ == 1:
            result = self.vector * other.vector
        else:
            result = self.vector * other
            if rank(result) == 1:  # external polymorphism
                result = type(self).fromvector(result)
        return result

    @debug
    def __rmul__(self, other):
        return type(self).fromvector(operator.mul(self.vector,
                                                  other))

    ## Tricky business here.
    # \param other who knows what this will be
    # \param _hash choose __mul__ or __rmul__ based on other.
    # \returns __(r)mul__(self, other) IF not rank(other).
    def __imul__(self, other, _hash=(__rmul__, __mul__)):
        return _hash[int(bool(rank(other)))](self, other)

    @debug
    def __div__(self, other):
        return type(self).fromvector(operator.div(self.vector,
                                                  other))

    @debug
    def __xor__(self, other):
        return type(self).fromvector(operator.xor(self.vector,
                                                  other.vector))

    @debug
    def __rshift__(self, other):
        return type(self).fromvector(operator.rshift(self.vector,
                                                     other.vector))

    @debug
    def __rrshift__(self, other):
        return type(other).fromvector(operator.rshift(other.vector,
                                                      self.vector))

    @debug
    def __lshift__(self, other):
        return type(self).fromvector(operator.lshift(self.vector,
                                                     other.vector))

    @debug
    def __rlshift__(self, other):
        return type(other).fromvector(operator.lshift(other.vector,
                                                      self.vector))

    @debug
    def __or__(self, other):
        return type(self).fromvector(operator.or_(self.vector, other.vector))

    @debug
    def __ror__(self, other):
        return type(other).fromvector(operator.or_(other.vector, self.vector))

    @debug
    def __and__(self, other):
        return operator.and_(self.vector, other.vector)

    @debug
    def __rand__(self, other):
        return operator.and_(self.vector, other.vector)

    @debug
    def __pow__(self, other):
        return operator.pow(self.vector, other)

    @debug
    def __getitem__(self, index):
        return self.fromvector(self.vector[index])

    @debug
    def __setitem__(self, index, value):
        raise TypeError

    @debug
    def __len__(self):
        return len(self.vector)

    @debug
    def __abs__(self):
        return abs(self.vector)

    @debug
    def __float__(self):
        return float(self.vector)

    @debug
    def __complex__(self):
        return complex(self.vector)

    def __str__(self):
        return ", ".join("{}={}".format(*item) for item in
                         zip(self.components, map(self.__getattribute__,
                                                  self.components)))
    def __repr__(self):
        return "{}({}, {}, {})".format(type(self).__name__, *self.iter())

    ## Highly Experimental
    # \returns An attribute or a closure around a method (disabled)?
    def __getattr__(self, attr):
#        print "DD K0", attr
        result = getattr(self.vector, attr)
        if callable(result):  # deeply confusion external polymorphism
 #           print "DD K1"
            if isinstance(result, type(self.vector)):  # type check
                print "DD K20"
                result = self.fromvector(result)
            elif isinstance(result, type(self)):
                print "DD K30"
            else:
                print "DD K40"
                return result
        return result

    ## Basis vector component
    # \param i
    # \returns \f$ v_i \f$
    def e(self, i):
        """Get basis element"""
        return getattr(self, self.components[i])

    ## Basis vector.Vector
    # \param i
    # \returns \f$ {\bf \hat{e}_i} \f$
    def ehat(self, i):
        """Unit vector along e-th basis?"""
        return self._hats[i](self)

    def local_basis(self):
        """Local basis Tensor (Matrix, really)"""
        from .... import ziprows
        return ziprows(*map(self.ehat, (0, 1, 2)))
