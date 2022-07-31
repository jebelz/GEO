"""This module defines interfaces so that clients can use numpy
without their clients knowning about implementation.

It has high cyclomatic complexity in some methods-- numpy is likely to grab
array structure when you don't want it to-- hence you must HIDE the fact
that your tensors look like arrays. Moreover, operator overloading is not
trivial:

You cannot let numpy run the show!

In spite of the long methods and hiddeous If-Then depth, there are some
winners (of the beauty contest-- it's all congenial):

enlister
dishcharger
matrix_wrap

are (nested) decorators to handle methods that cannot work on arraylike
data and singleton data with the same code (well, now they can).
"""
#pylint: disable=E1101

## \namespace geo.utils.arraylike Array_like interfaces
import abc
import functools
import itertools
import operator

from . import trig
## For debugging the getitem, getattr, getatribute complexities.
VERB = 0


## Got length?
def got_length(arg):
    """arg is a container-- old school check."""
    try:
        len(arg)
    except TypeError:
        return False
    return True


## operator.itemgetter, with Error protection.
def itemgetter(index):
    """itergetter w/ error protection."""

    def func(item):
        """get item or None"""
        try:
            item = item[index]
        except TypeError:
            pass
        return item
    return func


## For method that cannot act "array_like", array_like results are put
# into a list
def enlister(method):
    """Enlisted method is run on items in an iterable (Geometric),
    and the results are returned in a list."""

    method.__doc__ = method.__doc__ + " Or as list thereof."

    @functools.wraps(method)
    def wrapped_method(self, *args, **kwargs):
        """array_like's method call will be mapped to a list"""
        if self._got_length():  # external poly (array/singleton)
            result = [method(item, *args, **kwargs) for item in self]
        else:
            result = method(self, *args, **kwargs)
        return result
    return wrapped_method


## For enlister wrapped method who's results can be caste array_like:
# this decorator does it.
def discharger(method):
    """Method results--If enlisted-- are repackaged into correct
    Geometric Animal."""

    @functools.wraps(method)
    def wrapped_method(self, *args, **kwargs):
        """Result be de-enlisted, If needed."""
        # just call the method
        result = method(self, *args, **kwargs)
        # If it returned a list, then the results must be reconstituted
        if isinstance(result, list):  # external poly
            # chose type for where else? The list.
            if trig.NUMPY:
                from numpy import array
            else:
                array = list
            return type(result[0])(
                *itertools.imap(array,
                                itertools.izip(*[item.iter()
                                                 for item in result]))
            )
        return result
    return wrapped_method


## Decorator to wrap result in a method.
def matrix_wrap(method):
    """Wrap discharged method in a matrix."""

    @functools.wraps(method)
    @discharger
    @enlister
    def wrapped_method(self, *args, **kwargs):
        """discharger(enlister(method(*args, **kwargs)))"""
        return self.frommatrix(method(self)(self.asmatrix(), *args, **kwargs))
    return wrapped_method


## This base class controls "__getitem__" behavior and provides a component
#  iterator.
class ABCNumpy(object):
    """This class is for classes that may have singleton on numpy array
    attributes (Vectors, Coordinates, ....)."""

    __metaclass__ = abc.ABCMeta

    ## B/c they need to know
    _got_length = got_length

    ## Object is iterable ONLY If its attributes are-- this method
    #  works of tensors and coordinates, so there is some peg-point
    #  problems--the peg has to be set in an extended method--which
    #  may violate LSP.
    #  \param index a general index
    #  \retval instance of type(self)
    #  If self isn't indexable at run-time.
    def __getitem__(self, index):
        if len(self):  # external poly
            return type(self)(*itertools.imap(itemgetter(index), self.iter()))

    ## reversed(Tensor).
    def __reversed__(self):
        return type(self)(*map(list, map(reversed, self.iter())))

    ## \param index a general index
    #  \param value an instance type(self).
    #  \retval None this is a procedure
    #  \throws exceptions.NonGeometricTypeError on silly value
    #  \throws exceptions.PolymorphismError on going array_like
    #  when self doesn't play that.
    def __setitem__(self, index, value):
        try:
            values = list(value.iter())
            items = list(self.iter())
            if len(items) != len(values):  # guard
                raise TypeError
        except (AttributeError, TypeError):
            from .exceptions import NonGeometricTypeError
            msg = (type(value).__name__ +
                   " does not constitute a " +
                   type(self).__name__)
            raise NonGeometricTypeError(msg)
        else:
            try:
                for item, value in itertools.izip(items, values):
                    operator.setitem(item, index, value)
                    continue
            except TypeError:
                from .exceptions import PolymorphismError
                msg = (
                    type(self).__name__ +
                    " instance isn't built for item assignment"
                    )
                raise PolymorphismError(msg)

    def __delitem__(self, index):
        for item in self.iter():
            operator.delitem(item, index)

    ## len is the len of the components-If they all have the same length
    # \throws ValueError on illdefined length
    def __len__(self):
        lens = []
        for item in self.iter():
            try:
                item_length = len(item)
            except TypeError:
                continue
            else:
                lens.append(item_length)
            continue
        if lens:  # raise
            if lens.count(lens[0]) == len(lens):  # guard
                result = lens[0]
            else:
                # non-unique len--> check it before you wreck it.
                raise ValueError("object length is bad: %s" % str(lens))
            return result
        raise TypeError("{} instance isn't iterable".format(self.__class__))

    ## For container attributes: search, otherwise return NotImplemented so
    #  that the object doesn't look like a container when introspected.
    #  \param other
    #  \retval bool IFF other is in self
    #  \throws exceptions.PolymorphismError IFF self doesn't play arrays
    def __contains__(self, other):
        """any(map(other.__eq__, self))"""
        try:
            return any(itertools.imap(other.__eq__, self))
        except TypeError:
            from .exceptions import PolymorphismError
            msg = (
                type(self).__name__ +
                " instance isn't built like a container"
                )
            raise PolymorphismError(msg)

    ## Container like iter
    #  \retval generator
    #  \throws exceptions.PolymorphismError is iter() was uncalled for.
    def __iter__(self):
        try:
            return itertools.starmap(type(self),
                                     itertools.izip(*self.iter()))
        except TypeError:
            from .exceptions import PolymorphismError
            msg = (
                type(self).__name__ +
                " instance isn't built like a sequence."
                )
            raise PolymorphismError(msg)

    ## Component iter.
    # \yields Geometric Components
    def iter(self):
        """return a generator that generates components"""
        return itertools.imap(self.__getattribute__, self.components)

    ## List of Components
    # \returns list of components
    def tolist(self):
        """return a list of #components"""
        return list(self.iter())

    ## Recast data type
    # \param dtype a numpy dtype
    # \return geometric object, recast
    def asdtype(self, dtype):
        """Cast component arrays astype)"""
        return self.numpy_method('astype', dtype)

    ##  This allows you to broadcast functions (numpy functions) to the
    # attributes and rebuild a class
    # \param func \f$ f \f$
    # \param args arguments to f
    # \param kwargs keywords to f
    # \returns \f$ T'_{ijk} = f(T_{ijk}, *args, **kwargs) \f$
    def broadcast(self, func, *args, **kwargs):
        """vector.broadcast(func, *args, **kwargs) -->
        Vector(*map(func, vector.iter()))

        That is: apply func component wise, and return a new Vector
        """
        f_ak = functools.partial(func, *args, **kwargs)
        return type(self)(*itertools.imap(f_ak, self.iter()))

    ## Mean
    # \returns \f$ T_{i...k} \rightarrow
    # \frac{1}{n}\sum_0^{n-1}{T_{i...k}[n]} \f$
    def mean(self, *args, **kwargs):
        """broadcast np.mean (see broadcast.__doc__)"""
        from numpy import mean
        return self.broadcast(mean, *args, **kwargs)

    ## Sum
    # \returns \f$ T_{i...k} \rightarrow \sum_0^{n-1}{T_{i...k}[n]} \f$
    def sum(self, *args, **kwargs):
        """broadcast np.sum (see broadcast.__doc__)"""
        from numpy import sum as np_sum
        return self.broadcast(np_sum, *args, **kwargs)

    ## Cumulative Sum
    # \returns \f$ T_{i...k}[n] \rightarrow \sum_0^{n-1}{T_{i...k}[n]} \f$
    def cumsum(self, *args, **kwargs):
        """broadcast np.cumsum (see broadcast.__doc__)"""
        from numpy import cumsum
        return self.broadcast(cumsum, *args, **kwargs)

    ## Apply np.ravel to components
    def ravel(self, order=0):
        """broadcast numpy.ravel"""
        from numpy import ravel
        return self.broadcast(ravel, order=order)

    ## Conjugate via components' conjugate method
    # \returns \f$ T_{ij\ldots} \rightarrow T^{*}_{ij\ldots} \f$
    def conjugate(self, method_name='conjugate'):
        """broadcast conjugate methodcaller"""
        return self.numpy_method(method_name)

    ## Call a numpy method by name
    def numpy_method(self, method_name, *args, **kwargs):
        """np.method(method_name, *args, **kwargs)

        broadcast np.ndarray method (see broadcast.__doc__)"""
        return self.broadcast(operator.methodcaller(method_name,
                                                    *args, **kwargs))

    ## For a tensor, chart, or coordinate made of array_like objects:\n
    #  append a like-wise object onto the end
    def append(self, other):
        """For array_like attributes, append and object into the end."""
        #pylint: disable=E0203
        from numpy import append
        result = type(self)(*itertools.starmap(
            append,
            itertools.izip(self.iter(), other.iter())))
        ## TODO: get rid of this attrocity, pegged thingies need to
        # extend this method with the peg adder- otherwise it is burden
        # for all other objects in geo's universe.
        try:
            peg = self.peg
        except AttributeError:
            pass
        else:
            result.peg = peg
        return result

    ## The idea is to check equality for all tensor/coordinate types and for
    #  singleton and numpy arrays, \n so this is pretty concise-- the previous
    #  version was 10 If-then blocks deep.
    # \param other anything
    # \returns bool or ndarray of dtype=bool_.
    def __eq__(self, other):
        """Check tensor equality, for each component.

        Returns BOOL or INT (0,1) for:

        (1) inter-type comparison
        (2) non array-like comparison

        Othersie, returns numpy.ndarray(dtype=bool) for:

        (1) Any broadcastable comparison (intra-type).
        """
        try:
            iter_o = other.iter()
        except AttributeError:
            # or raise an Error?
            pass
        else:
            #pylint: disable=E0602
            for count, (left, right) in enumerate(itertools.izip(self.iter(),
                                                                 iter_o)):

                ## Set the value OR array of values, who knows?
                value = (left == right)
                if not count:  # initializer
                    # on 1st pass: set values
                    result = value
                else:
                    # on 2+ pass, update value
                    result *= value
                    # go to next case
                    continue
                # can't break on False, b/c of array_like behavior
                continue
            # result may be array_like
            return result
        # result is straight False after TypeError in comparison
        return False

    ## not __eq__,  while per preserving array_like behavior- not easy --
    #  note that function/statement calling enforces type checking
    #  (no array allowed),
    #  \n while method calling does not.
    def __ne__(self, other):
        """See __eq__"""
        inv = self.__eq__(other)
        try:
            result = operator.not_(inv)
        except ValueError:
            #pylint: disable=E1103
            result = (1-inv).astype(bool)
        return result

    ## This covers <   <= >  >= and raises a RuntimeError, unless numpy
    #  calls it, and then that issue is addressed (and it's complicated)
    def __cmp__(self, other):
        """This method is called in 2 cases:

        Vector Comparison: If you compare (<, >, <= , >= ) too non-scalar
        tensors, or rotations--that makes no sense and you get a TypeError.

        Left-Multiply with Numpy: this is a little more subtle. If you are
        working with array_like tensors, and say, do

        >>> (v**2).w*v   #instead of:
        >>> (v**2)*v     #(which works), numpy will take over and call
        __cmp__ in order to figure how to do linear algebra on
        the array_like components of your tensor.
        The whole point of the Scalar class is to avoid this
        pitfall.


        Side Note: any vector operation should be manIfestly co variant--
        that is-- you don't need to access the "w"-- so you should not get
        this error.

        But you might.
        """
        raise TypeError(
            """comparison operation not permitted on %s\ncheck for left
            "mul by a numpy array.\n Right operand  is %s""" % (
                type(self).__name__,
                type(other).__name__))

    ## In principle, we always want true division: this is here in case a
    #  client calls it
    def __truediv__(self, other):
        return self.__div__(other)

    ## In principle, we always want true division: this is here in case a
    #  client calls it
    def __rtruediv__(self, other):
        return self.__rdiv__(other)

    ## This is risky and pointless-- who calls float?
#    def __float__(self):
#        return float(abs(self).w)

    ## +T <a href="http://docs.python.org/reference/expressions.html#is">is</a>
    #  T   --- is as in python\n This method is used in coordinate transforms
    #  when an identity transform is required, dummy keyword allows
    #  polymorphism.
    def __pos__(self, peg=None):
        """+obj is obj"""
        #pylint: disable=W0613
        return self

    ## Get some numpy properties from the instance- numpy may request this.
    # \param attrgetter_ ?
    # \returns shapes of something
    # \throws TypeError for numpy.
    def _array_property(self, attrgetter_):
        shapes = []
        for item in self.iter():
            try:
                shapes.append(attrgetter_(item))
            except AttributeError:
                pass
            continue
        if shapes:  # external poly
            if shapes.count(shapes[0]) == len(shapes):  # raise
                return shapes[0]
        raise TypeError

    @property
    def shape(self):
        """Numpy shape of elements."""
        return self._array_property(operator.attrgetter('shape'))

    @property
    def ndim(self):
        """Numpy ndim of elements."""
        return self._array_property(operator.attrgetter('ndim'))

    @property
    def size(self, x=0):
        """Numpy size."""
        for count, item in enumerate(self.iter()):
            x = item.size
            if count:  # initializer
                if x != result[-1]:  # raise
                    #pylint: disable=E0601
                    msg = "Ill formed array structures in {}"
                    raise TypeError(msg.format(self))
                result.append(x)
            else:
                result = [x]
        return x

    ## Broadcast trig.Re
    # \f$ \Re{({\bf\vec{v}})} =
    # \Re{(v_x)}{\bf\hat{e}}_x +
    # \Re{(v_y)}{\bf\hat{e}}_y +\Re{(v_z)}{\bf\hat{e}}_z \f$ \n
    # and likewise for other ranks.
    @property
    def real(self):
        """v.real"""
        return self.broadcast(trig.Re)

    ## Broadcast trig.Im
    # \f$ \Im{{\bf\vec{v}}} = \Im{v_x}{\bf\hat{e}}_x +
    # \Im{v_y}{\bf\hat{e}}_y +\Im{v_z}{\bf\hat{e}}_z \f$ \n
    # and likewise for other ranks.
    @property
    def imag(self):
        """v.imag"""
        return self.broadcast(trig.Im)
