"""The purpose of the module is 2-fold:

(1) Help the user understand exatcly where they went wrong when using a
geomteric object.

(2) Make the user feel inqdequate--even Nilpotent, and that end, the
exceptions are giving absurdly wonky names."""

## \namespace geo::utils::exceptions
# <a href="http://docs.python.org/2/library/exceptions.html">Exceptions</a>


## This function should make a generic error message for bad operations
def error_message(op, left, right):
    """error_message(op, left, right)

    for class problems with <left><op><right>
    """
    return "%s(%s, %s)" % (op.__name__,
                           type(left).__name__,
                           type(right).__name__)


## Raised when you request a nonsense tensor index
class DimensionError(UserWarning, ValueError):
    """Thrown when a bad dimension is requested"""


## Tensor Bads.
class TensorError(UserWarning):
    """Some kind of Tensor Problem"""


## Mistreating __R3__  and __M4__ scalars like lowely __Z__ integers.
class ScalarAreNotNumbersError(TensorError, TypeError):
    """Scalars are not integers, or even numbers. They are tensors that
    share a lot of behaviors with them, but moduo arithmetic is not
    one of them. It simply does not have a geometric interpretation,
    and this exception is here to remind you, dear User, of that fact."""

    def __init__(self, *args):
        super(ScalarAreNotNumbersError, self).__init__(
            "I'm, not a number... I..  AM..  A TENSOR! {}".format(args)
            )


## Tensor index's value does not make sense
class TensorIndexValueError(TensorError, ValueError):
    """Tensor Index Problems."""
    
    def __init__(self, n, comment=""):
        msg = "Index {} is invalid".format(n)
        if comment:  # kwd
            msg += ": {}".format(comment)
        super(TensorIndexValueError, self).__init__(msg)


## Something is wrong with a set of indicies.
class IndicesError(TensorError, ValueError):
    """Indices are not appropriate for operations (wrong number,
    or duplicates)."""
        


## Raise for wrong rank, without reference to indices.
class RankError(TensorError, TypeError):
    """Bad Rank in an OP."""

    def __init__(self, n, m):
        super(RankError, self).__init__(
            "Operation required rank {} tensor, got {}".format(n, m)
        )



## Base class for geometric errors
class GeometricException(UserWarning, TypeError):
    """Base for Geometry Problems"""


## A reminder to treat geometric objects properly.
class NonCovariantOperation(GeometricException):
    """Raise when you do something that is silly[1], like adding
    a Scalar to a Vector.
    [1]Silly: (adj.) syn: non-covariant"""


## A reminder that Affine space are affine, and vector spaces are not.
class AffineSpaceError(GeometricException):
    """Raised when you forget the points in an affine space are
    not vectors in a vector space, and visa versa"""


## A catch-all for overloaded operations getting nonsense.
class UndefinedGeometricOperation(GeometricException):
    """This will raised If you do an operation that has been defined for
    a Tensor/Affine/Coordinate argument, but you just have a non-sense
    combination, like vector**vector.
    """


## Clifford Algebra is not Implemented for general cases.
class CliffordAlgebraError(UserWarning, TypeError):
    """Error raised when you do an operation that would be kosher
    in a Cl(3, 0) implementation of this package"""


## A senisible operator that has not been implemented
class OverloadNotImplementedError(UserWarning, NotImplementedError):
    """Raised when a sensible overload is not implemented,
    such as v/T being v * (T**-1)"""


## Inverting a non invertable affine transform
class DegenerateInversionError(UserWarning, ZeroDivisionError):
    """Attempt to invert and non-invertable affine transformation"""


## Raised when inverting a projection
class ProjectionInversionError(DegenerateInversionError):
    """Projection Inversion"""


## Raised when inverting a rejection
class RejectionInversionError(DegenerateInversionError):
    """Rejection Inversion"""


## Raise when you set a geometric object to something it can't be
class NonGeometricTypeError(GeometricException, TypeError):
    """Value does not make a geometric object"""


## Thrown when mixing singleton and array_like concretions badly
class PolymorphismError(UserWarning, TypeError):
    """Geo object is not array_like, and you treat it as such (or visa
    versa)"""


## E.g., parallel transport in non-cartesian frame.
class KoszulConnectionError(AffineSpaceError):
    """trying to do something without a fibre bundle, when you need one"""


## Creating non-sense SU(2) representions
class LieAlgebraError(UserWarning):
    """Nonsense parameters for a Lie Algebra."""


## Accessing a spherical tensor out-of bounds
class RepresentationError(LieAlgebraError, IndexError):
    """|n> when n is not right"""


## Representation is not (1/2) integer
class CasimirInvariantError(LieAlgebraError, ValueError):
    """Degree is not integer (or 1/2 integer)"""


## Nilpotentcy
class NilpotencyError(LieAlgebraError, ValueError):
    """|order| is too big"""


## Sliceing and indexable
class KetError(LieAlgebraError, TypeError):
    """Nonsense indices for ket container."""


## HomomorphismError
class HomomorphismError(UserWarning):
    """Raised when a homomorphism doesn't work."""
