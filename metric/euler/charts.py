"""Charts on SO(3) are the various way to express rotations, and these charts
are just some of the angle based representations.

Euler Angles (Active reps, intrinsic):

e1 * e2 * e1'

>>>Euler(alpha, beta, gamma)

Tait Bryan Angles (Active rep, intrinsic):

e1 * e2 * e3

>>>YPR(alpha, beta, gamma)
>>>RPY(alpha, beta, gamma)

Of course, you can mix and match rotations as needed, e.g:

>>>YPR(145,0,20) * Euler(180,20,35)

is the identity rotation. They also play with Versors and Tensors.

"*" Composes operators
"**" interpolates (from the identity)

pow(dest, fraction, start=Identity)

does SLERP from start, fractionally to dest.
"""
## \namespace geo.metric.euler.charts
# <a href="http://en.wikipedia.org/wiki/Charts_on_SO(3)">Charts in SO(3)</a>
# for rotations.

import abc
import functools

from ...metric.euclid import euclid
from ...metric.euclid import vector

#pylint: disable=R0904

__all__ = ('YPR', 'RPY', 'roll', 'pitch', 'yaw', 'Euler', 'euler_maker')


DEGREE = 360.


## Decorator converts numeric axes to intrinsic vector.Vector instances.
def intrinsic_(cls):
    """1, 2, 3 <--> x, y, z"""
    cls.AXES = tuple([vector.BASIS[n-1] for n in cls.AXES])
    return cls


## Decorator converts numeric axes to extrinsic vector.Vector instances.
def extrinsic_(cls):
    """1, 2, 3 <--> x, y, z"""
    cls.AXES = tuple([vector.BASIS[n-1] for n in reversed(cls.AXES)])
    return cls


## A Base class for <a href="http://en.wikipedia.org/wiki/Euler_angles">
# Euler Angles</a>: it defines operations, but does not define axis order or
# units.
# @image html euler.gif
# @image latex euler.gif
class ElementalRotationTriplet(euclid.Geometric, euclid.LinearMap):
    """ABC(alpha, beta, gamma)
    are the 3 Euler's angles"""

    __metaclass__ = abc.ABCMeta

    ## \f$ (\alpha, \beta, \gamma) \f$
    components = ('alpha', 'beta', 'gamma')

    ## Define angular unit (as degrees)
    circumference = DEGREE

    ## Abstract Property
    @abc.abstractproperty
    def AXES(self):
        """Subs need a tuple of axes of rotation as a 3-tuple
        of axes of rotation, e.g: (3, 1, 3)"""

    ## \f$ (\alpha, \beta, \gamma) \f$ -- units are unknown
    # \param alpha 1st rotation angle
    # \param beta 2nd rotation angle
    # \param gamma 3rd rotation angle
    def __init__(self, alpha, beta, gamma):
        ## \f$ \alpha \f$, 1st rotation
        self.alpha = alpha
        ## \f$ \beta \f$, 2nd rotation
        self.beta = beta
        ## \f$ \gamma \f$, 3rd rotation
        self.gamma = gamma

    ## Generic Inversion via hamilton.Versor
    # \return multiplicative-inverse.
    def __invert__(self):
        return (~(self.versor())).aschart(type(self))

    ## Generic mul via hamilton.Versor
    def __mul__(self, other):
        return (self.versor() * other.versor()).aschart(type(self))

    ## Generic or via hamilton.Versor
    def __or__(self, other):
        return self.versor() | other.versor()

    ## Generic lshift via hamilton.Versor
    def __lshift__(self, other):
        return (self.versor() << other.versor()).aschart(type(self))

    ## Generic rshift via hamilton.Versor
    def __rshift__(self, other):
        return (self.versor() >> other.versor()).aschart(type(self))

    ## \returns geo.desic.hamilton.quaternion.Quaternion
    def quaternion(self):
        """As a quaternion."""
        return self.versor().quaternion()

    ## SLERP is used!
    def __pow__(self, other, modulo=None):
        """3-arg pow is SLERP, see Versor.__pow__.__doc__"""
        return pow(self, other, p=modulo).aschart(type(self))

    @property
    def I(self):
        """inverse"""
        return ~self

    ## Transpose is Inverse
    def transpose(self):
        """inverse."""
        return ~self

    ## Alias transformation of argument, using versor(), but effectively:\n
    ## \f$  {\bf \vec{v}'} = {\bf \vec{v} \cdot M} \f$
    # \param other Any transformable
    # \returns aliased transformable
    def alibi_transform(self, other):
        """Apply transformation to argument"""
        return self.versor()(other)

    ## rotation, counted from 0.
    # \param n index of rotation wrt AXES
    # \returns geo.metric.euclid.versor.versor.Versor repping the rotation
    def _rotation(self, n):
        """for the 0, 1, 2 rotations, perform rotation 'n'"""
        return self.AXES[n].versor(getattr(self, self.components[n]),
                                   circumference=self.circumference)

    ## get 1st, 2nd, or 3 rotation Versor
    # \param n ordinal count of Euler rotation
    # \returns geo.metric.euclid.versor.versor.Versor repping the rotation
    def rotation(self, n):
        """versor = rotation(n) for n = 1, 2, 3

        gets the Versor representing the n-th rotation.
        """
        return self._rotation(n-1)

    ## Compose the 3 rotations using chain()
    # \returns geo.metric.euclid.versor.versor.Versor that is the rotation
    def versor(self):
        """Compute the equivalent Versor for all three rotations."""
        # Note the order was correctly set in 9/1/15, after the transition
        # from Alias to Alibi default.
        return self.chain(*map(self._rotation, xrange(euclid.DIMENSIONS)))

    @property
    def vector(self):
        """vector version of tensor is axis of rotation."""
        print "future warning... property of method?"
        return self.tensor().vector()

    ## Rotation Angle
    # \returns \f$ \theta \f$
    def angle_of_rotation(self):
        """Angle of rotation."""
        return self.tensor().angle_of_rotation()

    ## Rotation Axis
    # \return \f$ {\bf \hat{n}} \f$
    def axis_of_rotation(self):
        """Axis of rotation."""
        return self.tensor().axis_of_rotation()

    ## Convert to YPR
    # \returns YPR
    def ypr(self):
        """As a YPR"""
        return self.versor().ypr()

    ## Convert to RPY
    # \returns RPY
    def rpy(self):
        """As RPY"""
        return self.versor().rpy()

    ## Convert to <a href=http://en.wikipedia.org/wiki/Direction_cosine">
    # DCM
    # \returns euclid.tensor.Tensor
    def tensor(self):
        """As a rotation matrix (Tensor)."""
        return self.versor().tensor()

    ## Because it really is just a matrix full of scalars.
    dcm = tensor

    ## Convert to Euler angles
    # \returns Euler
    def euler(self):
        """As euler angles"""
        return self.versor().euler()


## Traditional Euler Angles: \f$(z, x', z'')\f$
# @image html Euler.png
# @image latex Euler.png
@intrinsic_
class Euler(ElementalRotationTriplet):
    """euler = Euler(alpha, beta, gamma)

    For intrinsic rotations about:
    Z (alpha) Precession
    X' (beta) Nutation
    Z" (gamma) Spin

    with angles in degrees.
    """
    AXES = (3, 1, 3)

    ## Identity
    euler = euclid.Geometric.__pos__

    ## And this method explains Euler's angles utility in a time before CPU.
    def __invert__(self):
        return type(self)(*reversed(list((-self).iter())))

    ## Invert sign of \f$\beta\f$, keeping same rotation.
    def flip_beta(self):
        """flip beta (rotation unchanged)"""
        return type(self)(
            (self.alpha + self.circumference/2.) % self.circumference,
            -self.beta,
            (self.gamma + self.circumference/2.) % self.circumference)


## Base class for Tait Bryan (3-spearate axis) rotations.
class TaitBryan(ElementalRotationTriplet):
    """Base class for TaitBryan: designed for 3 separate axes
    of rotation."""


## Yaw -> Pitch -> Roll
# @image html ypr.png
# @image latex ypr.png
@intrinsic_
class YPR(TaitBryan):
    """YPR(yaw, pitch, roll) --all in degrees
    and in that order, polymorphic with Versors and rotation matrices.
    """
    ## Yaw, Pitch, and *then * Roll
    AXES = (3, 2, 1)

    ## Identity
    ypr = euclid.Geometric.__pos__

    ## \f$ \alpha \f$
    @property
    def yaw(self):
        """yaw-->alpha"""
        return self.alpha

    ## \f$ \beta \f$
    @property
    def pitch(self):
        """pitch-->beta"""
        return self.beta

    ## \f$ \gamma \f$
    @property
    def roll(self):
        """roll-->gamma"""
        return self.gamma

    ## Negate the reversed inverse.
    def rpy(self):
        """As an RPY."""
        from operator import neg
        return RPY(*reversed(map(neg, (~self).iter())))


## Roll -> Pitch -> Yaw
@intrinsic_
class RPY(TaitBryan):
    """RPY(roll, pitch, yaw) --all in degrees
    and in that order
    """
    ## Yaw, Pitch, and *then * Roll
    AXES = (1, 2, 3)

    ## Identity
    rpy = euclid.Geometric.__pos__

    ## \f$ \gamma \f$
    @property
    def yaw(self):
        """yaw-->gamma"""
        return self.gamma

    ## \f$ \beta \f$
    @property
    def pitch(self):
        """pitch-->beta"""
        return self.beta

    ## \f$ \alpha \f$
    @property
    def roll(self):
        """roll-->alpha"""
        return self.alpha


## Roll coordinate transformation (in degrees) \n
# (http://en.wikipedia.org/wiki/Flight_dynamics)
# @image html roll.gif
# @image latex roll.gif
roll = functools.partial(vector.X.versor, circumference=DEGREE)

## Pitch coordinate transformation (in degrees) \n
# (http://en.wikipedia.org/wiki/Flight_dynamics)
# @image html pitch.gif
# @image latex pitch.gif
pitch = functools.partial(vector.Y.versor, circumference=DEGREE)

## Yaw coordinate transformation (in degrees) \n
# (http://en.wikipedia.org/wiki/Flight_dynamics)
# @image html yaw.gif
# @image latex yaw.gif
yaw = functools.partial(vector.Z.versor, circumference=DEGREE)


## Euler Angle Maker: chose your own axes.
class euler_maker(type):
    """Class = euler_maker(axis1, axis2, axis3 [, Intrinsic=True])

    For example:

    Euler = euler_maker(3, 1, 3, Intrinsic=True)
    """
    ## \param axis1 number
    # \param axis2 number
    # \param axis3 number
    # \kwd Intrinsic = True intrinsic() or extrinsic() angles.
    # \return cls Euler (TaitBryan) class is axis1 == (!=) axis3
    def __new__(mcs, axis1, axis2, axis3, Intrinsic=True):
        if axis2 in (axis1, axis3):  # raise
            raise ValueError("axis2 is Degenerate")
        axes = (axis1, axis2, axis3)
        return (intrinsic if Intrinsic else extrinsic)(  # kwd
            type('Euler%i%i%i' % axes,
                 ((Euler if axis1 == axis3 else TaitBryan),),
                 {'AXES': axes}))


## Euler Angle Maker: chose your own axes.
class euler_maker(type):
    """cls = euler_maker(axis1, axis2, axis3 [, Intrinsic=True])

    For example:

    Euler = euler_maker(3, 1, 3, intrinsic=True)"""
    
    def __new__(mcs, axis1, axis2, axis3, intrinsic=True):
        if axis2 in (axis1, axis3):  # raise
            raise ValueError("axis2 is Degenerate")
        decorator = intrinsic_ if intrinsic else extrinsic_  # kwd 
        base = Euler if axis1 == axis3 else TaitBryan  
        name = '{}{}{}{}'.format(base.__name__, axis1, axis2, axis3)
        axes = dict(AXES=(axis1, axis2, axis3))
        cls = type(name, (base,), axes)
        return decorator(cls)
