"""This module is full of abstract base classes for coordinates. They need an
ellipsoid to give them concrete meaning.

To use: 1st create the ellipsoid.Ellipsoid, and then use its methods:

ECEF, LLH, LTP, SCH

to make coordinates.

You'll may use:

origin.PegPoint  (or a tuple, If you must), to defines peg points.
-------------------------------------------------------------------------
The classes are completely polymorphic-- you don't need to know which you
have to work with an instance:

If you have an instance p:

p.ecef()
p.llh()
p.ltp(peg=None)
p.sch(peg=None)

which give the same point in a new coordinate system. Doesn't matter what p is
to start with. Likewise, of you need to change an ellipsoid, say, to airy1830,
just do:

p_new = p | airy1830

Doesn't matter what p is, p_new will be the same type of coordinate. Don't
operate on coordinates from dIfferent ellipses. [The feature that prevented
this--dynamic class creation--was deprecated so that objects could be pickled]

Note on Affine spaces: you can dIfference points and get vectors. You can add
vector to points. But you can't add two points.

Note: all coordinates have a cartesian counter part, which for a cartesian
system has the same orientation and origin:

LLH.otf --> ECEF
SCH.otf --> LTP   (with the same PegPoint)

With that, subtraction is defined for all coordinates classes:

v = p2 -p1

will give a vector in p2's cartesian counter part, and the inverse operation:

p2 = p1 +v

will give a new point p2 in the same type a p1. (So If p1 is in SCH, then v
will be interpreted as a vector in the congruent LTP). It's all about
polymorphism. If you don't like it, then write your code explicitly.


Methods:
--------
Coordinates have all the generalized methods described in euclid.__doc__, so
see that If you need to take an average, make an iterator, a list, or what not.

p1.bearing(p2)
p1.distance(p2)   compute bearing and distance bewteen to points on the
                  ellipsoid.

p.latlon()        get the latitude longitude tuple.
p.vector()        make a euclid.Vector
p.spacecurve()   makes a motion.SpaceCurve


Transformation functions are computed on the fly (& memoized),
and are available via:


ECEF._f_ecef2lla
ECEF.affine2ltp

LLH._f_lla2ecef

LTP.f_ltp2sch
LTP.affine2ecef

SCH.f_sch2ltp


these are either affine.Affine transformations, or nested dynamic functions--
either way, they do not exist unit requested, and at that point, the get made,
called, and forgotten.

Note on architecural dissatisfaction:

coordiantes are 3 things.

The addition of orgin + ellipsoid (a datum) complete i

"""
#pylint: disable=E1101,E1102

## \namespace geo::detic::coordinates Abstract Coordinates on Ellipsoids.
import abc
import functools
import operator

from ..utils import arraylike
from ..utils.exceptions import AffineSpaceError

from . import origins
from . import christoffel


## Create a single WGS84 SCH point as lat, lon, hdg:
# \param lat \f$ \phi \f$
# \param lon \f$ \lambda \f$
# \param hdg \f$ \psi \f$
# \returns SCH on WGS84 with s,c,h=0.
def sch_origin_at(lat, lon, hdg):
    """Return an SCH point at peg."""
    from .. import WGS84
    return WGS84.LLH(lat, lon, 0.).sch(origins.Peg(lat, lon, hdg))


## A decorator to add peg point to the results of extened methods
# \param method A method for pegged coordinates.
# \returns method with peg setting.
def pegged(method):
    """decorator for methods that take a peg"""

    @functools.wraps(method)
    def pegged_method(self, *args, **kwargs):
        """Call method with args, and then add pegpoint post facto"""
        # call nominal method
        result = method(self, *args, **kwargs)
        # add self's peg point to result
        result.peg = self.peg
        return result

    return pegged_method


## a coordinate abstract base class.
class _Coordinate(arraylike.ABCNumpy):
    """Abstract Base Class for a Coordinate triplets."""

    __metaclass__ = abc.ABCMeta

    ## Init 3 coordinates from class's cls.coordinates
    # the peg point keyword is under investigation
    # \param coordinate1 1st coordinate
    # \param coordinate2 2nd coordinate
    # \param coordinate3 3rd coordinate
    # \kwd ellipsoid=None Don't touch this.
    def __init__(self, coordinate1, coordinate2, coordinate3, ellipsoid=None):
        # Set triplet of coordinate attributes
        for name, value in zip(self.coordinates,
                               (coordinate1, coordinate2, coordinate3)):
            setattr(self, name, value)
        from .. import WGS84
        self.ellipsoid = ellipsoid or WGS84

    ## This is a protected class thing
    # \returns _Coordinates._ellipsoid
    @property
    def ellipsoid(self):
        """Ellipsoid is a class attribute"""
        from .newton.wgs84 import WGS84
        return self._ellipsoid or WGS84


    ## No Touchy
    # \throws TypeError
    @ellipsoid.setter
    def ellipsoid(self, e):
        """Ellipsoid is a class attribute, and you can't set it"""
        self._ellipsoid = e

    ## This method allows polymorphic constructor calls for Pegged and
    # Fixed classes.
    @abc.abstractmethod
    def kwiter(self):
        """Keyword iterator for subs"""

    ## Coordinate attribute names
    @abc.abstractproperty
    def coordinates(self):
        """Names of the coordinates, in order."""

    ## Put other's in my coordinate system
    @abc.abstractmethod
    def _like_me(self, other):
        """Makes others like self."""
        return other | self.ellipsoid #?

    ## Convert to ECEF
    @abc.abstractmethod
    def ecef(self):
        """Convert to ECEF."""

    ## Convert to LLH
    @abc.abstractmethod
    def llh(self):
        """Convert to Geodetic."""

    ## Convert to LTP (aka pegpoint xyz)
    # \kwd peg (optional)
    @abc.abstractmethod
    def ltp(self, peg=None):
        """Convert to Local Tangent Plane at peg."""

    ## Convert to SCH
    # \kwd peg (optional)
    @abc.abstractmethod
    def sch(self, peg=None):
        """Convert to Local Tangent Sphere at peg."""

    ## replace coordinates with components.
    # \returns coidinates().
    @property
    def components(self):
        """components are coordinates."""
        return self.coordinates

    ## Copy self.
    def __copy__(self):
        return type(self)(*self.iter(), **self.kwiter())

    ## \f$ -p \rightarrow vp+ (-2*v) \f$
    # \throws utils.exceptions.AffineSpaceError
    def __neg__(self):
        msg = ("Cannot negate an affine space point (p). Be Explicit: \n" +
               " >> >p + 2*(-p.vector())")
        raise AffineSpaceError(msg)

    ## Get n-th coordinates in a name independent manner.
    # \param n Integer in (0, 1, 2)
    # \returns Attribute that is the n-th component
    # \throws utils.exceptions.DimensionError
    def e(self, n):
        """P.e(n) --> n-th component of P"""
        try:
            result = getattr(self, self.components[int(n)])
        except (TypeError, IndexError):
            from ..utils.exceptions import DimensionError
            raise DimensionError("Invalid dimension %s" % str(n))
        return result

    ## Convert result to a Vector relative to origin
    # \retval vector geo.metric.euclid.vector.Vector from origin.
    def vector(self):
        """vector([peg=None])

        will convert self.x, self.y, self.z into a euclid.Vector If peg
        is None, otherwise, it will transform to the LTP defined by the
        peg and the call vector(None).
        """
        from ..metric.euclid import vector
        return vector.Vector(*self.otf().iter())

    ## Make into a SpaceCurve().
    # \returns geo.metric.frenet_seret.motion.SpaceCurve
    def spacecurve(self):
        """make a SpaceCurve"""
        return self.vector().spacecurve()

    ## string
    def __str__(self):
        return "".join(map("{}: {}\n".format, self.components, self.iter()))

    ## repr() and str() are now the same.
    __repr__ = __str__

    ## p1 << p2 : recaste p2 into p1's frame
    # \throws TypeError
    def __lshift__(self, other):
        ## todo: push error up.
        try:
            return self._like_me(other)
        except AttributeError:
            raise TypeError("excepted a Coordinate, got: %s" %
                            type(other).__name__)

    ## Recaste self into other's coordinates
    def __rshift__(self, other):
        return other << self

    ## crd | ellipsoid put coordinates on an ellipsoid
    # \param target A geo.detic.newton.ellipsoid.Ellipsoid
    # \retval _Coordinate self referenced to the new ellipsoid.
    def __or__(self, target):
        """change ellipsoids...."""
        return self.change_ellipsoid(target)

    ## Get like-named class on another ellipsoid.Ellipsoid by using
    # methodcaller and the class's __name__; this avoids eval or exec calls.
    # \param target A geo.detic.newton.ellipsoid.Ellipsoid
    # \retval there_ecef self on target's ECEF (that's 1.2 the problem, subs
    # do the rest).
    def change_ellipsoid(self, target):
        """Move coordinates to a new ellipsoid's ECEF--subs will take
        it from there."""
        x, y, z = self.ecef().iter()
        try:
            there_ecef = target.ECEF(x, y, z)
        except AttributeError as err:
            if isinstance(target, type(self.ellipsoid)):  # raise
                raise err
            raise TypeError(
                "change_ellipsoid expects an Ellipsoid, got: %s" %
                type(target).__name__
            )
        return there_ecef

    ## get EGM08 WGS84 ellipsoid relative to coordinate's ellipsoid
    # \param geoid_ A geoid
    # \retval array_like Geoid height at points (a call to geoid_)
    # \throws ValueError If ellipsoid is NOT .ellipsoid.wgs84.WGS84.
    def _geoid(self, geoid_):
        from .newton.wgs84 import WGS84
        if self.ellipsoid == WGS84:  # guard (only WGS84 works).
            llh_wgs84 = self.llh()
            dh = 0.
        else:
            raise ValueError("Geoid is for WGS84, not %s" %
                             self.ellipsoid.model)
        return geoid_(llh_wgs84.lat, llh_wgs84.lon) + dh

    ## Return egm08 height corrections
    def egm08(self):
        """EGM96 heights"""
        from ..ids.egm08 import geoid
        return self._geoid(geoid)

    ## Return egm96 height corrections
    def egm96(self):
        """EGM08 heights"""
        from ..ids.egm96 import geoid
        return self._geoid(geoid)

    ## Convert to LLH and get the bearing bewteen two coordinates, see module
    # function bearing() for more.
    # \param other Another Coordinate
    # \retval bearing The bearing from self to other
    def bearing(self, other):
        """psi = p1.bearing(p2) is the heading from p1 to p2."""
        return origins.bearing(*(self.llh().latlon() + other.llh().latlon()))

    ## Get distance between nadir points on the Ellipsoid, self
    # ellipsoid.Ellipsoid.distance()
    # \param other Another Coordinate
    # \retval bearing The approx. great-circle distance from self to other
    def distance_spherical(self, other):
        """calls self.ellipsoid.distance_spherical"""
        return self.ellipsoid.distance_spherical(
            *(self.llh().latlon()+other.llh().latlon())
            )

    ## \param other Another Coordinate
    # \retval bearing The exact great-circle distance from self to other
    def distance_true(self, other):
        """calls self.ellipsoid.distance_true"""
        return self.ellipsoid.distance_true(
            *(self.llh().latlon() + other.llh().latlon())
            )

    ## Choose distance method
    distance = distance_true

    def great_circle(self, other):
        """see ellipsoid.great_circle."""
        return self.ellipsoid.great_circle(
            *list(self.latlon() + other.latlon()))

    ## Get the origin.Graticule at self.PegPoint
    # \retval graticlue origins.Graticlue (lat, lon).
    def latlon(self):
        """makes a (lat, lon) namedtuple (Graticule)."""
        return origins.Graticule(*(self.llh().tolist()[:-1]))

    ## Subtraction ALWAYS does vector conversion in the left arguments
    # cartesian frame.
    # \param other A coordinate object
    # \retval vector self-other A geo.metric.euclid.vector.Vector from self
    # to other
    # \throws geo.utils.exceptions.AttributeError
    # \throws AttributeError
    def __sub__(self, other):
        """point1 - point2 gives the vector pointing from point1 to point2,
        in point1's cartesian coordinate system."""
        try:
            return self.vector() - (other.otf() >>
                                    self.otf()).vector()
        except AttributeError as err:
            from .. import Vector
            if isinstance(other, Vector):  # raise
                raise AffineSpaceError(
                    "can't subtract vector from a point. Try pint + (-vector)"
                    )
            raise err

    ## You can only add a vector, in the coordaintes cartesian frame, and you
    # get back an affine point in the same coordiantes--
    # \param vector_ A geo.metric.euclid.vector.Vector
    # \retval _Coordinate A new coordinate
    # \throws geo.utils.exceptions.AttributeError
    def __add__(self, vector_):
        """point1 + vector gives point2 in point1's coordinate system, with
        the vector interpreted in point1's cartesian frame"""
        # Do not allow adding of affine points (it is meaningless,
        # though it is computable, hence this guard statement is justified.
        if isinstance(vector_, _Coordinate):  # raise
            msg = ("Cannot add 2 affine points-- add a vector")
            raise AffineSpaceError(msg)
        # Add vector to self's point in otf, unpack coordinates. All that in
        # one line- for all four coordinate types. wow.
        return type(self.otf())(
            *(self.otf().vector() + vector_).iter(), **self.kwiter()
        ) >> self


## "Private" mix-in for _Coordinate with a fixed origin
class _Fixed(_Coordinate):
    """_Fixed is a base class for coordinates with a fixed origin"""

    __metaclass__ = abc.ABCMeta

    ## These have no orgins.PegPoint, explicitly.
    peg = None

    ## Keyword ITERator
    # \returns dict {'peg': self.peg}
    def kwiter(self):
        """{} is a keyword iterator-- it adds peg=peg when needed."""
        return dict(ellipsoid=self.ellipsoid)

    def pkwiter(self):
        """stub"""
        return dict()

    ## Use _Coordinate.change_ellipsoid() and call the result's method
    # to self's coordinates (without peg keyword)
    # \param target A geo.detic.ellipsoid.ellipsoid.Ellipsoid
    # \retval _Coordinate self referenced to the new ellipsoid
    def change_ellipsoid(self, target):
        """change ellipsoid to target, and re-compute coordinate"""
        result = super(_Fixed, self).change_ellipsoid(target)
        result = _Coordinate.change_ellipsoid(self, target)


        method = operator.methodcaller(type(self).__name__.lower(),
                                       **self.pkwiter())
        return method(result)

    ## ENU..
    # \param lat Latitude of origin
    # \param lon Longitude of origine
    # \retval LTP self in ENU at (lat, lon)
    def enu(self, lat, lon):
        """
        p.enu(lat, lon) ltp with an ENU peg at (lat, lon)
        """
        return self.ltp(origins.PegPoint(lat, lon, origins.EAST))

    #pylint: disable=W0221
    ## Check equality in affine space, not computer space.
    def __eq__(self, other):
        """p == p' means that the 2 point are the same on Earth, regarless
        of their coordinate representation-- within a 1mm tolerance. Should
        you need to change that, then use:


        p.__eq__(p', tolerance=tolerance)

        This returns a bool, or an array_like result."""
        return abs(self._like_me(other) - self) < self.tolerance


#pylint: enable=W0221
## "Private" mix-in for Coordinate with a run0time decided origin, thus
# it needs to extend method with a peg key-word (violating LSP)
class _Pegged(_Coordinate):
    """_Pegged is a base class for coordinates with a variable origin """

    __metaclass__ = abc.ABCMeta

    ## \param coordinate1 1st coordinate
    # \param coordinate2 2nd coordinate
    # \param coordinate3 3rd coordinate
    # \param peg=None A semi-MANDATORY origins.PegPoint
    # \param ellipsoid=None A DUMMY, the
    # geo.detic.newton.ellipsoid.Ellipsoid is internal
    def __init__(self, coordinate1, coordinate2, coordinate3,
                 peg=None, ellipsoid=None):
        super(_Pegged, self).__init__(coordinate1,
                                      coordinate2,
                                      coordinate3,
                                      ellipsoid=ellipsoid)
        ## These classes NEED a peg, but it can be set after __init__
        # via the pegged() decorator, or with the new kwiter iterator.
        self._peg = peg

    ## Magnitude wrt to Origin.
    # \returns geo.metric.euclid.scalar.Scalar
    def __abs__(self):
        return abs(self.vector())

    ## \returns {'peg': self.peg}
    def kwiter(self):
        """{'peg': self.peg}"""
        return dict(peg=self.peg, ellipsoid=self.ellipsoid)

    ## PegPoint getter
    @property
    def peg(self):
        """peg point"""
        return self._peg

    ## PegPoint setter
    @peg.setter
    def peg(self, pp):
        """peg point"""
        self._peg = pp

    ## getitem is decorated with pegged, as kwiter is difficult to apply.
    @pegged
    def __getitem__(self, index):
        return super(_Pegged, self).__getitem__(index)

    ## Call super and add the peg after the fact (this is when to use
    # lazy instatiation, need to check LSP).
    # \param func The func to boradcast
    # \param *args func's arguments
    # \kwd *kwargs func's keywords
    # \retval _Coordinate func(self) as type(self)
    @pegged
    def broadcast(self, func, *args, **kwargs):
        """broadcast(func, *args, **kwargs) will broadcast
        to coordinate components"""
        return super(_Pegged, self).broadcast(func, *args, **kwargs)

    ## Use _Coordinate.change_ellipsoid() and call the result's method
    # to self's coordinates (including peg keyword)
    # \param target A geo.detic.ellipsoid.ellipsoid.Ellipsoid
    # \retval _Coordinate self referenced to the new ellipsoid
    def change_ellipsoid(self, target):
        """change ellipsoid and re-compute coordinates"""
        result = super(_Pegged, self).change_ellipsoid(target)
#        result = _Coordinate.change_ellipsoid(self, target)
        # this is an antipattern, as subs will fail on method name?
        method = operator.methodcaller(type(self).__name__.lower(),
                                       peg=self.peg)
        return method(result)

    ## str: call super and added a peg point.
    def __str__(self):
        return "{}{}".format(super(_Pegged, self).__str__(), str(self.peg))

    ## convert ENU at current peg point:
    # \retval LTP in ENU.
    def enu(self):
        """Convert to East North Up at Peg."""
        return self.ltp(peg=self.peg.enu())

    ## Parallel Transport to another Peg Point or ECEF
    def parallel_transport(self, peg=None):
        """v' = parallel_transport([peg=None])

        express non-affine vector (e.g., velocity) in:

        ECEF  or peg frame
        """
        # affine to ecef
        self_to_ecef = self.ellipsoid.ltp2ecef_affine(self.peg)
        # rotation part
        rotation = self_to_ecef.rotation
        if peg:  # kwd
            # affine to other (from ecef)
            ecef_to_other = self.ellipsoid.ecef2ltp_affine(
                peg
                )
            # update rotation
            rotation *= ecef_to_other.rotation
        # rotate vector's coordinates to target frame.
        return rotation(self.vector())

    ## Compute local radius of curvature along peg heading
    # \returns \f$ R(\phi, \psi) \f$
    def local_radius_of_curvature(self):
        """Local Radius of Curvature for peg (and ellipsoid)."""
        return self.ellipsoid.local_radius_of_curvature(self.peg[0],
                                                        self.peg[2])

    ## pre memoizer.
    _R = None

    ## Memoized.
    # \returns local_radius_of_curvature().
    @property
    def R(self):
        """R is the local radius of curvature (memoized)."""
        if self._R is None:  # memoizer
            self._R = self.local_radius_of_curvature()
        return self._R


## A "private" mixin for Cartesian _Coordinate
class _Cartesian(object):
    """A mix-in for cartesian coordinates"""

    __metaclass__ = abc.ABCMeta

    ## The default is always "x" "y" "z"
    coordinates = ("x", "y", "z")

    ## This method is trival, and allow polymorphic calls with _NonCartesian
    # instances
    def otf(self):
        """self.otf() --> self"""
        return self

    ## Move (non-affine) vectors from self's frame to other's frame (via
    # rotation2ecef()).
    # \param other A _Cartesian() coordinate
    # \return geo.metric.euler.charts.versor.Versor
    def rotator(self, other):
        """Usage:
        >>>vector' = self.rotatorr(other)(vector)

        That is, it returns a versor that projects tangent vector to
        another affine space's orientation."""

        return self.rotation2ecef() * (~(other.rotation2ecef()))


## A "private" mixin for non-cartesian _Coordinate
class _NonCartesian(object):
    """A mixin for noncartesion classes, brings in the "vector" method."""
    __metaclass__ = abc.ABCMeta

    ## Sub's need names for their coordinates.
    @abc.abstractproperty
    def coordinates(self):
        """subs define a coordinate names"""

    ## Noncartesian conjugate is always the orthognal_tangent_frame().
    def __invert__(self):
        """see otf()."""
        return self.otf()

    ## Can't rotate from non-cartesian, yet.
    # \throws geo.utils.exception.KoszulConnectionError
    def rotation2ecef(self):
        """Basis rotations are only defined for global coordinates.

        If you need this, you'll have to work it out using the
        stuff in christoffel.py"""
        from ..utils.exceptions import KoszulConnectionError as Err
        raise Err(
            "NonCartesian Frame {} do not support Global Lifting".format(
                type(self)
                )
            )


## <a href="http://en.wikipedia.org/wiki/ECEF">Earth
# Centered Earth Fixed</a> coordinates.
# @image html ECEF.svg.png
# @image latex ECEF.svg.png
class ECEF(_Cartesian, _Fixed, christoffel.ECEFMixin):
    #pylint: disable=R0904
    """ECEF(x, y, z)

    Earth Centered Earth Fixed Coordinates.

    The x axis goes from the center to (lat=0, lon=0)
    The y axis goes from the center to (lat=0, lon=90)
    The z axis goes from the center to (lat=90)

    Methods to tranaform coordinates are:

    ecef()
    llh()
    ltp(peg)
    sch(peg)

    Other methods are:

    vector()            convert to a Vector object

    The Coordinates are defined on the Ellipsoid:
    {}
    """
    ## ECEF --> ECEF is trivial
    ecef = _Cartesian.otf

    def __invert__(self):
        return self.llh()

    ## ECEF --> LLH  via f_ecef2lla().
    # \param peg=None is NOT used
    # \retval LLH self in LLH
    def llh(self):  #, peg=None):
        """ecef.llh() puts coordinates in LLH"""
        return self.ellipsoid.LLH(*(self.f_ecef2lla(*self.iter())))

    ## ECEF --> LTP via affine2ltp() .
    # \kwd peg origins.PegPoint for new representation
    # \retval LTP self in LTP at peg.
    def ltp(self, peg=None):
        """ecef.ltp(peg)  put coordinates in LTP at peg."""
        assert peg is not None
        return self.ellipsoid.LTP(
            *(self.affine2ltp(peg)(self.vector()).iter()),
            peg=peg)

    ## ECEF --> LTP --> SCH   (derived)
    # \kwd peg origins.PegPoint for new representation
    # \retval SCH self in SCH at peg.
    def sch(self, peg=None):
        """ecef.sch(peg) puts coordinates in SCH at peg."""
        assert peg is not None
        return self.ltp(peg).sch()

    ## This read-only attribute is a function, f,  that does:\n (lat, lon, hgt)
    # = f(x, y, z) \n for a fixed Ellipsoid using
    # ellipsoid.Ellipsoid.XYZ2LatLonHgt .
    @property
    def f_ecef2lla(self):
        """Read-only function--> you access this, and the ellipsoid
        'compiles' a function in ellipsoid.XYZ2LatLonHgt."""
        return self.ellipsoid.XYZ2LatLonHgt

    ## This method returns the euclid::affine::Affine transformation from the
    # ECEF frame to LTP at a::PegPoint using
    # ellipsoid.Ellipsoid.ecef2ltp_affine
    def affine2ltp(self, peg):
        """given a peg, return an affine transform to the LTP of that peg, see:
        ellipsoid.ecef2ltp_affine."""
        try:
            # ellipsoid builds the function
            result = self.ellipsoid.ecef2ltp_affine(peg[0],
                                                    peg[1],
                                                    peg[2])
        except AttributeError as err:
            # Dx the exception
            if peg is None:  # raise
                msg = ("Attempt a coordinate conversion to an affine" +
                       "space with NO ORIGIN: peg is None")
                raise AffineSpaceError(msg)
            raise err
        return result

    ## We'll do it my way.
    # \param other A Coordinate
    # \retval ecef other iin ECEF
    def _like_me(self, other):
        """make other in my coordinates"""
        return other.ecef()

    ## Rotateing rotation to ECEF (Idenity)
    # \return geo.metric.euler.charts.versor.Versor
    @staticmethod
    def rotation2ecef():
        """Returns a rotation that rotate local tangent vector to ECEF"""
        from ..metric.euclid.tensor import DELTA
        return DELTA.versor()

    ## Coriouls for Calculator:
    # \param velocity \f$  \vec{v} \f$
    # \kwd omega [=None] \f$ \vec{\omega} \f$,
    # see ::Ellipsoid.angular_velocity()
    # \return acceleration \f$ -2 \vec{\Omega} \times \vec{v} \f$
    # @image html crls1.gif
    # @image latex crls1.gif
    def coriolus(self, velocity, omega=None):
        """Earth Only? Coriolus force for velocity ar point."""
        if omega is None:  # kwd
            omega = self.ellipsoid.angular_velocity()
        return -2 * (omega ^ velocity)


## <a href="http://en.wikipedia.org/wiki/Geodetic_coordinates#Coordinates">
# Geodetic Coordinates</a>: <a href="http://en.wikipedia.org/wiki/Latitude">
# Latitue</a>, <a href="http://en.wikipedia.org/wiki/Longitude">Longitude</a>,
# <a href="http://en.wikipedia.org/wiki/Elevation">Height</a>/
# @image html longlat2.gif
# @image latex longlat2.gif
class LLH(_NonCartesian, _Fixed, christoffel.LLHMixin):
    """point = LLH(lat, lon, hgt) or
    points = LLH(lat, lon, hgt)

    Geodetic Coordinates: geodetic latitude.

    lat-itude in degrees (NEVER RADIANS, EVER)
    lon-gitude in degrees (NEVER RADIANS, EVER)
    hgt elevation, height, in meters


    Methods are:

    point.ecef()
    point.llh()
    point.ltp(peg=peg)
    point.sch(peg=peg)

    The Coordinates are defined on the Ellipsoid:
    {}
    """

    ## Geodetic coordinate names
    coordinates = ("lat", "lon", "hgt")

    ## Units are as is
    UNITS = ("deg", "deg", "m")

    ## LLH --> ECEF via _f_lla2ecef()
    # \returns instance in ECEF
    def ecef(self):
        """llh.ecef() puts coordinates in ECEF."""
        return self.ellipsoid.ECEF(*self._f_lla2ecef(*(self.iter())))

    otf = ecef

    ## Trivial
    llh = arraylike.ABCNumpy.__pos__

    ## LLH --> ECEF --> LTP
    # \param peg origins.PegPoint for target coordinates
    # \returns instance in LTP at peg.
    def ltp(self, peg=None):
        """llh.ltp(peg) puts coordinates in LTP at peg."""
        assert peg is not None
        return self.ecef().ltp(peg)

    ## LLH --> ECEF --> LTP --> SCH
    # \kwd peg origins.PegPoint for target coordinates
    # \returns instance in SCH at peg.
    def sch(self, peg=None):
        """llh.sch(peg) puts coordinates in SCH at peg."""
        assert peg is not None
        return self.ecef().sch(peg)

    ## This read-only attribute is a function, f,  that does:\n (x, y, z) =
    # f(lat, lon, hgt) \n using ellipsoid.Ellipsoid.LatLonHgt2XYZ
    @property
    def _f_lla2ecef(self):
        """Read-only function--> ellipsoid.LatLonHgt2XYZ."""
        return self.ellipsoid.LatLonHgt2XYZ

    ## h is hgt
    @property
    def h(self):
        """lazy: h is hgt"""
        return self.hgt

    ## hardocde conversion to radians.
    @staticmethod
    def radians(x):
        """degs--> radians"""
        return 0.017453292519943295*x

    ## \f$ {\bf n}^e = \left[\begin{array}{c} \cos{\phi}\cos{\lambda} \\
    # \cos{\phi}\sin{\lambda} \\ \sin{\phi}  \end{array} \right] \f$ \n is
    # the <a href="http://en.wikipedia.org/wiki/N-vector">N-vector</a>.
    def n_vector(self):
        """Compute a Vector instance representing the N-Vector"""
        from ..utils.trig import sind, cosd
        from ..metric.euclid import vector
        return vector.Vector(cosd(self.lat)*cosd(self.lon),
                             cosd(self.lat)*sind(self.lon),
                             sind(self.lat))

    ## Make a peg (for a signelton)
    # \param heading Is a heading in degrees
    # \retval peg origins.PegPoint here pointing along heading
    def make_peg_here(self, heading):
        """make a peg point at this coordinate (SINGULAR)"""
        return origins.PegPoint(self.lat, self.lon, heading)

    ## \retval peg origins.PegPoint here pointing along origins.EAST
    def make_enu_peg_here(self):
        """make an ENU peg point at this coordinate (SINGULAR)"""
        return self.make_peg_here(origins.EAST)

    ## We'll my way.
    # \param other A Coordinate
    # \retval ecef other in LLH
    def _like_me(self, other):
        """put other in my coordinates"""
        return (other | self.ellipsoid).llh()


## Local Tangent Plane Cartesian coordinates
# @image html enu.png
# @image latex enu.png
class LTP(_Cartesian, _Pegged, christoffel.LTPMixin):
    #pylint: disable=R0904
    """point(s) = LTP(x, y, z, peg=<peg>)

    A point or points in a cartesian versions of SCH coordinates\n

    ARGS:
    ____
    x   along track cartesian coordinate in meters
    y   cross track cartesian coordinate in meters
    z   height above the ellipsoid in meters, at origin.


    KWARGS:
    ______
    peg      A PegPoint instance defining the coordinate system.\n


    Methods are:

    point.ecef()
    point.llh()
    point.ltp([peg=None])
    point.sch([peg=None])

    The Coordinates are defined on the Ellipsoid:
    {}
    """
    def __invert__(self):
        return self.sch()

    ## LTP --> ECEF via ellipsoid.Ellipsoid.ltp2ecef_affine .
    # \returns instance in ECEF.
    def ecef(self):
        """ltp.ecef([peg=None]) put it into ECEF.
        Keyword is ignored."""
        return self.ellipsoid.ECEF(*(self.affine2ecef(self.vector()).iter()))

    ## LTP --> ECEF --> LLH
    # \returns instance in LLH.
    def llh(self):
        """ltp.llh([peg=None]) put it into LLH
        Keyword is ignored."""
        return self.ecef().llh()

    ## Trivial OR LTP --> ECEF --> LTP'
    # \kwd peg=None None: use peg or self.peg
    # \returns instance in LTP at peg.
    def ltp(self, peg=None):
        """ltp.ltp(peg) transforms to a new peg"""
        return self if peg is None else self.ecef().ltp(peg)  # kwd

    ## LTP --> SCH via ellipsoid.Ellipsoid.TangentPlane2TangentSphere
    # OR LTP --> LTP' --> SCH
    # \kwd peg=None None: use peg or self.peg
    # \returns instance in SCH at peg.
    def sch(self, peg=None):
        """ltp.ltp(peg=None) transforms to SCH, possibly in a
        different peg"""
        return (self.ellipsoid.SCH(*self.to_sch_tuple, peg=self.peg)
                if peg is None else  # kwd
                self.ltp(peg).sch(None))

    ## Memoizer initialization
    _f_ltp2sch = None

    ## This read-only memoized attribute is a function, f,  that does:
    # (s, c, h) = f(x, y, z) \n for this coordinate's origin.PegPoint
    # \returns \f$ _{\rm sch}f_{\rm xyz}(x, y, z; R=R) \f$, with R
    # partially evaluated (once) at Peg.
    @property
    def f_ltp2sch(self):
        """Read-only function--> you access this, and the ellipsoid
        'compiles' a function in ellipsoid.TangentPlane2TangentSphere"""
        if self._f_ltp2sch is None:  # memoizer
            self._f_ltp2sch = self.ellipsoid.TangentPlane2TangentSphere(
                *self.peg)
        return self._f_ltp2sch

    @property
    def to_sch_tuple(self):
        """(s, c, h) tuple, via: f_ltp2sch()"""
        return self.f_ltp2sch(*self.iter())

    ## Memoizer initialization
    _affine2ecef = None

    ## return the euclid::affine::Affine transformation to ECEF coordinates
    @property
    def affine2ecef(self):
        """Let ellipsoid get the affine transform via its
        ltp2ecef_affine
        method."""
        if self._affine2ecef is None:  # memoizer
            self._affine2ecef = self.ellipsoid.ltp2ecef_affine(
                *self.peg)
        return self._affine2ecef

    ## Get the rotation between basis triplets from here to ECEF.
    def rotation2ecef(self):
        """Rotation to ECEF's global basis oritentation."""
        return self.affine2ecef.rotation

    ## \param other A Coordinate
    # \retval ecef other in LTP
    def _like_me(self, other):
        """put other in my coodinates"""
        return (other | self.ellipsoid).ltp(peg=self.peg)

    ## Rotate to ECEF, compute, and rotate back to LTP.
    # \returns \f$ \vec{a}_c = -2\vec{\omega} \times \vec{v} \f$
    def coriolus(self, velocity):
        """Coriolus acceleration for velcoity
        (all in LTP tangent vector space)"""
        rot = self.rotation2ecef()
        return (~rot)(self.ecef().coriolus(rot(velocity)))


## Local (2nd Order) Tangent Sphere coordinates
# @image html sch.jpeg
# @image latex sch.jpeg
class SCH(_NonCartesian, _Pegged, christoffel.SCHMixin):
    #pylint: disable=R0904
    """point.SCH(s, c, h, peg=<peg>)

    A point or points in the SCH coordinate system:

    ARGS:
    ____
    s   along track polar coordinate in meters
    c   cross track polar coordinate in meters
    h   height above the ellipsoid in meters (at origin).


    KWARGS:
    ______
    peg      A PegPoint instance defining the coordinate system.


    Methods are:

    point.ecef()
    point.llh()
    point.ltp([peg=None])
    point.sch([peg=None])


    Bonus Methods:
    ==============
    shat(), chat(), hhat() give unit vectors in the OTP--and local_frame
    combines them into an SO(3) matrix.


    The Coordinates are defined on the Ellipsoid:
    {}
    """

    ## These non-cartesian coordinates are called S, C and H.
    coordinates = ("s", "c", "h")

    ## SCH --> LTP --> ECEF
    # \returns instance in ECEF coordiantes
    def ecef(self):
        """ecef([peg=None]) converts to (same) LTP, and then
        to ECEF.  Keyword argument is ignored."""
        return self.otf().ecef()

    ## SCH --> LTP --> LLH
    # \returns instance in LLH coordiantes
    def llh(self):
        """llh([peg=None]) converts to (same) LTP, and then
        to LLH.  Keyword argument is ignored."""
        return self.otf().llh()

    ## SCH --> LTP using _to_ltp_tuple() If peg is None, otherwise:
    # SCH --> ECEF --> LTP'
    # \kwd peg=None None: use peg or self.peg
    # \returns instance in LTP at peg.
    def ltp(self, peg=None):
        """ltp([peg=None]) will convert to LTP coordinates
        with the same peg (If None), or to LTP in peg
        """
        return (
            self.ellipsoid.LTP(*self._to_ltp_tuple, peg=self.peg)
            if peg is None else  # kwd
            self.ecef().ltp(peg=peg)
            )

    ## Trivial If peg is None, otherwise: SCH --> ECEF --> SCH'
    # \kwd peg=None None: use peg or self.peg
    # \returns instance in SCH at peg.
    def sch(self, peg=None):
        """sch([peg=None]) is trivial If peg is None,
        otherwise:
        convert to ECEF and then to SCH at the new peg."""
        # self.ecef.sch(peg or self.peg) # 3000 times slower.
        return self if peg is None else self.ecef().sch(peg=peg)  # kwd

    otf = ltp

    ## Memoizer initialization
    _f_sch2ltp = None

    ## This read-only memoized attribute is a function, f, that does:
    # (x, y, z) = f(s, c, h) \n for this coordinate's origin.PegPoint
    # \returns \f$ _{\rm xyz}f_{\rm sch}(s, c, h; R=R) \f$, with R
    # partially evaluated (once) at Peg.
    @property
    def f_sch2ltp(self):
        """read-only property returns a 'compiled' function from
        ellipsoid.TangentSphere2TangentPlane"""
        if self._f_sch2ltp is None:  # memoizer
            self._f_sch2ltp = (
                self.ellipsoid.TangentSphere2TangentPlane(*self.peg)
            )
        return self._f_sch2ltp

    ## This read-only attribute is  a tuple of (x, y, z) computed from
    # f_sch2ltp.
    @property
    def _to_ltp_tuple(self):
        """(x, y, z) tuple via _f_sch2ltp()"""
        return self.f_sch2ltp(*self.iter())

    ## Put another coordinate in my frame (and ellipsoid).
    # \param other A Coordinate
    # \retval ecef other in SCH
    def _like_me(self, other):
        """put other in my coodinates"""
        return (other | self.ellipsoid).sch(peg=self.peg)


## Special case is a subclass- Liskov safe?
class ENU(LTP):
    """The circle ellipse problem with ENU and LTP."""

    # experimental class -not quire right
    def __init__(self, *args, **kwargs):
        super(ENU, self).__init__(*args, **kwargs)
        ## peg must point to origins.EAST
        self.peg.hdg = origins.EAST

    ## You need this for change_ellipsoid's methodcaller.
    enu = LTP.ltp


## Available Coordinates
COORDINATES = (ECEF, LLH, LTP, SCH)
