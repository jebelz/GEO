"""The wgs84 module is for beginners-- and that is not a put down.

The 4 coordinate classes (ECEF, LLH, LTP, SCH) provide objects for geodesy,
based, of course, on the WGS84 ellipsoid.

The 16 module functions <in>2<out>:

ecef2ecef, ..., sch2sch2

are implemented old school (almost): They take a triple input (with possible
peg or ellipsoid keyword) and return... a triplet? No -- that would be
old school-- they return the triplet wrapped in the appropriate class.

Nevertheless, should you require purely procedural functionality, just call the
object's tolist() method. Thus, for example, given a triplet in an SCH
coordinate system (at 'peg'):

(s, c, h)

you would get the local tangent plane coordinates via:

>>>x, y, z = sch2ltp(s, c, h, peg=peg)

Of course, (s, c, h) could be numbers or numpy arrays-- (x, y, z) will be
likewise.

Finally: there is a distance function:

>>>distance(lat1, lon1, lat2, lon2)

"""
## \namespace geo.detic.newton.wgs84.wgs84 The
## <a href="http://en.wikipedia.org/wiki/World_Geodetic_System">
## World Geodetic System 1984</a>.

import operator
from ..almanac import WGS84

__all__ = ('WGS84', 'distance', 'ECEF', 'LLH', 'LTP', 'SCH',
           'ecef2ecef', 'ecef2llh', 'ecef2ltp', 'ecef2sch',
           'llh2ecef', 'llh2llh', 'llh2ltp', 'llh2sch',
           'ltp2ecef', 'ltp2llh', 'ltp2ltp', 'ltp2sch',
           'sch2ecef', 'sch2llh', 'sch2ltp', 'sch2sch')


## The distance function on WGS84
distance = operator.attrgetter(WGS84)


ECEF = WGS84.ECEF
LLH = WGS84.LLH
LTP = WGS84.LTP
SCH = WGS84.SCH

ecef2ecef = WGS84.ecef2ecef
ecef2llh = WGS84.ecef2llh
ecef2ltp = WGS84.ecef2ltp
ecef2sch = WGS84.ecef2sch

llh2ecef = WGS84.llh2ecef
llh2llh = WGS84.llh2llh
llh2ltp = WGS84.llh2ltp
llh2sch = WGS84.llh2sch

ltp2ecef = WGS84.ltp2ecef
ltp2llh = WGS84.ltp2llh
ltp2ltp = WGS84.ltp2ltp
ltp2sch = WGS84.ltp2sch

sch2ecef = WGS84.sch2ecef
sch2llh = WGS84.sch2llh
sch2ltp = WGS84.sch2ltp
sch2sch = WGS84.sch2sch
