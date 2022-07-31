"""The ellispoid sub-package. The only concrete class is the Ellipsoid:
it is HUGE. This cannot be avoided, as it is and does so much.

It is a spheroid (_spheroid.py): while it only required 2 parameters,
there are so many combinations that can define it.  It also has functions
to convert bewtween the common lattitude and the myriad of others.

It is a place where datum have meaning (_datum.py). This provides the
functions to convert bewtween coordinate systems.
There are combined in ellipsoid.py, where the Ellipsoid class makes
abstract coordinates concrete.

The ellipsoids of yore are stored in alamanac.py, while the WGS84 ellipsoid
has its own sub-package  (wgs84/).

Finally, should you work on Mars or elsewhere, planets.py has some
well-known ellipsoids from around the Solar System.


PS: It was Newton (Sir Issac) who proved a rotating, self-gravitating,
body takes an oblate ellipsoidal shape [Reference 0: Principa].
"""

## \namespace geo.detic.newton  Ellipsoid Abstractions and Implementations
from .ellipsoid import *
from .almanac import *
from .wgs84 import *

## Density of the Earth (kg/m$^3$)
RHO = 5510.

## Mass of the Earth (kg)
M_E = 5.9722e24

## Uncertainty of M_E (kg)
dM_E = 0.0006e24
