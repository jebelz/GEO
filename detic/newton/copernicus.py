"""Ellipsoids from the solar system. The do need angular velocities."""
## \namespace geo.detic.newton.copernicus Planets, and other biaxial
# ellipsoids in the solar system.
from .ellipsoid import Ellipsoid

HOUR = 3600
DAY = 24 * HOUR


## Neither function nor decorator, this makes Earth ellipsoids work
# for celestiral bodies.
def planetoid(a, b=None, model="", day=0):
    result = Ellipsoid(a, b or a, model=model)
    # This is a major violation of all python sensibilites, as we are
    # mutating and instance attribute from a way over here--thus is the
    # nature of adding planets, and it's all for one method that is
    # never called in the package, but users may want.
    result.day = day
    return result


## <a href="http://www.jpl.nasa.gov/spaceimages/search_grid.php?category=sun">
# The Sun </a>
# @image html sun.jpg
# @image latex sun.jpg
SUN = Ellipsoid(696342e3, 111111.11111111111111, model="Sun")

## <a href="https://pds.jpl.nasa.gov/planets/choices/mercury1.htm">Heremes </a>
# @image html p1.jpeg
# @image latex p1.jpeg
MERCURY = planetoid(2439.7e3, model="Mercury", day=58.646*DAY)

## <a href="https://pds.jpl.nasa.gov/planets/choices/venus1.htm">Aphrodite</a>
# @image html p2.jpg
# @image latex p2.jpg
VENUS = planetoid(6051.8e3, model="Venus", day=-243.025*DAY)

## <a href="http://photojournal.jpl.nasa.gov/target/Moon">The Moon</a>
# @image html p3.jpeg
# @image latex p3.jpeg
MOON = planetoid(1738.14e3, 1735.97e3, model="Moon")

## <a href="https://pds.jpl.nasa.gov/planets/choices/mars1.htm">Aries</a>
# @image html p4.jpeg
# @image latex p4.jpeg
MARS = planetoid(3396.2e3, 3376.2e3, model="Mars", day=1.025957*DAY)

## <a href="http://dawn.jpl.nasa.gov">Demeter</a> Not a planet.
# @image html pc.jpg
# @image latex pc.jpg
CERES = planetoid(481.5e3, 455.5e3, model="Ceres", day=0.3781*DAY)

## <a href="https://pds.jpl.nasa.gov/planets/choices/jupiter1.htm">Zeus</a>
# @image html p5.jpeg
# @image latex p5.jpeg
JUPITER = planetoid(71492.e3, 66854.e3, model="Jupiter", day=9.925*HOUR)

## <a href="https://pds.jpl.nasa.gov/planets/choices/saturn1.htm">Cronos</a>
# @image html p6.jpeg
# @image latex p6.jpeg
SATURN = planetoid(60268.e3, 54364.e3, model="Saturn", day=10.55*HOUR)

## <a href="http://saturn.jpl.nasa.gov/science/index.cfm?SciencePageID=73">
# Titan</a>
# @image html p6a.jpeg
# @image latex p6a.jpeg
TITAN = planetoid(2575e3, model="Titan", day=15.945*DAY)

## <a href="https://pds.jpl.nasa.gov/planets/choices/neptune1.htm">Caelus</a>
# @image html p7.jpeg
# @image latex p7.jpeg
URANUS = planetoid(25559e3, 24973e3, model="Uranus", day=-0.71833*DAY)

## <a href="https://pds.jpl.nasa.gov/planets/choices/neptune1.htm">Poseidon</a>
# @image html p8.jpeg
# @image latex p8.jpeg
NEPTUNE = planetoid(24764e3, 24341e3, model="Neptune", day=0.6713*DAY)

## <a href="https://pds.jpl.nasa.gov/planets/choices/pluto1.htm">Dis Pater</a>
# a classic Planet.
# @image html p9.jpg
# @image latex p9.jpg
PLUTO = planetoid(1188.3e3, model="Pluto", day=6.387230*DAY)

## <a href="https://pds.jpl.nasa.gov/planets/choices/charon1.htm">Hades</a>
# a classic Planet.
# @image html charon.jpg
# @image latex charon.jpg
CHARON = planetoid(603.5e3, model="Charon", day=PLUTO.day)


PSR_J1748 = planetoid(17e3, model="Pulsar", day=0.00139595482)
