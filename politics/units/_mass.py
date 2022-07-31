"""mass"""
## \namespace geo.politics.units._mass
# <a href="http://en.wikipedia.org/wiki/Mass">Mass</a>
from ._bases import *
from .length import m2ft

__all__ = ('NperKg2slug', 'lbs2kg', 'kg2lbs', 'M_solar', 'M_earth', 'M_lunar', 'atmosphere_mass', 'ocean_mass', 'human_mass', 'planck_m', 'BMI', 'BMI_Dx')


class Mass(Unit):
    """Convert from/to Kilograms:

    Kilograms = 1 * unit  (convert 1 unit to Kilograms)
    units = 1 / unit   (convert 1 Kilogram to units)
    units = unit(1)    (convert 1 Kilogram to units)
    """


##<a href="http://en.wikipedia.org/wiki/Slug_(mass)">Slug</a>.
NperKg2slug = SI(float(m2ft), 'meters_per_second**2', 'slug')

## <a href="http://en.wikipedia.org/wiki/Pound_(mass)">Pound</a> It.
lbs2kg = SI(constants.pound, 'lbs', 'kg')
kg2lbs = ~lbs2kg


stone2lbs = SI(14, 'st', 'lbs')
quarter2stone = SI(2, 'quarter', 'st')

## <a href="http://en.wikipedia.org/wiki/Solar_mass">Solar Mass</a>
M_solar = 1.98892e30

## <a href="http://en.wikipedia.org/wiki/Earth_mass">Earth Mass</a>
M_earth = 5.9722e24

## Lunar Mass.
M_lunar = 7.3477e22


atmosphere_mass =(5.27e18)
ocean_mass = 1.42e21
human_mass = (3.e11)


planck_m = Mass(constants.hbar * constants.c / constants.G) ** 0.5


## <a href="http://en.wikipedia.org/wiki/Body_mass_index">Body Mass
# Index in AMERICAN Units</a>\n
# \f$ BMI \propto  m/h^2 \f$.
def BMI(pounds, feet, inch=0., Dx=False):
    """BMI(pounds, feet {, inch=0, Dx=False})

    if Dx, then pipe it BMI_Dx()
    """
    from .length import ft2m
    result = lbs2kg(pounds)/ ft2m(feet + inch/12.)**2

    if Dx:
        result = BMI_Dx(result)
    return result


## A function that makes a dictionary of {iz.Intervals: "Dx strings"}
def BMI_Dx(bmi=None):
    """BMI_Dx(bmi=None):

    return dictionary of weight intervals: BMI Dx's if bmi is None, else:
    """
    try:
        from my.math import iz
    except ImportError:
        raise NotImplementedError("Requires: belz/python::my.math.iz")
    Dx = {
        iz.LeftClosedRightOpenInterval(0, 15): "Very severly underweight",
        iz.LeftOpenRightClosedInterval(15, 16): "Severly underweight",
        iz.LeftOpenRightClosedInterval(16, 18.5): "Underweight",
        iz.LeftOpenRightClosedInterval(18.5, 25): "Normal (healthy weight)",
        iz.LeftOpenRightClosedInterval(25, 30): "Overweight",
        iz.LeftOpenRightClosedInterval(30, 35): "Moderately obese",
        iz.LeftOpenRightClosedInterval(35, 40): "Severley obese",
        iz.LeftOpenRightOpenInterval(40, iz.inf): "Very Severley obese",
        iz.LeftOpenRightClosedInterval(iz.inf, iz.inf): "Singularly obese",
        }

    return (
        Dx if bmi is None else
        filter(lambda a: bmi in a[0], Dx.iteritems())[0][-1]
        )
