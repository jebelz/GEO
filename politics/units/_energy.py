"""Energy"""
## \namespace geo.politics.units._energy
# <a href="http://en.wikipedia.org/wiki/Energy">Energy</a>.
from ._bases import *

__all__ = ('kWhours2J', 'J2kWhours', 'BTU2J', 'J2BTU', 'eV2J', 'J2eV',
           'foe2J', 'photons_per_joule', 'foe', 'eV', 'Ry', 'planck_E')


class Energy(Unit):
    """Convert from/to Joules:

    Joules = 1 * unit  (convert 1 unit to Joules)
    units = 1 / unit   (convert 1 Joule to units)
    units = unit(1)    (convert 1 Joule to units)
    """

## Energy
kWhours2J = SI(constants.hour*constants.kilo, 'kW-hour', 'J')
J2kWhours = ~kWhours2J
BTU2J = SI(1055., "BTU", "J")
J2BTU = ~BTU2J


## <a href="http://en.wikipedia.org/wiki/Electronvolt">eV</a>
eV2J = SI(constants.e, 'J', 'eV')

## <a href="http://en.wikipedia.org/wiki/Joule">Joule to eV</a>
J2eV = ~eV2J

## <a href="http://en.wikipedia.org/wiki/Foe_(unit)">\f$10^{51} \f$ ergs.</a>
foe2J = SI(1e44, 'foe', 'J')


## \f$ n = (h\nu N_A)^{-1} \f$ \n Moles per Joule at frequency \f$nu\f$.
def photons_per_joule(nu):
    """Get photons (moles) per joule at frequency (cps)"""
    return 1/(constants.h*nu)/ constants.N_A

foe = foe2J.unit(Energy)
eV = eV2J.unit(Energy)

## <a href="http://en.wikipedia.org/wiki/Rydberg_constant">Rydberg Energy</a>.
Ry = Energy(constants.Rydberg * constants.speed_of_light * constants.h)


planck_E = Energy(constants.hbar * constants.c**5 / constants.G) ** 0.5




