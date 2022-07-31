"""Pressure"""
## \namespace geo.politics.units._pressure
# <a href="http://en.wikipedia.org/wiki/Pressure">Pressure</a>
from ._bases import *

__all__ = ('bar2Pa', 'Pa2bar', 'atm2Pa', 'Pa2atm', 'bar2atm', 'atm2bar', 'torr2Pa', 'Pa2torr', 'Pa2mb', 'mb2Pa', 'Ba2Pa', 'Pa2Ba', 'pz2Pa', 'Pa2pz', 'atm')

## <a href="http://en.wikipedia.org/wiki/Bar_(unit)">Bar/a>
bar2Pa = SI(constants.bar, 'bar', 'Pa')
Pa2bar = ~bar2Pa

## <a href="http://en.wikipedia.org/wiki/Atmosphere_(unit)">Atmopshere</a>
atm2Pa = SI(constants.atm, 'atm', 'Pa')
Pa2atm = ~atm2Pa

bar2atm = bar2Pa*Pa2atm
atm2bar = ~bar2atm
## <a href="http://en.wikipedia.org/wiki/Torr_(unit)">Torr</a>.
torr2Pa = SI(constants.torr, 'torr', 'Pa')
Pa2torr = ~torr2Pa

Pa2mb = (Pa2bar*SI(1/constants.milli, From="bar", To="mbar"))
mb2Pa = ~Pa2mb


## <a href="http://en.wikipedia.org/wiki/Barye">Barye</a> to Pa.
Ba2Pa = SI(constants.kilo*constants.centi**2)
Pa2Ba = ~Ba2Pa

## <a href="http://en.wikipedia.org/wiki/Pieze">Pieze</a> to Pa
pz2Pa = SI(constants.kilo)
Pa2pz = ~pz2Pa

atm = atm2Pa.unit()
