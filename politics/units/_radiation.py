"""Radiation"""
## \namespace geo.politics.units._radiation
#<a href="http://en.wikipedia.org/wiki/Radiation">Radiation</a>.
from ._bases import *

__all__ = ('Ci2Bq', 'Gy2rad', 'Sv2rem', 'BED2Sv')

## <a href="http://en.wikipedia.org/wiki/Curie">Curiousity killed  Marie</a>.
Ci2Bq = SI(3.7e10, 'Ci', 'Bq')

##  <a href="http://en.wikipedia.org/wiki/Absorbed_dose">Absorbed Dose</a>
# and the <a href="http://en.wikipedia.org/wiki/Rad_(unit)">rad</a>.
Gy2rad = SI(100., 'Gy', 'JperKg')

## <a href="http://en.wikipedia.org/wiki/Roentgen_equivalent_man">REM</a>.
Sv2rem = SI(100., 'Sv', 'rem')

## <a href="http://en.wikipedia.org/wiki/Banana_equivalent_dose">BED</a>
BED2Sv = SI(36.e-6/365., 'BED', 'Sv')
