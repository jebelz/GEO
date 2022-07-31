"""Acceleration"""
## \namespace geo.politics.units._acceleration Acceleration
from ._bases import *
from .times import hour2second
from .length import mile2ft, ft2m, nauticalmile2meter


__all__ = ('NperKg2g', 'NperKg2gal', 'g', 'ftpersqsec2g')

## ( <a href="http://en.wikipedia.org/wiki/Standard_gravity">g</a>
NperKg2g = SI(1./constants.g, 'acceleration', """g's""")


## Guys and <a href="http://en.wikipedia.org/wiki/Gal_(unit)">Gals</a>.
NperKg2gal = SI(constants.hecto, 'NperKg', 'Gal')


ftpersqsec2g = NotImplemented

g = (~NperKg2g).unit()
