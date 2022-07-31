"""Volume"""
## \namespace geo.politics.units._volume
# <a href="http://en.wikipedia.org/wiki/Volume">Volume</a>.
from ._bases import *
from .length import yd2ft, inch2ft, ft2m
from .area import acre2sqft

cuyd2cuft = yd2ft**3

acreft2cuft = SI(float(acre2sqft), 'acre*ft', 'ft**3')

cubicinch2cc = (inch2ft * ft2m/constants.centi)**3
cc2cubicinch = ~cubicinch2cc
cubicinch2liter = cubicinch2cc/constants.kilo

minum2drop = SI(1, "minum", "drop")
dram2minum = SI(60, "dram", "minum")
tsp2min = SI(80, "tsp", "minum")
Tbsp2tsp = SI(3, "Tbsp", "tsp")
oz2Tbsp = SI(2, "oz", "Tbsp")
jig2Tbsp = SI(3, "jig", "Tbsp")
jig2oz = (jig2Tbsp * (~oz2Tbsp))
gil2oz = SI(4, "gil", "oz")
cup2gil = SI(2, "cup", "gil")
cup2oz = (cup2gil * gil2oz)
pint2cup = SI(2, "pint", "cup")
quart2pint = SI(2, "quart", "pint")
gallon2quart = SI(4, "gal", "quart")
barrel2gallon = SI(31.5, "barrel", "gal")
bbl2gallon = SI(42, "bbl", "gal")
hogshead2gal = SI(63, "hogshead", "gal")
hogshead2cuft = SI(8.421875, "hogshead", "cuft")
cuft2gal = SI(fractions.Fraction(576, 77))


peck2gallon = SI(2, "peck", "gal")
bushel2peck = SI(4, "bu", "peck")
