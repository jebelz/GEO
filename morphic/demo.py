"""Module demonstates usage of stuff"""
import numpy as np
import geo

N = 5

# 1 arc second spacing (latitude)
lat1 = np.linspace(33, 33 + N/3600., N)

# 1 arc second spacing (longitude)
lon1 = np.linspace(-118, -118 + N/3600., N)


## 2d grids of lat, lon 
lon, lat = np.meshgrid(lon1, lat1)

## make up a 0 dem
dem = 0*lat

## Put in an affine space: Geodetic Coordinates:
points = geo.LLH(lat, lon, dem)


peg = geo.Peg(lat1[N//2], lon1[N//2], 90)

## 
enu = points.sch(peg)

# now add a known dem after the fact:
enu.h = enu.s + 2 * enu.c

from geo.morphic import monge


patch = monge.MongePatch(enu.s, enu.c, enu.h)



q = patch.quadric()

print q.dip(), q.aspect()

