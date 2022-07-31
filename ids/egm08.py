#!/usr/bin/env python

"""EGM2008 1 or 2.5 minute posting bilinear interpolator:

geoid(lat, lon)
geoid[lat1: lat2: dlat, lon1:lon2:dlon]

egm08 lat1 lat2 ... latN lonN > height

See module global POSTING for posting resolution. All heights are relative
to the WGS-84 ellipsoid.
"""
## \namespace geo.ids.egm08 File dependent
# <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/">EGM08</a>
# linear interpolator.

#pylint: disable=E1102
import os

import numpy as np


## Full path to data file from:\n
#  http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
#  You must pick the file with the posting you want (1 or 2.5 minutes),
#  and the endianness of your machine.
EGM_FILE = os.path.join(os.getenv('HOME'), 'utils', 'egm08.dem')

## Posting, in arc-minutes (1.0, or 2.5).
POSTING = 1.0

## Posting, in degrees.
DELTA = POSTING/60.

## Number of longitude postings
NLON = 2 + int(360/DELTA)

## Number of latitude postings
NLAT = 1 + int(180/DELTA)

## Data shape
SHAPE = (NLAT, NLON)


## Latitude post --> index + fractional part
lat2i = lambda x: (90-x)/DELTA

## Longitude post --> index + fractional part
lon2j = lambda x: (x % 360)/DELTA


## covert (lat, lon) to nearest (i, j)
def latlon2index(lat, lon):
    """latlin2index(lat, lon) --> (i, j) tuple of indices"""
    return ([int(lat2i(lat)+0.5)], [int(lon2j(lon)+0.5)])


## Nearest Neighbor posting
def nearest_neighbor(lat, lon):
    """geoid(lat, lon)

    returns a memmap.

    Use float(geoid(lat, lon) for singletons
        array(....) for array input
    """
    try:
        return egm_map[latlon2index(lat, lon)]
    except IndexError as err:
        if lat < -90:  # guard
            raise ValueError("Latitude is not valid: %s" % str(lat))
        raise err


## Bilinear interpolator around a single lat, lon point.
def bilinear(lat, lon):
    """bilinear(lat, lon) for singleton lat, lon pairs"""
    r_i = lat2i(lat)
    r_j = lon2j(lon)

    i = int(r_i)
    j = int(r_j)

    fQ = egm_map[i:i+2, j:j+2]

    x1 = r_i - i
    y1 = r_j - j

    x2 = POSTING - x1
    y2 = POSTING - y1

    try:
        f11, f12, f21, f22 = map(float, fQ.ravel())
    except ValueError:
        ## To close to edge to work
        if lat == -90:  # semi guard
            lon = 90.
        return float(nearest_neighbor(lat, lon))

    return (f11 * x2 * y2 +
            f21 * x1 * y2 +
            f12 * x2 * y1 +
            f22 * x1 * y1)/POSTING**2


## Convert slice object to arange arrays
# \param slice_ A slice object
# \return array spanning slice's start, stop, stridding by step.
def slice2array(slice_):
    """Convert slice to array"""
    return np.arange(slice_.start, slice_.stop, slice_.step or 1)


## Convert get item slice objects in 1-d arrays, into 2 arrays, and
#  pass to geoid() as arguments
def getitem(index):
    """getitem(index) take a __getitem__ style argument,
    and then makes index-like floats with which to interpolate geoid."""
    lat, lon = map(slice2array, index)
    lat2 = lat.repeat(len(lon)).reshape(len(lat), len(lon))
    lon2 = lon.repeat(len(lat)).reshape(len(lon), len(lat)).T
    return geoid(lat2, lon2)


## A 1-off class with call, getitem, and a docstring, for function
#  emulation and geoid slicing.
geoid = type('Geoid', (object,),
             {'__doc__': "geoid(lat, lon) or geoid[...:...]",
              '__call__': lambda self, lat, lon: biv(lat, lon),
              '__getitem__': lambda self, index: getitem(index)})()


# NOW ADD NUMPY DEPENDENT CODE
DTYPE = np.float32

## A memmap to the data
egm_map = np.memmap(EGM_FILE, dtype=DTYPE, shape=SHAPE)

## Vectorized geoid interpolator
biv = np.vectorize(bilinear)


__all__ = ('geoid',)


if __name__ == '__main__':
    import sys
    ## parse args 1, 3, 5, ...
    lat_ = np.array(sys.argv[1::2], dtype=float)
    ## parse args 2, 5, 6, ...
    lon_ = np.array(sys.argv[2::2], dtype=float)
    for height in geoid(lat_, lon_):
        print height
