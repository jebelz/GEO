#!/usr/bin/env python

## \namespace geo.detic.newton.wgs84.ecef2llh Self-referential
## shell script for coordinate conversion
"""Script to convert input arguments to output stream:

ecef2llh x1 y1 z1 ... xn xn zn

Note: function used is derived from module name and imported from wgs84.
stdout is redirected until final output stream starts.
"""
## \brief{Generic script that changes its function based on its name}


null = type('Null', (object,), {'write': lambda self, *args: None})()
import sys
temp = sys.stdout
sys.stdout = null
from . import wgs84

funcname = str(sys.modules['__main__']).split()[-1][3:-5]
func = getattr(wgs84, funcname)

if __name__ == '__main__':
    sys.stdout = temp
    while True:
        try:
            x, y, z = map(float, sys.stdin.readline().split())
        except ValueError:
            break
        result = func(x, y, z)
        print result[0], result[1], result[2]
