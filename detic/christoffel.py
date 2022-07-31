"""LTP <--> SCH co/contra variant derivatives, Jacobians, tangent
bundles, and proto Christofel symbols...

Derivative functions used SYMPY string to compute them symbolically (which
you load as follows:

>>>from geo.detic.christoffel import *
>>>exec(SYMPY_SCH)
>>>exec(SYMPY_LLH)
>>>diff(X, lat)    # for example

Of course, you don't need simpy, as this has already been done--the functions
are in the module.

"""
#pylint: disable=C0301,W0613,C0103

# We need to true division for sympy's exponents.
from __future__ import division
import abc
import functools
import itertools
import operator

from ..utils.trig import pi, sqrt, sin, cos

## \namespace geo.detic.christoffel Diffeomorphisms and Tangent Bundles
# @image html christoffel.JPG .
# @image latex christoffel.JPG .


## Helper to take inner product of basis iterales and make a tensor:
# \param h1 \f$ {\bf \hat{h}}^{(1)}_i \f$ A tuple of basis vectors
# \param h2 \f$ {\bf \hat{h}}^{(2)}_j \f$ A tuple of basis vectors
# \returns \f$ g_{ij} = {\bf \hat{h}}^{(1)}_i \cdot {\bf \hat{h}}^{(2)}_j \f$
def _metric_tensor(h1, h2):
    from ...metric.euclid.tensor import Tensor
    return Tensor(*[(a*b).w for a, b in itertools.product(h1, h2)])


## Dual for any linear independent triplet:
# \param e1 \f$ {\bf \hat{e}}_1 \f$
# \param e2 \f$ {\bf \hat{e}}_2 \f$
# \param e3 \f$ {\bf \hat{e}}_2 \f$
# \returns \f$ e^1 = \frac{e_2 \times e_3}{e_1 \cdot (e_2 \times e_3)} \f$
def dual(e1, e2, e3):
    """(e1 ^ e2) / e1*(e1 ^ e2)-- is dual to e1, e2"""
    e23 = e2 ^ e3
    return e23 / (e1 * e23)


## Duals for any linear independent triplet:
# \param e1 \f$ {\bf \hat{e}}_1 \f$
# \param e2 \f$ {\bf \hat{e}}_2 \f$
# \param e3 \f$ {\bf \hat{e}}_2 \f$
# \returns list \f$ (e_1, e_2, e_3) \rightarrow (e^1, e^2, e^3) \f$
def duals(e1, e2, e3):
    """Get dual(e1, e2, e3) for all + permutations of arguments."""
    return [dual(x, y, z) for (x, y, z) in [(e1, e2, e3),
                                            (e3, e1, e2),
                                            (e2, e3, e1)]]


## >>>exex(SYMPY_SCH) will run these for you
SYMPY_SCH = """
from sympy import *

x, y, z = symbols('x y z')
s, c, h = symbols('s c h')
R= symbols('R')

S = R*atan2(x, R+z)
C = R*asin(y/sqrt(x**2 + y**2 + (R+z)**2))
H = sqrt(x**2 + y**2 + (R+z)**2) - R

ds_dx = diff(S, x)
ds_dy = diff(S, y)
ds_dz = diff(S, z)

dc_dx = diff(C, x)
dc_dy = diff(C, y)
dc_dz = diff(C, z)


dh_dx = diff(H, x)
dh_dy = diff(H, y)
dh_dz = diff(H, z)


## REMMBER: LTP's xyz is swapped from nominal spherical coordinate inversions.
X = (h + R) * sin(pi/2 - c/R) * sin(s/R)
Y = (h + R) * cos(pi/2 - c/R)
Z = (h + R) * sin(pi/2 - c/R) * cos(s/R) - R

dx_ds = diff(X, s)
dx_dc = diff(X, c)
dx_dh = diff(X, h)

dy_ds = diff(Y, s)
dy_dc = diff(Y, c)
dy_dh = diff(Y, h)

dz_ds = diff(Z, s)
dz_dc = diff(Z, c)
dz_dh = diff(Z, h)
"""


## >>>exex(SYMPY_LLH) will run these for you
SYMPY_LLH = """
from sympy import *

x, y, z = symbols('x y z')
lat, lon, hgt = symbols('lat lon hgt')
a, b = symbols('a b')

N = (a**2 / sqrt( (a * cos(pi*lat/180))**2 + (b * sin(pi*lat/180))**2 ))

M = (a*b)**2 / sqrt( (a * cos(pi*lat/180))**2 + (b * sin(pi*lat/180))**2 )**3

epsilon2 = (1-(b/a)**2)

X = cos(pi*lat/180) * cos(pi*lon/180) * ((a**2 / sqrt((a * cos(pi*lat/180))**2 + (b * sin(pi*lat/180))**2)) + hgt)
Y = cos(pi*lat/180) * sin(pi*lon/180) * ((a**2 / sqrt((a * cos(pi*lat/180))**2 + (b * sin(pi*lat/180))**2)) + hgt)
Z = sin(pi*lat/180) * ( (b/a)**2 *  (a**2 / sqrt((a * cos(pi*lat/180))**2 + (b * sin(pi*lat/180))**2)) + hgt)

dX_dlat = diff(X, lat)
dX_dlon = diff(X, lon)
dX_dhgt = diff(X, hgt)

dY_dlat = diff(Y, lat)
dY_dlon = diff(Y, lon)
dY_dhgt = diff(Y, hgt)

dZ_dlat = diff(Z, lat)
dZ_dlon = diff(Z, lon)
dZ_dhgt = diff(Z, hgt)


hdg = symbols('hdg')

LRC = 1/(cos(pi*hdg/180)**2/M + sin(pi*hdg/180)**2/N)

## lower metric
g_lon2 = diff(1/X, lon)**2 + diff(1/Y, lon)**2 + diff(1/Z, lon)**2
g_lat2 = diff(1/X, lat)**2 + diff(1/Y, lat)**2 + diff(1/Z, lat)**2
g_hgt2 = diff(1/X, hgt)**2 + diff(1/Y, hgt)**2 + diff(1/Z, hgt)**2

## Proto Christofel symbol
G_latlatlat = -0.5*(diff(g_lat2, lat))
G_latlatlon = -0.5*(diff(g_lat2, lon))
G_latlathgt = -0.5*(diff(g_lat2, hgt))
"""


## Computed Christoffel Symbol
# \returns \f$ \Gamma^m_{ij} = {\bf \epsilon}^m \cdot
#  \frac{\partial {\bf \epsilon}_i}{\partial q} \f$
def christoffel(e_i, q_j, e_m):
    """christoffel symbol:

    christoffel(e_i, q_j, e_m) --> e_m * (q_j >> e_i)
    """
    return (e_m * (q_j >> e_i)).w


## \f$ \frac{\partial s}{\partial x} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
# \returns  \f$ \frac{R(R+z)}{x^2 + (R+z)^2} \f$
def dsdx(R, x, y, z):
    """ds/dx"""
    return (
        R*(R + z)/(x**2 + (R + z)**2)
    )


## \f$ \frac{\partial s}{\partial y} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
# \returns \f$ 0 \f$
def dsdy(R, x, y, z):
    """ds/dy"""
    return (
        0*x
        )


## \f$ \frac{\partial s}{\partial z} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
# \returns \f$ \frac{-Rx}{x^2 + (R+z)^2} \f$
def dsdz(R, x, y, z):
    """ds/dz"""
    return (
        -R*x/(x**2 + (R + z)**2)
        )


## \f$ \frac{\partial c}{\partial x} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
# \returns \f$ \frac{-Rxy}{\sqrt{1+\frac{-y^2}{x^2+y^2+(R+z)^2}}(x^2+y^2+(R+z)^2)^{\frac{3}{2}}} \f$
def dcdx(R, x, y, z):
    """dc/dx"""
    return (
        -R*x*y/(sqrt(-y**2/(x**2 + y**2 + (R + z)**2) + 1)*(x**2 + y**2 + (R + z)**2)**(3/2.))
        )


## \f$ \frac{\partial c}{\partial y} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
# \returns \f$
#  R(-y^2/(x^2 + y^2 + (R + z)^2)^{\frac{3}{2}} + 1/\sqrt{x^2 + y^2 +
# (R + z)^2})/\sqrt{-y^2/(x^2 + y^2 + (R + z)^2) + 1} \f$
def dcdy(R, x, y, z):
    """dc/dy"""
    return (
        R*(-(y**2)/(x**2 + y**2 + (R + z)**2)**(3/2.) +
           1/sqrt(x**2 + y**2 + (R + z)**2))/sqrt(
               -y**2/(x**2 + y**2 + (R + z)**2) + 1)
        )


## \f$ \frac{\partial c}{\partial z} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
# \returns \f$ \frac{-Ry(R+z)}{\sqrt{1+\frac{-(y^2)}{x^2+y^2+(R+z)^2}}
# (x^2+y^2+(R+z)^2)^{\frac{3}{2}}} \f$
def dcdz(R, x, y, z):
    """dc/dz"""
    return (
        -R*y*(R+z)/(sqrt(-(y**2)/(x**2 + y**2 +
                                  (R + z)**2) + 1)*(x**2 + y**2 +
                                                    (R + z)**2)**(3/2.)))


## \f$ \frac{\partial h}{\partial x} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$ \frac{x}{\sqrt{x^2 +y^2 + (R+z)^2}}  \f$
def dhdx(R, x, y, z):
    """dh/dx"""
    return x/sqrt(x**2 + y**2 + (R + z)**2)


## \f$ \frac{\partial h}{\partial y} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$  \frac{y}{\sqrt{x^2 +y^2 + (R+z)^2}}  \f$
def dhdy(R, x, y, z):
    """dh/dy"""
    return y/sqrt(x**2 + y**2 + (R + z)**2)


## \f$ \frac{\partial h}{\partial z} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$ \frac{R+z}{\sqrt{x^2 +y^2 + (R+z)^2}}  \f$
def dhdz(R, x, y, z):
    """dh/dz"""
    return (
        (R + z)/sqrt(x**2 + y**2 + (R + z)**2)
        )


## \f$ \frac{\partial x}{\partial s} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$  \frac{(R+h)}{R}\cos{\frac{c}{R}}\cos{\frac{s}{R}}\f$
def dxds(R, s, c, h):
    """dx/ds"""
    return (R + h)*cos(c/R)*cos(s/R)/R


## \f$ \frac{\partial x}{\partial c} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$  \frac{-(R+h)}{R}\sin{\frac{c}{R}}\sin{\frac{s}{R}}\f$
def dxdc(R, s, c, h):
    """dx/dc"""
    return -(R + h)*sin(c/R)*sin(s/R)/R



## \f$ \frac{\partial x}{\partial s} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$  \sin{\frac{s}{R}}\cos{\frac{c}{R}}\f$
def dxdh(R, s, c, h):
    """dx/dh"""
    return sin(s/R)*cos(c/R)


## \f$ \frac{\partial y}{\partial s}  \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$ \frac{-(R+h)}{R}\sin{\frac{s}{R}} \f$
def dyds(R, s, c, h):
    """dy/ds = 0 (like s)"""
    return s*0


## \f$ \frac{\partial y}{\partial c} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$ \frac{R+h}{R}\cos{\frac{c}{R}}\f$
def dydc(R, s, c, h):
    """dy/dc"""
    return (R + h)*cos(c/R)/R


## \f$ \frac{\partial y}{\partial z} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$ \cos{\frac{s}{R}}\f$
def dydh(R, s, c, h):
    """dy/dh"""
    return sin(c/R)


## \f$ \frac{\partial z}{\partial s} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$ \frac{R+h}{R} \sin{\frac{c}{R}} \sin{\frac{s}{R}} \f$
def dzds(R, s, c, h):
    """dz/ds"""
    return -(R + h)*sin(s/R)*cos(c/R)/R



## \f$ \frac{\partial z}{\partial s} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$ -\frac{R+h}{R}\sin{\frac{c}{R}}\cos{\frac{s}{R}} \f$
def dzdc(R, s, c, h):
    """dz/dc"""
    return -(R + h)*sin(c/R)*cos(s/R)/R



## \f$ \frac{\partial z}{\partial s} \f$
# \param R local radius of curvature along x-direction
# \param x x-coordinate
# \param y y-coordinate
# \param z z-coordinate
#  \returns \f$ \sin{\frac{c}{R}}\sin{\frac{s}{R}}  \f$
def dzdh(R, s, c, h):
    """dz/dh"""
    return cos(c/R)*cos(s/R)



## \f$ \frac{\partial x}{\partial \phi} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
# \returns \f$ a^2(a^2\sin{\phi}\cos{\phi} - b^2\sin{\phi}\cos{\phi})\frac{\cos{\phi}\cos{\lambda}}{(a^2\cos^2{\phi}+b^2\sin^2{\phi})^{\frac{3}{2}}} - (h+\frac{a^2}{\sqrt{a^2\cos^2{\phi}+b^2\sin^2{\phi}}})\sin{\phi}\cos{\lambda} \f$
def dXdlat(a, b, lat, lon, hgt):
    """dX/dlat"""
    return (
        a**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*cos(pi*lat/180)*cos(pi*lon/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2.) - pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*sin(pi*lat/180)*cos(pi*lon/180)/180
        )

## \f$ \frac{\partial x}{\partial \lambda} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
# \returns \f$ -\frac{\pi}{180}(\frac{a^2}{\sqrt{a^2\cos^2{\phi}+b^2
# \sin^2{\phi}}}+h)\sin{\lambda}\cos{\lambda} \f$
def dXdlon(a, b, lat, lon, hgt):
    """dX/dlon"""
    return (
        -pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*sin(pi*lon/180)*cos(pi*lat/180)/180
    )


## \f$ \frac{\partial x}{\partial h} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
#  \returns \f$ \cos{\phi}\cos{\lambda} \f$
def dXdhgt(a, b, lat, lon, hgt):
    """dX/dhgt"""
    return cos(pi*lat/180)*cos(pi*lon/180)


## \f$ \frac{\partial y}{\partial \phi} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
# \returns \f$ a^2(a^2\sin{\phi}\cos{\phi} - b^2\sin{\phi}\cos{\phi})\frac{\cos{\phi?}\sin{\lambda}}{(a^2\cos^2{\phi}+b^2\sin^2{\phi})^{\frac{3}{2}}} - (h+\frac{a^2}{\sqrt{a^2\cos^2{\phi}+b^2\sin^2{\phi}}})\sin{\phi}\sin{\lambda} \f$
def dYdlat(a, b, lat, lon, hgt):
    """dY/dlat"""
    return (
        a**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*sin(pi*lon/180)*cos(pi*lat/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2.) - pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*sin(pi*lat/180)*sin(pi*lon/180)/180
        )


## \f$ \frac{\partial y}{\partial \lambda} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
# \returns \f$ +\frac{\pi}{180}(\frac{a^2}{\sqrt{a^2\cos^2{\phi}+b^2\sin^2{\phi}}}+h)\cos{\lambda}\cos{\lambda} \f$
def dYdlon(a, b, lat, lon, hgt):
    """dY/dlon"""
    return (
        pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*cos(pi*lat/180)*cos(pi*lon/180)/180
        )


## \f$ \frac{\partial y}{\partial h} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
#  \returns \f$ \cos{\phi}\sin{\lambda} \f$
def dYdhgt(a, b, lat, lon, hgt):
    """dY/dhgt"""
    return sin(pi*lon/180)*cos(pi*lat/180)



## \f$ \frac{\partial z}{\partial \phi}  \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
def dZdlat(a, b, lat, lon, hgt):
    """dZ/dlat"""
    return (
        b**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*sin(pi*lat/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2) + pi*(b**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*cos(pi*lat/180)/180
        )


## \f$ \frac{\partial z}{\partial h} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
#  \returns \f$  0 \f$
def dZdlon(a, b, lat, lon, hgt):
    """dZ/dlon"""
    return 0*lat



## \f$ \frac{\partial z}{\partial h} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
#  \returns \f$ \sin{\phi} \f$
def dZdhgt(a, b, lat, lon, hgt):
    """dZ/dhgt"""
    return sin(pi*lat/180)



## \f$ \frac{dN}{d\phi} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
def dNdlat(a, b, lat):
    """dN/dlat"""
    return (
        a**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2)
        )


## \f$ \frac{dM}{d\phi} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
def dMdlat(a, b, lat):
    """dM/dlat"""
    return (
        a**2*b**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/60 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/60)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(5/2)
        )


## \f$ \frac{\partial R}{\partial \phi} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
def dRdlat(a, b, lat, hdg):
    """dR/dlat"""
    return (
        (-(-pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 + pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*sin(pi*hdg/180)**2/(a**2*sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)) - sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)*(-pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/60 + pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/60)*cos(pi*hdg/180)**2/(a**2*b**2))/(sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)*sin(pi*hdg/180)**2/a**2 + (a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2)*cos(pi*hdg/180)**2/(a**2*b**2))**2
        )


## \f$ \frac{\partial R}{\partial \psi} \f$
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
def dRdhdg(a, b, lat, hdg):
    """dR/dhdg"""
    return (
        (-pi*sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)*sin(pi*hdg/180)*cos(pi*hdg/180)/(90*a**2) + pi*(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2)*sin(pi*hdg/180)*cos(pi*hdg/180)/(90*a**2*b**2))/(sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)*sin(pi*hdg/180)**2/a**2 + (a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2)*cos(pi*hdg/180)**2/(a**2*b**2))**2
        )


## Experiment christoffel symbols from the metric...
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
def g_lat_lon(a, b, lat, lon, hgt):
    """metric tensor component's derivaive"""
    return (
        pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*(a**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*sin(pi*lon/180)*cos(pi*lat/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2) - pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*sin(pi*lat/180)*sin(pi*lon/180)/180)*cos(pi*lat/180)*cos(pi*lon/180)/180 - pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*(a**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*cos(pi*lat/180)*cos(pi*lon/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2) - pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*sin(pi*lat/180)*cos(pi*lon/180)/180)*sin(pi*lon/180)*cos(pi*lat/180)/180
    )


## diff(X, lat)*diff(X, lat) + diff(Y, lat)*diff(Y, lat) + diff(Z, lat)*diff(Z, lat)
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
def g_lat_lat(a, b, lat, lon, hgt):
    """metric tensor component's derivaive"""
    return (
        (b**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*sin(pi*lat/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2) + pi*(b**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*cos(pi*lat/180)/180)**2 + (a**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*sin(pi*lon/180)*cos(pi*lat/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2) - pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*sin(pi*lat/180)*sin(pi*lon/180)/180)**2 + (a**2*(pi*a**2*sin(pi*lat/180)*cos(pi*lat/180)/180 - pi*b**2*sin(pi*lat/180)*cos(pi*lat/180)/180)*cos(pi*lat/180)*cos(pi*lon/180)/(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2)**(3/2) - pi*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)*sin(pi*lat/180)*cos(pi*lon/180)/180)**2
        )


## =          diff(X, lon)**2 + diff(Y, lon)**2 + diff(Z, lon)**2
# \param a Ellipsoid semi-major axis
# \param b Ellipsoid semi-minor axis
# \param lat Latitude, \f$ \phi \f$ (degrees)
# \param lon Longitude, \f$ \lambda \f$ (degrees)
# \param h Height relative to Ellipsoid (m)
def g_lon_lon(a, b, lat, lon, hgt):
    """metric tensor component's derivaive"""
    return (
        pi**2*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)**2*sin(pi*lon/180)**2*cos(pi*lat/180)**2/32400 + pi**2*(a**2/sqrt(a**2*cos(pi*lat/180)**2 + b**2*sin(pi*lat/180)**2) + hgt)**2*cos(pi*lat/180)**2*cos(pi*lon/180)**2/32400
        )


## Decorator to convert metric tensors into Levi Civita tensors.
def _E(method):
    """Decorator."""
    from ..metric.euclid.three import EPSILON

    @functools.wraps(method)
    def E(self):
        """Convert g to epsilon."""
        return abs(method(self).det().w) * EPSILON
    return E


## Base class for Frame differntiating mixins.
class _Frame(object):
    """ABC for Frame Mixins. Subs may be local or global frames."""

    __metaclass__ = abc.ABCMeta

    ## Concrete classes provide the Jacobian between conjugate coordinates:
    # \n ECEF <--> LLH
    # \n XYZ <--> SCH
    @abc.abstractmethod
    def jacobian(self):
        """Jacobian"""

    ## Covariant Basis Vectors
    # @image html covar.png
    # @image latex covar.png
    @abc.abstractmethod
    def basis_covariant(self):
        """Etc.."""

    ## Contravariant Basic Vectors
    # @image html contra.png
    # @image latex contra.png
    @abc.abstractmethod
    def basis_contravariant(self):
        """Etc.."""

    ## Concrete classes shall name this property, and it will point to
    # basis_contravariant() or basis_covariant().
    @abc.abstractproperty
    def _local_basis(self):
        """Etc.."""

    ## See jacobian().
    @property
    def J(self):
        """J is jacobian()."""
        return self.jacobian()

    ## Iterate sub's method named _local_basis
    def _iter_basis(self):
        return operator.methodcaller(self._local_basis)(self)

    ## The local basis is the normalized transform to the cannonicol
    # cartesian frame
    # \returns geo.metric.euclid.tensor.Tensor that is the basis unit vectors
    # in the right index.
    def local_basis(self):
        """Local basis vectors (zipped)."""
        from ..metric.euclid.tensor import ziprows
        return ziprows(*[item.hat() for item in self._iter_basis()])

    ## Convariant Metric Tensor
    # \returns \f$ g^{ij} \f$
    def metric_tensor_covariant(self):
        """Covariant metric tensor."""
        return _metric_tensor(*itertools.repeat(self.basis_covariant(), 2))

    ## Contravariant Metric Tensor
    # \returns \f$ g_{ij} \f$
    def metric_tensor_contravariant(self):
        """Contravariant metric tensor."""
        return _metric_tensor(*itertools.repeat(self.basis_contravariant(), 2))

    ## Mixed Metric Tensor
    # \returns \f$ g^i_j = \delta_{ij} \f$ (up to round-off error)
    def metric_tensor_mixed(self):
        """Mixed metric tensor (should be delta_ij)"""
        return _metric_tensor(self.basis_covariant(),
                              self.basis_contravariant())

    ## Covariant Alternating Tensor
    # \returns \f$ E^{ijk} = \epsilon^{ijk} \sqrt{||\det{g}||} \f$
    @_E
    def levi_civita_covariant(self):
        """Covariant Levi-Civita Tensor."""
        return self.metric_tensor_covariant()

    ## Contravariant Alternating Tensor
    # \returns \f$ E_{ijk} = \epsilon_{ijk} / \sqrt{||\det{g}||} \f$
    @_E
    def levi_civita_contravariant(self):
        """Contrvariant Levi-Civita Tensor."""
        return self.metric_tensor_contravariant()

    ## EXPERIMENTAL implementation of Chirstofel connection
    # \returns \f$ \Gamma^m_{ij} \f$ (maybe)- in a rank-3 Tensor body.
    def Gamma2(self):
        """Gamma 2?"""
        from ..metric.euclid.three import Three
        from ..metric.euclid.vector import BASIS
        return Three(*itertools.starmap(
            christoffel,
            itertools.product(BASIS,
                              self.basis_covariant(),
                              BASIS)))


## Mixin for coordinate with local basis vectors.
class _Local(_Frame):
    """Local frames have position dependent basis vectors."""

    __metaclass__ = abc.ABCMeta

    ## This is a tuple of method which compute tangent vectors.
    @abc.abstractproperty
    def _tan(self):
        """A tuple of tangent methods."""

    ## Local frames exprese their tangent vectors in their OTF's basis
    _local_basis = 'basis_covariant'

    ## Call methods in _tan() static attribute
    # \returns list of tangent vectors.
    def basis_covariant(self):
        """list of tangent (_tan) vectors in cartessian coordinates"""
        return [method(self) for method in self._tan]

    ## Call conjugate class's method
    # \returns list of cotangent vectors.
    def basis_contravariant(self):
        """list of cotangent basis vectors, from Cannonical
        Conjugate class (see __invert__.__doc__)"""
        return (~self).basis_contravariant()


## Mixin for coordinate with global basis vectors.
class _Global(_Frame):
    """Local frames have position independent basis vectors."""

    __metaclass__ = abc.ABCMeta

    ## Concretes have a tuple of gradient methods
    @abc.abstractproperty
    def _grad(self):
        """Tuple of gradient methods."""

    ## Global frames express the gradients of their local complements.
    _local_basis = 'basis_contravariant'

    ## Call methods in _grad static attribute
    # \returns list of coordinate gradients
    def basis_contravariant(self):
        """list of cotangent (_grad) vectors in cartessian coordinates:
        list(starmap(apply, izip(self._grad, repeat((self,), 3))))"""
        return [method(self) for method in self._grad]

    ## Call conjugate class's method
    # \returns list of basis tangent vectors
    def basis_covariant(self):
        """list of cotangent basis vectors, from Cannonical
        Conjugate class (see __invert__.__doc__)"""
        return (~self).basis_covariant()

    ## Cotangent synonym: it is the gradient.
    @property
    def _cotan(self):
        """_grad"""
        return self._grad


## Mixin for LTP's Tangent Bundle
class LTPMixin(_Global):
    """Mixin for the Contravariant SCH Basis Vectors."""

    ## Gradient of "s" coordinate.
    # \returns \f$ \vec{\nabla} s(x, y, z) \f$
    def grad_s(self):
        """Contravariant derivative of constant-s plane."""
        from ..metric.euclid.vector import Vector
        return Vector(dsdx(self.R, self.x, self.y, self.z),
                      dsdy(self.R, self.x, self.y, self.z),
                      dsdz(self.R, self.x, self.y, self.z))

    ## Gradient of "c" coordinate.
    # \returns \f$ \vec{\nabla} c(x, y, z) \f$
    def grad_c(self):
        """Contravariant derivative of constant-c plane."""
        from ..metric.euclid.vector import Vector
        return Vector(dcdx(self.R, self.x, self.y, self.z),
                      dcdy(self.R, self.x, self.y, self.z),
                      dcdz(self.R, self.x, self.y, self.z))

    ## Gradient of "h" coordinate.
    ## \returns \f$ \vec{\nabla} h(x, y, z) \f$
    def grad_h(self):
        """Contravariant derivative of constant-h plane."""
        from ..metric.euclid.vector import Vector
        return Vector(dhdx(self.R, self.x, self.y, self.z),
                      dhdy(self.R, self.x, self.y, self.z),
                      dhdz(self.R, self.x, self.y, self.z))

    ## The grad is for s, c, h
    _grad = (grad_s, grad_c, grad_h)

    ## This method is here to supply a docstring (and call super)
    def local_basis(self):
        """Normalized Tensor of:
        | s_x s_y s_z |
        | c_x c_y c_z |
        | h_x h_y h_z |
        from basis_contravariant."""
        return super(LTPMixin, self).local_basis()

    ## Jacobian from LTP (xyz) to SCH.
    #\f$ \bar{J}_{ij} = \frac{\partial q_i}{\partial {\bf x}_j}(x, y, z) \f$
    # \return \f$ \left[ \begin{array}{ccc}
    # \frac{\partial s}{\partial x} &  \frac{\partial s}{\partial y} &
    # \frac{\partial s}{\partial z} \\
    # \frac{\partial c}{\partial x} &  \frac{\partial c}{\partial y} &
    # \frac{\partial c}{\partial z} \\
    # \frac{\partial h}{\partial x} &  \frac{\partial h}{\partial y} &
    # \frac{\partial h}{\partial z} \end{array} \right]\f$
    def jacobian(self):
        """Jacobian from LTP:xyz to LTS:sch

        | ds/dx  ds/dy  ds/dz |
        |                     |
        | dc/dx  dc/dy  dc/dz |
        |                     |
        | dh/dx  dh/dy  dh/dz |
        """
        from ..metric.euclid.tensor import Tensor
        return Tensor(dsdx(self.R, self.x, self.y, self.z),
                      dsdy(self.R, self.x, self.y, self.z),
                      dsdz(self.R, self.x, self.y, self.z),

                      dcdx(self.R, self.x, self.y, self.z),
                      dcdy(self.R, self.x, self.y, self.z),
                      dcdz(self.R, self.x, self.y, self.z),

                      dhdx(self.R, self.x, self.y, self.z),
                      dhdy(self.R, self.x, self.y, self.z),
                      dhdz(self.R, self.x, self.y, self.z))


## Mixin for SCH's Tangent Bundle
class SCHMixin(_Local):
    """Mixin for the Covariant SCH Basis Vectors."""

    ## \f$s\f$-tangent vector
    # \returns \f$ \frac{\partial x}{\partial s} {\bf \hat{x}} +
    #  \frac{\partial y}{\partial s} {\bf \hat{y}} +
    #  \frac{\partial z}{\partial s} {\bf \hat{z}} \f$
    def s_tan(self):
        """Covariant tangent vector parallel to 's' at point."""
        from ..metric.euclid.vector import Vector
        return Vector(dxds(self.R, self.s, self.c, self.h),
                      dyds(self.R, self.s, self.c, self.h),
                      dzds(self.R, self.s, self.c, self.h))

    ## \f$c\f$-tangent vector
    # \returns \f$ \frac{\partial x}{\partial c} {\bf \hat{x}} +
    #  \frac{\partial y}{\partial c} {\bf \hat{y}} +
    #  \frac{\partial z}{\partial c} {\bf \hat{z}} \f$
    def c_tan(self):
        """Covariant tangent vector parallel to 'c' at point."""
        from ..metric.euclid.vector import Vector
        return Vector(dxdc(self.R, self.s, self.c, self.h),
                      dydc(self.R, self.s, self.c, self.h),
                      dzdc(self.R, self.s, self.c, self.h))

    ## \f$h\f$-tangent vector
    # \returns \f$ \frac{\partial x}{\partial h} {\bf \hat{x}} +
    #  \frac{\partial y}{\partial h} {\bf \hat{y}} +
    #  \frac{\partial z}{\partial h} {\bf \hat{z}} \f$
    def h_tan(self):
        """Covariant tangent vector parallel to 'h' at point."""
        from ..metric.euclid.vector import Vector
        return Vector(dxdh(self.R, self.s, self.c, self.h),
                      dydh(self.R, self.s, self.c, self.h),
                      dzdh(self.R, self.s, self.c, self.h))

    ## Tangents are s, c, h
    _tan = (s_tan, c_tan, h_tan)

    ## This method is here to supply a docstring (and call super).
    def local_basis(self):
        """Normalized Tensor of:
        | s_x s_y s_z |
        | c_x c_y c_z |
        | h_x h_y h_z |
        from basis_contravariant."""
        return super(SCHMixin, self).local_basis()

    ## Jacobian from SCH to LTP (xyz).
    # \f$ J_{ij} = \frac{\partial {\bf x}_i}{\partial q_j}(s, c, h) \f$
    # \return \f$ \left[ \begin{array}{ccc}
    # \frac{\partial x}{\partial s} &  \frac{\partial x}{\partial s} &
    # \frac{\partial x}{\partial h} \\
    # \frac{\partial y}{\partial s} &  \frac{\partial y}{\partial c} &
    # \frac{\partial y}{\partial h} \\
    # \frac{\partial z}{\partial s} &  \frac{\partial z}{\partial h} &
    # \frac{\partial z}{\partial h} \end{array} \right]\f$
    def jacobian(self):
        """Jacobian from sch to LTP: xyz:

        | dx/ds  dx/dc  dx/dh |
        |                     |
        | dy/ds  dy/dc  dy/dh |
        |                     |
        | dz/ds  dz/dc  dz/dh |
        """
        from ..metric.euclid.tensor import Tensor
        return Tensor(dxds(self.R, self.s, self.c, self.h),
                      dxdc(self.R, self.s, self.c, self.h),
                      dxdh(self.R, self.s, self.c, self.h),

                      dyds(self.R, self.s, self.c, self.h),
                      dydc(self.R, self.s, self.c, self.h),
                      dydh(self.R, self.s, self.c, self.h),

                      dzds(self.R, self.s, self.c, self.h),
                      dzdc(self.R, self.s, self.c, self.h),
                      dzdh(self.R, self.s, self.c, self.h))


## LLH Mixin
class LLHMixin(_Local):
    """LLH vs ECEF mixin:
    Covariance (tangent) of lat, lon, hgt wrt x, y z
    """

    ## All hail Demeter: this class's math is ellipsoid dependent.
    @property
    def _a(self):
        """ellipsoid's semi-major-axis"""
        return self.ellipsoid.a

    ## IBID.
    @property
    def _b(self):
        """ellipsoid's semi-minor-axis"""
        return self.ellipsoid.b

    ## Latitude Tangent
    # \returns \f$ \frac{\partial x}{\partial \phi}{\bf \hat{x}} +
    # \frac{\partial y}{\partial \phi}{\bf \hat{y}} +
    # \frac{\partial z}{\partial \phi}{\bf \hat{z}}  \f$
    def lat_tan(self):
        """Covariant tangent vector parallel to 's' at point."""
        from ..metric.euclid.vector import Vector
        return Vector(dXdlat(self._a, self._b, self.lat, self.lon, self.hgt),
                      dYdlat(self._a, self._b, self.lat, self.lon, self.hgt),
                      dZdlat(self._a, self._b, self.lat, self.lon, self.hgt))

    ## Longitude Tangent
    # \returns \f$ \frac{\partial x}{\partial \lambda}{\bf \hat{x}} +
    # \frac{\partial y}{\partial \lambda}{\bf \hat{y}} +
    # \frac{\partial z}{\partial \lambda}{\bf \hat{z}}  \f$
    def lon_tan(self):
        """Covariant tangent vector parallel to 'c' at point."""
        from ..metric.euclid.vector import Vector
        return Vector(dXdlon(self._a, self._b, self.lat, self.lon, self.hgt),
                      dYdlon(self._a, self._b, self.lat, self.lon, self.hgt),
                      dZdlon(self._a, self._b, self.lat, self.lon, self.hgt))

    ## Height Tangent
    # \returns  \f$ \frac{\partial x}{\partial h}{\bf \hat{x}} +
    # \frac{\partial y}{\partial h}{\bf \hat{y}} +
    # \frac{\partial z}{\partial h}{\bf \hat{z}}  \f$
    def hgt_tan(self):
        """Covariant tangent vector parallel to 'h' at point."""
        from ..metric.euclid.vector import Vector
        return Vector(dXdhgt(self._a, self._b, self.lat, self.lon, self.hgt),
                      dYdhgt(self._a, self._b, self.lat, self.lon, self.hgt),
                      dZdhgt(self._a, self._b, self.lat, self.lon, self.hgt))

    ## Tangents are lat, lon, hgt
    _tan = (lat_tan, lon_tan, hgt_tan)

    ## Jacobina from LLH to ECEF:
    #  \f$ J_{ij} = \frac{\partial {\bf x}_i}{\partial q_j}(s, c, h) \f$
    # \return \f$ \left[ \begin{array}{ccc}
    # \frac{\partial x}{\partial \phi} &
    # \frac{\partial x}{\partial \lambda} &
    # \frac{\partial x}{\partial h} \\
    # \frac{\partial y}{\partial \phi} &
    # \frac{\partial y}{\partial \lambda} &
    # \frac{\partial y}{\partial h} \\
    # \frac{\partial z}{\partial \phi} &
    # \frac{\partial z}{\partial \lambda} &
    # \frac{\partial z}{\partial h} \end{array} \right]\f$
    def jacobian(self):
        """Jacobian from LLH to ECEF:

        | dX/dlat  dX/dlon  dX/dhgt |
        |                     |
        | dY/dlat  dY/dlon  dY/dhgt |
        |                     |
        | dZ/dlat  dZ/dlon  dZ/dhgt |
        """
        from ..metric.euclid.tensor import Tensor
        return Tensor(dXdlat(self._a, self._b, self.lat, self.lon, self.hgt),
                      dXdlon(self._a, self._b, self.lat, self.lon, self.hgt),
                      dXdhgt(self._a, self._b, self.lat, self.lon, self.hgt),

                      dYdlat(self._a, self._b, self.lat, self.lon, self.hgt),
                      dYdlon(self._a, self._b, self.lat, self.lon, self.hgt),
                      dYdhgt(self._a, self._b, self.lat, self.lon, self.hgt),

                      dZdlat(self._a, self._b, self.lat, self.lon, self.hgt),
                      dZdlon(self._a, self._b, self.lat, self.lon, self.hgt),
                      dZdhgt(self._a, self._b, self.lat, self.lon, self.hgt))


## ECEF mixin
class ECEFMixin(_Global):
    """ECEF Mixin wrt to LLH:
    contravariance (gradient) of LLH w.r.t (x, y, z)"""

    ## Gradient of Latitude
    # \returns \f$ \vec{\nabla}\phi(x, y, z) \f$
    def grad_lat(self):
        """Gradient of 'lat'"""
        return ~(self.llh().lat_tan())

    ## Gradient of Longitude
    # \f$ \vec{\nabla}\lambda(x, y, z) \f$
    def grad_lon(self):
        """Gradient of 'lon'"""
        return ~(self.llh().lon_tan())

    ## Gradient of Height
    # \f$ \vec{\nabla}h(x, y, z) \f$
    def grad_hgt(self):
        """Gradient of 'h'"""
        return ~(self.llh().hgt_tan())

    ## Ordered tuple of gradient computing methods.
    _grad = (grad_lat, grad_lon, grad_hgt)

    ## This method is here to supply a docstring (and call super)
    def local_basis(self):
        """Normalized Tensor of:
        | lat_,x lat_,y lat_,z |
        | lon_,x lon_,y lon_,z |
        | hgt_,x hgt_,y hgt_,z |
        from basis_covariant."""
        return super(ECEFMixin, self).local_basis()

    ## The Jacobian of an iterative transformation is not so analytic,
    # hence, this inverts the conjugate-coordinates' Jacobian:
    # \f$ \bar{J}_{ij} =
    # \frac{\partial {\bf q}_i}{\partial x_j}(\bf{\vec{x}}) \f$
    # \return \f$ \left[ \begin{array}{ccc}
    # \frac{\partial\phi}{\partial x} &
    # \frac{\partial\phi}{\partial y} &
    # \frac{\partial\phi}{\partial z} \\
    # \frac{\partial\lambda}{\partial x} &
    # \frac{\partial\lambda}{\partial y} &
    # \frac{\partial\lambda}{\partial z} \\
    # \frac{\partial h}{\partial x} &
    # \frac{\partial h}{\partial y} &
    # \frac{\partial h}{\partial z} \end{array} \right]\f$
    def jacobian(self):
        """Jacobian from ECEF to LLH:

        | dlat/dX  dlat/dY  dlat/dZ |
        |                           |
        | dlon/dX  dlon/dY  dlon/dZ |
        |                           |
        | dhgt/dX  dhgt/dY  dhgt/dZ |"""
        return (~self).J.I
