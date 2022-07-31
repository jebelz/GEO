"""Straight Geomorphic Curvature Functions for Earth scientist."""
## \namespace geo.morphic.gauss From slopes to Gaussian curvature
import functools

import numpy as np


## \f$ S = \sqrt{f_x^2+f_y^2} \f$ \n via np.hypot
S = np.hypot


## Normalized a function by \f$ \frac{1}{S^n} \f$
# \param n Power of slope in denominator
# \returns function that normalized functions
def _slope(power):
    """Decorate Factory with a power of slope for normalization."""

    def decoration(func):
        """Decorator"""

        # Note: you need to wrap, since the docstring (antipattern)
        # has info-- I should use func.im_name or something....
        @functools.wraps(func)
        def normalized(f_xx, f_xy, f_yy, f_x, f_y):
            """func() / slope ** power"""
            return -(
                func(f_xx, f_xy, f_yy, f_x, f_y)/S(f_x, f_y)**power
                )
        return normalized
    return decoration


## \f$ \kappa_p =
# \frac{-(f_x^2f_xx + 2 f_xf_yf_xy + f_y^2f_yy)}{S^2} \f$\n
# The Profile is __very__ similar to CRA
@_slope(2)
def profile(f_xx, f_xy, f_yy, f_x, f_y):
    """Profile Curvature: The curvate of a cut aling the downhill slope."""
    return f_x**2*f_xx + 2*f_x*f_y*f_xy + f_y**2*f_yy


## Plan (contour) curvature: \n
# The rate of change of slope along a contour\n
# \f$ \kappa_c =
# \frac{-(f_y^2f_{xx} - 2 f_xf_yf_{xy} + f_x^2f_{yy})} {S^3} \f$
@_slope(3)
def plan(f_xx, f_xy, f_yy, f_x, f_y):
    """Plan Curvature: The Curvature of a contour line."""
    return f_y**2*f_xx - 2*f_x*f_y*f_xy + f_x**2*f_yy


## \f$ \kappa_s =
# \frac{(f_x f_y(f_{xx}-f_{yy}) + (f_y^2-f_x^2)f_{xy})}{S^3} \f$\n
# Streamline is similar to plan()
@_slope(3)
def streamline(f_xx, f_xy, f_yy, f_x, f_y):
    """Streamlime Curvature."""
    return f_x*f_y*(f_xx - f_yy) + (f_y**2-f_x**2)*f_xy


## \f$ \kappa_t =
# \frac{ f_y^2f_{xx} -2f_xf_yf_{xy} + f_x^2f_{yy}}
# {S^2\sqrt{1+S^2}} \f$\n
# Tangential is similar to CRX.
@_slope(2)
def tangential(f_xx, f_xy, f_yy, f_x, f_y):
    """Tangential Curvature."""
    s = S(f_x, f_y)
    return (f_y**2*f_xx - 2*f_x*f_y*f_xy + f_x**2*f_yy)/np.hypot(s, 1)


## \f$ T = f^2_{xx} + 2f_{xy} + f^2_{yy}\f$
# The OTHER Total Curvature
def total(f_xx, f_xy, f_yy, f_x, f_y):
    """Total Curvature No. 2"""
    return f_xx + 2*f_xy + f_yy


## \f$ K = (f^2_{xx} + 2f_{xy} + f^2_{yy})/S^2\f$
@_slope(2)
def gaussian(f_xx, f_xy, f_yy, f_x, f_y):
    """Gaussian Curvature."""
    return total(f_xx, f_xy, f_yy, f_x, f_y)


## \f$ H = \frac{1}{2}\f$
@_slope(3)
def mean(f_xx, f_xy, f_yy, f_x, f_y, f_z=-1):
    """Mean Curvature."""
    return ((1+f_x**2)*f_yy - 2*f_y*f_z*f_xy + (1+f_y**2)*f_xx)/2.


## __T__ errain __R__ ugeddness __I__ ndex (Riley, 1999).
# \param box A size=(3,3) array.
# \return \f$ R_I = \sum_{i=-1,1}\sum_{j=-1,1}{{{(h_{ij}-h_{00})^2}}}\f$
def tri(box):
    """
    TRI Categories
    LEVEL = 0-80m
    NEARLY LEVEL = 81-116m
    SLIGHTLY RUGGED = 117-161m
    INTERMEDIATELY RUGGED = 162-239m
    MODERATELY RUGGED = 240-497m
    HIGHLY RUGGED = 498-958m
    EXTREMELY RUGGED = 959-4367m
    """
    return ((box-box[1, 1])**2).sum()**0.5


## <a href="http://gis4geomorphology.com/roughness-topographic-position/">
# Surface Roughness Factor</a>.
def srf():
    """
    Hobson (1972) introduced the Surface Roughness Factor (SRF).
    Xi, Yi, and Zi are vectors normal to the land surface.

    SRF = sqrt (sum:ni=1 Xi)2 + (sum:ni=1 Yi)2 + (sum:ni=1 Zi)2 / n

    Xi = sin(s) * cos(a)
    Yi = sin(s) * sin(a)
    Zi = cos(s)
    n = number of cells in the analysis window (sample size)
    """
    return NotImplemented
