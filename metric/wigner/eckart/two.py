"""Pure rank-2 Conversion Formulae.


3 x 3  breaks up as follow:

Weyl Module  dimW Irreps (j; q)         T_lambda (all tensors are
Tableau                                            symmetrized)
=====================================================================
[0][1]         6     1   (0; 0)              s = T.ii
                     5   (2; 0)           T.S() - s * DELTA / 3

[0]            3     3   (1; 0)            v.i = (EPSILON & T).ijkjk / 2
[1]


The (forward) functions defined within are designed to operatate on a
Symmetric, Trace Free (aka, Natural Form), rank-2 Cartesian Tensor. 

To select the stength of m-eigenstates from a natural form Cartesian tensor:

m    function
--------------
+2    t22
+1    t21
+0    t20
-1    t2m1
-2    t2m2


To get the Cartesian form of the j-eigenstates (for all m) from a spherical
tensor:


j    function
--------------
0    T0
1    T1
2    T2
"""
## \namespace geo.metric.wigner.eckart.two Pure Rank-2 Formulae
from .one import ROOT2, ROOT3, SCALAR_NORM



## Doc string writer (strickly lazy meta code).
def doc(func):
    j = func.func_name[1]
    m = func.func_name[-1]
    if 'm' in func.func_name:
        m = '-' + m
    doc = "T[{}, {}] = {}(S).".format(j, m, func.func_name)
    doc += "\nS is a SYMMETRIC, TRACE-FREE Cartesian Rank 2 Tensor"
    func.__doc__ = doc
    return func
    
    
## \f$ \langle 2, 0|S|2, 0 \rangle \f$
# \param S A rank-2 natural form Cartesian tensor
# \returns \f$ \sqrt{\frac{3}{2}}S_{zz} \f$
@doc
def t20(S):
    return  ROOT3 * S.zz / ROOT2


## \f$ \langle 2, 1|S|2, 1 \rangle \f$
# \param S A rank-2 natural form Cartesian tensor
# \returns \f$-(S_{xz} - iS_{yz}) \f$
@doc
def t21(S):
    return -(S.zx - 1j*S.zy)


## \f$ \langle 2, -1|S|2, -1 \rangle \f$
# \param S A rank-2 natural form Cartesian tensor
# \returns \f$(S_{xz} + iS_{yz}) \f$
@doc
def t2m1(S):
    return S.zx + 1j*S.zy


## \f$ \langle 2, 2|S|2, 2 \rangle \f$
# \param S A rank-2 natural form Cartesian tensor
# \returns \f$\frac{1}{2}[(S_{xx}-S_{yy}) - 2iS_{xy}] \f$
@doc
def t22(S):
    return ((S.xx - S.yy) - 2*1j*S.xy) / 2.


## \f$ \langle 2, -2|S|2, -2 \rangle \f$
# \param S A rank-2 natural form Cartesian tensor
# \returns \f$\frac{1}{2}[(S_{xx}-S_{yy}) + 2iS_{xy}] \f$
@doc
def t2m2(S):
    return ((S.xx - S.yy) + 2*1j*S.xy) / 2.


## Hash Table of Natural Form Projections
PROJECTIONS = {2: t22,
               1: t21,
               0: t20,
               -1: t2m1,
               -2: t2m2}


## |2, m>, Cartesian
# \returns
# \f$ T^2 = \frac{1}{2} \left[\begin{array}{ccc}
# T^2_2 + T^2_{-2} - \sqrt{\frac{2}{3}} T^2_0 &
# -i(T^2_2 - T^2_{-2}) &
# -T^2_1 + T^2_{-1} \\
# -i(T^2_2 - T^2_{-2}) &
# -T^2_2 - T^2_{-2} - \sqrt{\frac{2}{3}} T^2_0 &
# i(T^2_1 + T^2_{-1}) \\
# -T^2_1 + T^2_{1} &
# i(T^2_1 + T^2_{-1}) &
# \sqrt{\frac{8}{3}}T^2_0
# \end{array}\right] \f$
# @image html T2m.tiff
# @image latex T2m.tiff
def T2(t):
    """Symmetric Trace-free part, as a Cartesian Tensor."""
    from ...euclid.tensor import from_natural_form

    xx_minus_yy = t[2, 2] + t[2, -2]

    zz = ROOT2 * t[2, 0] / ROOT3

    # 2i(ixy+ixy)/2 = 2ixy 
    xy = -(t[2, 2] - t[2, -2]) / 2. / 1j

    # -(zx+izy) - (zx-izy) = -zx-izy-zx+izy = -2zx
    xz = (t[2, 1] - t[2, -1]) / (-2.)
    # -(zx+izy) + (zx-izy) = -zx-izy+zx-izy = -2izy
    yz = -(t[2, 1] + t[2, -1]) / (-2.) / 1j

    xx_plus_yy = -zz

    xx = (xx_plus_yy + xx_minus_yy) / 2.
    yy = (xx_plus_yy - xx_minus_yy) / 2.

    return from_natural_form(xx, xy, xz, yy, yz)


## |1, m>, Cartesian
# \returns
# \f$ T^{(1)} = \sqrt{\frac{1}{3}} \left[
# \begin{array}{ccc} 0 &  \sqrt{2}T_0^1 & -i(T^{1}_1+ T^{1}_{-1}) \\
#  -\sqrt{2}T_0^{1} & 0  & -(T^{1}_1 - T^{1}_{-1})  \\
# i(T^{1}_ 1+T^{1}_{-1}) & (T^{1}_1 - T^{1}_{-1})&0\end{array} \right] \f$
# @image html T1m.tiff
# @image latex T1m.tiff
def T1(self):
    """Vector Components, as a 3 x 3 (Antisymmetric) Cartesian Tensor."""
    from ...euclid import vector
    from .one import ROOT2
    x = -(self[1, 1] - self[1, -1])/ROOT2
    y = -1j*(self[1, 1] + self[1, -1])/ROOT2
    z = self[1, 0]
    return vector.Vector(x, y, z).dual() *1j/ROOT2 
    
    
## |0, 0>, Cartesian.
# \returns
# \f$ T^{(0)} = \frac{1}{3}\left[
# \begin{array}{ccc}
# T_0^{0} &  0 & 0  \\
#  0 & T_0^{0} &  0   \\
#  0 & 0 & T_0^{0} \end{array} \right] \f$
# @image html T00.tiff
# @image latex T00.tiff
def T0(self):
    """Scalar Components, as a 3 x 3 Cartesian Tensor ~ Identity."""
    from ...euclid import tensor

    value = self[0][0]
    try:
        len(value)
    except TypeError:
        pass
    else:
        print "NEEDED 3 gets:", type(self)
        value = value[0]
    
    return tensor.DELTA * value / SCALAR_NORM

from .one import EPLUS, EMINUS, EZERO
ep, em, ez = EPLUS, EMINUS, EZERO

BASIS = {(2, 2): EPLUS & EPLUS,
         (2, 2): EMINUS & EMINUS,
         (0, 0): (-((ep&em) + (em&ep) - (ez&ez))/pow(3,0.5))}
