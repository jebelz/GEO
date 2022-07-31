"""We have Four Coordinates:


                                            Transform   Ortho-  Units
Class  Cartesian  Origin  Conjugate  OTF     Method     Normal
-------------------------------------------------------------------------
ECEF |   Yes    | Fixed  |   LLH   | ECEF | p.ecef()  |  Yes | m, m, m   |
LLH  |    No    | Fixed  |   ECEF  | ECEF | p.llh()   |No Way|deg, deg, m|
LTP  |   Yes    | Pegged |   SCH   | LTP  | p.ltp(peg)|  Yes | m, m, m,  |
SCH  |    No    | Pegged |   LTP   | LTP  | p.sch(peg)|   No | m, m, m   |
--------------------------------------------------------------------------

OTF = Orthogonal Tangent Frame (all automated vector ops occur in this frame)

For any coordinate "p", the transform method returns the SAME AFFINE point(s)
in the new class, ON THE SAME ELLIPSOID. Note the Ellipsoid is None, and
there is no method for you to supply it. You must start out with an
ellipsoid. So,

ECEF = wgs84.WG484.ECEF
LLH = wgs84.WG484.LLH

by default, but you can make you own ellipsoid, which will, boom, have its
own coordinates:

>>>ellipsoid.Ellipsoid(6000e3, 300.).ECEF   # etc, etc.

Let's assume WGS84 for now:

>>>p = LLH(34.201694, -118.171667, 345.)   # Boom, JPL.

>>> print p.ecef()
x:-2493238.95577
y:-4655397.3425
z:3565166.28058

Now make a Peg point (at the corner of JPL's SRTM cell):

>>>peg = Peg(34, -119, 90)  # 90 means East North Up

>>>pxyz =  p.ltp(peg=peg)  # Local Tangent Plane at SRTM Corner
>>>print pxyz
x:76345.4923882
y:22682.5221492
z:-151.908686075

>>>psch =  pxyz.sch()  # Local Tangent Sphere at SRTM Corner
>>>print psch
s:76343.6703336
c:22681.344927
h:344.818865087

Question: how far is JPL from the corner?

>>>print abs(psch)   # same as abs(pxyz)
geo.metric.euclid.scalar.Scalar(79643.9206424)

How about from the center of the Earth:
>>>print abs(p)
geo.metric.euclid.scalar.Scalar(6371763.92495)

THE POINT: Vector is taken relative to the origin.

Can we add a vector to JPL? Yes we can:

print psch + 1000*Y   # move 1km in the peg's Y direction.
s:76343.6703336
c:23681.28405
h:348.449547232
(34, -119, 90)

Can we subtract a vector? No, we cannot. But we can subtract another point,
here andother corner:

>>>print psch - LLH(33, -118, 0)
[-17102.9778509, 133134.261659, 1492.03312464]

The result is in the left operands OTF. So in ECEF:

 print p - LLH(33,-118,0)
[20534.6903323, 72323.2816167, 111207.639445]


NOTE ON PEGS: If no PEG is supplied, the existing (self's) peg is used. So
to from SCH to SCH' (here NWU):

print psch.sch(peg=Peg(34, -118, 0))
s:22386.0707412
c:15821.825881
h:345.091452248
(34, -118, 0)


SO: All coordinates behave the same and are simply transformed based on the same
method calls (it knows). Any vector operations are done in the OTF.

So what If you want to transform a vector in the tangent space (e.g. velocity),
or have something represented in non-linear co/contra-variant basis?
Well, this is a Cartesian package, there is no s-hat, c-hat, h-hat vector
with an operator overload contractin on the correct metric tensor.
You need to do it yourself, but here are the tools:


What are the vectors tangent to constant s, c, h:
>>>print psch.basis_covariant()
Vector(0.999976207848, 0.0, -0.0119573385009)
Vector(-4.24772387383e-05, 1.00004769597, -0.00355231459829)
Vector(0.0119566927675, 0.0035523767029, 0.999922206033)

What are the cotangent (gradients) s, c, h:
>>>print pxyz.basis_covariant()
(etc)

How about the normalized transformation of s, c, h to x, y, z:

>>>print psch.local_basis()
[0.999928515292, 0.0, -0.0119567682113]
[-4.24749448357e-05, 0.99999369029, -0.00355212276229]
[0.0119566927675, 0.0035523767029, 0.999922206033]

Since we used a pegged coordinate, that was relative to the peg. To get
it relative to ECEF, us the LLH instance:

>>>print p.basis_covariant()
Vector(29439.0914536, 54968.9261839, 91747.9047636)
Vector(81252.0116152, -43515.2288172, 0.0)
Vector(-0.390469215704, -0.729087496772, 0.562107830969)

The 1st 2 tangents are big, as they represent a 1 degree change.

There are 3 metric tensors:

print p.metric_tensor_covariant()
[12305920979.9, -4.76837158203e-07, 0.0]
[-4.76837158203e-07, 8495464530.52, 0.0]
[0.0, 0.0, 1.0]

In [47]: print p.metric_tensor_contravariant()
[8.1261695214e-11, -3.23117426779e-27, 0.0]
[-3.23117426779e-27, 1.17709867001e-10, 4.23516473627e-22]
[0.0, 4.23516473627e-22, 1.0]

--with metric_tensor_mixed being the kornecker delta (computed).

There is also a jacobian to the conjugate frame:

In [48]: print psch.J
[0.999976207848, -4.24772387383e-05, 0.0119566927675]
[0.0, 1.00004769597, 0.0035523767029]
[-0.0119573385009, -0.00355231459829, 0.999922206033]

In [49]: print pxyz.J
[0.999880825011, 0.0, -0.0119561979489]
[-4.24726510571e-05, 0.999939687531, -0.00355193093665]
[0.0119566927675, 0.0035523767029, 0.999922206033]

You can also get Christoffel symbols, If you need them.
"""

## \namespace geo::detic  Geodesy
from .origins import *
from .newton import *
from .dms import *
