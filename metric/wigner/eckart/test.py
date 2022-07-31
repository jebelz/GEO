import operator
import geo
from geo.metric.wigner.eckart import one
from geo.metric.wigner.eckart import two
from geo.metric.wigner.eckart import three

from geo.metric.wigner import racah


n, p, m = racah.VECTOR

w = operator.methodcaller('spherical')
c = operator.methodcaller('tocart')

wc = lambda dum: w(c(dum))
cw = lambda dum: c(w(dum))

X, Y, Z = geo.BASIS
x, y, z = map(w, geo.BASIS)

XX = X & X
YY = Y & Y
ZZ = Z & Z

xx = x & x
yy = y & y
zz = z & z



z0, p0, m0 = one.BASIS
z1, p1, m1 = map(w, one.BASIS)


u0, v0 = geo.Vector.random(2)
u1, v1 = [item.spherical() for item in (u0, v0)]


u_dot_v0 = u0 * v0
u_dot_v1 = u1 * v1


u_x_v0 = u0 ^ v0
u_x_v1 = u1 ^ v1

uv0 = u0 & v0
uv1 = u1 & v1


zz0 = (z0 & z0)
pp0 = (p0 & p0)
mm0 = (m0 & m0)

zz1 = z1 & z1
pp1 = p1 & p1
mm1 = m1 & m1


def test1():
    """Create random cartesian vectors and test all ops"""
    
    c1, c2 = geo.Vector.random(2)
    s1, s2 = map(w, (c1, c2))
    d1, d2 = map(c, (s1, s2))
    print "cart->sph-->cart:\n\t", abs(d1-c1), abs(d2-c2)


    print "dot:\t", complex(c1*c2) - complex(s1*s2)

    print "cross:\t", abs((c1^c2) - c(s1^s2))
    
    tc = c1 & c2
    sc = s1 & s2

    print "Outer:", abs(tc-c(sc)), abs(w(tc)-sc)


def test2():
    t1 = geo.Tensor.random(1)[0]
    s1 = t1.spherical()
    print (t1-s1.tocart()).clean()
    
