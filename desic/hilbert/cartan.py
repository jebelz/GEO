"""Spinors: what's a spinor. It's something you have to get used to. It
is the 'i' (or 1j in python) to Vector's -1, to Tensor's 1.

Like a tensor is the even outer product of 2 (odd--under reflection) vectors,
so the vector is the odd outer product of 2 odd spinors.


As Cartesian objects, they don't really have an intuitive description.
If a vector is something with "a magnitude and a direction", what's the square
root of that? It's something with a root-magnitude, a direction, and a phase.
Still: not intuitive.

Well it's better to think of a vector as something defined by
how it transforms-- and that's even easier in the spherical basis.

Recall (eckart) vectors in the spherical basis |j=1, m>

|1, -1> is |-1> is |->
|1, 0>  is |0>
|1, 1>  is |1> is |+>

where the |m> are eigenstates of Z-rotations (by choice of an arbitrary axis)
with eigenvalue exp(1j * m * theta).

Likewise we can write down the spinor basis states (j=1/2):

|+1/2>   or  "u"
|-1/2>   or  "d"


They too are eigenstates of rotations with eigenvalue exp(1j * m * theta).
Note that "u", or up, is not "-d", the negative of down. In fact, the
eigenvalue equation shows that rotating d by 360 degrees results in -d: you
have to rotate them 720 degrees to remain invariant.

So it is established that they are not very intuitive.

They are transformed by 2 x 2 complex matrices (SU(2), aka the Pauli
Matrices)-- which is equivalent to multiplication by a quaternion, Q, so
that:

Q*u

rotates u by theta, about some axis n, where

Q = {cos(theta/2), n * sin(theta/2)}

This is why it takes 2 quaternions to rotate a vector, e.g.:

uu --> |1>  (racah.py)

so:

QuQu --> Quu(~Q) --> Q|+>(~Q).

And hence the quaternion works.


Recall:

|+> = -(x + iy)/sqrt(2) = uu
|-> =  (x - iy)/sqrt(2) = dd

so that u *is* the square root of |+>, and likewise for d and |->.

The other normalized symmetric combination of base spinors is:

(ud + du)/sqrt(2) = |1, 0> = z (plane old z-hat contains both spinors).


Finally, the antisymmetric combo is:

(ud - du)/sqrt(2) = w... it's a scalar. Even "1" has spinors in it.


So, the tensor product of 2 spinors spaces is:

2 x 2 = 3 + 1, and that is a quaternion:

w     ud - du
x     uu - dd
y    (uu + dd)/i
z     ud + du

So that's a core dump in spinors in R3. They exist, they're real (c.f.
electrons).

The real kicker is that they're Lorentz invariant, and can be kicked up to
Minkowski space, just like that. So in a 3D Newtonian world that is
also a 3 + 1 Minkowski spacetime, well, they're fundamental."""

## \namespace geo.desic.hilbert.cartan 2-componet Spinors.
import fractions
import operator

#from geo.utils.trig import exp, pi
from geo.metric.wigner import racah


## Complexifier:
# \f$ z = a + ib\f$
def C(a, b):
    """Complexify: C(a, b) --> z"""
    return a + 1j*b


## Quaternion-Spinor Rotation (with a humorous Matlab-like name).
# \param q Quaternion
# \param s Spinor
# \returns Spinor
def qsrot(q, s):
    """qsrot: right mul of a spinor by a quaternion"""
    u, d = s.iter()
    w, x, y, z = q.tolist()
    return type(s)(C(w, -z)*u - C(y, x)*d,
                   C(y, -x)*u + C(w, z)*d)


## Spinor Inner Product
# \param s1 Spinor \f$ s_1 \f$
# \param s2 Spinor \f$ s_2 \f$
# \returns Scalar \f$ \langle s_1 | s_2 \rangle \f$
def inner(s1, s2):
    """<s1|s2>"""
    from ...metric.euclid.scalar import Scalar
    return Scalar(
        reduce(operator.add,
               map(operator.mul,
                   s1.C.iter(),
                   reversed(list(s2.iter())))))


dot = inner


## Spinor Outer Product
# \param s1 Spinor \f$ s_1 \f$
# \param s2 Spinor \f$ s_2 \f$
# \returns Quaternion \f$ |s_1 \rangle \langle s_2 | \f$
def outer(s1, s2):
    """quaternion = |s1><s2|"""
    from ..hamilton import Quaternion
    m = [[s1.u * s2.u, s1.u * s2.d],
         [s1.d * s2.u, s1.d * s2.d]]
    # this piece of code is WET: see versor.frompsinor.
    w = (m[0][0] + m[1][1]) / 2.
    z = (-m[0][0] + m[1][1]) / 2j

    x = -(m[0][1] + m[1][0]) / 2j
    y = (-m[0][1] + m[1][0]) / 2.
    return Quaternion.fromwxyz(w, x, y, z)


## Dilation
def dilation(xi, alpha):
    """dilation(xi, alpha)"""
    return type(xi)(*[item*alpha for item in xi.iter()])


## \f$ {\bf \xi} = \xi e^{-\alpha/2} \left[ \begin{array}{c}
# \cos{\frac{\theta}{2}}e^{-i\frac{\phi}{2}} \\
# \sin{\frac{\theta}{2}}e^{i\frac{\phi}{2}}
# \end{array}  \right] \f$
class Spinor(object):
    """(Weyl) Spinor(u, d)"""

    rank = fractions.Fraction(1, 2)

    ## Construct from a vector plus a phase
    # \param v A vector (Polar or not)
    # \param alpha A phase
    # \returns \f$ \left[\begin{array}{c}
    # \sqrt{r}\cos{\frac{\theta}{2}}e^{i(-\alpha-\phi)/2} \\
    # \sqrt{r}\sin{\frac{\theta}{2}}e^{i(-\alpha+\phi)/2} \end{array}\right]f$
    @classmethod
    def fromvector(cls, v, alpha=0):
        """spinor = Spinor.fromvector(v [, alpha=0])"""
        from geo.utils.trig import exp, cos, sin, sqrt
        p = v.polar()
        r = sqrt(p.radius)
        a = r * cos(p.theta/2) * exp(1j*(-alpha - p.phi)/2.)
        b = r * sin(p.theta/2) * exp(1j*(-alpha + p.phi)/2.)
        return cls(a, b)

    ## TBD
    def tocart(self):
        """Not working"""
        from ..hamilton import Quaternion
        return Quaternion.fromwxyz(self.t, self.x, self.y, self.z)


    ## Two component constructor.
    # \param u \f$ \psi^+ \f$
    # \param d \f$ \psi^- \f$
    def __init__(self, u, d):
        ## \f$ \left[\begin{array}{c}1\\0\end{array}\right] \f$ component.
        self.u = u
        ## \f$ \left[\begin{array}{c}0\\1\end{array}\right] \f$ component.
        self.d = d


    ## \returns \f$ u\bar{d} + d\bar{u} \f$
    @property
    def x(self):
        """x projection."""
        return self.u*self.d.conjugate() + self.d*self.u.conjugate()

    ## \returns \f$ i(u\bar{d} - d\bar{u}) \f$
    @property
    def y(self):
        """y projection."""
        return 1j*(self.u*self.d.conjugate() - self.d*self.u.conjugate())

    ## \returns \f$ ||u||^2 - ||d||^2 \f$
    @property
    def z(self):
        """z projection."""
        return abs(self.u)**2 - abs(self.d)**2

    ## \returns \f$ ||u||^2 + ||d||^2 \f$
    @property
    def t(self):
        """t projetion"""
        return abs(self.u)**2 + abs(self.d)**2

    def __str__(self):
        return " | %s | \n | %s | " % (str(self.u), str(self.d))

    __repr__ = __str__

    ## Obvioulsy I don't understand the difference between
    # \f$ {\bf 2} \f$ and \f$ {\bf \bar{2}} \f$.
    def __invert__(self):
        return type(self)(self.d.conjugate(), self.u.conjugate())
    #Spinor(self.u.conjugate(), self.d.conjugate())

    ## \returns \f$ \bar{xi} \f$
    @property
    def C(self):
        """~"""
        return ~self

    ## \returns \f$ \sqrt{||u||^2 + ||d||^2} \f$
    def __float__(self):
        return (abs(self.u)**2 + abs(self.d)**2) ** 0.5

    ## \returns \f$ \sqrt{|u|^2 + |d|^2} \f$
    def __complex__(self):
        return (self.u**2 + self.d**2) ** 0.5

    ## \returns \f$ +\xi =
    #  \left[\begin{array}{c}\xi^+ \\ \xi^-\end{array}\right] \f$
    def __pos__(self):
        return type(self)(*map(operator.pos, self.iter()))

    ## \returns \f$ =\xi =
    #  \left[\begin{array}{c}-\xi^+ \\ -\xi^-\end{array}\right] \f$
    def __neg__(self):
        return type(self)(*map(operator.neg, self.iter()))

    ## \returns \f$ \sqrt{\xi\bar{\xi}} \f$
    def __abs__(self):
        return abs(self*self)**0.5

    ## \returns \f$ \psi + \xi =
    #  \left[\begin{array}{c}\psi^+ + \xi^+ \\
    #   \psi^- +  \xi^-\end{array}\right] \f$
    def __add__(self, other):
        return type(self)(*[operator.add(left, right) for left, right in
                            zip(self.iter(), other.iter())])

    ## \returns \f$ \psi - \xi =
    #  \left[\begin{array}{c}\psi^+ - \xi^+ \\
    #   \psi^- -  \xi^-\end{array}\right] \f$
    def __sub__(self, other):
        return self + (-other)

    ## \returns \f$ \alpha\xi =
    #  \left[\begin{array}{c}\alpha\xi^+  \\
    #  \alpha\psi^- \end{array}\right] \f$
    def __rmul__(self, alpha):
        return dilation(self, alpha)

    def __div__(self, alpha):
        return dilation(self, 1./alpha)

    ## This is confusing, as the outer product could be a |ket><bra|,
    # which makes a 2 x 2 pauli like object, or could be a
    # |ket 1> x |ket 2> --> |1,2> which might be a vector: this is
    # the latter.
    # \return Quaternion Spinor outer product has a Scalar and Vector part.
    def __and__(self, other):
        return self.ket()*other.ket()

    ## Ket version of a spinor
    # \returns metric.euler.wigner.racah.Ket
    def ket(self):
        """Convert to a racah.Ket"""
        return racah.ITS([a*k for a, k in zip(self.iter(), racah.SPINOR) if a])


    def __mul__(self, other):
        from ..hamilton import Quaternion
        if isinstance(other, type(self)):  # todo: fix
            return dot(self, other)
        elif Quaternion:
            print "this total screws u quaternion__mul__..."

    ## Iter components
    # \yields u and d
    def iter(self):
        """yield u and d."""
        yield self.u
        yield self.d


U = Spinor(1, 0)
D = Spinor(0, 1)

#uu = outer(~U, U)
#dd = outer(~D, D)

#ud = outer(~U, D)
#du = outer(~D, U)

## Noraml
#x = (ud + du)/1j
#y = (du - ud)
#z = (uu - dd)/1j

#w = (uu + dd)

## Conjuagte
#y = -(ud + du)
#x = -(du - ud)/1j
#w = (uu - dd)
#z = (uu + dd)/1j
