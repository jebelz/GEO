"""Mobius Transformations"""
## \namespace geo.desic.plucker.mobius <a
# href="http://en.wikipedia.org/wiki/Mobius_transformation">Mobius Transform
# </a> GL(2, C) or PGL(2, C).


## Mobius Transform: PGL(2, C).
# @image html mobius.jpeg
# @image latex mobius.jpeg
class Mobius(object):
    """f(z) = Mobius(a, b, c ,d)

              a*z + b
      f(z) = ---------
              c*z + d

      (f * g)(z) = f(g(z))

      (f * ~f)(z) = z
    """

    ## Contruct from
    # \f$ \left[
    # \begin{array}{cc} a & b \\ c & d \end{array}
    # \right]
    # \f$
    # \param m A 2 x 2 number matrix
    # \returns Mobius instance
    @classmethod
    def frommatrix(cls, m):
        """Build from a numpy 2 x 2:
        [[a, b],
        [c, d]]"""
        return cls(m[0, 0], m[0, 1], m[1, 0], m[1, 1])

    ## Construct from 3 fixed points  via \f$ \left[\begin{array}{cc} z_2-z_3 & -z_1(z_2-z_3) \\z_2-z_1 & -z_2(z_2-z_1) \end{array} \right] \f$
    # \param z1 \f$ z_1\f$
    # \param z2 \f$ z_2\f$
    # \param z3 \f$ z_3\f$
    # \returns Mobius \f$ \left[\begin{array}{cc} z_2-z_3 & -z_1(z_2-z_3) \\z_2-z_1 & -z_2(z_2-z_1) \end{array} \right] \f$
    @classmethod
    def from3(cls, z1, z2, z3):
        return cls(z2-z3, -z1(z2-z3), z2-z1, -z3(z2-z1))

    ## \f$ f(z) = \frac{az+b}{cz+d} \f$
    # \param a
    # \param b
    # \param c
    # \param d
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def __str__(self):
        """
                %s*z + %s
        f(z) =  %s
                %s*z + %s
        """
        a, b, c, d = map(str, (self.a, self.b, self.c, self.d))

        m = 5 + max(len(a) + len(b), len(c) + len(d))
        l = '-' * m
        return type(self).__str__.__doc__ % (a, b, l, c, d)

    ## Yield Serial Components: \n
    # \yields Translation: \f$ f_0(z) = z+\frac{d}{c} \f$ \n
    # \yields Inversion+Reflection: \f$ f_1(z) = \frac{1}{z} \f$ \n
    # \yields Homothety+Rotation: \f$ f_2(z) = \frac{ab-cd}{c^2}z  \f$  \n
    # \yields Translation: \f$ f_3(z) = z + \frac{a}{c} \f$
    def __iter__(self):
        # Translation
        yield type(self)(1, self.d/self.c, 0, 1)
        # Inversion and Reflection wrt real axis
        yield type(self)(0, 1, 1, 0)
        # Homothety and Rotation
        yield type(self)(-abs(self), 0, 0, self.c**2)
        # Translation
        yield type(self)(1, self.a/self.c, 0, 1)

    ## Det
    def __abs__(self):
        return self.det()

    ## Tr
    def __complex__(self):
        return self.trace()

    ## Determinant
    # \returns  \f$ ||f|| = ab-cd \f$
    def det(self):
        """Determinant."""
        return self.a * self.d - self.b * self.c

    ## Trace
    # \returns \f$ a+d \f$
    def trace(self):
        """Trace"""
        return self.a + self.d

    ## Not Singular
    # \returns bool(\f$ ||f(z)||\f$ )
    def __nonzero__(self):
        return bool(self.det())

    ## Inverse
    # \returns \f$ f^{-1}(z) = \frac{dz-b}{-cz+a} \f$
    def __invert__(self):
        return type(self)(self.d, -self.b, -self.c, self.a)

    ## Apply:
    # \return \f$  \frac{az+b}{cz+d} \f$
    def __call__(self, z):
        return (self.a * z + self.b)/(self.c * z + self.d)

    ## Matrix Form
    # \returns \f$ \left[\begin{array}{cc}a&b\\c&d\end{array}\right]\f$
    def __array__(self):
        from numpy import matrix
        return matrix([[self.a, self.b], [self.c, self.d]])

    def __div__(self, g):
        return type(self)(self.a/g, self.b/g, self.c/g, self.d/g)

    def __rmul__(self, g):
        return self.__div__(1./g)

    ## Composition:
    # \returns \f$ h(z) = f(z)g(z) = \frac{(aa'+bc')z+(ab'+bd')}{(ca'+dc')z+
    # (cb'+dd')} \f$
    def __mul__(self, other):
        return type(self)(
            self.a * other.a + self.b * other.c,
            self.a * other.b + self.b * other.d,
            self.c * other.a + self.d * other.c,
            self.c * other.b + self.d * other.d
            )
