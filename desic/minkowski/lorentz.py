"""Lorentz Transformations"""

## \namespace geo.desic.minkowski.lorentz Lorentz transformations

## THE <a href="http://en.wikipedia.org/wiki/Minkowski_space#Standard_basis">
# Metric Tensor </a>
# \returns g_uv \f$ g_{\mu\nu} \equiv -\eta_{\mu\nu}\f$ (a Lorentz)
def g_uv():
    """Compute (1, -1, -1, -1) eta-metric"""
    from ...metric.euclid.scalar import ONE
    from ...metric.euclid.vector import NULL
    from ...metric.euclid.tensor import DELTA
    return Lorentz(ONE, NULL, NULL, -DELTA)


## A <a href="http://en.wikipedia.org/wiki/Four-tensor">4-tensor</a>
class Lorentz(object):
    """Lorentz(tt, ts, st, ss)

    t --> time-like    (3-scalar)
    s --> space-like   (3-vector)
    """
    components = ("tt", "ts", "st", "ss")

    ## \f$ \Lambda_{\mu\nu}\f$
    # \param tt metric.euclid.scalar.Scalar time-time component
    # \param ts metric.euclid.vector.Vector time-space components
    # \param st metric.euclid.vector.Vector space-time components
    # \param ss metric.euclid.tensor.Tensor space-space components
    def __init__(self, tt, ts, st, ss):
        ## \f$ \Lambda_{00} \f$
        self.tt = tt
        ## \f$ \Lambda_{0j} \f$
        self.ts = ts
        ## \f$ \Lambda_{i0} \f$
        self.st = st
        ## \f$ \Lambda_{ij} \f$
        self.ss = ss

    def __str__(self, root=" | {} | %s %s %s | "):
        r = [root.format(self.tt.w) % tuple(map(str, self.ts.iter()))]
        for i in range(3):
            r.append(
                root.format(self.st.e(i)) % tuple(map(str, self.ss.e(i).iter()))
            )
            continue
        return "\n".join(r)

    ## left * right --> right.rmul(left)
    # \param self A Lorentz object
    # \param other A Lorentz or minkowski.FourVector
    # \retval self*other type(other) via other.rmul
    def __mul__(self, other):
        return other.rmul(self)

    ## \f$ \Lambda''_{\mu\nu} =\Lambda'_{\mu\sigma}\Lambda_{\sigma\lambda} \f$
    # \param self A Lorentz object
    # \param lam Another Lorentz object
    # \retval lam*self A 3rd Lorentz object.
    def rmul(self, lam):
        """
        |TT' TS'|    |TT   TS| |tt  ts|    |TT*tt + TS*st    TT*ts + TS*ss|
        |       | =  |       | |      | = |                               |
        |ST' SS'|    |ST   SS| |st  ss|    |ST*tt + SS*st    ST&TS + SS*ss|
        """
        return type(self)(lam.tt*self.tt + lam.tt*self.st,
                          lam.tt*self.ts + lam.ts*self.ss,
                          lam.st*self.tt + lam.ss*self.st,
                          (lam.st & lam.ts) + lam.ss*self.ss)

    ## NotImplemented
    def __invert__(self):
        from ...metric.euclid.scalar import Scalar
        from ...metric.euclid.vector import Vector
        from ...metric.euclid.tensor import Tensor


        if 1:
            raise NotImplementedError("Sorry...")

        m00 = self.tt
        m01, m02, m03 = self.ts.iter()
        m10, m20, m30 = self.st.iter()
        m11, m12, m13, m21, m22, m23, m31, m32, m33 = self.ss.iter()

        t00 = (m12*m23*m31 - m13*m22*m31 + m13*m21*m32 -
               m11*m23*m32 - m12*m21*m33 + m11*m22*m33)
        t01 = (m03*m22*m31 - m02*m23*m31 - m03*m21*m32 +
               m01*m23*m32 + m02*m21*m33 - m01*m22*m33)
        t02 = (m02*m13*m31 - m03*m12*m31 + m03*m11*m32 -
               m01*m13*m32 - m02*m11*m33 + m01*m12*m33)
        t03 = (m03*m12*m21 - m02*m13*m21 - m03*m11*m22 +
               m01*m13*m22 + m02*m11*m23 - m01*m12*m23)
        t10 = (m13*m22*m30 - m12*m23*m30 - m13*m20*m32 +
               m10*m23*m32 + m12*m20*m33 - m10*m22*m33)
        t11 = (m02*m23*m30 - m03*m22*m30 + m03*m20*m32 -
               m00*m23*m32 - m02*m20*m33 + m00*m22*m33)
        t12 = (m03*m12*m30 - m02*m13*m30 - m03*m10*m32 +
               m00*m13*m32 + m02*m10*m33 - m00*m12*m33)
        t13 = (m02*m13*m20 - m03*m12*m20 + m03*m10*m22 -
               m00*m13*m22 - m02*m10*m23 + m00*m12*m23)
        t20 = (m11*m23*m30 - m13*m21*m30 + m13*m20*m31 -
               m10*m23*m31 - m11*m20*m33 + m10*m21*m33)
        t21 = (m03*m21*m30 - m01*m23*m30 - m03*m20*m31 +
               m00*m23*m31 + m01*m20*m33 - m00*m21*m33)
        t22 = (m01*m13*m30 - m03*m11*m30 + m03*m10*m31 -
               m00*m13*m31 - m01*m10*m33 + m00*m11*m33)
        t23 = (m03*m11*m20 - m01*m13*m20 - m03*m10*m21 +
               m00*m13*m21 + m01*m10*m23 - m00*m11*m23)
        t30 = (m12*m21*m30 - m11*m22*m30 - m12*m20*m31 +
               m10*m22*m31 + m11*m20*m32 - m10*m21*m32)
        t31 = (m01*m22*m30 - m02*m21*m30 + m02*m20*m31 -
               m00*m22*m31 - m01*m20*m32 + m00*m21*m32)
        t32 = (m02*m11*m30 - m01*m12*m30 - m02*m10*m31 +
               m00*m12*m31 + m01*m10*m32 - m00*m11*m32)
        t33 = (m01*m12*m20 - m02*m11*m20 + m02*m10*m21 -
               m00*m12*m21 - m01*m10*m22 + m00*m11*m22)

        return Lorentz(Scalar(t00),
                       Vector(t01, t02, t03),
                       Vector(t10, t20, t30),
                       Tensor(t11, t12, t13,
                              t21, t22, t23,
                              t31, t32, t33))

    #pylint: enable=R0914
    ## Transpose
    @property
    def T(self):
        """Transpose"""
        return type(self)(self.tt, self.st, self.ts, self.ss.T)


    ## The <a href="http://en.wikipedia.org/wiki/Lorentz_transformation">
    # Lorentz Transformation</a> of four-vectors.
    # \param v \f$ v_{\mu} \f$
    # \returns  \f$ v'_{\mu} = \Lambda_{\mu}^{\nu} \ v_{\nu} \f$
    def __call__(self, v):
        """        |TT   TS||t|
        (t', s') = |       || |
                   |ST   SS||s|
        """
        return type(v)(self.tt*v.t - self.ts*v.s,
                       self.st*v.t + self.ss*v.s)

    ## Dual
    # \return \f$ \epsilon_{\mu\nu\sigma\lambda}\Lambda^{\sigma\lambda} \f$
    def dual(self):
        """Dual"""
        return type(self)(self.tt, -self.st, -self.ts, self.ss.T)
