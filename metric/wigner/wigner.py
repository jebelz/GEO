"""Wigner d, D matrices:

dj_mm
Dj_mm

are curried function of angles.

D is the (2j+1)**2 matrix of fucntions, which supports tensor products
and sums:

D(j1) * D(j2) --> D(j1+j2) + D(|j1+j2|)


"""
## \namespace geo.metric.wigner.wigner Wigner Matrices
# \f$ D^j_{m'm}(\alpha, \beta, \gamma)\f$
import operator


from ...utils.trig import cos, sin
from ..schur_weyl.pascal import Binomial


## \f$ d^j_{m'm}(\beta)\f$
class dj_mm(object):
    """dj_mm(j, m1, m)(beta)"""

    ## \param j
    # \param m1
    # \param m
    def __init__(self, j, m1, m):
        ## \f$ j\f$
        self.j = j
        ## \f$ m'\f$
        self.m1 = m1
        ## \f$ m\f$
        self.m = m
        ## Solve for a, l.
        self.a, self.l = self._k_hash()[self.k]
        ## \f$ b = 2j - 2k -a \f$
        self.b = 2*self.j - 2*self.k - self.a
        ## Jacobi polynomial
        self.jacobi = Binomial(self.k).jacobi(self.a, self.b)


    ## \f$ k = j + min(m, -m, m', -m') \f$
    @property
    def k(self):
        """k = j + min(m, -m, m1, -m1)"""
        return self.j + min(self.m, -self.m, self.m1, -self.m1)

    ## Find boundaries.
    def _k_hash(self):
        m1_m = self.m1 - self.m
        return {self.j + self.m: (m1_m, m1_m),
                self.j - self.m: (-m1_m, 0),
                self.j + self.m1: (-m1_m, 0),
                self.j - self.m1: (m1_m, m1_m)}

    ## \f$ d^j_{m'm}(\beta) = (-1)^{\lambda}
    # {2j-k \choose k+a}^{+\frac{1}{2}}
    # {k+b \choose b}^{-\frac{1}{2}}
    # \sin^a{\frac{\beta}{2}}
    # \cos^b{\frac{\beta}{2}}
    # P_k^{a, b}(\cos{\beta}) \f$
    def __call__(self, beta):
        """(beta)"""
        return (
            operator.pow(-1, self.l) *
            operator.pow(Binomial(2*self.j-self.k)(self.k+self.a), 0.5) *
            operator.pow(Binomial(self.k+self.b)(self.b), -0.5) *
            operator.pow(sin(beta/2), self.a) *
            operator.pow(cos(beta/2), self.b) *
            self.jacobi(cos(beta))
            )


## \f$ D^j_{m'm}(\alpha, \beta, \gamma)\f$
class Dj_mm(dj_mm):
    """Dj_mm(j, m1, m)(alpha, beta, gamma)"""

    ## \f$ D^j_{m'm}(\alpha, \beta, \gamma) = e^{-im'\alpha}d^j_{m'm}(\beta)
    # e^{-im\gamma} \f$
    def __call__(self, alpha, beta, gamma):
        """(alpha, beta, gamma) --dIffers from super."""
        #pylint: disable=W0221
        import scipy as sp
        return (
            sp.exp(-1j * self.m1 * alpha) *
            super(Dj_mm, self).__call__(beta) *
            sp.exp(-1j * self.m * gamma)
            )


## \f$ {\bf \hat{D}}^j \f$
class D(object):
    """Irreducible Wigner Big D Matrix"""

    def __init__(self, j):
        self.j = j

    ## \f$ {\bf \hat{D}}^{j_1} \otimes  {\bf \hat{D}}^{j_2} =
    #  {\bf \hat{D}}^{j_1+  j_2} \oplus {\bf \hat{D}}^{j_1 + j_2- 1} \oplus
    # \cdots  {\bf \hat{D}}^{|j_1 - j_2|} \f$
    def __mul__(self, other):
        return reduce(operator.add, ())

    ## \f$ {\bf \hat{D}}^{j_1} \oplus  {\bf \hat{D}}^{j_2} \f$
    def __add__(self, other):
        return DSum(self, other)

    def __array__(self):
        """To: solve index convention, and object mode"""
        return NotImplemented


## Tensor Sum of Irreducible Representations (TBD)
class DSum(object):
    """DSum is tensor sum rep- i.e, block diagonal."""

    def __init__(self, *js):
        self.js = list(reversed(sorted(js)))

    def __array__(self):
        return "Kroneckor Product of Matricies/or block diagonal sum: TND"
