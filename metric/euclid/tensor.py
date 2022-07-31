"""The Rank-2 Tensor, and its functions:

Tr
det
adj
inv
norm

Some helpers:
ziprows, zipcols

and constants:

DELTA
DYADS = (XX, XY, XZ, ..., ZZ)


Note:
Tensors are both tensors and functions from:

R3 --> R3
R3 x R3 -->R

and rotation matrices.
"""
## \namespace geo::metric::euclid::tensor Rank-2 Tensor classes and functions.

import operator
import itertools

from ...utils.arraylike import enlister, discharger, matrix_wrap
from ...utils.trig import degrees, arctan2, arccos, arcsin
from . import euclid
from .vector import BASIS as _vector_basis

__all__ = ('Tr', 'det', 'adj', 'inv', 'cof', 'hyd', 'dev', 'Tensor',
           'skew', 'symm', 'posterior_product', 'from_natural_form',
           'ziprows', 'zipcols', 'DELTA', 'DYADS', 'XX', 'XY', 'XZ',
           'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ')


## Build a Symmetric Tensor (ala Poinsot's Ellipsoid)
# \param xx \f$ S_{xx} \f$
# \param xx \f$ S_{yy} \f$
# \param xx \f$ S_{zz} \f$
# \param xx \f$ S_{yz} = S_{zy} \f$
# \param xx \f$ S_{zx} = S_{xz} \f$
# \param xx \f$ S_{xy} = S_{yx} \f$
# \param xx \f$ S_{yz} = S_{yz} \f$
# \returns \f$ \left[\begin{array}{ccc}
# S_{xx} & S_{xy} & S_{zx} \\
# S_{xx} & S_{xy} & S_{zx} \\
# S_{xx} & S_{xy} & S_{zx} \end{array} \right] \f$
def poinsot(xx, yy, zz, yz, zx, xy):
    """T = pointsot(xx, yy, zz, yz, zx, xy) is a symmetric tensor"""
    return Tensor(xx, xy, zx,
                  xy, yy, yz,
                  zx, yz, zz)
    

## Build from 5 unique components
# \param xx \f$ S_{xx} \f$
# \param xx \f$ S_{xy} = S_{yx} \f$
# \param xx \f$ S_{xz} = S_{zx} \f$
# \param xx \f$ S_{yy} \f$
# \param xx \f$ S_{yz} = S_{yz} \f$
# \returns \f$ S_{xx}\bf{\hat{x}\hat{x}} +
# S_{xy}(\bf{\hat{x}\hat{y}} + \bf{\hat{y}\hat{x}}) +
# S_{xz}(\bf{\hat{x}\hat{z}} + \bf{\hat{z}\hat{x}}) +
# S_{yy} \bf{\hat{y}\hat{y}} +
# S_{yz}(\bf{\hat{y}\hat{z}} + \bf{\hat{z}\hat{y}}) -
# (S_{xx} + S_{yy}) \bf{\hat{z}\hat{z}} \f$
def from_natural_form(xx, xy, xz, yy, yz):
    """T = from_natural_form(xx, xy, xz, yy, yz) returns a cartesian Tensor
    that is symmetric and trace free."""
    return poinsot(xx, yy, -(xx+yy), yz, xz, xy)


## __Tr__ calls the method (for Tensor or Versor).
# \param tensor_ A Tensor: \f$ T_{ij} \f$
# \retval result  Trace: \f$ T^i_i \f$
def Tr(self):
    """Scalar(reduce(add, self.diag()))"""
    return self.ZOO[0](self.xx + self.yy + self.zz)


## __det__(T) = __Tr__((T ^ T.T ^ T).contract(0, 2))/6
# \param tensor_ A Tensor: \f$ T_{ij} \f$
# \returns \f$ T_{ii'}T_{jj'}T_{kk}\epsilon_{ijk}\epsilon_{i'j'k'} \f$
# @image html det.png
# @image latex det.png
def det(tensor_):
    """scalar_ = det(tensor_) =
    t.xx*(t.yy*t.zz - t.zy*t.yz) -
    t.xy*(t.yx*t.zz - t.zx*t.yz) +
    t.xz*(t.yx*t.zy - t.zx*t.yy)"""
    from .vector import scalar_triple_product
    return scalar_triple_product(*tensor_.rows())


## <a href="http://en.wikipedia.org/wiki/Permanent">Permanent</a>
# \param tensor_ A Tensor: \f$ T_{ij} \f$
# \returns \f$ \sum_{\sigma \in S_3}{\prod_{i=0}^{n-1}{T_{i, {\sigma(i)}}}}\f$
def perm(t):
    """Permanent of a Tensor (is not an invariant, other than
    interchange of axes). It has no simple geometric interpretation,
    and thus does not deserved method status."""
    from ..schur_weyl.monte import Sym
    #pylint: disable=W0212
    return reduce(operator.add,
                  (reduce(operator.mul,
                          (t.e(i, sigma[i]) for i in t._list_index()))
                   for sigma in Sym(euclid.DIMENSIONS)))


## Antisymmetric Part
# \param T
# \returns \f$ \frac{1}{2}(T_{ij} - T_{ji}) \f$
def skew(T):
    """Antisymmetric part of a rank-2 Tensor."""
    return (T - T.T)/2


## Symmetric Part
# \param T
# \returns  \f$ \frac{1}{2}(T_{ij} + T_{ji}) \f$
def symm(T):
    """Symmetric part of a rank-2 Tensor."""
    return (T + T.T)/2


## <a href="http://en.wikipedia.org/wiki/Adjugate_matrix">
# Adjugate Matrix</a>.
# \f$ (t_l b_r - b_l t_r ) (-1)^p \f$ Is a private combinator.
# \returns: \f$ {\rm adj}A =  \left(\begin{array}{ccc} +
# \left|\begin{array}{cc}A_{yy}&A_{yz}\\A_{zy}&A_{zz}\end{array}\right|&
# -\left|\begin{array}{cc}A_{xy}&A_{xz}\\A_{zy}&A_{zz}\end{array}\right|&
# +\left|\begin{array}{cc}A_{xy}&A_{xz}\\A_{yy}&A_{yz}\end{array}\right|
# \\
# -\left|\begin{array}{cc}A_{xy}&A_{xz}\\A_{yy}&A_{yz}\end{array}\right|&
# +\left|\begin{array}{cc}A_{yx}&A_{yz}\\A_{zx}&A_{zz}\end{array}\right|&
# -\left|\begin{array}{cc}A_{xx}&A_{xz}\\A_{yx}&A_{yz}\end{array}\right|
# \\
# +\left|\begin{array}{cc}A_{yx}&A_{yy}\\A_{zx}&A_{zy}\end{array}\right|&
# -\left|\begin{array}{cc}A_{xx}&A_{xy}\\A_{zx}&A_{zy}\end{array}\right|&
# +\left| \begin{array}{cc}A_{xx}&A_{xy}\\A_{yx}&A_{yy}\end{array}\right|
# \end{array}\right)   \f$
def adj1(tensor_):
    """tensor_ = adj(tensor_) is the transpose of the cofactors:

    e.inner(((e*t).transpose(0, 2, 1)*t).transpose(2,1,0), 2)/2
    """
    return zipcols(*[vector_*((-1)**(1+parity))
                     for parity, vector_ in
                     enumerate(
                         [a ^ b for a, b in itertools.combinations(
                             reversed(list(tensor_.rows())), 2)])])


## \f$ \epsilon_{ijk}T_{kl}T_{jm} \f$
def adj2(t):
    """21 times slower than adj()"""
    from .three import EPSILON as e
    return e.inner(((e*t).transpose(0, 2, 1)*t).transpose(2, 1, 0), 2)/2


adj = adj2


## Cofactor tensor
def cof(tensor_):
    """cof(A) --> adj(A).T"""
    return adj(tensor_).T


## Inverse.
# \param tensor_ A Tensor: \f$ T_{ij} \f$
# \returns Inverse \f$ (T^{-1})_{ij}T_{jk} = {\bf \delta}_{ik} \f$
# \throws geo.utils.exceptions.DegenerateInversionError
def inv(tensor_):
    """T' = inv(T) = adj(t) / det(t)

    T'T = I

    (Note: it's not == DELTA, b/c
    >>>inv(T)*T
    will retain the array_like structure of T."""
    try:
        return adj(tensor_) / det(tensor_)
    except ZeroDivisionError:
        from ...utils.exceptions import DegenerateInversionError
        raise DegenerateInversionError(
            "Can't invert: \n{}".format(str(tensor_)))


## Explicit Dual: \f$ {\bf \tilde{T}} \equiv \frac{1}{2} {\bf{ T:\epsilon}} \f$
# \param tensor_ \f$ T_{ij} \f$
# \returns vector.Vector \f$ v_i = \frac{1}{2} \epsilon_{ijk} T_{jk} \f$
def dual(tensor_):
    u"""vector = dual(tensor_).

    v_i = \u00bd\u2211_jk[T_jk] = (t.yz - t.zy, t.zx - t.xz, t.xy - t.yz) / 2

    Inverts vector.dual(v)."""
    return euclid.ZOO[1](tensor_.yz - tensor_.zy,
                         tensor_.zx - tensor_.xz,
                         tensor_.xy - tensor_.yx) / 2


## Deviatoric 
# \param T
# \returns \f$ {\bf T} - \frac{1}{3}(tr{\bf T}){\bf I} \f$
def dev(T):
    """T.dev() --> the deviatoric."""
    return T - hyd(T)


## Hydrostatic (or spherical, or isotropic) Part.
# \param T
# \return \f$ Tr{T}\delta_{ij}/3\f$
def hyd(tensor_):
    """Hydrostatic part."""
    return tensor_.fromdiag(*itertools.repeat(
        Tr(tensor_).w/float(euclid.DIMENSIONS), euclid.DIMENSIONS))


## A synonym in the continuum community
isotropic = hyd


## \f$ {\bf \vec{v'}}={\bf M\vec{v}} \f$
# \param m A Tensor: \f$ M_{ij} \f$
# \param u vector.Vector \f$ u_i \f$
# \returns vector.Vector (\f$ v_i = m_{ij}u_j \f$)
def posterior_product(m, u):
    """v_i = m_ij u_j"""
    return type(u)(m.xx*u.x + m.xy*u.y + m.xz*u.z,
                   m.yx*u.x + m.yy*u.y + m.yz*u.z,
                   m.zx*u.x + m.zy*u.y + m.zz*u.z)


## \f$ {\bf M'}={\bf MT} \f$
# \param m A Tensor: \f$ M_{ij} \f$
# \param t A Tensor: \f$ T_{ij} \f$
# \returns Tensor \f$ m_{ij} = m'_{ik}m''_{kj} \f$
def _air_quote_matrix_mul(m, t):
    """m'_ik = m_ij t_jk (explicitly)."""
    return Tensor(
        m.xx*t.xx + m.xy*t.yx + m.xz*t.zx,
        m.xx*t.xy + m.xy*t.yy + m.xz*t.zy,
        m.xx*t.xz + m.xy*t.yz + m.xz*t.zz,

        m.yx*t.xx + m.yy*t.yx + m.yz*t.zx,
        m.yx*t.xy + m.yy*t.yy + m.yz*t.zy,
        m.yx*t.xz + m.yy*t.yz + m.yz*t.zz,

        m.zx*t.xx + m.zy*t.yx + m.zz*t.zx,
        m.zx*t.xy + m.zy*t.yy + m.zz*t.zy,
        m.zx*t.xz + m.zy*t.yz + m.zz*t.zz)


## Dilation from the right: these matter when using internal DoF.
_tensor_dilation = euclid.Tensor_.right_dilation


## Multiply a Tensor and a Scalar to get a Tensor .
# \param m A Tensor: \f$ M_{ij} \f$
# \param s scalar.Scalar
# \returns Tensor
def _tensor_times_scalar(m, s):
    """m' = _tensor_times_scalar(m, s)

    s       is a Scalar
    m, m'   is a Tensor
    """
    return _tensor_dilation(m, s.w)


## This is here because this package needs multimethods, and the
# simple workaround got out of hand, fast.
def _spacecurve(t, sc):
    return (t * sc.vector).spacecurve(sc.t)


## Mixin for scalar invariants
class _Invariant2Mixin(object):
    """Invariants of the Rank-2 Tensor."""

    ## \f$ J_1(T) \equiv \frac{1}{\vec{x}\cdot(\vec{y} \times \vec{z}}
    # \vec{x}_i T_{ij} (\vec{y} \times \vec{z})_j
    # + \mathrm{ cycl}\ldots \rightarrow T_{ii} \f$
    # \return Scalar trace (with non-Cartesian calculation).
    def J1(self):
        """1st (linear) invariant (is the trace):

        ([T.xi, Y, Z] + [X, T.yi, Z] + [X, Y, T.zi]) / [X, Y, Z]
        """
        from .vector import scalar_triple_product as A # associator
        from .vector import X, Y, Z
        return (A(self.xi, Y, Z) +
                A(X, self.yi, Z) +
                A(X, Y, self.zi)) / A(X, Y, Z)

    ## The sum of the principal minors
    # \returns II_A = \f$ \frac{1}{2}\big[
    # Tr^2(A) - Tr(A^2)\big] \f$
    def J2(self):
        """2nd (quadratic) Invariant:

        [Tr(A)**2 - Tr(A**2)]/2
        """
        #pylint: disable=W0212
        return (self.J1()**2 - (self**2).J1())/2

    ## Alternate DET
    # \returns III_A = \f$ \frac{1}{6}\big[
    # Tr^3(A) - 3Tr(A^2) Tr^2(A) + 2 Tr^3(A)\big] \f$
    def J3(self):
        """3rd (cubic) invariant:

        [Tr(A)**3 - 3*Tr(a)*Tr(A**2) + 2*Tr(A**3)]/6
        """
        t = self.J1()
        #pylint: disable=W0212
        t2 = (self**2).J1()
        return (t**3 - 3*t*t2 + 2*(self**3).J1())/6

    ## Invariants:
    # \param n \f$ n \in (1, 2, 3) \f$
    # \return \f$ J_n(T) \f$
    # \throws ValueError
    def J(self, n):
        """J(n) for n=1, 2, 3 is the n-th invariant (the coefficients
        in Tensor.poly().

        See: J1, J2, J3."""
        if int(n) not in (1, 2, 3):  # raise
            raise ValueError("{} is not 1, 2, or 3.".format(n))
        return [self.J1, self.J2, self.J3][int(n)-1]()


## Developmental mix in
class _DevMixIn(object):
    """Methods not totally tested or necessary"""

    def _generalized_power(self, n):
        #pylint: disable=E1121
        u = self.diagonalizer()
        t = u(self)
        return u.I(self.fromdiag(*[complex(item) ** n for item in t.diag()]))

    ## Rotation angle (for a rotation)
    # \returns array_like \f$ \theta = \cos^{-1}{\frac{Tr(M)-1}{2}} \f$
    def angle_of_rotation(self):
        """for orthogonal tensors..."""
        return arccos((self.trace().w-1)/2)

    ## Rotation axis (for a rotation)
    # \return Vector \f$ -\hat{\tilde{\bf T}} \f$
    def axis_of_rotation(self):
        """axis of rotation"""
        return -self.dual().hat()

    ## Combine axis_of_rotation() and angle_of_rotation() via Vector.versor().
    # \returns Versor
    def possible_versor(self):
        """Is this the versor?"""
        return self.axis_of_rotation().versor(self.angle_of_rotation())

    ## Tensor Logarithm
    # \return \f$ T' = \log{T} \f$
    @matrix_wrap
    def log(self):
        """Only for diagonalizable"""
        #pylint: disable=R0201
        from scipy.linalg import logm
        return logm

    ## Numerical check on diagonalness
    # \kwd tol=1.e-14
    # \returns bool If it's diagonal.
    def is_diagonal(self, tol=1.e-14):
        """True IFF the tensor is close to diagonal."""
        return self._how_diagonal().w < tol

    ## Some metric where "0" is diagonal.
    def _how_diagonal(self):
        D = self.fromdiag(*self.diag())  # diagonal elements
        return (self-D).dot(C=True) / D.dot(C=True).w

    ## Find Diagonalizing Unitary Transformation:
    # \param \f$ T \f$ (implicit)
    # \returns \f$U \rightarrow U(T) = U^{-1}TU \f$ is diagonal.
    @discharger
    @enlister
    def diagonalizer(self):
        """Diagonalizing Tensor, U:

        >>>U.I * T * U  # is diagonal."""
        return zipcols(*(self.eig()[-1]))

    ## Diagonalize
    def diagonalized(self):
        """U.I*T*U for U = T.diagonalizer())"""
        U = self.diagonalizer()
        return U.I * self * U

    ## Tensor square root:
    # \returns \f$ {\bf{T}'} = {\bf{T}}^{\frac{1}{2}} \f$
    @matrix_wrap
    def sqrt(self):
        """Sqrt, via scipy.linalg"""
        #pylint: disable=R0201
        from scipy.linalg import sqrtm
        return sqrtm


##  <a href="http://en.wikipedia.org/wiki/Tensor">Rank-2 Tensor</a>
# class transforms as \f$ T_{ij}' = M_{ik}M_{jl}T_{kl} \f$
# @image html tensor.png
# @image latex tensor.png
class Tensor(euclid.Tensor_, euclid.LinearMap, _Invariant2Mixin, _DevMixIn):
    """T = Tensor(xx, xy, xz, yx, yy, yz, zx, zy, zz)

    Is a Cartesian tensor, and it is a function from E3-->E3
    (a rotation matrix).


    TRANSFORMING VECTORS:

    >>>v_prime = T(v)

    That is, the __call__ method does it for you. Use it-- do not multiply.
    You can multiply or do an explicit transformation with any of the
    following:

    T.alias_transformation(v) --> T*v
    T.alibi_transformation(v) --> T.T*v = v*T

    You can convert to other charts on SO(3) with:

    T.versor()
    T.ypr()

    Or get individual angles:

    T.yaw, T.pitch, T.roll

    Other matrix/tensor methods are:

    T.T       (same as T.transpose())
    ~T        (same as T.I inverts it)
    abs(T)    ...det?
    T.trace()    trace (a Scalar)
    T.L2norm()   L2-norm (a Scalar)
    T.A()          antisymmetric part
    T.S()          symmetric part
    T.tensor(2)    symmetric trace free part
    T.dual()   Antisymmetric part (contracted with Levi-Civita tensor)

    Finally:

    T.row(n), T.col(n), T.rows(), T.cols(), T.tolist() are all

    pretty simple... except: rotation matrices are not tensors-but I
    use them as such. The tensor is just a map from R3->R3.

    Basis Dyads [e_i][e_j] are available:
    Tensor.baseunit(i, j)
    """
    from ...utils.exceptions import OverloadNotImplementedError as DivisionError

    __metaclass__ = euclid.ranked

    ## The rank is 2.
    rank = 2

    ## The hash table assigns multiplication
    _dispatch_mul = {None: _tensor_dilation,
                     0: _tensor_times_scalar,
                     1: posterior_product,
                     2: _air_quote_matrix_mul,
                     'spacecurve': _spacecurve}

    ## vector helper.
    # \param v0 \f$ T_0j \f$
    # \param v1 \f$ T_1j \f$
    # \param v2 \f$ T_2j \f$
    # \returns \f$ T_{ij} = (\vec{v}_i)_j \f$
    @classmethod
    def ziprows(cls, v0, v1, v2):
        """T = ziprows(v0, v1, v2) --> T.row(n) = v<n>

        (zips 3 Vectors into a Tensor)."""
        return cls(*itertools.chain(*[item.iter() for item in (v0, v1, v2)]))

    ## vector helper.
    # \param v0 \f$ T_i0 \f$
    # \param v1 \f$ T_i1 \f$
    # \param v2 \f$ T_i2 \f$
    # \returns \f$ T_{ij} = (\vec{v}_i)_j \f$
    @classmethod
    def zipcols(cls, v0, v1, v2):
        """T = zipcols(v0, v1, v2) --> T.col(n) = v<n>

        (zips 3 Vectors into a Tensor).

        Hence, If you know (x, y, z) --> (x', y', z'), then the
        transformation that does that is:

        T = zipcols(x', y', z')
        """
        return cls.ziprows(v0, v1, v2).T

    ## Construct from a 3 x 3 numpy matrix
    # \param matrix_ A numpy 3 x 3 matrix
    # \returns \f$ T_{ij}=M[i][j]\f$
    @classmethod
    def frommatrix(cls, matrix_):
        """Contruct from a 3x3 or 9 element numpy matrix."""
        from numpy import array
        return cls.fromarray(array(matrix_).ravel())

    ## Construct from diagonal elements
    # \param xx \f$ T_{xx} \f$
    # \param yy \f$ T_{yy} \f$
    # \param zz \f$ T_{zz} \f$
    # \returns \f$ {\bf T} =
    # \left[ \begin{array}{ccc} T_{xx} & 0 & 0 \\
    # 0 & T_{yy} & 0 \\
    # 0 & 0 & T_{zz} \end{array} \right] \f$
    @classmethod
    def fromdiag(cls, xx, yy, zz):
        """fromdiag(xx, yy, zz)
        Contructs Tensor from diagonal elements:
        >>>print Tensor.fromdiag('xx', 'yy', 'zz')
        [xx, 0.0, 0.0]
        [0.0, yy, 0.0]
        [0.0, 0.0, zz]"""
        return cls(xx, 0., 0.,
                   0., yy, 0.,
                   0., 0., zz)

    ## explicit 9 argument init.
    # \param xx \f$T_{xx}\f$ component
    # \param xy \f$T_{xy}\f$ component
    # \param xz \f$T_{xz}\f$ component
    # \param yx \f$T_{yx}\f$ component
    # \param yy \f$T_{yy}\f$ component
    # \param yz \f$T_{yz}\f$ component
    # \param zx \f$T_{zx}\f$ component
    # \param zy \f$T_{zy}\f$ component
    # \param zz \f$T_{zz}\f$ component
    # \returns \f$ {\bf T} =
    # \left[ \begin{array}{ccc} T_{xx} & T_{xy} & T_{xz} \\
    # T_{yx} & T_{yy} & T_yz{} \\
    # T_{zx} & T_{zy} & T_{zz} \end{array} \right] \f$
    def __init__(self, xx, xy, xz, yx, yy, yz, zx, zy, zz):
        ## <a href="http://en.wikipedia.org/wiki/Tensor">Cartesian
        # Component:</a> \f$ m_{xx} = {\bf{T}}^{({\bf e_x})}_x \f$
        self.xx = xx
        ## <a href="http://en.wikipedia.org/wiki/Tensor">
        # Cartesian Component:</a> \f$ m_{xy} = {\bf{T}}^{({\bf e_y})}_x \f$
        self.xy = xy
        ## <a href="http://en.wikipedia.org/wiki/Tensor">
        # Cartesian Component:</a> \f$ m_{xz} = {\bf{T}}^{({\bf e_z})}_x \f$
        self.xz = xz
        ## <a href="http://en.wikipedia.org/wiki/Tensor">
        # Cartesian Component:</a> \f$ m_{yx} = {\bf{T}}^{({\bf e_x})}_y \f$
        self.yx = yx
        ## <a href="http://en.wikipedia.org/wiki/Tensor">
        # Cartesian Component:</a> \f$ m_{yy} = {\bf{T}}^{({\bf e_y})}_y \f$
        self.yy = yy
        ## <a href="http://en.wikipedia.org/wiki/Tensor">
        # Cartesian Component:</a> \f$ m_{yz} = {\bf{T}}^{({\bf e_z})}_y \f$
        self.yz = yz
        ## <a href="http://en.wikipedia.org/wiki/Tensor">
        # Cartesian Component:</a> \f$ m_{zx} = {\bf{T}}^{({\bf e_x})}_z \f$
        self.zx = zx
        ## <a href="http://en.wikipedia.org/wiki/Tensor">
        # Cartesian Component:</a> \f$ m_{zy} = {\bf{T}}^{({\bf e_y})}_z \f$
        self.zy = zy
        ## <a href="http://en.wikipedia.org/wiki/Tensor">
        # Cartesian Component:</a> \f$ m_{zz} = {\bf{T}}^{({\bf e_z})}_z \f$
        self.zz = zz

    ## Tensor as a linear function from R3-->R3 or (R3 x R3) --> R.
    #  \param v1 \f$ \rightarrow {\bf T(\vec{v}_1)} \f$ (call super).
    #  \param v2 \f$ \rightarrow {\bf T}({\bf \vec{v}_1, \vec{v}_2}) \f$
    # \return \f$ T_{ij}v^{(1)}_j, T_{ij}v^{(1)}_j v^{(2)}_i \f$
    def __call__(self, v1, v2=None):
        """The Tensor as a function of 1 or 2 vectors:
        T(a)    --> T*a
        T(a, b) --> T(a)*b"""
        if v1 is None:  # kwd
            return self if v2 is None else self(v2)  # kwd
        if v2 is None:  # kwd
            result = super(Tensor, self).__call__(v1)
        else:
            print "Check ORder of T(u, v)"
            result = super(Tensor, self).__call__(v2)(v1)
        return result
            

    ## Trace
    trace = Tr

    ## Deviatoric
    dev = dev

    ## Hydrostatics
    hyd = hyd

    ## Isotropicerical part.
    isotropic = isotropic
    
    ## Determinant
    det = det

    ## Adjugate: \f$ {\mathrm{adj}}({\bf A})_{jk} = \frac{1}{2}[(A_{ii}^2 -
    # (A^2)_ii)\delta_{jk} - A_{jk}A_{ii} + (A^2)_{jk}] \f$
    adj = adj

    ## Cofactor
    cof = cof

    ## Dual
    dual = dual

    ## Skew Symmetric:  30 Times faster than S
    skew = skew

    ## Symmetric part: 30 Times faster than A
    symm = symm

    ## Inverse
    I = property(inv)

    ## Generalized Scalar Trace (see J1).
    # \returns \f$ Tr(T) \equiv \frac{1}{\vec{x}\cdot(\vec{y} \times \vec{z})}
    # \vec{x}_i T_{ij} (\vec{y} \times \vec{z})_j
    # + \mathrm{ cycl}\ldots \rightarrow T_{ii} \f$
    def generalized_trace(self):
        """Scalar Trace (non cartesian version)


        ([Tx, y, z] + [x, Ty, z] + [x, y, Tz]) / [x, y, z]
        """
        cycle = (itertools.imap(euclid.DIMENSIONS.__rmod__,
                                itertools.imap(n.__add__,
                                               self._list_index()))
                 for n in self._list_index())
        return sum(e_i*self*(e_j ^ e_k) for e_i, e_j, e_k in
                   (itertools.imap(DELTA.e, ind) for ind in cycle))

    ## Spherical (Fundamental, complex) Representation:
    # \returns \f$T\f$ as a geo.metric.wigner.eckart.Tensor
    # @image html wig1.png
    # @image latex wig1.png
    def spherical(self):
        """Spherical representation."""
        from ..wigner import eckart
        return eckart.Tensor.fromcart(self)

    ## Coordinate Change
    # \param scalar_ a scalar.Scalar or number
    # \returns \f$ S' = S ||T|| \f$
    # \bug The determinant scales the operation! Be Careful, also it abs
    # because there is no pseudo scalar.
    # @image html Y00_comp.gif
    # @image latex Y00_comp.gif
    def _alibi_transform_scalar(self, scalar_):
        """"scalar transform is trivial -unless the determinant is
        not one-- then you got trouble."""
        return abs(self.det()) * scalar_

    ## Coordinate Change
    # \param vector_ a vector.Vector
    # \returns \f$ V_j' = T_{ij}V_i \f$ \n (note summation index)
    def _alibi_transform_vector(self, vector_):
        """vector transform is left mul"""
        return self * vector_

    ## Coordinate Change
    # \param tensor_ a tensor.Tensor
    # \returns \f$ T_{ij}' = T_{kl}M_{ik}M_{jl} \f$
    # @image html trot.png
    # @image latex trot.png
    def _alibi_transform_tensor(self, tensor_):
        """tensor transform is a similarity transformation"""
        #pylint: disable=E1101
        return self.sandwich(tensor_)

    ## Coordinate Change (ad hoc).
    # \param q A Versor or Quaternion
    # \returns \f$ (q'; q_i') = (q; T_{ij}q_j) \f$
    def _alibi_transform_quaternion(self, q):
        """farm out quaternion parts to scalar & vector."""
        return type(q)(*map(self, q.iter()))

    ## Extract a "row-vector" (run on last index).
    # \param n dimension in left tensor index
    # \returns vector.Vector from \f$ v_i = T_{ni} \f$
    def row(self, n):
        """row is 'e'"""
        # why is this a call super?
        return super(Tensor, self).e(n)

    ## Extract a "column-vector" (run on first index).
    # \param n dimension in right tensor index
    # \returns vector.Vector from \f$ v_i = T_{in} \f$
    def col(self, n):
        """Run on 1st index. See row.__doc__ """
        return self.e(Ellipsis, n)

    ## Iterate row() calls.
    # \yields rows
    def rows(self):
        """iterator over row(n) for n in 0, 1, 2 """
        return itertools.imap(self.row, self._list_index())

    ## Iterate col() calls.
    # \yields columns
    def cols(self):
        """iterator over col(n) for n in 0, 1, 2 """
        return itertools.imap(self.col, self._list_index())

    ## (Explicit) transpose indices-- just for rank 2 (for 8x speed).
    # \returns  \f$ [{\bf T}^T]_{ij} = T_{ji} \f$
    @property
    def T(self):
        """Explicit transpose(1, 0)."""
        return type(self)(self.xx, self.yx, self.zx,
                          self.xy, self.yy, self.zy,
                          self.xz, self.yz, self.zz)

    ## The is here to maintain the transpose with arguments protocol, even
    # though there is only 1 sensible value: args=(1, 0).
    # \param args (1, 0) makes sense
    # \returns \f$T_{ji}\f$ for args=(j=1, i=0).
    def transpose(self, *args):
        """tensor.transpose(*args):

        args must be in -1, 0, 1.
        """
        return super(Tensor, self).transpose(*args)

    ## left permutation of indices (trivial for rank 2).
    _left = transpose

    ## right permutation of indices (trivial for rank 2).
    _right = transpose

    ## reversed permutation of indices
    _reversed_indices = transpose

    ## Mutliplication Inverse:
    # \returns \f$ {\bf T}^{-1} \ | \  {\bf T}^{-1} {\bf T} = {\bf I} \f$
    def __invert__(self):
        """~T --> T.I"""
        return self.I

    ### Axial Vector part of Tensor
    ## \returns \f$ v_k \rightarrow  m_ij \epsilon_{ijk} \f$ \n
    #def vector(self):
    #    """Convert to a vector w/o scaling"""
    #    raise DeprecationWarning("Use 2*dual()")

    ## not quite right note the clean call....
    def __str__(self):
        return "\n".join(map(str, self.rows()))

    ## Integer and 1/2-integer exponents.
    # \returns \f$ T^{\alpha} \f$ for \f$ |\alpha| \in
    # (0, \frac{1}{2}, 1, \frac{3}{2}, 2, \frac{5}{2}, \ldots) \f$
    def __pow__(self, n):
        if n < 0:  # power branch
            result = self.I.__pow__(-n)
        elif n == 0:
            result = self * self.I
        elif n == 1:
            result = self
        elif n == 0.5:
            result = self.sqrt()
        else:
            i, f = divmod(n, 1)
            result = reduce(operator.mul,
                            itertools.repeat(self, i)) if i else self ** 0
            if f:  # todo verify, array v singleton
                from numpy import log, exp
                result *= exp(log(self) * f)

        return result

    ## DivisionError()
    # \throws DivisionError()
    def __rdiv__(self, other):
        msg = "/tensor => * tensor ** (-1) overload not implemented"
        raise self.DivisionError(msg)

    ## geo.metric.euler.charts.YPR pitch
    @property
    def pitch(self):
        """YPR pitch (deg)"""
        return degrees(arcsin(self.xz))

    ## geo.metric.euler.charts.YPR roll
    @property
    def roll(self):
        """YPR roll (deg)"""
        return degrees(arctan2(-self.yz, self.zz))

    ## geo.metric.euler.charts.YPR yaw
    @property
    def yaw(self):
        """YPR yaw (deg)"""
        return degrees(arctan2(-self.xy, self.xx))

    ## Convert to a tuple of (yaw(), pitch(), roll())
    # \returns \f$(\alpha, \beta, \gamma)\f$ for a
    # geo.metric.euler.charts.YPR instance.
    def ypr_tuple(self):
        """compute to angle triplet"""
        roll = self.roll
        pitch = self.pitch
        return (self.yaw, pitch, roll)

    ## geo.metric.euler.charts.YPR
    # \returns YPR see ypr_tuple()
    def ypr(self):
        """convert to YPR class"""
        from ..euler.charts import YPR
        return YPR(*self.ypr_tuple())

    ## geo.metric.euler.charts.RPY
    def rpy(self):
        """Convert to RPY class"""
        return self.ypr().rpy()

    ## geo.metric.euler.charts.Euler
    # \returns \f$ \left( \begin{array}{c}
    # \pi - \tan^{-1}{T_{xz} / T_{yz}  } \\
    # \cos^{-1}{T_{zz}} \\
    # -\tan^{-1}{-T_{zz} / T_{zy} } \end{array} \right) \f$
    def euler(self):
        """Convert to Euler class"""
        from numpy import pi
        from ..euler.charts import Euler
        try:
            choice = bool(self.zz < 0.99999999)
        except ValueError:
            # handle array_like (could be a decorator)
            for count, item in enumerate(self):
                eul = item.euler()
                if count:  # init
                    result = result.append(eul)
                else:
                    result = eul
            return result

        if choice:  # degeneracy branch
            theta = arccos(self.zz)  # these are passive euler angles?
            phi = -arctan2(-self.zx, self.zy)
            psi = pi-arctan2(self.xz, self.yz)
        else:
            phi = 0.
            theta = 0.
            psi = arctan2(self.yx, self.xx)
        return Euler(*[(degrees(item) % 360) for item in (psi, theta, phi)])

    ## Convert to a rotation versor via YPR()
    def versor(self):
        """Convert to a rotation versor--w/o checking for
        viability."""
        return self.ypr().versor()

    ## [A, B] = AB - BA
    # \returns \f$ {\bf [A, B] = AB - BA } \f$
    def commutator(self, other):
        """commutator(other) --> [self, other]"""
        return (self * other) - (other * self)

    ## {A, B} = AB + BA
    # \returns \f$ {\bf \{A, B\} = AB + BA } \f$
    def anticommutator(self, other):
        """commutator(other) --> [self, other]"""
        return (self * other) + (other * self)

    ## \param self (implicit) \f$ A \f$
    # \param other A Tensor (probably), \f$ B \f$
    # \returns \f$ A B A^T \f$
    # todo: complexify (or not?)
    def sandwich(self, other):
        """Sandwich product"""
        return self * other * self.T

    ## \param self (implicit) \f$ A \f$
    # \param other A Tensor (probably), \f$ B \f$
    # \returns \f$ A^{-1} B A \f$
    # todo: complexify (or not?)
    def similarity(self, other):
        """T' = U.similarity(T) = U.I * T * U"""
        return self.I * other * self

    ## Characteristic Invariants
    def characteristics(self):
        """returns a list of characteristics:
        [-Tr(A), II_A, ||A||]"""
        return [-self.trace(), self.J(2), self.det()]

    ## Characteristic Polynomial by Direct Computation
    # \returns \f$ p(x) = ||{\bf T}-x{\bf I}|| \f$
    def poly(self):
        """Characteristic polynomial (as a numpy.poly1d)."""
        from numpy import poly1d
        x = poly1d([1, 0])
        # need to convert self.ij into poly1d[self.ij]
        return (self.broadcast(poly1d) - DELTA * x).det().w

    ## Convert to a numpy matrix.
    # \retval matrix or a list of them.
    @enlister
    def asmatrix(self):
        """As a numpy matrix, or a list of them"""
        from numpy import matrix
        return matrix([item.tolist() for item in self.rows()])

    ## geo wrapped <a href=
    # "http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.eig.html">linalg.eig.</a>.
    # \retval tuple (Scalar eigenvalues, Vector eigenvectors), or a list
    # of them.
    @enlister
    def eig(self):
        """returns scalar eigenvalue and vector eigenvector,
        or list of them:

        >>>(s1, s2, s3), (v1, v2, v3) = t.eig()
        """
        from numpy import linalg
        from numpy import array
        lams, vecs = linalg.eig(self.asmatrix())
        return (map(self.ZOO[0], lams),
                list(itertools.starmap(self.ZOO[1], array(vecs.T))))

    ## Why? Not an invariant.
    def diag(self):
        """Tuple of diagonals? See svd()"""
        return self.xx, self.yy, self.zz

    ## Dyadic Decomposition
    @enlister
    def svd(self):
        """(u1, v1), .., (uN, vN) = T.svd() implies:

        T = (u1 & v1) + ... + (uN & vN)

        and N is as small as possible...(1, 2, or 3).
        """
        from numpy.linalg import svd
        u, sigma, v = svd(self.asmatrix())
        return [(s_i*u_i, v_i) for u_i, s_i, v_i in
                zip(self.frommatrix(u).cols(),
                    sigma,
                    self.frommatrix(v).rows())
                if s_i]  # filter

    ## Dot product (contracts tensor from right with other from left).
    # \param other
    # \returns  \f$ T_{ij}U_{ij} \f$
    # \throws TypeError for non-tensor arguments.
    def dot(self, other):
        """T.dot(U) in the inner product with the max possible contracitons."""
        try:
            return self.inner(other, min(self.rank, other.rank))
        except (AttributeError, TypeError):
            raise TypeError("inner expected a tensor, not: {}".format(
                type(other)))
        
    ## Double dot product
    # \returns \f$ (A:B) = A_{ij}B_{ji} \f$
    # \throws ValueError for non rank-2 args.
    def double_dot(self, other):
        """A.double_dot(B) --> A:B"""
        if euclid.rank(other) == 2:  # raise
            return self.inner(other, 2)
        raise ValueError(
            "Double Dot is for rank 2, not: {}".format(euclid.rank(other)))
            
    ## <a href="http://en.wikipedia.org/wiki/Dyadics#Product_of_dyadic_and_dyadic">Dot-Cross Product</a>
    # \param self \f$ [AB]_{ij} \f$
    # \param other \f$ [CD]_{lm} \f$
    # \returns \f$ v_k = [AB]_{ij}\epsilon_{jkl}[CD]_{li} \f$
    # @image html dot_cross.png
    # @image latex dot_cross.png
    def dot_cross(self, other):
        u"""dot cross product of 2 tensors:
        (A\u00b7\u2a2fB)_k = \u2211_ijl[A_ij*\u025b_lkl*B_li]"""
        return (self ^ other.T).contract(0, 2)

    ## <a href="http://en.wikipedia.org/wiki/Dyadics#Product_of_dyadic_and_dyadic">Cross-Dot Product</a>
    # \param self \f$ [AB]_{ij} \f$
    # \param other \f$ [CD]_{lm} \f$
    # \returns \f$ v_k = [AB]_{ji}\epsilon_{ikl}[CD]_{lj} \f$
    # @image html cross_dot.png
    # @image latex cross_dot.png
    def cross_dot(self, other):
        u"""cross dot product of 2 tensors:
        (A\u2a2f\u00b7B)_k = \u2211_ijl[A_ij*\u025b_ikl*B_lj]"""
        return (self.T ^ other).contract(0, 2)

    ## <a href="http://en.wikipedia.org/wiki/Dyadics#Product_of_dyadic_and_dyadic">Double Cross Product</a>
    # \param self \f$ [AB]_{ij} \f$
    # \param other \f$ [CD]_{lm} \f$
    # \returns \f$ T_{kn} = [AB]_{ji}\epsilon_{ikl}[CD]_{lm}\epsilon_{jnm} \f$
    # @image html double_cross.png
    # @image latex double_cross.png
    def double_cross(self, other):
        u"""Double cross product of 2 tensors:
        (A \u2a37 B)_kn = \u2211_ijlm[A_ij*\u025b_ijk*B_lm*\u025b_jnm]"""
        from .three import EPSILON
        return (self.T ^ other).transpose(1, 0, 2).inner(EPSILON, 2)

    ## Nearest Orthogonal Matrix:
    # \returns \f$ {\bf T} \rightarrow {\bf T}\sqrt{{\bf T}^T{\bf T}}^{-1} \f$
    def nom(self):
        """Nearest Orthogonal Matrix?"""
        return self * (self.T * self).sqrt().I

    ## \f$ e^T \f$
    # \kwd n_terms = 20
    # \returns \f$ \sum_{k=0}^{n-1}{\frac{1}{k!}T^k} \f$
    def exp(self, n_terms=20):
        """Taylor series computation of exp(T)."""
        result = DELTA
        m = DELTA
        f = 1.
        for i in range(1, n_terms):
            f *= i
            m *= self
            result += m/f
        return result

#    def exp2(self, n_terms=20):
#        return reduce(operator.add, (t/f

    ## Once a good idea, is now bogged down by generality.
    _dispatch_call = {None: _alibi_transform_scalar,
                      0: _alibi_transform_scalar,
                      1: _alibi_transform_vector,
                      2: _alibi_transform_tensor,
                      1j: _alibi_transform_quaternion,
                      'spacecurve': _alibi_transform_vector}

    ## Alias transforms are from the right \n
    # \param other Anything to rotate
    # \returns other in rotated coordiated
    def alibi_transform(self, other):
        """TODO:
        [0] --> T(S) --> S
        [1] --> T(V) --> T*V
        [2] --> T(M) --> T*M*T.T
        """
        from ...metric.euclid.euclid import rank
        try:
            func = self._dispatch_call[euclid.rank(other)]
        except KeyError:
            ## higher order? use this antipattern
            from .three import transformer as func
        return func(self, other)

    @enlister
    ## For Tensor Rotations, convert to a 9 x 1 matrix
    # \returns numpy.matrix 9 x 1 column matrix.
    def nineXone(self):
        """Convert to a numpy matrix"""
        from numpy import matrix
        return matrix(list(
            itertools.starmap(self.e, itertools.product(
                range(euclid.DIMENSIONS), repeat=2)))).T

    ## Kronecker Product of self and inverse.
    # \returns \f$ {\bf T \bigotimes T^I} \f$
    def nineXnine(self):
        """Kronecker product with inverse."""
        return self.kronecker(self.I)

    ## For Tensor Rotations, convert to a 9 x 9 matrix
    # \returns \f$ {\bf T \bigotimes T} \f$
    @enlister
    def kronecker(self, other):
        """T" = T.kronecker(T') is the kronecker product
        (as a 9 x 9 number matrix)."""
        from numpy import kron
        return kron(self.asmatrix(), other.asmatrix())

    ## The Natural Form (Symmetric and Trace Free in All indices)
    # \return \f$ T_{(2)} = T_{\{ij\}} - \frac{1}T_{ii}\delta_{ij} \f$
    def natural_form(self):
        """Symmetric, trace-free"""
        ## this is lame, as it should be able to express this
        # as a totally symmetic module minus projections onto
        # isotropic tensors.
        return self.S() - self.isotropic()

    def det2(self):
        """test implementation of determinant."""
        from ..schur_weyl.monte import Sym
        result = 0.
        for sig in Sym(3):
            result += sig.sign()*reduce(
                operator.mul,
                [self.e(sig[i], i) for i in range(3)])
        return result

    ## Project Irreducible Representations
    # \param j Weight
    # \param m Azimuth projection
    # \returns \f$ T^{(j; m)} \f$
    def irrep(self, j, m):
        """irrep(j, m) --> |j, m> part"""
        return super(Tensor, self).irrep(j, 0, m)


## Project out scalar tensor.
def p00(t):
    """Scalar module projector"""
    return DELTA * Tr(t).w

## Tensor constructor alias (rows are typographic rows according to __str__)
ziprows = Tensor.ziprows

## Tensor constructor alias (cols are row's compliment).
zipcols = Tensor.zipcols

## \f$ \delta_{ij} \f$\n
# <a href="http://en.wikipedia.org/wiki/Kronecker_delta">Identity</a> Tensor
DELTA = ziprows(*_vector_basis)


from collections import namedtuple

## \f$ \left[\begin{array}{ccc}
# {\bf \hat{x}\hat{x}} &  {\bf \hat{x}\hat{y}} &  {\bf \hat{x}\hat{z}} \\
# {\bf \hat{y}\hat{x}} &  {\bf \hat{y}\hat{y}} &  {\bf \hat{y}\hat{z}} \\
# {\bf \hat{z}\hat{x}} &  {\bf \hat{z}\hat{y}} &  {\bf \hat{z}\hat{z}}
# \end{array} \right] \f$
BASIS = namedtuple("Basis", sorted(DELTA.__dict__.keys()))(
    *Tensor.basis())


DYADS = BASIS

# unpack Dyads
XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ = DYADS
