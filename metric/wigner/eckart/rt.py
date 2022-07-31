"""This Module defines the conversions from Cartesian to Spherical bases,
and back."""

## \namespace geo.metric.wigner.eckart.rt Coordinate Conversions
import itertools
import operator



## Cartesians <--> Spherical Transformation.
class SphericalTensor(object):
    """Generic spherical tensor conversions"""

    ## Get heighest weight components from a NF-tensor
    # \param S A natural form rank-L tensor 
    # \returns tuple 2L+1 components
    @classmethod
    def m_project(cls, S):
        """(0, ...L, -L, ... 1) = m_project(T.natural_form())

        argument should be a symmetric trace-free tensor, result is
        a (2L+1)-tuple (iterator) of spherical components."""
        return (cls.projectors[m](S) for m in cls.miter())


## From __1__ to __1__ L Trivial
class Scalar(SphericalTensor):
    """Scalar(w)"""

    ## Spherical Scalar from a Cartesian Scalar
    # \param s A Cartesian scalar
    # \returns eckart.Scalar A spherical Scalar
    @classmethod
    def fromcart(cls, s):
        """s = Scalar.fromcart(s')
        s' is a cartesian scalar,
        s is the spherical (trivial) SO(3) representaion."""
        # Yes: a lot of code for a trivial (literally) operator.
        return cls(*cls.m_project(s))

    ## __1__
    def tocart(self):
        """Convert to 1, the Cartesian Represenation."""
        from ...euclid.scalar import Scalar as ScalarCart
        return ScalarCart(self[0])


## Fundamental/Adjoint Rep conversions.
class Vector(SphericalTensor):
    """Vector(m=0, m=1, m=-1)"""

    ## <a href="https://en.wikipedia.org/wiki/Tensor_operator#Spherical_vector_operators">
    #  Construct </a> from cartesian Vector (NB: phase convention)
    # \param v \f$ {\bf \vec{v}} \f$ geo.metric.euclid.vector.Vector
    # \returns \f$ \begin{array}{ccc}
    #  & J=0 & J=1 \\
    # M=1 &-& -\sqrt{\frac{1}{2}}(v_x + i v_y) \\
    # M=0 &-& v_z \\
    # M=-1 &-& +\sqrt{\frac{1}{2}}(v_x - i v_y) \end{array} \f$
    @classmethod
    def fromcart(cls, v):
        """spherical = Vector.fromcart(<Cartsian Vector>):

        |+> = (-x - Iy) / 2 ** 0.5
        |0> = z
        |-> = (x - Iy) / 2 ** 0.5

        This is the Fundamental vs. the  Adjoint represenation."""
        return cls(*cls.m_project(v))

    ## Convert to cartesian: \n
    # \returns \f$ \left[ \begin{array}{c}
    # -\sqrt{\frac{1}{2}}(v^+-v^-) \\
    # i\sqrt{\frac{1}{2}}(v^++v^-) \\
    # v^0 \end{array} \right] \f$
    def tocart(self):
        """Convert to 3, the Cartesian Represenation."""
        from ...euclid.vector import Vector as VectorCart
        from .one import Vx, Vy, Vz
        return VectorCart(Vx(self), Vy(self), Vz(self))

    ## vector is tocart()-- this is ill-advised polymorphism
    vector = property(tocart)


## \f$ {\bf 3 \otimes 3 = 5 \oplus 3 \oplus 1} \f$
class Tensor(SphericalTensor):
    """Tensor(scalar, [vector x 3], [tensor x 5])"""

    # This is a direct import of a method: makes natural form Cartesian.
    from .two import T2, T1, T0


    ## Constructor from Cartesian Tensor (we can't us ::U) since
    # it would give 3 x 3, and we want 5 + 3 + 1.\n
    # \f$ T_0^{(0)} = \frac{1}{3}T_{ii} \f$ \n
    # \f$ T^{(1)}_m \leftarrow {\bf U}(\epsilon_{ijk}[T_{jk}-T_{kj}]/2) \f$ \n
    # \f$ T_{\pm 2}^{(2)}  = \frac{1}{2}(S_{xx} - S_{yy} \pm 2iS_{xy}) \f$
    # \f$ T_{\pm 1}^{(2)}  = S_{xz} \pm iS_{yz} \f$
    # \f$ T_0^{(2)} = \sqrt{\frac{3}{3}}S_{zz} \f$ \n
    # with \f$ {\bf S} = \frac{1}{2}[{\bf T} + {\bf T}^T] \f$
    @classmethod
    def fromcart(cls, t):
        """T = Tensor.fromcart(T')
        T' is a 3 x 3 Cartesian tensor,
        T  is 5 + 3 + 1 spherical tensor."""
        from .one import ROOT2, SCALAR_NORM
        #  The scalar part is computed as follows:
        # (1) Take the trace, and convert to spherical-- this lets ec.py
        # handle the list-form correctly.
        # (2) Use the SCALAR_NORM defined in one.py: this makes the transform
        # consistent with the dyadic of 2 spherical vectors the same as the
        # fromcart(u&v) for 2 cartesian vectors.
        scalar_ = t.trace().spherical()[0] / SCALAR_NORM
        
        # Tensor's dual goes straight to |1, m>, with some scale factors
        # that can be derived using racah.py. Calling dual() converts the
        # antisymmetric part into a vector, who's spherical components are
        # extracted using Vector.fromcart (by calling spherical()).
        A = t.dual().spherical() * ROOT2 / 1j

        # Note: write T = (v1)(v2), then
        # <1, m1, 1, m2| 2, m> = ClebshGordan(1, m1, 1, m2, 2, m) gives the
        # coefficient of T[2, m] in terms of v1[m1]v2[m2], which are then
        # converted to Cartesian (this should be automated for higher rank,
        # but it's not simple: 6j, 9j symbols...
        return cls(scalar_,
                   list(A),
                   list(cls.m_project(t.natural_form())))

    ## Cartesian Representation:
    # \return T0() + T1() + T2()
    def tocart(self):
        """5 + 3 + 1 --> 3 x 3"""
        # old school: apply all methods and sum.
        return reduce(operator.add, (map(apply, [self.T0, self.T1, self.T2])))


## \f$ {\bf 3 \otimes 3 \otimes 3 = 7 \oplus 3 \oplus }\ \ 2
# \cdot {\bf (5 \oplus 3) \oplus 1} \f$
class Three(SphericalTensor):
    """Three([[scalar]],
    [vector_0, vector_1, vector-2],
    [tensor_0, tensor_1],
    [m=0, m=1, m=2, m=3, m=-3, m=-2, m=-1])
    """

    ## Construct from a Cartesian tensor
    # \param t A Cartesian tensor
    # \returns t as a Spherical tensor
    @classmethod
    def fromcart(cls, t):
        """Convert Cartesian Tensor"""
        ## Get the Weyl Modules.
        from . import three
        # use the JQM class to organize all the pieces.
        return cls.fromlist(list(three.JQM(t)))

    ## Convert whole or parts to Cartesian
    # \param j None Weight
    # \param q None Seniority Index
    # \param m None Azimuth Order
    # \returns geo.metric.euclid.three.Three
    # \throws ValueError for garbage j, q, m.
    def tocart(self, j=None, q=None, m=None):
        """cartesian = spherical.tocart(j=None, q=None, m=None)
        
        converts to Cartesian, unless a keyword is set, then it
        picks out only the weight j, azimuth order m, and seniority
        index q (counting from zero)."""
        from .three2 import THREE
        if j is None:  # keyword
            return reduce(operator.add,
                          itertools.starmap(self.tocart, self.ijqm()))
        try:
            return THREE[(j, q, m)] * self[j, q, m]
        except KeyError:
            raise ValueError
        except TypeError as err:
            if THREE[(j, q, m)] is NotImplemented:  # raise
                raise NotImplementedError
            raise err


## \f$ {\bf 3 \otimes 3 \otimes 3 \otimes 3 = (9 \oplus 5 \oplus 1 ) \oplus
# } \ \ 3\cdot {\bf (7 \oplus  5 \oplus 3) \oplus} \ \ 2 \cdot
# {\bf (5 \oplus 1) \oplus} \ \ 2\cdot{\bf 3} \f$
class Four(SphericalTensor):
    """Not Coded."""


    def tocart():
        return super(Four, self).tocart()



    

