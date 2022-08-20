THE Package for Geometry on Earth.

Intro Introduction
The point of this package is to allow you to do manifestly covariant
     vector operations (and transformations), both abstractly and with
     respect to geodesy via affine spaces on Earth. Manipulation of
     Earth-bound coordinates and abstract vectors can be done without
     regard to their coordinate system. Moreover, the operations do
     not depend on any array-like properties of your objects: operations
     on single vectors or coordinates and the same as operations on multi-
     dimensional arrays of coordinates: _The_ _code_ _you_
     _write_ _is_ _the_ _same_ _for_ _either_ _case_ .

     Furthermore, all transforming objects are (polymorphic)
     
     [function emulators](https://docs.python.org/2/reference/datamodel.html#emulating-callable-objects), so that
     the code you write depends neither on the type of the object doing
     the transformation (rotation matrix, quaternion, Euler angles, etc.)
     nor on the object being transformed (Vector, Quaternion, Tensor, Scalar,
     etc.):

     \verbatim
>>>x2 = T(x1)

     \endverbatim
     
     will rotate _x1_ to _x2_, regardless of how __T__ represents
     the rotation for  _x1_ a tensor of any rank, or quaternion. Both
     __T__ and _x1_ can be singleton or array-like, and _x2_ will have the
     approriate array structure (or not) according to numpy's broadcasting
     rules.

     @subsection mc Manifest Covariance
     <a href="https://en.wikipedia.org/wiki/Manifest_covariance">
     Manifest Covariance</a> (aka
     <a href="https://en.wikipedia.org/wiki/Coordinate-free">
     "coordinate free"</a>) vector operations means you
     can write equations and manipulate vectors without regards to their
     coordinate representation. This package is all about that. Scalars,
     Vectors, and Tensors must be created with a specific coordinate system
     in mind, but they need not be manipulated with regards to their
     components: operations are coordinate free.

     When applied to Geodesy, the idea of a coordinate independent coordinate
     system may seem nonsensical, but it is not. Instances of geodesy
     classes represent a point, or an array of points on Earth, and while
     they must be created in a coordinate system (say, ECEF, geodetic, or
     some local tangent system), transformation to other coordinate systems
     do not require explicit reference to their internal representation.

     Needless-to-say, the concept of coordinate-independent objects conforms
     well with
     <a href="https://en.wikipedia.org/wiki/Object-oriented_programming">
     object-oriented programming</a>. While the
     <a href="https://en.wikipedia.org/wiki/Procedural_programming">
     procedural</a> and/or
     <a href="https://en.wikipedia.org/wiki/Imperative_programming">
     imperative</a> programmer may
     represent
     <a href="https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory">
     JPL's</a> coordinates as:

     \verbatim jpl = [[34, 12, 6.1], [-118, 10, 18], 250] \endverbatim
     he has to remember that address 0 contains 2 integers (latitude degrees
     and minutes) and a float (seconds) and so on, along with which 
     ellipsoid we're using, the OO user just says:
     
     \verbatim 
     
     jpl = WGS84.LLH(DMS(34, 12, 6.1), DMS(-118, 10, 18), 210)
     
     \endverbatim
     
     and the object knows, and remembers, exactly what that all means-- and
     is much more than just records in a structure, because when you write 

     \verbatim  
     
     d = jpl - gsfc
     
     \endverbatim
     
     it knows that the difference between 2 points in an affine space is
     a vector, and it knows how to compute it, regardless of how the
     coordinates of NASA Goddard were specified. That is, it might be that:

	\verbatim

	gsfc = almanac.ALMANAC["AIRY 1830"].LLH(38.9915337403, -76.85, 548)

	\endverbatim

which is a different coordinate system with respect to a different
      ellipsoid -- but the meaning of the point does not depend on its
     representation. Hence:

	\verbatim

>>>print repr(jpl - gsfc)

geo.metric.euclid.vector.Vector(-3622444.0, 177651.0, -426628.0)


	\endverbatim

     is the ECEF vector connecting them.
     Furthermore, if you add a vector to a point:

     \verbatim p = jpl + v
     \endverbatim
     it doesn't matter if _v_ is represented in Cartesian, polar, or
     cylindrical coordinates- because that doesn't change the meaning of _v_,
     just the representation. Likewise, the affine point could be in ECEF,
     and the results are the same:

     \verbatim  >>>jpl + v == jpl.ecef() + v.polar() 
     	True \endverbatim

     (It is the duty of the OO-programmer to implement these choices
     without _if_ - _then_ clauses: unlike functions acting on arrays,
     objects know themselves and thus
     <a href="https://en.wikipedia.org/wiki/Dynamic_dispatch">
     dispatch</a> the right code
     <a href="https://en.wikipedia.org/wiki/Late_binding">when
     called upon</a> to do so.)
     
     Operations are independent of representations. (NB: that's, little _r_
     representations. geo/ allows various group-theoretic big _R_
     Representations, and those you should not mix--but who would do that
     anyway?).

     @section xvec Geo-Metric:
     There are 2 main and 2 (or more) ancillary sub-packages:

     @subsection euclid Euclid and Euler

     @subsubsection svt Cartesian Vectors in R3.
     The geo.metric package (euclid.py) defines the basics for rank 0
     (scalar.py), 1 (vector.py),
     2 (tensor.py), and higher objects in real Euclidean
     \f$\mathbb E^3\f$ space, aka \f$\mathbb R^3\f$.

     @subsubsection svt2 Transformations in SO(3)
     Matrix __SO(3)__  (via tensor.py), quaternion __SU(2)__  (euler/hamilton.py),
     transformations are implemented.
     Furthermore, various
     representation of the charts on __SO(3)__ (charts.py) allow Euler-angle
     and Tait-Bryant angle based representations. Needless-to-say,
     transformation of vectors is independent of the representation
     of the transform, _T_:

     \verbatim v' = T(v)
     \endverbatim

     for all classes of _T_. That is, T could be a "tensor", euler angle,
     versor, etc. Moreover, their argument doesn't need to be a vector:

     \verbatim
s' = T(v*v)  # == T(v) * T(v) == v*v          Rotates a scalar
t' = T(v&v)  # == T(V) & T(v) == T*(v&v)*T.T  Rotates a rank-2 tensor.
     \endverbatim

     It is an object's job, as a polymorphic function emulator, to use the
     correct method to transform its argument.

     @subsubsection svt3 Transformations GL(3, R) (aka Aff(3))
     In the context of geo, an
     <a href="https://en.wikipedia.org/wiki/Affine_transformation">affine
     transformation</a> is a general linear transformations from
     \f$\mathbb R^3\f$ to \f$\mathbb R^3\f$:
     \n\n
     \f$ A(\vec x) = {\bf M}\vec x + \vec b \f$
     \n\n
     and are supported in affine.py. While __M__ can be any object that
     transforms vectors linearly, we are primarily concerned with rotations.
     The translation __b__ is always from the target space (which is trivial
     in geo's implementation). With general forms of __M__, one can
     implement
     <a href="https://en.wikipedia.org/wiki/Scaling_(geometry)">scalings</a>,
     <a href="https://en.wikipedia.org/wiki/Similarity_(geometry)">similiarity
     transformations</a>, reflections,
     <a href="https://en.wikipedia.org/wiki/Shear_mapping">shear
     mappings</a>, pure translations,
     <a href="https://en.wikipedia.org/wiki/Homothetic_transformation">
     homotheties</a>, and rotations. The 3 pure
     infitesimal translations and their duals (3 infitesimal rotations)
     together form the so-called generators of rigid Euclidean space.

     @subsubsection ncc Non-Cartesian Coordinates.
     Support for non-Cartesian vector representations
     is provided. Spherical (aka polar)
     (geo.metric.euclid.kant.spherical.Polar), cylindrical
     (geo.metric.euclid.kant.cylindrical.Cylinder), and parabolic
     (geo.metric.euclid.kant.parabolic.Parabola) coordinates
     are implemented. Through the wonders or __OOP__ and
     <a href="https://en.wikipedia.org/wiki/First-class_function">1st-class
     functions</a>
     you can dynamically create fully operational non-Cartesian 
     coordinate representations via the metaclass
     geo.metric.euclid.kant.CoordinateFactory. All you need to do
     is provide a function and its inverse (that is, to and from
     Cartesian coordinates).

     So the Cartessian basis vector can be expressed in other coordinate
     systems:
     \verbatim
>>>print X.polar()	
radius=1.0, theta=1.57079632679, phi=0.0

>>>print Y.cylinder()
rho=1.0, theta=1.57079632679, z=0.0

>>>print Z.parabola()
u=1.41421356237, v=0.0, theta=0.0
     \endverbatim

     These vector classes are
     <a href="https://en.wikipedia.org/wiki/Adapter_pattern">
     adapter pattern</a>
     interfaces to the Cartesian Vector: that is, when used in operations,
     they convert themselves to Cartesian, do the operation, and then
     convert back to their original coordinate system. That allows
     you to mix coordinates systems seemlessly, because a vector is a
     vector, regardless of the coordinate system in which it is expressed.
     Hence \f$ \hat x \times \hat y = \hat z \f$ in any coordainte system:
     \verbatim
>>>print X.polar() ^ Y.cylinder()
radius=1.0, theta=0.0, phi=3.14159265359

>>>print (X.polar() ^ Y.cylinder()) * Z
w=1.0
     \endverbatim

     @subsubsection reiman Non-Euclidean Geometry
     While it is tempting to go fully Non-euclidean, the
     computational and architectural hit is too high for every day purposes.
     (That is: every dot product becomes a double contraction with the
     metric tensor; architecurally, every vector needs it's own set of
     metric co/contra/mixed metric tensors.)
     Nevertheless,
     their is an experimental implementation (riemann.py), where Vectors
     are attached to a metric space.

     @subsubsection ssc Space Curves
     Space curves and the Frenet-Serret formalism is available in
     (motion.py). The space curve comprises an (array_like) Vector and a
     one dimensional time array, which serves as a parameter in the
     parameterization of the 3D curve.

     Take an array like Vector representing a helical path:
     \verbatim
     >>>t = linspace(0, 2*pi, 100000)  # parameter
     >>>v = X*sin(3*t) + Y*cos(3*t) + Z*t/1000
     >>>print type(v)
     <class 'geo.metric.euclid.vector.Vector'>
     \endverbatim
     That array-like Vector can be made into a spacecurve by formalizing
     the parameter $t$:
     \verbatim
     >>>v = v.spacecurve(t)
     >>>print type(v)
     <class 'geo.metric.frenet_serret.motion.SpaceCurve'>
     \endverbatim
     The SpaceCurve now has veclocity, acceleration (linear and angular),
     torsion and curvature (and radii thereof), along with the Tangent,
     Normal, and Binormal coordinate system at every point, e.g.:
     \verbatim
     >>>print v.TNB()[100000//6].broadcast(round)
[-1.0, -0.0, 0.0]
[-0.0, 1.0, 0.0]
[-0.0, -0.0, -1.0]
     \endverbatim     			
     where the result has been rounded for clarity.
     
     @subsubsection complex Complexification
     The bulk of geo is all about real vectors spaces: while vectors
     can have complex components, they still live in \f$\mathbb R^3\f$,
     and as such,
     their bilinear inner product is not positive definite:
     \verbatim
     >>>v = X + 1j*Y
     >>>print v*v
     w=0j
     \endverbatim
     In gibbs.py, the
     real vector.Vector objects are placed into
     \f$\mathbb C^3\f$: they are the
     same vectors, but their sesqui-linear inner product is positive
     definite. The outer product is also complexified, making them
     the ideal tool for representing polarization density matrices
     (faraday.py).
     \verbatim
     >>>v = v.gibbs()  # same vector, now lives in a new space
     >>>print type(v)
     <class 'geo.metric.euclid.gibbs.Gibbs'>

     >>>print v*v  # square is positive definite now
     w=(2+0j)

     >>> print v&v  #  outer product is Hermetian
     [(1+0j), -1j, 0j]
     [1j, (1+0j), 0j]
     [0j, 0j, 0j]
     \endverbatim

     @subsubsection ho Higher Rank Tensors
     3rd and 4th rank tensors, and a class-factory capable of making
     any order tensor are available in euclid/three.py and
     euclid/four.py, respectively.
     Inner, outer, and wedge products
     are supported, along with index contraction and other forms of
     index gymnastics.

     In practice, rank 2 and rank 4 tensors usually have some degree of
     symmetry;
     voigt.py implements a traditional method for dealing with
     major and or minor symmetric fourth rank tensors such as electric
     susceptibility and elasticity, .e.g:
     >>>print voigt.voigt(((X&X&X&X) +2*(Y&Y&Y&Y) + 3*(Z&Z&Z&Z) - 4*(X&X&Y&Y))
     ...: )
[[ 1. -4.  0.  0.  0.  0.]
 [ 0.  2.  0.  0.  0.  0.]
 [ 0.  0.  3.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.]
 [ 0.  0.  0.  0.  0.  0.]]
     \endverbatim()

     @subsubsection albert Einstein Summation Notation
     Through the wonders of python and
     <a href="https://docs.python.org/2/reference/datamodel.html#customizing-attribute-access">
     complete class customization</a>,
     <a href="https://en.wikipedia.org/wiki/Einstein_notation">Einstein
     summation notation</a> is _fully_ supported (einstein/albert.py).
     For example,
     a second rank Tensor only has 9 attributes: 'xx', 'xy', 'xz',
     'yx', 'yy', 'yz', 'zx', 'zy', and 'zz', and it has a trace
     method:

     \verbatim
     s = T.trace()
     \endverbatim

     which explicitly computes:
     
     \verbatim
     s = T.xx + T.yy + T.zz
     \endverbatim

     which is indeed the trace of the tensor. One may also code:
     
     \verbatim
     s = T.ii
     \endverbatim

     Of course, their is no 'ii' attribute (or property). Nevertheless,
     geo regonizes this as Einstein summation and goes ahead and computes
     the trace. Likewise for the transpose:

     \verbatim
     t_ji = T.ji  # same as T.T, or T.transpose(1,0)
     \endverbatim

     (The reversed alphabetical order indicates transposition). Rows and columns
     can also be extracted (as Vectors):

     \verbatim
     row0 = T.xi
     col2 = T.ix
     \endverbatim

     This is the preferred method of extracting "rows" and "columns" because
     Tensor don't have rows and columns, and the names are antithetical to
     geo's ethos.

     Of course, this seems rather cute for rank-2 Tensors (and down right
     silly for Vectors), it is quite powerful when dealing with high
     rank tensors: You really can write code that looks _exactly_ _like_
     the equations in your text book. For example:
     \n
     The explicit cross product (as the partial trace of a rank-5
     tensor):
     \n\n
     \f$ ({\bf{\vec a  \times \vec b}})_i \equiv \epsilon_{ijk}a_jb_k =
     ({\bf{\epsilon \vec a\vec b}})_{ijkjk} \f$
     \n\n
     can be coded as:

     \verbatim
     a_cross_b = (EPSILON & a & b).ijkjk
     \endverbatim

     Similarly, the determinant
     \n
     \f$ \det{\bf T} = \left| \begin{array}{ccc}
     T_{xx} & T_{xy} & T_{xz} \\
     T_{yx} & T_{yy} & T_{yz} \\
     T_{zx} & T_{zy} & T_{zz}
     \end{array} \right| =
     \epsilon_{ijk}T_{xi}T_{yj}T_{zk} \f$
     \n\n     
     can be computed via triple contraction of a rank-6
     tensor (with explicit labeling of the fixed indices):
     \verbatim
     det_t = (EPISLON & T.xi & T.yi & T.zi).ijkijk
     \endverbatim
     Another formulation as 
     the  6-fold contraction of a rank-12 tensor;
     \n\n
     \f$\det(T)=\frac{1}{6}\epsilon_{ijk}\epsilon_{lmn}T_{il}T_{jm}T_{kn}\f$
     \n\n
     is coded as:
     \verbatim
     det_t = (EPSILON & EPSILON & T & T & T).ijklmniljmkn / 6
     \endverbatim

     Note that the code matches math equations _exactly_. The latter
     implicitly runs a 12-deep loop to build and then contract the
     531,441 components of a rank-12 tensor--while it is not the most
     computationaly efficient algorithm, it does demonstrate the
     simplicity of Einstein summation notation. Moreover, it is not
     useless, as formulations with similar levels of complexity arise
     when considering the rotational symmetries  of rank-3+ tensors...


     @subsection irreps Irreducible Representations and Subspaces
     Here we delve into a rather involved topic: the irreducible
     representations and subspaces of the general Cartesian tensor.
     An irreducible subspace of a tensor is some linear combination of
     a tensor and its transposes that is closed under rotations. For example,
     an antisymmetric rank-2 tensor remains an antisymmetric rank-2 tensor
     under all rotations: it is a irreducible subspace with 3 non-zero
     components and 6 components that are zero (in any reference frame).
     As an irreducible representation, those 3 components rotate like
     a normal Cartesian vector. Likewise, a tensors trace is invariant under
     rotations: it behaves like a scalar.

     


     @subsubsection we Spherical Tensors: 3 x 3 = 5 + 3 + 1
     This does not refer to spherical coordinates, rather the decomposition
     of Vectors (wigner/eckart/one.py) and
     Tensors (wigner/eckart/two.py) into irreducible representations
     of __SO(3)__.
     Thus, the _x_ or _xy_-like
     attributes are replaced by
     orbital and
     azimuthal quantum numbers: \f$(l, m)\f$-- also knows as degree and order.
     The signifcance of these is that they are eigentensors of rotations
     about a fixed axis.

     For rank-1, the 3 spherical basis states are eigenvectors of z-rotations.
     This is similar to Cartesian vector addition being isomorphic to
     the translation operators acting on position eigenstates (which are plain
     old position vectors). Because translation commute, the construction
     of higher rank tensors is straightforward. A rank-N basis state is just
     the polyadic product of N basis vectors (e.g.
     \f$ {\bf \hat{x}\hat{x}\hat{y}}\f$). Rotations do not commute,
     so polyadic products of spherical basis vectors are not necessarily
     eigentensors of rotations. This is entirely analgous to the construction
     of
     <a href="https://en.wikipedia.org/wiki/Spherical_harmonics#Spherical_harmonics_in_Cartesian_form">spherical harmonics in Cartesian coordinates</a>,
     and it gets complicated, fast. Nevertheless, geo goes there--though
     it requires and entirely complex subpackage (geo.metric.wigner).

     In geo's implementation, the elements of spherical tensors are not
     accessed by
     attributes, rather, the instances are containers
     of irreducible
     representations and their azimuthal eigenstates.
     Methods allow conversions between
     the 2 representations (ec.py, ka.py, and rt.py).
     So for instance:
     \verbatim
     
     >>>M = Tensor(xx, xy, xz, ..., zz)
     >>>print M.xx, M.xy, M.xz, ..., M.zz
     >>>
     >>>T = M.spherical()
     >>>print T[0, 0], T[1, -1], T[2, -2], ..., T[2, 2]

     \endverbatim
     are represenations of the same geometric object, as very different
     python objects.

     Applications to polarization observables are
     demonstrated in faraday.py .

     The full transformation matrices associated with
     spherical representations of any dimensions 2j+1 for integer and
     1/2-integer representations (j)  are in wigner.py, which also
     implements the decomposition of tensor-product spaces into
     sums of irreducible spaces, with the help of clebsch_gordan.py. 
     Some assistance with decomposition of representations provided in
     racah.py, where products of kets (eigenstates of individual
     angular momentum/z-projection operators) are turned into sums
     over different kets (eigenstates of the total angular momentum
     and total z-projection). [This is truly dialed-in python, _Eds._]

     The preliminary implementation of rank 3 (wigner/eckart/three.py) and
     rank 4 (wigner/eckart/four.py) is underway; however, these are difficult
     (i.e., publishable) to implement.

     Spherical decomposition of complete Cartesian tensors requires and
     understanding of isotropic tensors, which are available up to 
     rank 6 in syzygy.py: these
     are tensors that are independent of coordinate system, which for
     rank 2 and 3 are the familiar Kroncker Delta
     (geo.metric.euclid.tensor.DELTA) and Levi-Civitta Tensor
     (geo.metric.euclid.three.EPSILON), repsectively.
     For higher ranks there are more than one isotropic tensor per rank,
     and they are not all linearaly independent.

     Additional tools for the alternating tensor, aka Levi-Civita (pseudo)-
     tensor/symbol are in levi_civita.py. The
     <a href="https://en.wikipedia.org/wiki/Kronecker_delta#Generalizations_of_the_Kronecker_delta">
     generalized Kronecker deltas</a>
     can be computed in kronecker.py.
     

     @subsubsection symrep Schur-Weyl Duality
     Rank-2 tensors have two subspaces closed under rotations which
     are defined by their symmetric and antisymmetric parts. These
     are constructed by symmetric and antisymmetric permutations of their
     two indices--which is rather simple-as there are only 2 permutations.
     Higher rank _n_ tensors with _n_ indices have more complicated
     subspaces,and they can be computed from the various _n!_ permutations
     of their indices. This is the domain of
     <a href="https://en.wikipedia.org/wiki/Schur–Weyl_duality">Schur-Weyl
     Duality</a>: 
     The symmetric group on _n_-letters,
     \f$ Sym(n) \f$ can be decomposed into irreducible representations based
     on the integer partitions of \f$ n \f$ (really?). So to that end,
     there is monte.py (rimshot), which provides classes for the
     <a href="https://en.wikipedia.org/wiki/Symmetric_group">
     symmetric group</a>,
     <a href="https://en.wikipedia.org/wiki/Permutation">permutations</a>,
     and
     <a href="https://en.wikipedia.org/wiki/Cyclic_permutation">cycles</a>.

     The
     <a href="https://en.wikipedia.org/wiki/Robinson–Schensted_correspondence">
     Robinson-Schensted Correspondence</a> then relates irreducible
     representations of the symmetric group to
     <a href="https://en.wikipedia.org/wiki/Young_tableau">Young Tableaux</a>
     (young.py). They provides a
     link to partitions of integers, which is taken up in pascal.py,
     which starts with Pochammer symbols and winds up with integer
     partitions. Additionaly, the number triangles are extended to
     many of the combinatoric polynomials, leading naturally to the
     computation of Jacobi
     polynomials, which are required for the Wigner D matrices.
     Also: since Young diagrams are like a multiset, a multiset object
     is available in de_bruijn.py, but it's really just a hashable
     <a href="https://docs.python.org/2/library/collections.html#collections.Counter">
     collections.Counter.</a>
     (The theory of conjugacy classes in finite
     groups necessitates multisets with keys that are themselves multisets).

     This subpackage also supports some results from the representation
     theory of finte groups, permutations, and combinatorics. The
     Young Tableaux formalism also has applications to the representation
     of Lie groups, which are commonly used to understand the quark
     (and separately the gluon) structure of matter
     (gelfand_tsetlin/__init__.py).

     @section yvec Geo-Detic:

     @subsection coo Coordinates and Ellipsoids
     The geo.detic package implements geo.metric objects on the Earth
     via the introduction of coordinates (coordinates.py), which can
     be represented in Earth Centered Earth Fixed
     (geo.detic.coordinates.ECEF), geodetic
     (geo.detic.coordinates.LLH), or local tangent plane
     (geo.detic.coordinates.LTP)
     or sphere (geo.detic.coordinates.SCH) relative to a
     geo.detic.origins.PegPoint.

     The coordinates are
     <a href="https://docs.python.org/2/library/abc.html">abstract</a>,
     and only have
     <a href="https://en.wikipedia.org/wiki/Class_(computer_programming)#Abstract_and_concrete">
     concrete</a> meaning with respect
     to an ellipsoid (ellipsoid.py). 
     Although preference is given to WGS84 (wgs84.py), a variety of
     historic ellipsoids are available (almanac.py), as well as user-defined
     ellipsoids and others in the Solar System (copernicus.py).
     
     Analytic computations of covariant and contravariant basis vectors
     linking LLH and SCH to their canonical Cartesian coordinates
     systems (ECEF and LTP, respectively) are included (christoffel.py).
     That is: you always have both the tangent basis vectors
     (_e.g._, \f$
     (\frac{\partial  \vec r}{\partial s},
     \frac{\partial  \vec r}{\partial c},
     \frac{\partial  \vec r}{\partial h})\f$
     which are not-necessarity
     normailzed vectors running along the direction of change of s, c, and h,
     and the cotangent basis vectors
     (_e.g._ \f$
     {\bf\vec{\nabla}}s,
     {\bf\vec{\nabla}}c,
     {\bf\vec{\nabla}}h)
     \f$), which are the gradients of the local coordinates (that is, they're
     normal to surfaces of constants s, c, and h--and likewise for latitude,
     longitude and height with resepct to ECEF).
     
     @subsection The Geoids     
     Support for interpolation
     on geoids is provided (egm08.py, egm96.py). While the former
     requires the
     <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/">
     1 minute grid</a>, the latter has the 1-degree grid self-contained:
     __EGM96__ is contained within this package.


     @subsection utm Universal Transverse Mercator
     UTM (utm.py, kruger.py) support is under development. A class for seamless
     use of DD MM SS.ssss format is provided (preliminarily) in dms.py

     
     @section four Geo-Desic
     Because it's there.
     
     The geo.desic sub-package goes beyond 6-DoF vector transformations in
     __E3__ with a fully developed Clifford algebra (clifford.py), the complete
     quaternion group __H__ (geo/desic/hamilton.py), and dynamic creation of
     Cayley-Dickson extension algebras (cayley_dickson.py).
     Preliminary implementation of 
     minkowski.py space(time) and lorentz.py transformations (with
     application to radar in de_broglie.py). Non-relativistic
     and relativistic spinor wrappers (pauli.py and dirac.py) available--for
     __all__ finite-dimensional irreducible
     representations (whaatttt?).
     Needless-to-say, a spinor representation is given
     in cartan.py.
  
     You'll also note SO3.py: it only has 2 lines of computation (2!),
     but what it does is amazing. You assign a vector to each of
     the generators of rotations in the fundamental rep (the complex
     2x2 matrices that rotate 2 component spinors), and then you write
     down the matrices that express the 2x2 matrices' commutations
     relations: boom, you have just derived the real adjoint
     representation, and those matrices ARE the generators of
     infinitesimal rotations-- plop them into a \f$ \exp{i\ldots}
     \f$, and you have rotation matrices. It's a little confusing,
     because the 2 component spinor are fundamental, while one would
     think the \f$2^2-1=3\f$ are-- just because their are 3 generators
     of the algebra.. it's the proliferation of "3"s that have
     totally different meanings that causes confusion (at least for me).

     Preliminary implementation of Lie algebras (lie.py), and SU(3)
     (gell_mann.py) are dubious, at best-- but you can convert the
     3-component fundamental representation (per the prior paragraph)
     into the \f$3^2-1=8\f$ adjoint representation
     (geo.desic.hilbert.gell_mann.EIGHTFOLD_WAY)
     and that is why
     3-color charges need 8 gluons, while \f$ 2 \f$ (weak) isospins
     need 3 pions (W bosons). (This package also makes clear the difference
     between fundamental and adjoint representations, which is a bit
     obtuse in the rotation group).

     @section comp Computational Implementation
     The point of this package is to provide support for using manifestly
     covariant formulations in __E3__ as applied to geodesy.

     The idea of manifest covariance is important: scalars and not numbers,
     vectors are not arrays, and tensors and not matrices. Each is a geometric
     object who's essence does not depend on any particular representation
     in a coordinate system; moreover, they have well defined transformational
     properties based on their rank.

     Similarly, when applied to geodesy: points on the Earth and their,
     relation to each other do not depend any particular coordinate system
     choosen to represent them.

     Thus: relations and equation should be expressible in term that do
     not depend on coordinates (manifest covariance); nevertheless, when
     implementing the equations, a representation must be chosen: the computer
     needs numbers.

     This package is the tool to do just that: one can create geometric
     objects in suitable coordinates, do covariant operations on them,
     and then express the results in a suitable coordinate system.
     
     ARRAYS vs. VECTORS:
     -------------------
     This package is a fully polymorphic object oriented vector and coordinate
     package, focused on geo-location on an arbitrary oblate ellipsoid of
     revolution.

     It is for use with large data sets represented by
     <a href="http://www.numpy.org">numpy.ndarrays</a>-
     or not- numpy is only required if you use it. (The numpy-interface
     is factored out into mix-in classes (arraylike.py), and functions
     chosen at import time (trig.py). (You do not need numpy installed
     to use the core package).)

     It is for use with manifestly covariant vector transformations: that is:
     the user does operations with vectors and tensors-- not their
     components. Yes, you start off with components, and you finish with
     components, but in between-- they don't really matter.
     
     RESULT: vectors have attributes (_x_, _y_, ...) as do
     coordinates (_lat_, _lon_, ...)- these are NOT, and NEVER WILL BE,
     represented as sequences
     ([0], [1], ...). That you would do it in lesser languages such as Fortran,
     or IDL, or Matlab(TM) is irrelevant--don't think about, even though you're
     natural urge, which evolved over years of primitive obsession, is to do
     so. Fight it.

     Array, sequence, container, and iterator properties are ALWAYS separated
     from vector, tensor, coordinate, Euler angle, quaternion properties: THEY
     NEVER MIX. EVER. EVER. EVER. EVER.

     __EVER__.

		 ....XXXXXXXX XX     XX XXXXXXXX XXXXXXXX  
     		 XX       XX     XX XX       XX     XX 
     		 XX       XX     XX XX       XX     XX 
     		 XXXXXX   XX     XX XXXXXX   XXXXXXXX  
     		 XX        XX   XX  XX       XX   XX   
     		 XX         XX XX   XX       XX    XX  
     		 XXXXXXXX    XXX    XXXXXXXX XX     XX 

   
     Needless-to-say, there is a caveat for spherical tensors--which have
     NO attributes, as they _are_ containers of eigenstates.
		 
     POLYMORPHISM:
     ------------
     This is the point of OO, and it is strongly supported. The point is, you
     DO NOT need to know which kind of coordinate you have when using it. For
     example if you have data, "platform" referenced to a coordinate system,
     and you NEED that data in geodetic coordinates (LLH), you just call the
     method:

     \verbatim	    >>>platform.llh()
     \endverbatim

     It does not matter if "platform" is in LLH already, or in ECEF, or in SCH,
     or in a tangent plane coordinate system. Moreover, if you have another
     point, "target"

     \verbatim	    >>>platform - target
     \endverbatim

     is the vector joining them: REGARDLESS OF THEIR COORDINATE SYSTEM-- you
     do not need to specify it-- the vector will be computed in the Cartesian
     coordinate system at platform's origin (even if the coordinate system is
     not Cartesian). This because coordinates are treated elements of
     an affine space--
     vector spaces without an origin.

     All coordinates are defined on an ellipsoid-- and they carry that
     information with them. You can translate to another ellipsoid, e.g.,
     AIR1830, via:

     \verbatim	    >>>platform | almanac.AIR1830
     \endverbatim

     Again-- you don't need to specify if "platform" is LLH, SCH, or the local
     tangent plane (LTP)-- "platform" knows that already and converts itself
     to the same sub-class of coordinates. That's the point of polymorphism.

     First Class Transformation:
     ---------------------------
     _Question_:

     Why are their so few functions for coordinate transformations? Why are
     there no 3x3 matrices filled with (bug prone) long combinations of trig
     operations?

     _Answer_:

     Transformation are 1st class object that permit arithmetic operations.
     Hence, you build a a full 3-DoF rotation as a product of 3 simple 1-DoF
     rotations. With the affine transformation, you can build any 6-DoF
     coordinate transformation out of its building blocks-- and compose
     them as they are needed.

     With that, most coordinate transformations are 1 to 3 lines of code,
     and can be easily understood. This avoids an important antipattern:
     "long-method" -- which claims any method over 7 lines is too long--
     it becomes hard to understand, hard to maintain, and hard to extend.
     
     Furthermore, the composition and inversion of transformations is
     handled with
     <a href="https://docs.python.org/2/reference/datamodel.html#emulating-numeric-types">operator overloading</a>
     of "*" ("__mul__")and "~"
     ("__invert__"),
     respectively, with some support for interpolating between transformations
     via overloading "**" ("__pow__") with non-integer exponents. (Python's
     built-in <a href="https://docs.python.org/2/library/functions.html#pow">
     pow</a> function has a 3 argument option, which implements
     <a href="http://en.wikipedia.org/wiki/Slerp">spherical linear
     interpolation</a>.)
     
     Array or Singleton:
     ====================
     In the above example, "platform" or "target" could be singletons or arrays
     --it doesn't matter. The resultant vector will match the inputs-- without
     you having to tell it. If the objects are made out of numpy-arrays, then
     the result will be too-- numpy does all the work-- neither the user nor
     the developer has to worry about loop indices when using
     _or_ writing code.
    
     With that in mind, here are the details:

     @section Introduction
     Vector exists so that coordinates can exist. The ellipsoid exists so that
     the coordinates make sense. The goal of this sub-package is to make ALL
     vector and coordinate operations simple overloaded math operations or
     method calls.

     Moreover, the objects are completely divorced from computer "arrays", and
     hence, you can use a single instance to represent a bunch of vectors, or
     points on a map, or orientations, etc., and ALL the broadcasting is done
     by numpy. This is a very important point to digest, for a vector, v, or
     tensor T:

\verbatim     	    >>> v[2], T[0, 1]
\endverbatim
       have NO MEANING with respect to vectors and tensors. They are not:

\verbatim     	    >>> v.z, T.xy
\endverbatim
     V[2] only has meaning if v.x[2], v.y[2], and v.z[2] have meaning, and then
     it should be obvious what it means. Same for T.xy[0,1] and the other 8
     rank-2 tensor attributes.

     What this means is that you are free to make Vector, Tensor,
     quaternion, Euler-angle triplets, and coordinate objects out of
     arbitrarily large numpy arrays of any dimension. That they are captured
     in a single Tensor or Coordinate object should not make you think you have
     a singleton. A Tensor object can contain an entire yaw, pitch,and roll
     history of a platform, while a single Vector can represent the look
     vector to every post in an 2-D Digital Elevation Model.

     This is an important consideration. Tensors are not defined by their
     number of elements-- they are defined by their transformation properties.
     That is, a Vector object can be one vector, or N vectors arranged in
     some one-or-more numpy.ndarray shape, or it can be a generator or
     iterator that produces vectors-- they ALL have the same transformation
     properties, and hence are all Vector objects with the same interface.
     
     (Note that a vector can be a single vector, and still have non-singleton
     elements. Take for example, geo.metric.wigner.eckart.J(), which is
     indeed a single vector, yet it's components are Tensors, such that
     \f$ e^{i{\bf\vec J\cdot{\hat n}}\phi} \f$ generates rotations
     about a unit vector \f$ {\bf \hat n} \f$. Even worse is
     geo.detic.pauli.J(). Here the vector's components are numpy matrices
     of dimension \f$2j+1\f$, represented internal-degrees of freedom.
     This is all implemented polymorphically without any conditionals
     (_if_ _then_); you simply cannot do this with your procedural Matlab
     arrays and matrices.)

     DO NOT CONFUSE TENSOR PROPERTIES WITH ARRAY PROPERTIES:

     "__getattr__"  gets tensor properties, e.g:
     \verbatim m.xy  \endverbatim\n
     "__getitem__"  gets array properties,  e.g:
     \verbatim v[10:1000:3, -30:]
     \endverbatim \n\n
     The payoff:
     
     The user/developer should NEVER have to do a loop to transform ENTIRE
     motion histories or 2-d images.
     
     Of course you may want a loop, even an implied loop. For those who like
     python's iterator-pattern and support for functional programming style,
     all vector/ tensor/coordinate objects can be created with lazy-evaluation
     attributes. That is, you can make them out of iterators, generators,
     generator-expression, arrays, lists, etc, and then, using the "__iter__"()
     and next() method (i.e., the iter() and next() built-ins) you can cycle
     through an object and stream it to various functions or transformations.
     A Vector's vector-like interface to the world is not affected by the
     computer-science-like interfaces to its components.
     

     @section euclid Euclid
     Vectors (and other Tensors) in 3-Dimensional Euclidean Space (__E3__)

     Vectors are important objects. We're talking about a physicist's vectors.
     That's a geometric object that lives in Euclidean 3-space (for us) and does
     NOT depend on the reference frame. That "vectorized" and "vector" have
     meaning in computer science is entirely unrelated and is just a point of
     confusion. Try to FORGET IT.

     There are several ways to consider vectors. One way is as a rank-1
     tensor, that is, as a linear function from
     \f$ {\mathbb R}^3 \rightarrow {\mathbb R}\f$--this is a very powerful view
     point for characterizing physical laws--and is supported in geo:
     \n\n
     \f$ {\bf \vec a(\vec b} \equiv {\bf \vec a\cdot \vec b} \f$\\n\n.

     Nevertheless, we'll define them by how they transform.
     Vectors are defined by their transformation properties: you have to rotate
     them by 360 degrees to remain unchanged. Well, that is entirely useless
     for computer programming. What is useful is that:

     (1) They transform as:\n\n
			\f$	v'_i = T_{ij}  v_j  \f$
			\n\n
     (2) In a given reference frame, they can be represented by 3 components:
                      \f$v_i\f$ for \f$i = 1, 2, 3.\f$
     So there you go. Vectors have 3 things, and that defines their 
     representation in a frame. That leads to an important pythonic note:

     Vector (Tensor, Coordinate, etc..) objects are entirely defined by their
     "__init__" arguments. To that end, a vector is given
     by its x, y, and z attributes:

\verbatim     	       	  >>>v = Vector(x, y, z)
\endverbatim
     There shall be no nullary instantiation followed by setter calls:

\verbatim
	>>>v = Vector()  # NO
	>>>v.setX(x)     # NO
	>>>v.setY(y)     # NO
	>>>v.setZ(z)     # NO NO NO
\endverbatim
     Fair enough, now we have a vector. What do you do with that? Well, you
     don't call a function. The vector should have everything you want to do
     with it defined by methods and operators (magic methods):
     
     Magic Methods:
     =============
     The following invoke magic methods:
     
\verbatim  +v = v.__pos__()
\endverbatim
     Doesn't do anything.It returns "v", though it could return a copy of v.

\verbatim     	     -v = v.__neg__()
\endverbatim
     Negation, implemented as negation of each component

\verbatim     	    v + v' = v.__add__(v')
\endverbatim
     vector addition [Note: no constants allowed]:
\verbatim     	    v + 7 = v.__add__(7)
\endverbatim
raises an
     geo.utils.exceptions.NonCovariantOperation error (exceptions.py)]

\verbatim           v - v' = v.__sub__(v')
\endverbatim
     vector subtraction 

\verbatim
v/4	  =  v.__div__(4)
v*0.25	  =  v.__mul__(0.25)    
0.25*v    =  v.__rmul__(0.25)
\endverbatim
are all dilations. Reflected multiplication is the simplest, as it can
_only_ be a dilations

\verbatim	     v * v' = v.__mul__(v')
\endverbatim
     Scalar product

\verbatim     	    v**2 = v.__pow__(2)
\endverbatim
    "__pow__" takes any positive integer, though "2" is the only defensible
    argument

\verbatim abs(v) = v.__abs__() -> ||v**2|| ->  v.L2norm() = v.norm(l=2)
\endverbatim
     all need to be defined, and some don't, like:

\verbatim     	 v ^ v' = v.__xor__(v')
\endverbatim
    cross (wedge) product

\verbatim     	  v & v' = v.__and__(v')
\endverbatim
    outer (dyad) product

\verbatim     	  ~v   =  v.__invert__() = v / v*v
\endverbatim
    is a special case for converting covariant and contravariant vectors
    in local orthonormal bases and their conjugate global basses.
    \n
    Projection, Rejections, Reflection:
\verbatim
a >> b = a.__rshift__(b) = a.projection(b) = Proj(b)(a)
a << b = a.__lshift__(b) = a.rejection(b) = Rej(b)(a)
a | b  = a.__or__(b)     = a.reflection(b) = Ref(b)(a)
\endverbatim
	Respectively. (Projection, rejection, and reflection operations
	for vectors, tensors, and quaternions are unified in cauchy.py.)
	The right column are tensor operators:\n
	geo.metric.euclid.vector.Proj'\n
        geo.metric.euclid.vector.Rej'\n
	geo.metric.euclid.vector.Ref.\n
\n
        Moreover, what about "__getitem__":

\verbatim
v[3]                = v.__getitem__(3)
v[-1]               = v.__getitem__(3)
v[start:stop:step]  = v.__getitem__(slice(start, stop, step))
\endverbatim
     Do those have meaning? Yes they do, but they HAVE NOTHING TO DO WITH
     VECTORS. That is:

\verbatim	v.x is not v[0]\endverbatim
\verbatim       v.y is not v[1]\endverbatim
\verbatim       v.z is not v[2]\endverbatim
     No.
\verbatim
     v[start:stop:step] is
     Vector(v.x[start:stop:step], v.y[start:stop:step], v.z[start:stop:step])
\endverbatim
     This is not just good, but great. Moreover, it is true for all tensors and
     all coordinate classes. It is the key for making all objects work with or
     without numpy in a manner that allows the user or developer to never have
     to consider it. You can make objects out of singletons or multi-
     dimensional arrays. Your operations will be broadcast properly, and the
     resultant Tensor objects will be singletons or numpy arrays, as needed-- 
     and
     of course, you can mix and match singleton and array arguments in
     expressions-- again it's all about Tensors controlling their own destiny
     while leaving unto numpy what is numpy's.

     Thus, if you iterate over a multidimensional (in the computer array sense)
     vector, your iterator will return vectors one by one (or arrays of vectors-
     it's all passed to numpy or whatever controls the attributes' iterators).

     At A Deeper Level
     =================
     That Vectors can be more than singletons is not just a computer science
     convenience-- it's a fact of life. Consider the Pauli matrices (pauli.py).
     They can be combined into a single \f$\mathbb R^3\f$ vector,
     \f$ {\bf \vec{\sigma} } \f$-- there is but one vector, but its components
     are 2 x 2 matrices, representing internal degrees of freedom-- the
     fundamental representation of SU(2). They could be spin up and the
     orthogonal
     spin down (note: up is NOT minus down!... mind-blown), or
     isospin "up flavor" and
     "down flavor", or polarization degrees of freedom. It's up to you.
     The point is to convince the dear reader beyond all objections that
     vectors are not 3 numbers in an array.

     Regular Methods:
     ===============
	
     The 1st regular method is:

\verbatim      	     v.iter()
\endverbatim
      it returns "x", "y" and "z" in order via an iterator.
      
\verbatim		v.tolist() --> [v.x, v.y, v.z] = list(v.iter())
\endverbatim
      does it all at once.

      Other methods, some called by magic methods, are straightforward:

\verbatim      	    v.dual()
\endverbatim
      contract with Levi-Civita symbol and make a Tensor (which
      will do the cross product via posterior multiplication)

\verbatim      	   v.dot(u)
\endverbatim
      Dot Product: contract with another Vector and make a Scalar

\verbatim      	  v.cross(u)
\endverbatim
      Cross Product: contract the dual with another Vector and make
      a Vector (There is no distinction between vectors and axial
      vectors)

\verbatim			v.outer(u)
\endverbatim
       Outer Product: with another Vector, make a dyad, that is a Matrix

\verbatim       	     v.dilation(alpha)
\endverbatim
       change length with multiplication or division by a Scalar or number.
       There are also right_dilation and left_dilation method, if your
       vector's components don't commute with the scalar (which is only
       a scalar in \f$\mathbb R^3\f$, and not in the 'other' space.)

\verbatim      	     v.dyad(u)
\endverbatim
      Is outer

\verbatim      	     v.hat()
\endverbatim
      Is dilation by 1/abs(v) --> returns a Vector's unit Vector.

\verbatim      	   v.L2norm()
\endverbatim
      L2-norm is just the norm.
      
      There are also bonus functions for not only computing projections,
      rejections, and reflections, but also functions for creating
      like-wise operators:
      \verbatim
v.projection(u)   v >> u    u.projector()(v)    Proj(u)(v)	
v.rejection(u)    v << u    u.rejector()(v)     Rej(u)(v)	
v.reflection(u)   v | u     u.reflector()(v)    Ref(u)(v)	
\endverbatim	
Rotations about arbitrary vectors:
\verbatim      	      v.versor(angle)
\endverbatim	
      are given by the resulting unit quaternion.

      Special Numpy Methods:
      ======================
      You can call numpy array methods like "mean", "sum", "cumsum" on tensors
      and get a tensor result--or on a coordinate to get a coordinate result.
      That is:

\verbatim      	   >>>platform.mean()
\endverbatim
      will be the average position of a platform in LLH, ECEF, SCH, LTP, etc--
      whatever it is in. While for a 2-d array like SCH coordinate, "dem ",
      representing a digital elevation model (DEM) (with geo-located pixels):

\verbatim      		   >>>dem.mean(axis=0).h
\endverbatim
 be an array of heights that represent the height averaged over cross-
      track, at fixed along track coordinate.
      
      There is also an append method, that allows you to stack vectors onto
      your object by applying numpy.append to each component.
      
      It's really quite slick. It's all done via the
      geo.utils.arraylike.ABCNumpy.broadcast
      method, which allows you to apply functions to elements and
      reconstruct tensors/coordinates from the results (and you can use
      numpy's "axis" keyword seamlessly). You could also broadcast numeric
      derivatives and what-not if you're dealing with a motion track, for
      instance.

      More sophisticated operations are quite simple, for instance, if "v" is a
      Vector made up of numpy arrays:

\verbatim      	     sigma = ((v-v.mean())&(v-v.mean())).mean()
\endverbatim
      is their correlation tensor. (Do you know how much code that would be with
      Fortran arrays?).

      Furthermore, you can project that correlation tensor into it's
      irreducible spherical components at will-- it's just a method call that
      executes 1 or 2 lines of code.

     ------------
     |Non-Vectors|
     -------------

     Scalar:
     ======
     Scalars are NOT PRIMITIVES! A single scalar can be represented (poorly)
     by a single floating point number-- but that is no reason to do it that
     way. 

     When you square a vector, you are going to get a scalar. Why is this not
     just a number? There are several reasons:

     (1) Scalars are rank-0 tensors, and hence, have tensor properties that need
     to be captured in a class. For example: they are invariant under rotations
     and transform according to:

\verbatim                                 s' = s
\endverbatim
	While this is trivial, it is formally a trivial
	representation of SO(3), it is
     by no means "just a number". Consider the would-be scalar:
     \n\n
     \f$ s = \bf{ \vec A\cdot}(\bf{\vec B \times \vec C}) \f$,
     \n\n
     which is really a pseudo-scalar. In the full Clifford algebra off
     Euclidean 3-space, it's a trivector, or 3-blade, and that cannot
     be confused with "a number".

     (2) On a more practical front: If you are dealing with numpy arrays and 
     don't have scalars, you will 
     get burned when doing reflected multiplication. Let's say you want to scale
     a 100,000 point motion history, v, with:

\verbatim 	    	     (1/(v**2))*v
\endverbatim
     Here v.x, v.y, v.z are 100,000 element numpy arrays. If v**2 is just a
     number (a 100,000 element array), the operator will use numpy's overloaded
     "*" are try to make a:

     100,000 x 100,1000 =10,000,000,000 element array.

     Try it-- you will not like it. 
     
     To that end, I have a geo.metric.euclid.scalar.Scalar() class, and it has
     one attribute "w", so said
     operation will be executed by geo.metric.euclid.scalar.Scalar.__mul__
     which yields a single
     geo.metric.euclid.vector.Vector who's attributes are 100,000 element
     numpy arrays.

     Many a python vector module has failed w.r.t. to numpy broadcasting
     because the author(s) failed to give the rank-0 tensor its full due.

     Rank 2 Tensor:
     ==============     

     There are (at least) 3 ways to look at rank-2 tensors:
     @subsection ten1 Component-wise
     (1) A Component-wise view, which can be either:\n
          A triplet of irreducible spherical parts:\n
     	   	   1 scalar component (proportion to the identity)\n
     		   3 rank 1 components (The antisymmetric part-- which
		     	    	       rotates like a vector, but doesn't
				       change sign under reflection)
		+  5 pure rank-2 dyads.  (The symmetric, traceless part)\n
	-------------------------------------------------------------------\n
	or	   9 Cartesian components: xx, xy, .., yz, zz (We'll use this\n
		     	       		      	      	     for attributes)\n
     @subsection ten2 As Transforming Objects
     (2) Geometric animals, U,  that transform as\n\n

     	 	   	    \f$U'_{ij} = T_{ik}  T_{jl}  U_{kl} \f$\n
     @subsection ten3 (Bi)Linear Maps
     Tensors can be defined as bilinear maps from:\n
     	 (3A) \f$ {\mathbb R}^3 \times {\mathbb R}^3 \rightarrow {\mathbb R} \f$ \n
	 (3B)  \f$ {\mathbb R}^3  \rightarrow {\mathbb R}^3 \f$, \n

      geo chooses to represent tensors 9-components, representing the
      weights of 9 Cartesian dyads--this is standard. Moreover, they
      are used to represent the so-called
     direction-cosine-matrix (DCM), even though a DCM is technically not
     a tensor-- the python Zen (import this) "practicality beats purity"
     is respected--the user should look upon their rotation matrices as just
     that, eventhough they live in a tensor class.

      Rows and Columns:
      -----------------
      Rows and Columns have no physical meaning and are hence, meaningless.
      They are at best a linear algebra concept: matrices have rows and
      columns, but tensors are not matrices. With respect to tensors,
      they are a typographical construct-- a relic of writing frame dependent
      representation of geometric objects on paper. I do not support them, and
      frankly, I don't know what "row major" or "column major" means. My tensors
      have indices, and each one is as "major" as the other.
      
      Nevertheless, there are methods:

\verbatim
tensor.cols() -->  T.ix, T.iy, T.iz
tensor.rows() -->  T.xi, T.yi, T.zi
\endverbatim
 iterate over 1st and 2nd tensor indices -- and
      indices (slots, really) have meaning.
      I just happen to call it rows (cols)--but they're NOT rows (cols).
      It's handy though, if
      you want to know one of the important tensor invariants:
      \n\n
      \f$ m_{1i}  m_{2j}  m_{3k}  \epsilon_{ijk} \f$
      \n\n
      which is the determinant, you can compute it by unpacking the rows in the
      geo.metric.euclid.vector.scalar_triple_product function:
      
\verbatim >>>scalar_triple_product(*tensor.rows())
\endverbatim
      and that does make a (pseudo) scalar our of it. The other scalars are the
      trace and the L2-Norm, and there are methods for those too. To summarize,
      for a matrix/tensor T (with functions from the tensor.py modules):
      \n\n      
      Rank-2 Tensor Invariants:
      ----------------------
      Some invariants:
\verbatim  abs(T) = T.L2norm() = T.norm([l=2])
\endverbatim
      \f$ \sqrt{\sum_{ij}T_{ij}^2}\f$
\verbatim Tr(T) = T.trace() = T.contract(0, 1) = T.ii
\endverbatim
      \f$ \sum_{ij}T_{ij}\delta_{ij}\equiv T_{ii}\f$.
\verbatim det(T) = T.det()
\endverbatim
      \n\n
      \f$ m_{1i}  m_{2j}  m_{3k}  \epsilon_{ijk} \f$
      \n\n
      The 3 invariants can also be expressed as coefficients of
      the charactertics polynomial:
      \verbatim
T.poly1d()              # a numpy.poly1d object
T.characteristics()     # a tuple of Scalars.
T.J(n); n = (1, 2, 3)   # The traditional n-th order invariants.
\endverbatim
\n\n
      Rank-2 Tensor's Vectors:
      ---------------------- 
      Vectors from Tensors:

\verbatim      	      T.dual()
\endverbatim
      \f$ \frac{1}{2}\epsilon_{ijk}T_{jk} \f$
 
\verbatim	      T.vector()
\endverbatim
      \f$ \epsilon_{ijk}T_{jk} \f$ 
\n\n
      Rank-2 Tensor's Rank-2 Tensor 
      ----------------------------
      Tensors from Tensors:

\verbatim
T.C                                    # Complex Conjugate
T.T = T.transpose() = T.ji             # Transpose
T.H = T.T.C                            # Hermitian conjugate
~T    --> T.__invert__() = T.I = T**-1 # Inverse
adj(T) = T.adj() = ~T/det(T)           # Adjugate tensor
cof(T) = T.cof() = adj(T).T            # Cofactor tensor
symm(T)= T.S()                         # Symmetric Part
       = (T + T.T)/2                   #   via direct computation
       = T.symmetrize([[0, 1]])        #   via Schur-Weyl Duality
       = getattr(T, "{ij}")            #   via Einstein summation notation
skew(T)=T.A()                          # Antisymmetric Part
       = (T - T.T)/2                   #   cia direct computation
       = T.symmetrize([[0], [1]])      #   via Schur-Weyl Duality       
       = getattr(T, "[ij]")            #   via Einstein summation notation
hyd(T) = T.hyd() = DELTA * Tr(T)       # Hydrostatic part (c.f. Stress tensor)
dev(T) = T.dev() = T - hyd(T)          # Deviatoric part (c.f. Stress tensor)
T.natural_form() = symm(T) - hyd(T)    # Symmetric, trace-free part (pure J=2).
U = T.diagonalizer()                   # Unitary transform that diagonalizes T
T.diagonalized() = U.I* T * U          # Diagonalized form.
T.irrep(j, m)                          # Cartesian part that is pure J=j, M=m
\endverbatim
      \n\n
      Other methods include:
      \n\n
\verbatim
T.spherical()            # Irreducible Spherical Representation
T.eig()                  # Eigenvalues and Eigenvectors
T.angle_of_rotation()    # Angle of rotations (for a DCM)
T.axis_of_rotation()     # Axis.....
T.sqrt()                 # Square Root
T.exp() = exp(T)         # 1 + T + T**2/2 + T**3/6 ...
T.log() = log(T)         # Inverse of exp(T)
T.svd()                  # Single Value Decomp into dyads.
\endverbatim
\n\n
      Some other operations are:
\verbatim
T.commutator(M) = T*M - M*T     # [T, M]
T.anticommutator(M) = T*M + M*T # {T, M}
T.sandwich(M) = T * M * T.T     # aka conjugate product
T.row(n) = T.e(n)                      # n-th "row"
T.col(n) = T.e(Ellipsis, n)            # n-th "col"
for p in T.polyadics()...              # Projections onto Cartesian basis
for ee in T.basis()...                 # Cartesian Basis: (X&X&, X&Y, ... Z&Z)
T.inner(U)                     # Inner product (any rank)
T.outer(U)                     # Outer product (any rank)
T.wedge(U)                     # Generalized Cross Product (any rank)
\endverbatim
The last 7 are accomplished for any rank via:
\verbatim
>>list(T.iter_indices())
>>[(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
\endverbatim
That allows index gymnastics without resorting to a _for_-loop.
     @subsection tdps Tensor (dyadic) Products
     The double/mixed dot and cross products are:\n
     geo.metric.euclid.tensor.Tensor.dot()   \n
     geo.metric.euclid.tensor.Tensor.double_dot()  (same as dot) \n
     geo.metric.euclid.tensor.Tensor.dot_cross()   \n
     geo.metric.euclid.tensor.Tensor.cross_dot()   \n
     geo.metric.euclid.tensor.Tensor.double_cross() \n
     \n
     @subsection tps Tensor Projection, Rejection and Reflection
     Dyadic projection, rejection and reflection generalize their
     vector counterparts via the double dot product:
\verbatim
T" = T >> T'  # project T onto T' dyads.
T" = T << T'  # reject T onto T' dyads.
T" = T | T'  # reflect T onto T' dyads.     
\endverbatim
      Simple enough. Now if T' is replace by a vector "v", then it gets
      tricky, because you have to choose and index on which to perform
      the operation, and hence, the operator overloads return functions
      that take the index as an argument:
\verbatim
T' = (T >> X)(0)  # project out X "row"
T' = (T >> X)(1)  # project out X "column"
T' = ((T >> X)(0) >> Y)(1) # only T.xy survives.
\endverbatim
      I have used basis vector to make the geometric concept understandable-
      it's not a common operation-- but you can replace "X" with _any_
      vector. (An application, for instance, would be with a stress tensor
      and some arbitrary vector to get the momentum flux in that direction).
      The rejections work like-wise on the perpendicular components.
      Reflection is straight forward (here with unit vectors):
\verbatim
T' = (T | X)(0)  # invert sign of X "row"
T' = (T | Y)(1)  # invert sign of Y "column"
\endverbatim
      These operations work on rank 3 and higher rank tensors, too. (See
      http://en.wikipedia.org/wiki/Cauchy_stress_tensor for apps).

The Dyadic Basis
================
      Finally, tensors are often represented as vectors (1-index obejcts)
      on a 9-dimensional vector space, with bases: xx, xy, ..., zy, zz.
      You can convert your tensors to this space (rather procedurally as
      numpy-arrays) via:
\verbatim

v = T.nineXone()      # convert to a numpy column.

\endverbatim
	A rotation matrix can then be converted, via a Kronecker product:
\verbatim

R = M.nineXnine()      # convert to a 9 x 9 numpy matrix

\endverbatim
      The tensor transformation to a different frame then becomes a
      a linear transformation:
\verbatim      
v' = R * v.
\endverbatim

Note: R can be block diagonalized in to 5, 3, and 1 sized blocks. These blocks
represent the rank 2, 1, and 0 parts of a general rank-2 tensor.

      @subsection rank34 Higher Rank Tensors
      Cartesian rank-3 and rank-4 tensors are fully supported. Higher rank
      tensors are created on-the-fly when requested (e.g. by the outer
      products of existing tensors).

      The Levi-Civita (Pseudo)Tensor is a package constant:
\verbatim

>>>print geo.EPSILON
    [0, 0, 0] [0, 0, -1] [0, 1, 0]
    [0, 0, 1] [0, 0, 0] [-1, 0, 0]
    [0, -1, 0] [1, 0, 0] [0, 0, 0]

\endverbatim
     A expected, the attributes are triplets fromed from 'x', 'y', and 'z'.
\verbatim

>>>print geo.EPSILON.xyz
1

\endverbatim
      Alegbraic operations "+", "-", "*", "&" are supported. Inner products
      and contractions can be done with the inner() and contract() methods,
      though for clarity, Einstein summation is preffered. For example:
      \n\n \f$ \epsilon_{ijk}\epsilon_{ijl} = 2 \delta_{kl} \f$ \n\n
      can be verfied via:
\verbatim      

>>>print (EPSILON & EPSILON).ijkijl
[2.0, 0.0, 0.0]
[0.0, 2.0, 0.0]
[0.0, 0.0, 2.0]

\endverbatim
      Higher rank tensors are displayed to the screen recursively, for
      example:
\verbatim

>>>print Four(*range(3**4))
======================================
[0, 1, 2] [9, 10, 11] [18, 19, 20]
[3, 4, 5] [12, 13, 14] [21, 22, 23]
[6, 7, 8] [15, 16, 17] [24, 25, 26]
--------------------------------------
[27, 28, 29] [36, 37, 38] [45, 46, 47]
[30, 31, 32] [39, 40, 41] [48, 49, 50]
[33, 34, 35] [42, 43, 44] [51, 52, 53]
--------------------------------------
[54, 55, 56] [63, 64, 65] [72, 73, 74]
[57, 58, 59] [66, 67, 68] [75, 76, 77]
[60, 61, 62] [69, 70, 71] [78, 79, 80]
======================================
\endverbatim
	Einstein summation notation is supported up to rank-15 (which requires
	over 16,000,000 lines to print). Let's hope your tensor problems aren't
	that deep.

      @section charts Charts on SO(3)
      The euclid.py module is all about Tensors-- and how you add, subtract,
      and multiply them. The rank-2 tensor also transforms vectors--but it is
      not alone. There are many ways to represent rotations in
      \f$\mathbb R^3\f$, and
      collectively, they are known as charts on __SO(3)__-- the rotation group.
      They live in charts.py. A nice introduction is provided here:

      http://en.wikipedia.org/wiki/Charts_on_SO(3)

      Versors
      =======
      Some people like rotation matrices. I don't. Not just because Euler
      angles always seem ambiguous-- it's because __SO(3)__ -- the group of
      rotation
      matrices--is not simply connected and you get problems interpolating
      rotations. Moreover, some representation are degenerate, which leads
      to the well known problem of gimbal lock.
      Fortunately, in any
      dimension "N", the rotation group, __SO(_N_)__ has a simply connected
      universal convering group, and the group is the so-called spin-group: 
      __Spin(N)__.
      But what is __Spin(3)__ and how do you represent it? I have no idea.
      Fortunately, it's isomorphic to the much more popular
      __SU(2)__....  but that
      uses complex 2x2 matrices--which you have to exponentiate, and has spinors
      that need to be rotated 720 degrees return to their original state
      (what's up with that?)- it's all simply too much. The good news is
      __H__:
      the quaternion group from the old days, aka the hyper-complex numbers,
      will do just as well:

      \f$ 1, i, j ,k  \f$  (called "w", "x", "y", "z" axis) \n\n
      with: \f$i^2 = j^2 = k^2 = ijk = -1 \f$ \n\n
      As a group, they're a simply connected representation of rotations, are
      numerically stable, can be spherical-linearly interpolated (slerp), are
      easy to represent on a computer, and in my opinion, easier to use than
      matrices. (They are also the standard for on-board real-time spacecraft
      control). Their only draw back is that they transform vectors
      quadratically (via a sandwich product)-- that fact right there should
      tell you vectors are not fundamental--that they're made out of
      "something"....but I digress...).

      Hold-up. We only care about unit quaternions, that is, quaternions with
      norm == 1. They are a subset of the quaternions call versors
      (euler/hamilton.py)
      that is isomorphic to the hyper-sphere __S4__ and closed
      under the Grassmann product.
      
      Still, there is always an choice in how to represent them. If you're doing
      math, and don't care about rotations, the Cayley-Dickson extension of the
      complex numbers is best, two complex numbers and another imaginary unit
      "j":\n
      \f$ z + jz' \f$\n
      Hence, with:\n
      \f$ q  \equiv  z + z'i \f$\n
      \f$  q = (a + ib) +  (a' + ib')j  \f$\n
      \f$  q = (a + ib +  ja' +  kb') \f$\n\n
      you get a quartet of reals:\n
\verbatim      	  q = (a, ib, ja', kb').
\endverbatim
      You would think that would end it, but it isn't.

      There is always a question, if you have a primitive obsession and are
      using arrays to represent quaternions, no one is certain--in spite of the
      unambiguous Cayley-Dickson construction--if you are doing:

\verbatim      		  	(w, x, y, z)
\endverbatim
      or
\verbatim      		  	(x, y, z, w)
\endverbatim
      Really. No, REALLY. People put "w" last, even though it destroys the
      alphabet, because it preserves the indexing
      into the vector, while adding on a scalar--Well, I DON'T INDEX VECTORS--
      they don't have items (unless they do-- A Vector of ducks quacks like 3
      ducks).
      
      I break them down into a scalar part and a vector part. The order doesn't
      really matter. Hence a Versor, q,  is:

      a Scalar (w)
      a Vector (x,y,z).

      q.w is q.scalar.w
      q.x is q.vector.x and so on for y and z, and THERE IS NO q[i]. 

      Hence:\n
\verbatim			>>>q = Versor(scalar, vector)
\endverbatim
      Your quaternions can be singletons or numpy arrays, just like Scalar and
      Vectors.

      They transform with their "__call__" method (NOTE:
      __everything__ transforms
      with its "__call__" method), which overloads the function syntax:

\verbatim      	       >>>q(v)
\endverbatim
      and you compose rotations with multiplication:

\verbatim      	      q(q'(v)) = (q*q')(v) or (q'*q)(v)
\endverbatim
      Which is it, do you left or right multiply? It depends. The Versor
      rotations can be converted into Matrix objects via:
\verbatim      q.AliasMatrix()    --> Alias transformation matrix\endverbatim
\verbatim      q.AlibiMatrix()    --> Alibi transformation matrix
\endverbatim
      Alias or Alibi?
      ===============
      Another point of confusion is "what is a transformation?". Well, alias
      transformations leave the vector unchanged and give its representation in
      a different coordinate system. (Think ALIAS: same animal, different name).
      \n
      Meanwhile the Alibi transformation leaves the coordinates fixed and
      transforms the vector. (Think ALIBI: I wasn't there, because I was here)
      \n
      What about pre or post multiplication? That is:
      \n\n
      \f$ v'_i = M_{ij}v_j \f$ or \f$ v'_j = v_i M_{ij} \f$?
      \n\n
      I have chosen the former for alibi (active) transformations. There
      are 2 reasons for this choice:
      \n\n
      1) Euler-angle and Tait-Bryan representations of rotations are
      traditionally active; moreover, the more familiar pre-multiplication
      leads to break down of Tait-Bryan aircraft coordinates: Yaw, Pitch,
      and Roll to be:
      \n\n
      \f$ {\bf M}_{ypr} = {\bf M}_{\rm yaw}{\bf M}_{\rm pitch}{\bf M}_{\rm roll} \f$.
      \n\n
      That leads to a great deal of simplification when converting between
      representations.
      \n\n
      2) When dealing with quantum operators, state-changing
      pre-multiplication is the most common.
      \n\n
      3) For non-unitary / non-orthogonal transforms, it makes sense.
      \n\n
      4) When transforming objects of equivalent rank (Rotation matrices
      operating on rank-2 tensors, or quaternions transforming vectors),
      one recovers the standard conjugate (or sandwich) product:
      \n\n
      \f$ a' = {\bf \hat O} a {\bf \hat O}^T \f$.
      \n\n
      For tensors, this is a direct result of tensors being composed
      of the outer products of pairs of vectors. Likewise for spinors
      composing vectors (and scalars, for that matter).

      (That is, a quaternion transforms a vector via:
      \n\n
      \f$ (0, \vec v') = {\bf q} (0, \vec v) {\bf \bar q} \f$.
      \n\n
      and likewise for scalars
      \n\n
      \f$ (s', \vec 0) = {\bf q} (s, \vec 0) {\bf \bar q} \f$.
      \n\n
      Thus, the transformation of vectors and scalars is unified via the
      quaternion. This is because quaternions act fundamentally on spinors:
      \n\n
      \f$ \xi'={\bf q}\xi \f$
      \n\n
      and both vector and scalar can be written as the outer product of
      2 spinors:
      \n\n
      \f$ v = \xi\xi \f$
      \n\n.
      In terms of representation theory, this can be expressed as the
      tensor product of the 2 state spinor space with itself:
      \n\n
      \f$ {\bf 2} \otimes {\bf 2} = {\bf 3} \oplus {\bf 1}\f$
      \n\n
      being the tensor sum of a vector-space and a scalar space.

      But: this is a bit much for the dear user to worry about, hence,
      ALL operations are overloaded via function emulation ("__call__");
      thus, all transformation representations act like functions that
      can transform whatever kind of object you give them.
      \n
      For example, for quaternions, one might code"

\verbatim
v' = M * v
v' = (q * v * (~q)).vector  # quat. mul does v --> (0, v).
\endverbatim
      and get the right answer, but would you rather use the same form
      for each transformation:
\verbatim
v' = M(v)
v' = q(v)
\endverbatim
      Likewise for vectors and tensors (and anything else). You transform
      your angular momentum vector, l, and your inertia tensor, I,  via:
\verbatim
l' = M * l
I' = M * I * M.T
I' = M * I * ~M
\endverbatim
      but isn't it sweeter to just use:
\verbatim
l' = M(l)
I' = M(I)
\endverbatim
       \n
       and let the object figure out the right answer?
       \n
      Euler Angle Classes
      ==================
      There is a base class, EulerAngleBase, for transformations represented as
      Euler angles. Like matrices and versors, you can compose and call them.
      Nevertheless, there are 2 points of confusion:
      \n\n
      (1) What are the axes
      \n
      (2) What are the units.
      \n\n
      The answer:
      \n\n
      I don't know. No, really, the EulerAngle class doesn't know. It uses its
      <a href="https://en.wikipedia.org/wiki/Class_variable">
      static class attributes</a>
      (<a href="https://en.wikipedia.org/wiki/Method_(computer_programming)#Abstract_methods">
      abstract</a> in the base class):
      \n\n
      geo.metric.euler.charts._ElementalRotationTriplet.AXES
      ( a length 3 tuple of vectors representing ordered intrinsic
      	      rotations), and
\n	      
      Circumference    ( a float that should be 2*pi or 360)
      \n
      to figure it out. So, to support common radar problems like platform
      motion, there is a subclass:
      \n\n
      geo.metric.euler.charts.YPR
      \n\n
      which has AXES set to (z, y, x)-- so that you get a yaw rotation followed
      by a pitch rotation followed by a roll rotation; Circumference=360, so
      that if you have an ASCII motion file from a gyroscope (in degrees),
      for instance, you can do:
      \verbatim
      >>>attitude = YPR(*numpy.loadtxt("gyro.txt")
      \endverbatim
      Yeah-- it's like that. Plus, there's a RPY class, because some
      people do the rotations backwards. (Imagine the FORTRAN code to do this:
      \n(1) get a logical unit
      \n(2) open the file
      \n(3) compute a mess of a 3 x 3 matrix to do the rotation
      \n(4) create an empty array in which to put the result
      \n(4) loop over the file
      \n(5) read a line of the file
      \n(6) do the matrix multiplication (in a 2-deep do loop)
      \n(7) assign the results correct elements of the output array
         \n(7.1) --(7.3) for x, for y, and the for z
      \n(8) end do
      \n(9) close the file
      \n(10) free the logical unit
      \n(11) compile it
      \n(12) fix the bugs \n(goto 11, unit it compiles)
      \n(13) run it
      \n(14) read stack trace, \n(goto 1 unit it runs)
      \n(15) write code to put the results into a file
      \n(16) write code to read the file.
      \n(17) goto 1 until there are no bugs
      \n(18) hope you didn't miss anything
      \n(19) Take a lunch break-- no, call it a day.
      )
      \n\n
      Rotation Summary:
      ================
      In the spirit of OO: you do need to know which object you have when
      performing transformation. If you have "T", then:

\verbatim      		 T(v)
\endverbatim
     Transforms v
\verbatim		  ~T
\endverbatim
     Inverts the transformation
\verbatim     	     T*T'
\endverbatim
 composes 2 transformations, so that

\verbatim     	      T*(~T)
\endverbatim
is the identity transformation (for any rotation object). Specific forms
of the rotation are as follows:

\verbatim     	  	T.dcm()
\endverbatim    
     returns the equivalent direction cosine matrix (as a Tensor)

\verbatim		 T.versor()
\endverbatim
     returns the equivalent versor

\verbatim		 T.ypr()
\endverbatim
     returns the equivalent YPR triplet 
	
\verbatim			T.rpy()
\endverbatim
     returns the equivalent RPY triplet .
     \n\n
     T can be a Matrix Versor YPR or RYP instance. And that is POLYMORPHISM in
     a nutshell.

      @subsection Q Quaternions
      The full quaternion algebra, __H__, is represented in desic/hamilton.py.
      While it has no use for transformations in __E3__, it is included
      for completeness. While versors only support the Grassmann product,
      quaternions allow for addition, subtraction, and dilation operations
      as well as inner, outer, odd, and even products.

      @section Affine
      The Affine transformation is the most general linear function of
      a vector:
      \n\n
      \f$ A(\vec x) = {\bf M}\vec x + \vec b\f$.
      \n\n
      For proper euclidean transformation (det(M) = 1), than any rotating
      object will do. Dilations, skews, reflections, etc.. will require a
      fully general matrix (rep'd as a rank-2 tensor).
      \n
      The magic methods support:

\verbatim
A(v)                    # function emulation
(A * A')(v) = A(A'(v))  # composition
~A(A(v)) = v            # inversion
A ** (m/n)              # interpolation
\endverbatim
      MODULE FUNCTIONS 
      =================
      First of, I believe all functions should be PURE functions: that is,
      they do not change state. Methods can change the state of an objects,
      but functions should not. This is an important practice for proper OO
      design.
      \n\n
      You will notice that the euclid module is full of functions for doing the
      vector operations. They are there so that methods can call them (i.e.,
      private, but that's not enforced). You can use them if you want.
      (An example
      is doing X^Y. The super says that is X.wedge(Y), which is laborious
      summation
      over the Levi-Civita tensor, and that is not acceptable for real
      computations.
      \n\n
      Hence, it is overridden to call vector._cross_product(X, Y)-- which does
      the usual antisymmetric combination over the non-zero elements of the
      generalized sum).
      \n\n
      Some
      functions are there for public use:
      \n\n
      geo.metric.euler.charts.roll  (aliased to geo.roll)\n
      geo.metric.euler.charts.pitch (aliased to geo.pitch)\n
      geo.metric.euler.charts.yaw   (aliased to geo.yaw)\n 
      \n\n
      These function are constructor-helpers for defining rotation
      objects.
      They are bumped up to the top namespace, which is supposed to
      just contain things you need access too- but that may need a
      little code review. The return Versors (unit Quaternions).
      
      (Note: the geo.metric.euclid.euclid.LinearMap.aschart() method allows
      the user to convert any SO(3) to any form via an name, instance, or
      Type).
      \n\n
      MODULE CONSTANTS
      =================
      The module constant "BASIS" is a collections.namedtuple with  "x", "y",
      and "z" elements that provide basis vectors. A similar namedtuple
      "DYADS" has the 9 tensor basis dyads ("xx", "xy", .., "zz").
      \n\n
      geo.metric.euclid.vector.X = \f$ \bf \hat x \f$ \n
      geo.metric.euclid.vector.Y = \f$ \bf \hat y \f$ \n
      geo.metric.euclid.vector.Z = \f$ \bf \hat z \f$ \n
      geo.metric.euclid.vector.BASIS = \f$ ({\bf \hat x,\hat y,\hat z})\f$\n
      geo.metric.euclid.vector.NULL \f$ \bf \vec 0\f$ \n \n
      geo.metric.euclid.scalar.ZERO = 0 \n
      geo.metric.euclid.scalar.ONE = 1 \n
      are basic Scalars. \n\n
      geo.metric.euclid.tensor.DELTA = \f$ \delta_{ij} =
      {\bf \hat x\hat x + \hat y\hat y + \hat z\hat z}\f$ \n is the 
      Idempotent Tensor (Kronecker delta).\n
      geo.metric.euclid.tensor.PARITY = \f${\bf P} =  -\delta_{ij}\f$ \n
      is the coordinate inversion operator.\n
      geo.metric.euclid.tensor.BASIS = \f$ ({\bf \hat{x}\hat{x},\hat{x}\hat{y},\ldots, \hat{z}\hat{z} })\f$ \n\n
      geo.metric.euclid.three.LEVI_CIVITA = \f$ \epsilon_{ijk} \f$ \n is the
      fully antisymmetric rank-3 tensor. \n\n
      The Versors are:\n
      geo.metric.euler.hamilton.W \f$ {\bf w} \f$ \n
      geo.metric.euler.hamilton.I \f$ {\bf i} \f$ \n
      geo.metric.euler.hamilton.J \f$ {\bf j} \f$ \n
      geo.metric.euler.hamilton.K \f$ {\bf k} \f$ \n
      geo.metric.euler.hamilton.BASIS \f$ ({\bf w, i, j, k}) \f$\n
      and they are trapped on the unit hyper-sphere and only
      support multiplication and slerp, while:\n\n
      the Quaternions are:\n
      geo.desic.hamilton.W \f$ {\bf w} \f$ \n
      geo.desic.hamilton.I \f$ {\bf i} \f$ \n
      geo.desic.hamilton.J \f$ {\bf j} \f$ \n
      geo.desic.hamilton.K \f$ {\bf k} \f$ \n
      geo.desic.hamilton.BASIS \f$ ({\bf w, i, j, k}) \f$\n
      are free to roam __H__ (actually __H__ x __C__), as they are
      complexified (Biquaternions).\n\n
      Ellipsoids:\n
      geo.detic.newton.wgs84.wgs84.WGS84 is __The__ ellipsoid.\n
      geo.detic.newton.almanac.ALMANAC are __all__ the ellipsoids.\n\n	

      @section Coordinates
      Coordinates are a lot like vectors, since they can be represented by 3
      frame-dependent components. Both Vector and the coordinate base classes
      inherit from the ABCNumpy interface. Hence, all coordinate
      instances can be singletons or numpy arrays, and they all have
      broadcasting, mean, "__getitem__", "__getslice__",
      "__iter__", tolist(), and
      iter() and next() behavior at hand.
      \n\n
      Moreover, they're all instantiated with 3 arguments and a possible peg
      point keyword or positional argument, but I'm getting ahead of myself.
      \n
      At this time, there are 4 (four) coordinate systems:
      \n\n
         2 have origins at the center    \
	\n\n   
      	 2 are tangent to the ellipsoid   \   and the inheritance diagram
	 \n\n
	 2 are Cartesian                 /   reflects these facts
	 \n\n
	 2 are not                       / 
	 \n\n
      They are:
      \verbatim
      ECEF             Earth Centered Earth Fixed Cartesian (x, y, z)

      LLH              Geodetic   (lat, lon, hgt)
      
      LTP              Tangent Plane Cartesian (x, y, z, peg_point=peg_point)
      
      SCH              offset Tangent Sphere   (s, c, h, peg_point=peg_point)
      \endverbatim
      A PegPoint is just a namedtuple consisting of (lat, lon, hdg). Of course,
      all arguments are in degrees-- radians are simply not for geolocation.
      Presumably, there are a bunch of ships off the coast of west Africa, full
      of people using radians for navigation.
      \n\n
      Transformations:
      ----------------
      The idea here is to be as polymorphic as possible, so for instance, if
      you have a coordinate instance, "r":
\verbatim
r.ecef()	
r.llh()	
r.ltp()  or  r.ltp(peg)
r.sch()  or  r.sch(peg)
\endverbatim
      will convert r to ECEF, LLH, LTP, SCH (or LTP, SCH with a different peg).
      Moreover, non-Cartesian coordinates can be readily converted to Cartesian
      coordinates, via their "cartesian_counter_part" method, so that:
\verbatim
llh.cartesian_counter_part() --> llh.ecef()
sch.cartesian_counter_part() --> llh.ltp()
\endverbatim
      Moreover, all classes have a "vector()" method that converts the points
      into a euclid.Vector instance relative to the origin of the coordinate
      system. Thus:
\verbatim      	      r.vector()
\endverbatim
      always works, yielding a Vector relative to the center of a coordinate
      system or relative to a peg point.
      \n\n   
      The point is that all instances have the same methods, so you don't need
      to know what you're dealing with to get the result you want. If I have a
      motion history, r,  in LTP coordinates and a target point on the ground,
      p, I get the vector in the peg pointing from the platform to the target
      via:

\verbatim		>>>d = p.vector() - r.vector()
\endverbatim
      Now if the whole thing is in SCH coordinates, the computation is:
\verbatim      	        >>>d = p.vector() - r.vector()
\endverbatim
      (because non-Cartesian coordinates call their Cartesian_counter_part()
      transformation via vector()).
      \n\n
      Now if you have the motion in SCH and the target in LTP, you DO NOT have
      to type-check, then the formula is 
      
\verbatim	        >>>d = p.vector() - r.vector()
\endverbatim
      Clearly, a pattern is forming-- but all good thing must come to an end.
      What if you motion is in LLH (or ECEF, or a different SCH or LTP system)?
      Then:
      
\verbatim	       >>>d = p.vector() - r.ltp(p.peg_point).vector()
\endverbatim
      The strategy to overload the "__sub__" operator as follows:

\verbatim      	       >>>d = p - r
\endverbatim
      will make a euclid.Vector instance relative to p's origin and Cartesian
      orientation, even p is not Cartesian. Meanwhile, r can be in any
      coordinate system.
      \n\n
      Further Note: numpy will handle all the broadcasting, you can have 1
      coordinate be a bunch of arrays, and the other be a like shaped array,
      or a singleton. It's going to work.
      \\n\n
      Of course, "__add__", inverts the whole thing: you can add a vector to a
      point:

\verbatim			>>>p' = p + v
\endverbatim
      and p' will be in p's coordinates, with v interpreted relative to p's
      Cartesian coordinates.
      
      @section TEP The Ellipsoid Problem
      The ellipsoid is defined by a semi-major axis, an inverse
      flattening, and and optional model name (ellipsoid.py)
      \n\n
      Its method/properties convert to many of the various metrics
      associated with an ellipsoid of revolution. It also converts
      between all the various latitude definitions.
      \b\b
      COORDINATES: Here is the key fact: coordinates only work on ellipsoids.
      So, the ellipsoid has methods that construct coordinate objects.
      \n\n
      Suppose you have the WGS-84 ellipsoid:
      \n\n
\verbatim      	      >>>wgs84 = ellipsoid.Ellipsoid(6378137.0,
		       	 		     298.25722293286969,
	       	 			          model="WGS84")
						  \endverbatim
      has four methods:
      \\n\n
\verbatim      	  wgs84.ECEF\endverbatim
\verbatim      	  wgs84.LLH	\endverbatim
\verbatim      	  wgs84.LTP	\endverbatim
\verbatim      	  wgs84.SCH	\endverbatim
      Transformation are then done using the coordinate methods, so to get
      JPL's coordinates in ECEF:
\verbatim      	    jpl = wgs84.LLH(34.197, -118.175, 345.).ecef()
\endverbatim
      So, this makes an LLH instance and then converts it to an ECEF instance.
      Moreover, the ellipsoid has functions for conversion:

\verbatim      		x, y, z = wgs84.llh2ecef(34.197, -118.175, 345.)
\endverbatim
      That is, a triplet goes in, a triplet comes out.
      \n\n
      WGS84:
      -----
      The wgs84.py module has the classes and functions for WGS84 as module
      constants and functions, thereby allowing an entire procedure usage.

      @section pm Platform Motion 
      After all that, you still don't have platform motion. Enter the motion.py
      It requires numpy, since you will have array_like attributes.
      The SpaceCurve class is basically Vector which takes that into
      consideration. 
      (Recall the fundamental theorem of Space Curves? That curvature
      and torsion define a unique space curve? Yes, that one-- well space
      curves define all that: velocity, normal, acceleration, angular velocity,
      yadda yadda. They key property is you can define a local tangent frame,
      with:
      \verbatim
      x   parallel to the curve's velocity
      y   = z ^ x
      z   is some "z" orthogonal to x. The default "z" is DOWN, but you can
      	     	      		        make it UP, or something else.
					\endverbatim
      Hence, given a 3-DoF motion history, you get the transformation from
      level Cartesian space to the tangent frame. Now if you throw in attitude,
      represented by any kind of rotation, boom, you have an affine
      transformation to body coordinates.
      \n\n
      But wait: these things all have array_like attributes, that means in
      one object, you have the transformation from a local pegged coordinate
      system the body frame AT EVERY POINT IN THE MOTION HISTORY.
      	     \n\n
      Now stop and think about that for a minute. IF you were still suffering
      primitive obsession, using arrays for vectors, and using stand-alone
      functions for coordinate transformations--you be in a real pickle.
      All these arrays, which are JUST NUMBERS and have NO intrinsic meaning-
      no you the developer has to keep it straight. Then, you have to pass
      them to functions, and then to other functions-- how you deal with the
      fact that the functions are not constant--- I don't know- but you do.
\n\n
      None of that. You got a GPS history and an attitude history:

      	 \verbatim  f = SpaceCurve(*gps).tld2body(imu)
	 \endverbatim

\verbatim      f(look\endverbatim), \verbatim f(target)\endverbatim,
etc...
      does the whole affine transformation at every point.
      
      @subsection Connections
      LLH and SCH are local coordinates: their (co)tangent basis vectors depend
      on position, and we may need to know them in terms of their canonical
      global bases: ECEF and LTP, respectively. Hence: christoffel.py.
      It adds mix-ins that computer the jacobians, covariant and
      	contravariant bases
      , and connection coefficients between the cannonical frames.

      
      @section Utilities

      @subsection trig Trigonometry
      The proliferation of "degrees" necessitates basic trig functions
      that take and return degree arguments. They are implemented in
      trig.py. Moreover, it is in this module where you try to use numpy,
      but get python's built in math library should numpy be unavailable.

      @subsection Err Exceptions
      Various errors specific to the package are defined in exceptions.py.
      They're basically wonky ways to catch nonsense operations.
      

      @section exp Experimental Sub Packages: __Geodesic__
 
      @subsection cliff Clifford Algebra
      One will note the Versor and Vector seems to be related, but their
      is something missing-- they don't line up in the sense of rank,
      or which product is antisymmetric, etc.. Moreover, there is no
      distinction between the pseudo-scalar from the triple scalar product
      and the true scalar from the dot product, and likewise with axial
      vectors produced in the cross product. Finally, the versors hint at
      the existence of spinors, but why- there is NO motivation. All
      these problems are resolved in Clifford Algebras, specifically
      __Cl__(3, R),
      the Clifford algebra of
      \f$\mathbb R^3\f$. It is implemented in blades.py
      
      @subsection quat Quaternions
      Versors are limited, in that they must have unit norm, thus only the
      Grassmann product is allowed. The space is not closed under the inner,
      outer, odd, and even products as well as addition and subtraction.
      Hence, for fun, they are implemented in
      geo/desic/hamilton.py.

      @subsection cd Cayley-Dickson Construction
      Once you have full blown quaternions, there really is no reason to
      not have octonions and sedinions, tesserines, bi-quaternions,
      co-quaternions, and their "split" varieties (including the split
      complex numbers) They can all be implemented via the Cayley Dickson
      construction in extension.py

      @subsection mink Minkowski Space
      The minkowski.py module implements 4-vectors as a composition
      of a Scalar and a Vector- so that the metric is implicitly coded.
      Explicit coding, then Zen of python after all, would require implementing
      a metric tensor, thereby dumping a lot of multiplications by and
      additions of 0 into the package, and that is not warranted.

      Additional
      modules are lorentz.py (Lorentz transformation and general 4-tensors,
      as composed of rank 0, 1, and 2 3-tensors). de_broglie.py implements
      Lorentz transformation of wave-vector-- thus allowing simple
      computations of Doppler shifts and relativistic aberrations for
      your interplanetary radars (remember: your radar transmits and
      receives null 4-vectors). The pauli.py module is TBD to implement a
      spinor representation of spacetime.

      @subsection pluck Projective Geometry
      In implementation of homogeneous coordinates is TBD
      (geo.desic.plucker).

      @subsection hot  Higher Order Tensors
      euclid/three.py and euclid/four.py implement rank-3 and rank-4 tensors
      explicitly,
      while 5th and above rank tensors are created dynamically as needed
      (inheriting their behavior from geo.metric.euclid.three.HigherOrder).
      Be advised that higher-rank objects are not fully mature.

      All Cartesian tensor classes live in the
      geo.metric.euclid.euclid.Tensor_.ZOO, so that tensors of
      any-rank know about other ranks, which may be created via
      inner/outer products or index contraction.

      There is a levi_civita.py to facilitate the wedge product 
      in arbitrary dimensions (see
      geo.metric.euclid.three.LEVI_CIVITA).

      The symmetries of high rank tensors is complicated. While the rank-2
      tensor has only 2 obvious symmetries: \f$ {\bf T} \pm {\bf T}^T \f$,
      which can be computed manually.
      For rank-N, the number of symmetries is equal to the number of
      integer partitions of N, and it gets completely intractable by hand.
      Tensor symmetries are computed using the Schur-Weyl formalism:
      N is partitioned into integers and written as a Young Diagram. The
      diagram is filled to form all possible standard Young Tableaux.
      The permutations on N-letters that leave the Tableau row or
      column equivalent are combined to create the tableau's symmetrizer.
      Those permutations are applied to the tensor's indices and all
      the combinations are summed-up (and weighted), and boom: done.
      It is a complex problem.


      Index Gymnastics
      ================

      The 'e' method allows you to select or run on indices, as follows
      for a rank 5 tensor \f$ T_{ijkmn} \f$:

\verbatim      	  >>>M = T.e(0, Ellipsis, 2, 1, Ellipsis)\endverbatim
      Or in other words:
      \n
      \f$ M_{jn} = T_{0j21n}\f$
      While a contraction to a rank 3 tensor:
      \n
\verbatim      	      >>>A = T.contraction(0, 3)\endverbatim
	\n    
      represents:
\n
      \f$ A_{jkm} = T_{ijkim} \f$.
\n
      The inner and outer products, respectively are, for example:
\n
\verbatim      	  	>>>C = A.inner(B)\endverbatim
\n	        
      means:
\n
      \f$ C_{ijm} = A_{ijk}B_{km} \f$
\n
      and
\n

\verbatim		>>>C = A.outer(B)\endverbatim
\n
      \f$ C_{ijkmn} = A_{ijk}B_{mn} \f$.
\n
      In the later case, the rank 5 tensor is created at run-time and added
      to the Zoo. An antisymmetric product is also possible:


\verbatim   	     >>>C = A ^ B\endverbatim
		     \n
      which becomes:
      \n
      \f$ C_{ijkln} = A_{ijk}\epsilon_{klm}B_{mn} \f$.
      \n
      (though this not a true wedge product-- use
      geo.desic.clifford for that).

      Of course, for computational speed, the normal Vector and rank-2
      products are overridden with specific implementations. Never-the-less,
      problems arise:

      While it is straightforward to implement specific formulae for
      "normal" rank-1,2 operations, and it is also straightforward to
      do rank-N >= 1 with slower generalized formulae  (that computer and
      expand index iterable through a generous sprinkling of itertools),
      it is difficult to mix-them without resorting to so-called external
      polymorphism (read '_if_-then'). But I try. Apologies for any less-than
      perfect code. A similar problem arises with rank-0 objects-- but that
      is fairly well handled.

      @section geomorphic Morphic
      The geo.morphic sub-package is for geomorphology: the shape of terrain.
      monge.py lets you define surfaces, as sampled by \f$z_i = f(x_i, y_i)\f$,
      whence you can fit them with qudarics (taylor.py) and compute their
      various forms of slope and curvatures (gauss.py). It's all in good fun.

      Basically, you start with a Vector instance, v, that represents a local
      neighborhood in a map (so v.x, v.y are from a let-lon mesh-grid, and
      v.z is the height) and the put it into a monge patch

\verbatim      	        >>>m = Monge.fromvector(v)\endverbatim
\verbatim      		    >>>m = Monge(v.x, v.y, v.z)\endverbatim

      and let geo.morphic.monge.Monge do the work:

\verbatim		     		      >>>quad = m.quadric()\endverbatim
					      \n
      In fact, if someone gives you lat, lon and ellipsoid height in an array
      p, you can get the geoid corrected principal curvatures at the i'th
      point with axes defined wrt to north:
     
\verbatim                >>>ltp = LLH(*p).ltp(peg=(p[i, 0].lat, p[i, 1].lon, 0))\endverbatim
\verbatim                >>>k1, k2 = Monge(ltp.x, ltp.y, ltp.z + ltp.egm96()).quadrc().kmaxkmin()\endverbatim
			 \n
      I mean wow. For most packages, that request is an entire module--even
      a "task" with a stuckee; here, it's 2 lines, and on a whim at that.

      @section Magnetic
      The geo.magnetic subpackage has some modules related to
      electro-magnetism: The aforementioned faraday.py for
      polarization, and legendre.py--which has vector spherical
      harmonics. (The ease with which they are computed should
      convince you of the wonders of geo).
      
      There is a preliminary module, voigt.py, for anisotropic materials
      in both elasticity \f$ \tilde{\epsilon} \f$ and susceptibility
      \f$ \chi \f$. They require rank-4 tensors; however, because of symmetry
      they are reduced from 81 DoF to 21 and 36 DoF, respectively, and hence
      they can be represented by a symmetric (or not) 6 by 6 matrix.

      @section Centric
      The idea here is orbits, but it will probably never happen. keppler.py
      is empty.

      @section Politics
      Geopolitics includes various modules that involve some sort of
      interational standards, such as: time, units, and spectrum.

      @subsection LS Leap Seconds
      leapseconds.py provides a datetime object
      that can be substituted in-place of traditional
      <a href="https://docs.python.org/2/library/datetime.html#datetime-objects">
      python datetime.datetime</a> instances-- but it is leap second aware and
      will compute correct time differences across leap-seconds. Of course,
      this module's geo.centric.leapseconds.LEAPS cannot be predicted beyond
      6 months.

      @subsection Si System International
      The si.py module provides objects with dimensions, both base and derived.

      @subsection spectra Electromagnetic Spectrum
      hertz.py provides a collection of dictionaries that represent the
      various divisions of the electromagnetic spectrum, from static DC fields
      to the
      <a href="http://en.wikipedia.org/wiki/Planck_units">Planck Scale</a>,
      and a function to identify where a given frequency lies.

      @section tutorialu User's Tutorial
      See tutorial.py. It does all the stuff that you want to do, and some
      that you don't.
      
      @section tutoriald Developer's Tutorial
      This is
      <a href="http://en.wikipedia.org/wiki/SOLID_(object-oriented_design)">
      SOLID</a> OO, pythonic code. It attempts to maximize its maintainability
      score as determined by
      <a href="http://www.pylint.org">pylint</a>, the current standard for
      coding python. Exorcising anti-patterns leads to ravioli-code, which
      can be hard to read (but __"__ easy __"__ to maintain). The difficulty in reading
      arises because the implementation is hidden, while the algorithm
      is exposed. We learn from primitive-obsessed procedural code that
      to understand the code, we must understand the implementation. Begin
      unlearing that __now__. It is not easy, but it is well worth the
      effort.

      Good Luck.

      The goal is highly cohesive code, with low coupling. The mathematical
      nature of the objects with-in, make low coupling difficult--coordinates
      will depend on vector: they will always be coupled; like-wise for
      everything else.

      @subsection abc Abstract Base Classes.
      That so many mathematical object share behavior has lead to base
      classes (many base classes), where behavior is implemented once,
      abstractly, and the made concrete in the leaf classes. The python
      <a href="https://docs.python.org/2/library/abc.html">abstract base
      class</a> library is used extensively to signal to the developer:
      this is an ABC. Moreover, the abstractproperty and abstractmethod
      decorators are used to indicate where a concrete class *must*
      implement a class attribute or method (for example "rank"). Tensors
      must have a rank, but algorithms that are shared amoung tensors
      of different rank would be put into a method in an ABC. That method
      may reference a non-existent "rank" that __must__ be implemented
      in concrete tensor classes. Hence: the ABC restricts usage in an
      attempt to increase clarity-- it adds _no_ functionality. You may
      completely ignore any references to the _abc_ library when trying
      to understand what python is doing, while you can use references
      to the _abc_ library to understand how geo does business.

      @subsection poly Polymorphism
      What is it? Well having to use "print", "printf", or "sprintf"
      is NOT IT.
      Polymorphism means "do the right thing without being told."
      Hence, if you code:

\verbatim      		   >>>c = a * b\endverbatim
			   \n
      "c" better be the product of a and b. But what are a and b?
      	  \n\n
      \f$ c_i = ab_i \f$
      \n\n
      differs from
\n\n
      \f$ c = a_ib_i \f$
\n\n
      differs from
\n\n
      \f$ c_i = a_{ij}b_j \f$
\n\n
      differs from
\n\n
      \f$ c_{ij} = a_{ijk}b_k \f$
\n\n
      and so on. In these cases, internal polymoprhism--the
      highest--form works:
      The class dispatches the correct "a.__mul__" method to
      deal with the "*" operand.
\n\n
      What about:
\n\n
      \f$ c_i = a_i b \f$
    \n\n  
      differs from
\n\n
      \f$ c = a_i b_j \f$
\n\n
      differs from
\n\n
      \f$ c_i = a_i b_{ij} \f$
\n\n
      differs from
\n\n
      \f$ c_{ij} = a_i b_{ijk} \f$
\n\n
      Here, each case calls
      geo.metric.euclid.vector.Vector.__mul__, which there is none,
      it's kicked up to the super:
      geo.metric.euclid.euclid.Tensor_.__mul__. Now what?
      
      So one solution is to use case-switch like _if_-_then_ blocks to
      call the right function. That's too primitive. Instead, a hash
      table bound to the left operand looks up the rank of the
      right operand to get the functions (e.g.,
      geo.metric.euclid.vector.Vector._dispatch_mul). [Aside:
      That shows that
      things are further complicated by the quaternion and the
      geo.metric.frenet_serret.SpaceCurve]. That is  external
      polymorphism, where explicit code decides the function to call
      based on the argument-- this can only be solved elegantly
      with a language that supports
      <a href="http://en.wikipedia.org/wiki/Multiple_dispatch">
      multiple dispatch</a>, aka: multimethods.

      Well, we don't have that in python, so this is what you get.
      The solution presented is an attempt to avoid branchy
      dynamic code.
      The complexity is represented in a complex __static__ data
      structure--with one point of evaluation. It takes some
      getting used to, but it the right thing to do. If you disagree-
      read the next section. 

      @subsection flat Flatness
      The Zen of python requires flatness. The prior paragraph
      may seem in total violation of that--but it is not. That
      geo has deep taxonomies of subpackages and modules is a
      result of the complexity of the problem. The code itself
      is flat-- that is: each object does not violate the
      Law of Demeter: access is (mostly) 1-dot-deep. Moreover,
      deep _if_-then blocks and _for_-loops are avoided. The former,
      for example, uses the hash table to avoid something like:

      \verbatim
def __mul__(self, other):
    if isinstance(other, Geometric):
        if isinstance(other, Tensor_):
            if isinstance(other, Scalar):
                return dilation(self, other)
            elif isinstance(other, Vector):
                return dot(self, other)	
            elif isinstance(other, Tensor):
                return anterior_product(self, other)	
	    else:
		return (other.transpose() * self).transpose()
	elif isinstance(other, NonCartesian):
	    if isinstance(other, Polar):
	        return self * Vector(other.r * sin(other.theta) * cos(other.phi),
		                     other.r * sin(other.theta) * sin(other.phi),
				     other.r * cos(other.theta))
            elif isinstance(other, Parabola):
	         <suite>        
        else:
	    if isinstance(other, Quaternion):
	         return Quaternion(ZERO, self) * other
	      elif isinstance(other, SpaceCurve):
	         return self * SpaceCurve.vector 
     else:
         return Vector(other*self.x, other*self.y, other*self.z)
      \endverbatim
      that's not pretty. Not only is the cyclomatic complexity too high,
      it requires
      Cartesian vectors know about other coordinate systems--a
      total disaster (read: _antipattern_)--and
      utterly NOT FLAT.
      \n\n
      _for_-loops are generally avoided
      by internal iteration and the wonders of the
      standard library's
      <a href="https://docs.python.org/2/library/itertools.html">itertools
      module</a> (mad props
      to Raymond Hettinger).\n\n
      There is also question (previously "alluded too") of singleton vs.
      array components. That is mostly handled via internal polymorphism
      built into numpy's array broadcasting--but there are corner
      cases that can cause grief-- so you if get an odd result, it's
      quite possible that the operator overload was caught by
      numpy and not geo. (Note: this is why there is no "__array__"
      overloading in any of geo, as numpy looks for this and considers
      it a green-light to grab control).

      @subsection wc Wildcard Imports
      geo is full of wildcard imports:
      \verbatim
      	     from .module import *
      \endverbatim
      we all know wildcard imports are bad-- very bad. So
      what's up?

      __No__ module uses objects from a wildcard import, rather top-level
      subpackages collect their sub-modules key objects into their
      namespace. In support of that idea, all wildcard imports are from
      modules using the
      <a href="https://docs.python.org/2/tutorial/modules.html?highlight=__all__#importing-from-a-package">
      "__all__"</a> protocol, or from a sub-package's "__init__.py" which
      also used "__all__".
      The point is: this brings objects up to the top namespace. Hence
      "Vector" lives in "geo.Vector" -- and you have access to it via
      a __flat__ look up. That it got there via:
1 
\verbatim geo: from .metric import *
	  metric: from .euclid import *
          euclid: from .vector import * \endverbatim

      which got it from vector.py:

\verbatim
__all__ = ('scalar_triple_product', 'vector_triple_product',
           'scalar_quadruple_product', 'vector_quadruple_product',
           'Proj', 'Rej', 'Ref',
           'Vector', 'NULL', 'BASIS', 'X', 'Y', 'Z')
\endverbatim
      is not required knowledge. It's flat to the user.


      @section young Tensor Symmetries Revisited
      We all know the symmetries of a rank 2 tensors:
      \n\n
      \f$ S_{ij} = \frac 1 2 [T_{ij} + T_{ji}] \f$
      \n\n
      \f$ S_{ij} = \frac 1 2 [T_{ij} - T_{ji}] \f$      
      \n\n
      are the symmetric and antisymmetric parts of the tensor, \f$ T_{ij}\f$.
      They are not just (anti)symmetric under interchange of the indicies,
      they are closed under rotations. So the question is: how is this
      extended
      to higher rank tensors? What is the math?

      Well, the math is a deep dive. The invariant subspaces,
      aka irreducible representations (heretofore: irreps), of a rank \f$ N\f$
      tensor are related to the symmetric group \f$ S_N\f$, which is the
      group of permutation on \f$ N\f$-letter, via Schur-Weyl Duality (SWD). These
      are computed from Young Tableuax with \f$ N\f$ boxes via the
      Robinson-Schensted Correspondence (RSC) which in turn are related
      to the integer partitions of \f$ N\f$. The following seciotn works through
      the lowest nontrivial example, rank 3.

      What are the integer paritions of \f$ N=3\f$? While you can work these out
      by hand, the pascal.py package will do it for you:
      \verbatim
       >>>p = pascal.P(3)

       >>>print p
       P(3, k)

       >>>print type(p).mro()
       [<class 'geo.metric.schur_weyl.pascal.P'>,
       <class 'geo.metric.schur_weyl.pascal._Partition'>,
       <class 'geo.metric.schur_weyl.pascal.Triangle'>,
       <type 'long'>,
       <type 'object'>]
       \endverbatim
       Note that the
       <a href="http://mathworld.wolfram.com/PartitionFunctionP.html">
       Parition Function \f$ P\f$</a> object is fact just an extended long
       integer.

       The geo.metric.schur_weyl.pascal.P.partitions() method generates
       the paritions:
       \verbatim
       >>>list(p.partitions())
       [[3], [2, 1], [1, 1, 1]]
       \endverbatim
       of which there are three.

       Each parition can be represented by a
       geo.metric.schur_weyl.young.Diagram, respectively as follows:
       \verbatim
       In [169]: for d in p.young():
     ...:     print unicode(d)
     ...:     print
     ...:     
☐☐☐

☐☐
☐

☐
☐
☐
	\endverbatim
	Let's look at the seconds one (\f$ 3 = 2 + 1\f$):
	\verbatim
d = young.Diagram(2, 1)
	\endverbatim
	To get to the permutation group, one get a Standard Tableau
	corresponding to the diagram. A standard tableau has numbers from
	\f$ 1, ..., N\f$ in the boxes, with each row and each column strictly
	increasing. There are 2 ways to fill the diagram that meet those
	criteria:
	\verbatim
>>>for t in d.standard_tableaux():
     ...:     print t
     ...:     print
     ...:     
[i][j]
[k]

[i][k]
[j]
	\endverbatim
Note that I have filled the boxes not with number, but with ordered tensor
indices _i_, _j_, and _k_.
	Now we have to pick one of those tableaux:
	\verbatim
>>>t = d.fill(0, 1 ,2)

>>>t
Tableau_[[0, 1], [2]]

>>>print t
[0][1]
[2]
	\endverbatim
	(Note that the numeric representation starts at \f$ 0\f$, not \f$ 1\f$. This
	just makes it easier to work in a computer language that starts
	its indexing at zero.)

	At this point, we move into Schur-Weyl duality: what does this
	tableau have to do with the permutation group? Let start with the
	symmetric group on 3 letters:
	\verbatim
	>>>S3 = monte.Sym(3)
	\endverbatim
	It has 6 elements:
	\verbatim
>>>:for count, perm in enumerate(S3):
     ...:     print count, repr(perm)
     ...:     
0 
1 (12)
2 (01)
3 (012)
4 (021)
5 (02)
	\endverbatim
	Each geo.metric.schur_weyl.monte.Perm permutation is represented
	by its cycle structure, which shows the orbit of an element. (Note:
	the indentity, aka the nuetral elememt, is an empty cycle) They
	can also be show in two-line-format. For the last permutation, that
	is:
	\verbatim
>>>perm
(02)

>>>print perm
 | (0, 1, 2) | 
 | (2, 1, 0) | 
       \endverbatim
       where the latter format shows the trajectory of each element.
       
       Note that the Young tableau represents a permutation:
       \verbatim
>>>print monte.Perm.fromcycles(*t.cycles())
 | (0, 1, 2) | 
 | (1, 0, 2) | 
 \endverbatim
 But there is more. One can also consider the set of permutations
 that don't mix the letters in different rows:
\verbatim
>>>for count, perm in enumerate(t.Row()):
...:     print count, repr(perm)
...:     
...:     
0 
1 (01)
\endverbatim
The Row method yields permutations.
Likewise, one can consider the set of permutations that don't mix the
colums; however, there is a twist: when collecting these, we must keep
track of the parity of the permutation:
\verbatim
>>>for count, perm in enumerate(t.Col()):
...:     print count, repr(perm)
...:     
...:     
0 (1, )
1 (-1, (02))
\endverbatim
The Col method yields 2-ples of (parity, permutation).
The parity is either +1 or -1, and depends on the parity of the number of
transpositions in the permutation.

These two sets are called the Row and Column symmetrizers, respectively.
The Young Symmetrizer for the tableeua is constructed my taking the
product of these 2 sets (along with the parity):
\verbatim
for count, (parity, perm) in enumerate(t.S()):
     ...:     print count, ":",  parity, repr(perm)
     ...:     
0: 1 
1: 1 (01)
2: -1 (02)
3: -1 (021)
\endverbatim
Finally: if you apply those permutations to the tensor indices and
add the results up, you get an irreducible subspace of the tensor, which
can be show lexigraphically as follows:
\verbatim

>>>t.lexigraphicS()
'(+T.ijk +T.jik -T.kji -T.kij) / 3'

\endverbatim
This represents a mixed symmetry iredducible subspace.
\n
There are several way to apply that symetrizer to a tensor. We can use
the tableua's "syymetriz" method as follows:
\n
Start with a rank three tensor:
\verbatim
>>>T = Three(*(3*arange(27.)**2))
>>>print T.broadcast(lambda x: str(int(x)).zfill(4)) # nice up for display
[0000, 0003, 0012] [0243, 0300, 0363] [0972, 1083, 1200]
[0027, 0048, 0075] [0432, 0507, 0588] [1323, 1452, 1587]
[0108, 0147, 0192] [0675, 0768, 0867] [1728, 1875, 2028]
\endverbatim
and use the method:
\verbatim
q = t.symmetrize(T)
>>>print q.broadcast(lambda x: str(int(x)).zfill(3))  # again, nicen up
[000, -88, -352] [088, 000, -264] [352, 264, 000]
[000, -172, -520] [172, 000, -348] [520, 348, 000]
[000, -256, -688] [256, 000, -432] [688, 432, 000]
\endverbtim
OK. So how do we know that that is rotationally closed? Let's rotate it
by an arbitrary axis:
\verbatim
R = roll(34)*pitch(12)*yaw(-56)
qprime = R(q)
\endverbatim
Note that the versor \f$ R\f$ knows exactly how to rotate a rank 3 tensor--you,
the user do not have be concerned with all the index gymnastics.
\verbatim
>>>print qprime.broadcast(lambda x: round(x, 3))
[0.0, 222.464, -877.141] [-222.464, 0.0, 281.535] [877.141, -281.535, 0.0]
[-0.0, -86.987, 352.2] [86.987, 0.0, -113.674] [-352.2, 113.674, -0.0]
[0.0, 70.665, -516.829] [-70.665, -0.0, 182.12] [516.829, -182.12, -0.0]
\endverbatim
So I guess that's closed? The zeros are in the same spot, and, e.g.
\f$ q'_{yzx} = -q'_{xzy}\f$, but I am not convinced. Let's take the tensor
\f$ q'\f$ $ and ask it to show us its non-zero modules:
\verbatim
>>>qprime.show_modules()
[i][j]
[k] 
w=1641.91108164 
\endverbatim
which matches the orginal \f$q\f$:
\verbatim
>>>q.show_modules()
[i][j]
[k] 
w=1641.91108164 
\endverbatim

Note that that procedure was deepy complicated. Neverthelss, it reduces
to the simple stuff that we understand from lower rank tensors:

Scalars are less than trivial, they have ther permutation group on nothing:
there
is simply no index to permute. Meanwhile, vectors have one index, and
the permutation group on 1 element is trivial.

\f$N=2\f$ has two paritiotns (\f$2=2\f$ and \f$2=1+1\f$),
corresponding to the following diagrams:
\verbatim
>>>p = pascal.P(2)
>>>for d in p.young():
     ...:     print unicode(d)
     ...:     print
☐☐

☐
☐
\endverbatim
It should be easy to see that the first (last) diagram's Young symmetrizer
is total (anti)symmetric, corresponding to the familiar symmetric and
antisymetric tensors.

Dimensionality
==============
But wait, there is more. The Young diagrams have a remarkable formula,
called the "Hook Length Formula". The hook length is computed from
the arm and leg lengths. See:
\n\n
geo.metric.schur_weyl.young.Diagram.arm()\n
geo.metric.schur_weyl.young.Diagram.leg()\n
geo.metric.schur_weyl.young.Diagram.hook()\n
\n\n
for more. These formula can calculate the dimensions of the closed
subspaces, for any rank of tensor, over _any_ field.

Let's start with the basis diagram:
\verbatim
>>>d = Diagram(1)
>>>print unicode(d)
☐

\endverbatim
It is just a box. We can assign it any dimension and meaning we want:
\n\n
\f$D = 2 \rightarrow \f$ spinor / spin 1/2 fermion \n
\f$D = 3 \rightarrow \f$ Cartesian Vector, SU(3)-Flavor quark triplet \n
\f$D = 4 \rightarrow \f$ Minkowski 4-Vector \n
\n\n
Now, here comes the amazing part. We can take the tensor product of boxes
(which is formed by combining them _graphically_), and compute the
dimensions of the irreducible subspaces:
\verbatim
>>>for item in d ^ d:
...:     print unicode(item)
...:     print "{}\t{}\t{}".format(*map(item.dimW, [2, 3, 4]))
...:     print

☐☐
3	6	10

☐
☐
1	3	6
\endverbatim
So for combinig 2 quantum spinors (\f$D=2\f$), we get a symmetric triplet and
an antisymmetric singlet, ala:
\n\n
\f$ |1, 1\rangle = \uparrow\uparrow \f$
\n\n
Meanwhile, \f$D=3\f$ tells us:
\n\n
\f$ {\bf 3} \otimes {\bf 3} = {\bf 6}_S \oplus {\bf 3}_A \f$
\n\n
which says that (anti)symmetric rank 2 tensors have 6 (3) components in three
dimensions. In special relativity:
\n\n
\f$ {\bf 4} \otimes {\bf 4} = {\bf 10}_S \oplus {\bf 6}_A \f$
\n\n
so that the antisymmetric electromagnetic field strength tensor,
\f$ F_{\mu\nu} = \partial_{\mu}A_{\nu} - \partial_{nu}A_{\mu} \f$
has six componets (3 for the Electric field, and 3 for the magnetic field).
The symmetric stress-energy tensor has 10 components (1 for mass/energy
density, 3 for energy flux or momentum time derivative, and 6 for
spatial stress).

Now we can take it to 3rd rank:
\verbatim
In [378]: for item in d ^ d ^ d:
     ...:     print unicode(item)
     ...:     print "{}\t{}\t{}".format(*map(item.dimW, [2, 3]))
     ...:     print

☐☐☐
4	10

☐☐
☐
2	8

☐☐
☐
2	8

☐
☐
☐
0	1
\endverbatim
From the \f$D=2\f$ columns we learn: there is no antisymmetric combination
of 3 electrons. There is one symmetric combination of 4 dimensions
(\f$J = 3/2 \f$), and 2 mixed symmetry \f$J=1/2\f$ combinations.

For \f$D=3\f$, there is one antisymmetric tensor, the Levi-Civita tensor.
While in the quark model, this irrep is the \f$ \Omega^-\f$n baryon,
famoulsy predicted by Gell-Mann.

The application of Young tableau and symmetry spans a phenomenal amount
of physics.

Quantum Spins
=============
Spin 1/2:
---------
The addition of quantum spins is a special case of __SU__(2) irreducible
representations. For instance, a spin 1/2 particle has 2 eigenstates of
"alignment": the spin along and arbitary axis (taken to $J_z$)
has 2 eigenvalues:
$m = \pm frac 1 2$. When the spins of 2 indentical particales are combined,
the eigenstates of total $J_z$ are not eigenstates the individual particles'
$J_z$.
\verbatim
>>>u = racah.Ket(0.5, 0.5)
>>>print unicode(u)
|½, ½⟩
 
>>>print unicode(d)
|½, -½⟩

>>>d = racah.Ket(0.5, -0.5)

>>>print unicode(u*u)
[ |1, 1⟩]

>>>print unicode(u*d)
[√（½) |0, 0⟩ +
√（½) |1, 0⟩]

>>>print unicode(d*u)
[-√（½) |0, 0⟩ +
√（½) |1, 0⟩]

>>>print unicode(d*d)
[ |1, -1⟩]
\endverbatim
So that now, by inspection, we can construct linear combinations that a
(normalized) eigenstates of $J^2$ and $J_z=0$:
\verbatim
>>>print (u*d + d*u)/sqrt(2)
[ |1, 0>]

>>>print (u*d - d*u)/sqrt(2)
[ |0, 0>]
\endverbatim

Spin 1
------
The racah.py has the basis kets for spin 1:
\verbatim
>>>Z, P, M = [Ket(1, m) for m in (0, 1, -1)]
\endverbatim
so that:
\verbatim
>>>print racah.Z, racah.P, racah.M
 |1, 0>  |1, 1>  |1, -1>
\endverbatim
Then we can take linear combination:

Spherical Vectors:
==================
With the above tools we can now delve in to spherical vectors.
Mathematically, the basis vectors transform as the fundemental represent




\copyright
JEB WORLD LTD.
*/

