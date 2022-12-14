# THE Package for Geometry on Earth.

## Intro Introduction

The point of this package is to allow you to do manifestly covariant
     vector operations (and transformations), both abstractly and with
     respect to geodesy via affine spaces on Earth. Manipulation of
     Earth-bound coordinates and abstract vectors can be done without
     regard to their coordinate system. Moreover, the operations do
     not depend on any array-like properties of your objects: operations
     on single vectors or coordinates are the same as operations on multi-
     dimensional arrays of coordinates: _The_ _code_ _you_
     _write_ _is_ _the_ _same_ _for_ _either_ _case_ .

Furthermore, all transforming objects are (polymorphic)
[function emulators](https://docs.python.org/2/reference/datamodel.html#emulating-callable-objects), so that
     the code you write depends neither on the type of the object doing
     the transformation (rotation matrix, quaternion, Euler angles, etc.)
     nor on the object being transformed (Vector, Quaternion, Tensor, Scalar,
     etc.):

   
  		>>>x2 = T(x1)

 
 will rotate _x1_ to _x2_, regardless of how __T__ represents
     the rotation for  _x1_ a tensor of any rank, or quaternion. Both
     __T__ and _x1_ can be singleton or array-like, and _x2_ will have the
     approriate array structure (or not) according to numpy's broadcasting
     rules.

### Manifest Covariance

[Manifest Covariance](https://en.wikipedia.org/wiki/Manifest_covariance),bb
[coordinate free](https://en.wikipedia.org/wiki/Coordinate-free")
vector operations means you
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
     well with [object-oriented programming](https://en.wikipedia.org/wiki/Object-oriented_programming)
     While the [procedural](https://en.wikipedia.org/wiki/Procedural_programming)
     [imperative](https://en.wikipedia.org/wiki/Imperative_programming)
     and/or programmer may [JPL's](https://en.wikipedia.org/wiki/Jet_Propulsion_Laboratory)
     represent 
     coordinates as:

     jpl = [[34, 12, 6.1], [-118, 10, 18], 250]
     
     
he/she has to remember that address 0 contains 2 integers (latitude degrees
     and minutes) and a float (seconds) and so on, along with which 
     ellipsoid we're using, the OO user just says:
     
     
     jpl = WGS84.LLH(DMS(34, 12, 6.1), DMS(-118, 10, 18), 210)
     

     
and the object knows, and remembers, exactly what that all means-- and
     is much more than just records in a structure, because when you write 
     
     d = jpl - gsfc
     
it knows that the difference between 2 points in an affine space is
     a vector, and it knows how to compute it, regardless of how the
     coordinates of NASA Goddard were specified. That is, it might be that:

	gsfc = almanac.ALMANAC["AIRY 1830"].LLH(38.9915337403, -76.85, 548)

which is a different coordinate system with respect to a different
      ellipsoid -- but the meaning of the point does not depend on its
     representation. Hence:


	>>>print repr(jpl - gsfc)
	geo.metric.euclid.vector.Vector(-3622444.0, 177651.0, -426628.0)

is the ECEF vector connecting them.
     Furthermore, if you add a vector to a point:

     p = jpl + v
     
it doesn't matter if _v_ is represented in Cartesian, polar, or
     cylindrical coordinates- because that doesn't change the meaning of _v_,
     just the representation. Likewise, the affine point could be in ECEF,
     and the results are the same:

     >>>jpl + v == jpl.ecef() + v.polar() 
     	True 
	
	
(It is the duty of the OO-programmer to implement these choices
     without _if_ - _then_ clauses: unlike functions acting on arrays,
     objects know themselves and thus [dispatch](https://en.wikipedia.org/wiki/Dynamic_dispatch)
     the right code [when called upon](https://en.wikipedia.org/wiki/Late_binding)
     
Operations are independent of representations. (NB: that's, little _r_
     representations. [geo/](https://github.com/jebelz/GEO/)
     allows various group-theoretic big _R_
     Representations, and those you should not mix--but who would do that
     anyway?).

## Geo-Metric:
There are 2 main and 2 (or more) ancillary sub-packages:

### Euclid and Euler

#### Cartesian Vectors in R3.
The geo.metric package
[euclid.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/euclid.py) defines the basics for rank 0
[scalar.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/scalar.py), 1,
[vector.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/vector.py), 2
[tensor.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/tensor.py), and higher objects in real Euclidean
$\mathbb E^3 $
space,
     aka
$\mathbb R^3$.

#### Transformations in SO(3)

Matrix __SO(3)__  (via tensor.py), quaternion __SU(2)__ 
[hamilton.py](https://github.com/jebelz/GEO/blob/main/metric/euler/hamilton.py),
     transformations are implemented.
     Furthermore, various
     representation of the charts on __SO(3)__ ([charts.py](https://github.com/jebelz/GEO/blob/main/metric/euler/charts.py) allow Euler-angle
     and Tait-Bryant angle based representations. Needless-to-say,
     transformation of vectors is independent of the representation
     of the transform, _T_:

     v' = T(v)
     
for all classes of _T_. That is, T could be a "tensor", euler angle,
     versor, etc. Moreover, their argument doesn't need to be a vector:

	s' = T(v*v)  # == T(v) * T(v) == v*v          Rotates a scalar
	t' = T(v&v)  # == T(V) & T(v) == T*(v&v)*T.T  Rotates a rank-2 tensor.

It is an object's job, as a polymorphic function emulator, to use the
correct method to transform its argument.

##### Transformations GL(3, R) (aka Aff(3))
In the context of geo, an
     [affine transformation](https://en.wikipedia.org/wiki/Affine_transformation)is a general linear transformations from
     ${\mathbb R}^3 \rightarrow {\mathbb R}^3$
and are supported in affine.py. While __M__ can be any object that
     transforms vectors linearly, we are primarily concerned with rotations.
     The translation __b__ is always from the target space (which is trivial
     in geo's implementation). With general forms of __M__, one can
     implement
     [scalings](https://en.wikipedia.org/wiki/Scaling_(geometry)">)
     [similiarity](https://en.wikipedia.org/wiki/Similarity_(geometry))
     transformations
     [reflections]
     <a href="https://en.wikipedia.org/wiki/Shear_mapping">shear
     mappings</a>, pure translations,
     <a href="https://en.wikipedia.org/wiki/Homothetic_transformation">
     homotheties</a>, and rotations. The 3 pure
     infitesimal translations and their duals (3 infitesimal rotations)
     together form the so-called generators of rigid Euclidean space.

### Non-Cartesian Coordinates.
     
Support for non-Cartesian vector representations
     is provided. Spherical (aka polar)
     (geo.metric.euclid.kant.spherical.Polar), cylindrical
     (geo.metric.euclid.kant.cylindrical.Cylinder), and parabolic
     (geo.metric.euclid.kant.parabolic.Parabola) coordinates
     are implemented. Through the wonders or __OOP__ and
     [1st-class functions](https://en.wikipedia.org/wiki/First-class_function)
     you can dynamically create fully operational non-Cartesian 
     coordinate representations via the metaclass
     geo.metric.euclid.kant.CoordinateFactory. All you need to do
     is provide a function and its inverse (that is, to and from
     Cartesian coordinates).

So the Cartessian basis vector can be expressed in other coordinate
     systems:

	>>>print X.polar()	
	radius=1.0, theta=1.57079632679, phi=0.0

	>>>print Y.cylinder()
	rho=1.0, theta=1.57079632679, z=0.0

	>>>print Z.parabola()
	u=1.41421356237, v=0.0, theta=0.0

These vector classes are [adapter pattern](https://en.wikipedia.org/wiki/Adapter_pattern)
interfaces to the Cartesian Vector: that is, when used in operations,
     they convert themselves to Cartesian, do the operation, and then
     convert back to their original coordinate system. That allows
     you to mix coordinates systems seemlessly, because a vector is a
     vector, regardless of the coordinate system in which it is expressed.
     Hence $ \hat x \times \hat y = \hat z $ in any coordainte system:
     
     
	>>>print X.polar() ^ Y.cylinder()
	radius=1.0, theta=0.0, phi=3.14159265359

	>>>print (X.polar() ^ Y.cylinder()) * Z
	w=1.0


### Space Curves
     
Space curves and the Frenet-Serret formalism is available in
     [motion.py](https://github.com/jebelz/GEO/blob/main/metric/frenet_serret/motion.py).
     The space curve comprises an (array_like) Vector and a
     one dimensional time array, which serves as a parameter in the
     parameterization of the 3D curve.

Take an array like Vector representing a helical path:
   
     >>>t = linspace(0, 2*pi, 100000)  # parameter
     >>>v = X*sin(3*t) + Y*cos(3*t) + Z*t/1000
     >>>print type(v)
     <class 'geo.metric.euclid.vector.Vector'>
     
That array-like Vector can be made into a spacecurve by formalizing
     the parameter $t$:
     
     >>>v = v.spacecurve(t)
     >>>print type(v)
     <class 'geo.metric.frenet_serret.motion.SpaceCurve'>
     
The SpaceCurve now has veclocity, acceleration (linear and angular),
     torsion and curvature (and radii thereof), along with the Tangent,
     Normal, and Binormal coordinate system at every point, e.g.:
  
     >>>print v.TNB()[100000//6].broadcast(round)
	[-1.0, -0.0, 0.0]
	[-0.0, 1.0, 0.0]
	[-0.0, -0.0, -1.0]
       			
where the result has been rounded for clarity.
     
#### Complexification

The bulk of geo is all about real vectors spaces: while vectors
     can have complex components, they still live in $\mathbb R^3$ ,
     and as such,
     their bilinear inner product is not positive definite:

     >>>v = X + 1j*Y
     >>>print v*v
     w=0j
     
In [gibbs.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/gibbs.py), the
real vector.Vector objects are placed into
$\mathbb C^3$ : they are the
     same vectors, but their sesqui-linear inner product is positive
     definite. The outer product is also complexified, making them
     the ideal tool for representing polarization density matrices
     ([faraday.py](https://github.com/jebelz/GEO/blob/main/magnetic/faraday.py)).
     
     >>>v = v.gibbs()  # same vector, now lives in a new space
     >>>print type(v)
     <class 'geo.metric.euclid.gibbs.Gibbs'>

     >>>print v*v  # square is positive definite now
     w=(2+0j)

     >>> print v&v  #  outer product is Hermetian
     [(1+0j), -1j, 0j]
     [1j, (1+0j), 0j]
     [0j, 0j, 0j]
 

#### Higher Rank Tensors
     
3rd and 4th rank tensors, and a class-factory capable of making
     any order tensor are available in
    [three.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/three.py) and
    [four.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/four.py), respectively.
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
   
   
###   Einstein Summation Notation
     
Through the wonders of python and [complete class customization](https://docs.python.org/2/reference/datamodel.html#customizing-attribute-access)
     [Einstein summation notation](https://en.wikipedia.org/wiki/Einstein_notation) is _fully_ supported
     ([einstein/albert.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/einstein/albert.py)).
     For example,
     a second rank Tensor only has 9 attributes: 'xx', 'xy', 'xz',
     'yx', 'yy', 'yz', 'zx', 'zy', and 'zz', and it has a trace
     method:


     s = T.trace()


which explicitly computes:
     
    
     s = T.xx + T.yy + T.zz
    

which is indeed the trace of the tensor. One may also code:
     
    
	s = T.ii
   
Of course, their is no 'ii' attribute (or property). Nevertheless,
     geo regonizes this as Einstein summation and goes ahead and computes
     the trace. Likewise for the transpose:

     
     t_ji = T.ji  # same as T.T, or T.transpose(1,0)
     
(The reversed alphabetical order indicates transposition). Rows and columns
     can also be extracted (as Vectors):

     
     row0 = T.xi
     col2 = T.ix
    

This is the preferred method of extracting "rows" and "columns" because
     Tensor don't have rows and columns, and the names are antithetical to
     geo's ethos.

Of course, this seems rather cute for rank-2 Tensors (and down right
     silly for Vectors), it is quite powerful when dealing with high
     rank tensors: You really can write code that looks _exactly_ _like_
     the equations in your text book. For example:
 
 
The explicit cross product (as the partial trace of a rank-5 tensor) is coded as:


     
     a_cross_b = (EPSILON & a & b).ijkjk
     

Similarly, the determinant
can be computed via triple contraction of a rank-6
     tensor (with explicit labeling of the fixed indices):
     
     det_t = (EPISLON & T.xi & T.yi & T.zi).ijkijk
    
Another formulation as 
     the  6-fold contraction of a rank-12 tensor;
     $$\det(T)=\frac{1}{6}\epsilon_{ijk}\epsilon_{lmn}T_{il}T_{jm}T_{kn}$$
     
is coded as:
    
     det_t = (EPSILON & EPSILON & T & T & T).ijklmniljmkn / 6
    

Note that the code matches math equations _exactly_. The latter
     implicitly runs a 12-deep loop to build and then contract the
     531,441 components of a rank-12 tensor--while it is not the most
     computationaly efficient algorithm, it does demonstrate the
     simplicity of Einstein summation notation. Moreover, it is not
     useless, as formulations with similar levels of complexity arise
     when considering the rotational symmetries  of rank-3+ tensors...


### Irreducible Representations and Subspaces
     
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

     


 #### Spherical Tensors: 3 x 3 = 5 + 3 + 1
     
 This does not refer to spherical coordinates, rather the decomposition
     of Vectors ([wigner/eckart/one.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/eckart/one.py)) and
     Tensors ([wigner/eckart/two.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/eckart/two.py)) into irreducible representations
     of __SO(3)__.
     Thus, the _x_ or _xy_-like
     attributes are replaced by
     orbital and
     azimuthal quantum numbers: $(l, m)$ -- also knows as degree and order.
     The signifcance of these is that they are eigentensors of rotations
     about a fixed axis.

 For rank-1, the 3 spherical basis states are eigenvectors of z-rotations.
     This is similar to Cartesian vector addition being isomorphic to
     the translation operators acting on position eigenstates (which are plain
     old position vectors). Because translation commute, the construction
     of higher rank tensors is straightforward. A rank-N basis state is just
     the polyadic product of N basis vectors (e.g.
     $ {\bf \hat{x}\hat{x}\hat{y}}$ ). Rotations do not commute,
     so polyadic products of spherical basis vectors are not necessarily
     eigentensors of rotations. This is entirely analgous to the construction
     of [spherical harmonics in Cartesian coordinates](https://en.wikipedia.org/wiki/Spherical_harmonics#Spherical_harmonics_in_Cartesian_form),
     and it gets complicated, fast. Nevertheless, geo goes there--though
     it requires and entirely complex subpackage
     ([geo.metric.wigner](https://github.com/jebelz/GEO/blob/main/metric/wigner/)).

In geo's implementation, the elements of spherical tensors are not
     accessed by
     attributes, rather, the instances are containers
     of irreducible
     representations and their azimuthal eigenstates.
     Methods allow conversions between
     the 2 representations (
     [ec.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/eckart/ec.py),
     [ka.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/eckart/ka.py), and
     [rt.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/eckart/rt.py)).
     So for instance:
     
     
     >>>M = Tensor(xx, xy, xz, ..., zz)
     >>>print M.xx, M.xy, M.xz, ..., M.zz
     >>>
     >>>T = M.spherical()
     >>>print T[0, 0], T[1, -1], T[2, -2], ..., T[2, 2]

     
are represenations of the same geometric object, as very different
     python objects.

Applications to polarization observables are
     demonstrated in faraday.py .

The full transformation matrices associated with
     spherical representation of any dimensions 2j+1 for integer and
     1/2-integer representations (j)  are in
     [wigner.py]https://github.com/jebelz/GEO/blob/main/metric/wigner/wigner.py), which also
     implements the decomposition of tensor-product spaces into
     sums of irreducible spaces, with the help of
     [clebsch_gordan.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/clebsch_gordan.py). 
     Some assistance with decomposition of representations provided in
     [racah.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/racah.py),
     where products of kets (eigenstates of individual
     angular momentum/z-projection operators) are turned into sums
     over different kets (eigenstates of the total angular momentum
     and total z-projection). (This is truly dialed-in python, _Eds._)

The preliminary implementation of rank 3 ([wigner/eckart/three.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/eckart/three.py)) and
     rank 4 ([wigner/eckart/four.py](https://github.com/jebelz/GEO/blob/main/metric/wigner/eckart/three.py)) is underway; however, these are difficult
     (i.e., publishable) to implement.

Spherical decomposition of complete Cartesian tensors requires and
     understanding of isotropic tensors, which are available up to 
     rank 6 in
     [syzygy.py](https://github.com/jebelz/GEO/blob/main/metric/euclid/syzygy.py): these
     are tensors that are independent of coordinate system, which for
     rank 2 and 3 are the familiar Kroncker Delta
     (geo.metric.euclid.tensor.DELTA) and Levi-Civitta Tensor
     (geo.metric.euclid.three.EPSILON), repsectively.
     For higher ranks there are more than one isotropic tensor per rank,
     and they are not all linearaly independent.

Additional tools for the alternating tensor, aka Levi-Civita (pseudo)-
     tensor/symbol are in
     [levi_civita.py](https://github.com/jebelz/GEO/blob/main/metric/levi_civita.py). The
     [generalized Kronecker deltas](https://en.wikipedia.org/wiki/Kronecker_delta#Generalizations_of_the_Kronecker_delta)
     can be computed in
     [kronecker.py](https://github.com/jebelz/GEO/blob/main/metric/kronecker.py).
     

## Schur-Weyl Duality
     
Rank-2 tensors have two subspaces closed under rotations which
     are defined by their symmetric and antisymmetric parts. These
     are constructed by symmetric and antisymmetric permutations of their
     two indices--which is rather simple-as there are only 2 permutations.
     Higher rank _n_ tensors with _n_ indices have more complicated
     subspaces,and they can be computed from the various _n!_ permutations
     of their indices. This is the domain of
     [Schur-Weyl Duality]("https://en.wikipedia.org/wiki/Schur???Weyl_duality")
     The symmetric group on _n_-letters,
     $ Sym(n) $ can be decomposed into irreducible representations based
     on the integer partitions of $ n $ (really?). So to that end,
     there is monte.py (rimshot), which provides classes for the
     [symmetric group](https://en.wikipedia.org/wiki/Symmetric_group),
     [permutations](https://en.wikipedia.org/wiki/Permutation),
     and [cycles](https://en.wikipedia.org/wiki/Cyclic_permutation)
     The [Robinson-Schensted Correspondence]
     (https://en.wikipedia.org/wiki/Robinson???Schensted_correspondence)
     then relates irreducible
     representations of the symmetric group to
     [Young Tableaux](https://en.wikipedia.org/wiki/Young_tableau)
     [young.py](https://github.com/jebelz/GEO/blob/main/metric/schur_weyl/young.py). They provides a
     link to partitions of integers, which is taken up in
     [pascal.py](https://github.com/jebelz/GEO/blob/main/metric/schur_weyl/pascal.py),
     which starts with Pochammer symbols and winds up with integer
     partitions. Additionaly, the number triangles are extended to
     many of the combinatoric polynomials, leading naturally to the
     computation of Jacobi
     polynomials, which are required for the Wigner D matrices.
     Also: since Young diagrams are like a multiset, a multiset object
     is available in
     [de_bruijn.py]((https://github.com/jebelz/GEO/blob/main/metric/schur_weyl/young.py)),
     but it's really just a hashable
     [collections.Counter.](https://docs.python.org/2/library/collections.html#collections.Counter)
     
The theory of conjugacy classes in finite
     groups necessitates multisets with keys that are themselves multisets.

This subpackage also supports some results from the representation
     theory of finte groups, permutations, and combinatorics. The
     Young Tableaux formalism also has applications to the representation
     of Lie groups, which are commonly used to understand the quark
     (and separately the gluon) structure of matter.
     [gelfand_tsetlin/__init__.py)](https://github.com/jebelz/GEO/blob/main/metric/wigner/gelfand_tsetlin/__init__.py)

## Geo-Detic:

#### Coordinates and Ellipsoids
 
The [geo.detic](https://github.com/jebelz/GEO/tree/main/detic)
package implements geo.metric objects on the Earth
     via the introduction of coordinates
([coordinates.py](https://github.com/jebelz/GEO/blob/main/detic/coordinates.py)), which can
     be represented in Earth Centered Earth Fixed
     (geo.detic.coordinates.ECEF), geodetic
     (geo.detic.coordinates.LLH), or local tangent plane
     (geo.detic.coordinates.LTP)
     or tangent sphere (geo.detic.coordinates.SCH) relative to a
     geo.detic.origins.PegPoint.

The coordinates are [abstract](https://docs.python.org/2/library/abc.html),
     and only have [concrete](https://en.wikipedia.org/wiki/Class_(computer_programming)#Abstract_and_concrete)
     meaning with respect
     to an ellipsoid (see: [ellipsoid.py](https://github.com/jebelz/GEO/blob/main/detic/newton/ellipsoid.py)). 
     Although preference is given to WGS84 (wgs84.py), a variety of
     historic ellipsoids are available ([almanac.py](https://github.com/jebelz/GEO/blob/main/detic/newton/almanac.py)), as well as user-defined
     ellipsoids and others in the Solar System
     ([copernicus.py](https://github.com/jebelz/GEO/blob/main/detic/newton/copernicus.py)).
     
Analytic computations of covariant and contravariant basis vectors
     linking LLH and SCH to their canonical Cartesian coordinates
     systems (ECEF and LTP, respectively) are included ([christoffel.py](https://github.com/jebelz/GEO/blob/main/detic/christoffel.py)).
     That is: you always have both the tangent basis vector which are not-necessarity
     normailzed vectors running along the direction of change of s, c, and h,
     and the cotangent basis vectors
     which are the gradients of the local coordinates (that is, they're
     normal to surfaces of constants s, c, and h--and likewise for latitude,
     longitude and height with resepct to ECEF).
     
## The Geoids     
     
Support for interpolation
     on geoids is provided (egm08.py, egm96.py). While the former
     requires the
     <a href="http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/">
     1 minute grid</a>, the latter has the 1-degree grid self-contained:
     __EGM96__ is contained within this package.


### Universal Transverse Mercator
UTM (utm.py, kruger.py) support is under development. A class for seamless
     use of DD MM SS.ssss format is provided (preliminarily) in dms.py

     
## Geo-Desic
     
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
     infinitesimal rotations-- plop them into a $ \exp{i\ldots}
     $, and you have rotation matrices. It's a little confusing,
     because the 2 component spinor are fundamental, while one would
     think the $2^2-1=3$ are-- just because their are 3 generators
     of the algebra.. it's the proliferation of "3"s that have
     totally different meanings that causes confusion (at least for me).

Preliminary implementation of Lie algebras (lie.py), and SU(3)
     (gell_mann.py) are dubious, at best-- but you can convert the
     3-component fundamental representation (per the prior paragraph)
     into the $3^2-1=8$ adjoint representation
     (geo.desic.hilbert.gell_mann.EIGHTFOLD_WAY)
     and that is why
     3-color charges need 8 gluons, while $ 2 $ (weak) isospins
     need 3 pions (W bosons). (This package also makes clear the difference
     between fundamental and adjoint representations, which is a bit
     obtuse in the rotation group).


# Computational Implementation
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
     
## ARRAYS vs. VECTORS:
   
This package is a fully polymorphic object oriented vector and coordinate
     package, focused on geo-location on an arbitrary oblate ellipsoid of
     revolution.

It is for use with large data sets represented by [numpy.ndarrays](http://www.numpy.org)-
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
		 
## POLYMORPHISM:
    
This is the point of OO, and it is strongly supported. The point is, you
     DO NOT need to know which kind of coordinate you have when using it. For
     example if you have data, "platform", referenced to a coordinate system,
     and you NEED that data in geodetic coordinates (LLH), you just call the
     method:

     >>>platform.llh()
     
It does not matter if "platform" is in LLH already, or in ECEF, or in SCH,
     or in a tangent plane coordinate system. Moreover, if you have another
     point, "target",

     >>>platform - target
     

is the vector joining them: REGARDLESS OF THEIR COORDINATE SYSTEM-- you
     do not need to specify it-- the vector will be computed in the Cartesian
     coordinate system at platform's origin (even if the coordinate system is
     not Cartesian). This is because coordinates are treated as elements of
     an affine space--
     vector spaces without an origin.

All coordinates are defined on an ellipsoid-- and they carry that
     information with them. You can translate to another ellipsoid, e.g.,
     AIR1830, via:

     >>>platform | almanac.AIR1830
     

Again-- you don't need to specify if "platform" is LLH, SCH, or the local
     tangent plane (LTP)-- "platform" knows that already and converts itself
     to the same sub-class of coordinates. That's the point of polymorphism.

#### First Class Transformation:
     
_Question_:

Why are their so few functions for coordinate transformations? Why are
     there no 3x3 matrices filled with (bug prone) long combinations of trig
     operations?

_Answer_:

Transformation are 1st class objects that permit composing operations.
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
     [operator overloading](https://docs.python.org/2/reference/datamodel.html#emulating-numeric-types)
     of "\*" ("\_\_mul\_\_")and "~"
     ("\_\_invert\_\_"),
     respectively, with some support for interpolating between transformations
     via overloading "\*\*" ("\_\_pow\_\_") with non-integer exponents. (Python's
     built-in [pow](https://docs.python.org/2/library/functions.html#pow)
     function has a 3 argument option, which implements
     [spherical linear](http://en.wikipedia.org/wiki/Slerp)
     interpolation.
     
### Array or Singleton:
     
In the above example, "platform" or "target" could be singletons or arrays
     --it doesn't matter. The resultant vector will match the inputs-- without
     you having to tell it. If the objects are made out of numpy-arrays, then
     the result will be too-- numpy does all the work-- neither the user nor
     the developer has to worry about loop indices when using
     _or_ writing code.
    
With that in mind, here are the details:

## Introduction
     
Vector exists so that coordinates can exist. The ellipsoid exists so that
     the coordinates make sense. The goal of this sub-package is to make ALL
     vector and coordinate operations simple overloaded math operations or
     method calls.

Moreover, the objects are completely divorced from computer "arrays", and
     hence, you can use a single instance to represent a bunch of vectors, or
     points on a map, or orientations, etc., and ALL the broadcasting is done
     by numpy. This is a very important point to digest, for a vector, v, or
     tensor T:

	>>> v[2], T[0, 1]

have NO MEANING with respect to vectors and tensors. They are not:

 	    >>> v.z, T.xy

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
     $e^{i{\bf\vec J\cdot{\hat n}}\phi}$ generates rotations
     about a unit vector ${\bf \hat n}$ . Even worse is
     geo.detic.pauli.J(). Here the vector's components are numpy matrices
     of dimension $2j+1$, represented internal-degrees of freedom.
     This is all implemented polymorphically without any conditionals
     (_if_ _then_); you simply cannot do this with your procedural Matlab
     arrays and matrices.)

DO NOT CONFUSE TENSOR PROPERTIES WITH ARRAY PROPERTIES:

 	"__getattr__" 
	
gets tensor properties, e.g:
     
     m.xy  
     
While:
     
     "__getitem__"
     
gets array properties,  e.g:
   
   	v[10:1000:3, -30:]
   
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
     

## Euclid
Vectors (and other Tensors) in 3-Dimensional Euclidean Space (__E3__)

Vectors are important objects. We're talking about a physicist's vectors.
     That's a geometric object that lives in Euclidean 3-space (for us) and does
     NOT depend on the reference frame. That "vectorized" and "vector" have
     meaning in computer science is entirely unrelated and is just a point of
     confusion. Try to FORGET IT.

There are several ways to consider vectors. One way is as a rank-1
     tensor, that is, as a linear function from
     $ {\mathbb R}^3 \rightarrow {\mathbb R}$--this is a very powerful view
     point for characterizing physical laws--and is supported in geo:
     
$$ {\bf \vec a(\vec b} \equiv {\bf \vec a\cdot \vec b} $$.

Nevertheless, we'll define them by how they transform.
     Vectors are defined by their transformation properties: you have to rotate
     them by 360 degrees to remain unchanged. Well, that is entirely useless
     for computer programming. What is useful is that:

(1) They transform as:
			$$	v'_i = T_{ij}  v_j  $$
		
(2) In a given reference frame, they can be represented by 3 components:
                      
	v_i  for i = 1, 2, 3 
	
	
So there you go. Vectors have 3 things, and that defines their 
     representation in a frame. That leads to an important pythonic note:

Vector (Tensor, Coordinate, etc..) objects are entirely defined by their
     "__init__" arguments. To that end, a vector is given
     by its x, y, and z attributes:

	       	  >>>v = Vector(x, y, z)

There shall be no nullary instantiation followed by setter calls:


	>>>v = Vector()  # NO
	>>>v.setX(x)     # NO
	>>>v.setY(y)     # NO
	>>>v.setZ(z)     # NO NO NO

Fair enough, now we have a vector. What do you do with that? Well, you
     don't call a function. The vector should have everything you want to do
     with it defined by methods and operators (magic methods):
     
## Magic Methods:
     
The following invoke magic methods:
     
   	+v = v.__pos__()

Doesn't do anything.It returns "v", though it could return a copy of v.

     	-v = v.__neg__()


Negation, implemented as negation of each component

   	v + v' = v.__add__(v')

vector addition [Note: no constants allowed]:

     	v + 7 = v.__add__(7)

#raises a
     geo.utils.exceptions.NonCovariantOperation error (exceptions.py)]

    	v - v' = v.__sub__(v')

vector subtraction 


		v/4	  =  v.__div__(4)
		v*0.25	  =  v.__mul__(0.25)    
		0.25*v    =  v.__rmul__(0.25)

are all dilations. Reflected multiplication is the simplest, as it can
_only_ be a dilations

     		v * v' = v.__mul__(v')

Scalar product

  		v**2 = v.__pow__(2)

spherical linear Interpolation:c"\_\_pow\_\_" takes any positive integer, though "2" is the only defensible
argument

		abs(v) = v.__abs__() -> ||v**2|| ->  v.L2norm() = v.norm(l=2)

all need to be defined, and some don't, like:

        	v ^ v' = v.__xor__(v')

cross (wedge) product

 	  	v & v' = v.__and__(v')

outer (dyad) product

		~v   =  v.__invert__() = v / v*v

is a special case for converting covariant and contravariant vectors
    in local orthonormal bases and their conjugate global basses.
    
### Projection, Rejections, Reflection:

	a >> b = a.__rshift__(b) = a.projection(b) = Proj(b)(a)
	a << b = a.__lshift__(b) = a.rejection(b) = Rej(b)(a)
	a | b  = a.__or__(b)     = a.reflection(b) = Ref(b)(a)


Respectively. (Projection, rejection, and reflection operations
	for vectors, tensors, and quaternions are unified in cauchy.py.)
	The right column are tensor operators:\n
	geo.metric.euclid.vector.Proj'\n
        geo.metric.euclid.vector.Rej'\n
	geo.metric.euclid.vector.Ref.\n

Moreover, what about "\_\_getitem\_\_":


	v[3]                = v.__getitem__(3)
	v[-1]               = v.__getitem__(3)
	v[start:stop:step]  = v.__getitem__(slice(start, stop, step))

Do those have meaning? Yes they do, but they HAVE NOTHING TO DO WITH
     VECTORS. That is:

	v.x is not v[0]
    	v.y is not v[1]
	v.z is not v[2]
     
No.
     

     v[start:stop:step] is
     Vector(v.x[start:stop:step], v.y[start:stop:step], v.z[start:stop:step])


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
     
### At A Deeper Level
     
That Vectors can be more than singletons is not just a computer science
     convenience-- it's a fact of life. Consider the Pauli matrices (pauli.py).
     They can be combined into a single $\mathbb R^3$ vector,
$ {\bf \vec{\sigma} } $ -- there is but one vector, but its components
     are 2 x 2 matrices, representing internal degrees of freedom-- the
     fundamental representation of SU(2). They could be spin up and the
     orthogonal
     spin down (note: up is NOT minus down!... mind-blown), or
     isospin "up flavor" and
     "down flavor", or polarization degrees of freedom. It's up to you.
     The point is to convince the dear reader beyond all objections that
     vectors are not 3 numbers in an array.

#### Regular Methods:
     
The 1st regular method is:

	 v.iter()

it returns "x", "y" and "z" in order via an iterator.
      
	v.tolist() --> [v.x, v.y, v.z] = list(v.iter())

does it all at once.

Other methods, some called by magic methods, are straightforward:

   	    v.dual()

contract with Levi-Civita symbol and make a Tensor (which
      will do the cross product via posterior multiplication)

      	   v.dot(u)

 Dot Product: contract with another Vector and make a Scalar

     	  v.cross(u)

Cross Product: contract the dual with another Vector and make
      a Vector (There is no distinction between vectors and axial
      vectors)

		v.outer(u)

Outer Product: with another Vector, make a dyad, that is a Matrix

    	        v.dilation(alpha)

change length with multiplication or division by a Scalar or number.
       There are also right_dilation and left_dilation method, if your
       vector's components don't commute with the scalar (which is only
       a scalar in $\mathbb R^3$, and not in the 'other' space.)

   v.dyad(u)

Is the outer product.

     	     v.hat()

Is dilation by 1/abs(v) --> returns a Vector's unit Vector.

      	   v.L2norm()

L2-norm is just the norm.
      
There are also bonus functions for not only computing projections,
      rejections, and reflections, but also functions for creating
      like-wise operators:
    
	v.projection(u)   v >> u    u.projector()(v)    Proj(u)(v)	
	v.rejection(u)    v << u    u.rejector()(v)     Rej(u)(v)	
	v.reflection(u)   v | u     u.reflector()(v)    Ref(u)(v)	
	
Rotations about arbitrary vectors:

	v.versor(angle)

are given by the resulting unit quaternion.

#### Special Numpy Methods:
      
You can call numpy array methods like "mean", "sum", "cumsum" on tensors
      and get a tensor result--or on a coordinate to get a coordinate result.
      That is:

	>>>platform.mean()

will be the average position of a platform in LLH, ECEF, SCH, LTP, etc--
      whatever it is in. While for a 2-d array like SCH coordinate, "dem ",
      representing a digital elevation model (DEM) (with geo-located pixels):

	   >>>dem.mean(axis=0).h

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

    	     sigma = ((v-v.mean())&(v-v.mean())).mean()

is their correlation tensor. (Do you know how much code that would be with
      Fortran arrays?).

Furthermore, you can project that correlation tensor into it's
      irreducible spherical components at will-- it's just a method call that
      executes 1 or 2 lines of code.


#### Scalar:
   
Scalars are NOT PRIMITIVES! A single scalar can be represented (poorly)
     by a single floating point number-- but that is no reason to do it that
     way. 

When you square a vector, you are going to get a scalar. Why is this not
     just a number? There are several reasons:

(1) Scalars are rank-0 tensors, and hence, have tensor properties that need
     to be captured in a class. For example: they are invariant under rotations
     and transform according to:

                               s' = s

While this is trivial, it is formally a trivial
	representation of SO(3), it is
by no means "just a number". Consider the would-be scalar:

$$ s = \bf{ \vec A\cdot}(\bf{\vec B \times \vec C}) $$
    
which is really a pseudo-scalar. In the full Clifford algebra off
     Euclidean 3-space, it's a trivector, or 3-blade, and that cannot
     be confused with "a number".

(2) On a more practical front: If you are dealing with numpy arrays and 
     don't have scalars, you will 
     get burned when doing reflected multiplication. Let's say you want to scale
     a 100,000 point motion history, v, with:

  	     (1/(v**2))*v

Here v.x, v.y, v.z are 100,000 element numpy arrays. If v\*\*2 is just a
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


##### Rank 2 Tensor:
    


There are (at least) 3 ways to look at rank-2 tensors:

##### Component-wise
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
##### As Transforming Objects

(2) Geometric animals, U,  that transform as

$$ U'_{ij} = T_{ik}  T_{jl}  U_{kl} $$

##### (Bi)Linear Maps
 Tensors can be defined as bilinear maps from:
 
$$ {\mathbb R}^3 \times {\mathbb R}^3 \rightarrow {\mathbb R} $$
	 

geo chooses to represent tensors 9-components, representing the
weights of 9 Cartesian dyads--this is standard. Moreover, they
      are used to represent the so-called
     direction-cosine-matrix (DCM), even though a DCM is technically not
     a tensor-- the python Zen (import this) "practicality beats purity"
     is respected--the user should look upon their rotation matrices as just
     that, eventhough they live in a tensor class.


###### Rows and Columns:

Rows and Columns have no physical meaning and are hence, meaningless.
      They are at best a linear algebra concept: matrices have rows and
      columns, but tensors are not matrices. With respect to tensors,
      they are a typographical construct-- a relic of writing frame dependent
      representation of geometric objects on paper. I do not support them, and
      frankly, I don't know what "row major" or "column major" means. My tensors
      have indices, and each one is as "major" as the other.

Nevertheless, there are methods:


	tensor.cols() -->  T.ix, T.iy, T.iz
	tensor.rows() -->  T.xi, T.yi, T.zi

iterate over 1st and 2nd tensor indices -- and
      indices (slots, really) have meaning.
      I just happen to call it rows (cols)--but they're NOT rows (cols).
      It's handy though, if
      you want to know one of the important tensor invariants:
  
$$ m_{1i}  m_{2j}  m_{3k}  \epsilon_{ijk} $$
      
which is the determinant, you can compute it by unpacking the rows in the
      geo.metric.euclid.vector.scalar_triple_product function:
      
	 >>>scalar_triple_product(*tensor.rows())

and that does make a (pseudo) scalar our of it. The other scalars are the
      trace and the L2-Norm, and there are methods for those too. To summarize,
      for a matrix/tensor T (with functions from the tensor.py modules):
       
## Rank-2 Tensor Invariants:
      
Some invariants:

	abs(T) = T.L2norm() = T.norm([l=2])

$$ \sqrt{\sum_{ij}T_{ij}^2} $$
 
	Tr(T) = T.trace() = T.contract(0, 1) = T.ii

$$ \sum_{ij}T_{ij}\delta_{ij}\equiv T_{ii} $$

	det(T) = T.det()

      
$$ m_{1i}  m_{2j}  m_{3k}  \epsilon_{ijk} $$
 
The 3 invariants can also be expressed as coefficients of
      the charactertics polynomial:

	T.poly1d()              # a numpy.poly1d object
	T.characteristics()     # a tuple of Scalars.
	T.J(n); n = (1, 2, 3)   # The traditional n-th order invariants.

###### Rank-2 Tensor's Vectors:
      
Vectors from Tensors:

   	      T.dual()

$$ \frac{1}{2}\epsilon_{ijk}T_{jk} $$
 
     	T.vector()

$$ \epsilon_{ijk}T_{jk} $$ 

#### Rank-2 Tensor 
      
Tensors from Tensors:


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


Other methods include:
     
	T.spherical()            # Irreducible Spherical Representation
	T.eig()                  # Eigenvalues and Eigenvectors
	T.angle_of_rotation()    # Angle of rotations (for a DCM)
	T.axis_of_rotation()     # Axis.....
	T.sqrt()                 # Square Root
	T.exp() = exp(T)         # 1 + T + T**2/2 + T**3/6 ...
	T.log() = log(T)         # Inverse of exp(T)
	T.svd()                  # Single Value Decomp into dyads.
	
Some other operations are:

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

### Higher Rank Tensors:

For example, a rank-4 tesnor is diplayed as follow:

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

Einstein summation notation is supported up to rank-15 (which requires
	over 16,000,000 lines to print). Let's hope your tensor problems aren't
	that deep.

### Charts on SO(3)
The euclid.py module is all about Tensors-- and how you add, subtract,
      and multiply them. The rank-2 tensor also transforms vectors--but it is
      not alone. There are many ways to represent rotations in
      $\mathbb R^3$, and
      collectively, they are known as charts on __SO(3)__-- the rotation group.
      They live in charts.py. A nice introduction is provided here:

      http://en.wikipedia.org/wiki/Charts_on_SO(3)

#### Versors
     
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

$$ 1, i, j ,k  $$

(called "w", "x", "y", "z" axis)
      with: $i^2 = j^2 = k^2 = ijk = -1 $ 
      
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


###### Alias or Alibi?
      
Another point of confusion is "what is a transformation?". Well, alias
      transformations leave the vector unchanged and give its representation in
      a different coordinate system. (Think ALIAS: same animal, different name).
      
Meanwhile the Alibi transformation leaves the coordinates fixed and
      transforms the vector. (Think ALIBI: I wasn't there, because I was here)
      
What about pre or post multiplication? That is:

$$ v'_i = M_{ij}v_j $$

or

$$ v'_j = v_i M_{ij} $$
      
I have chosen the former for alibi (active) transformations. There
      are 2 reasons for this choice:

1) Euler-angle and Tait-Bryan representations of rotations are
      traditionally active; moreover, the more familiar pre-multiplication
      leads to break down of Tait-Bryan aircraft coordinates: Yaw, Pitch,
      and Roll to be:
      
$$ {\bf M}_{ypr} = {\bf M}_{\rm yaw}{\bf M}_{\rm pitch}{\bf M}_{\rm roll} $$
      
That leads to a great deal of simplification when converting between
      representations.
     
2) When dealing with quantum operators, state-changing
      pre-multiplication is the most common.

3) For non-unitary / non-orthogonal transforms, it makes sense.
     
4) When transforming objects of equivalent rank (Rotation matrices
      operating on rank-2 tensors, or quaternions transforming vectors),
      one recovers the standard conjugate (or sandwich) product:
     
$$ a' = {\bf \hat O} a {\bf \hat O}^T $$
      
For tensors, this is a direct result of tensors being composed
      of the outer products of pairs of vectors. Likewise for spinors
      composing vectors (and scalars, for that matter).

(That is, a quaternion transforms a vector via:
     
$$ (0, \vec v') = {\bf q} (0, \vec v) {\bf \bar q} $$
  
and likewise for scalars
      
$$ (s', \vec 0) = {\bf q} (s, \vec 0) {\bf \bar q} $$
     
Thus, the transformation of vectors and scalars is unified via the
      quaternion. This is because quaternions act fundamentally on spinors:
      
$$ \xi'={\bf q}\xi $$

and both vector and scalar can be written as the outer product of
 2 spinors:
      
$$ v = \xi\xi $$
   
In terms of representation theory, this can be expressed as the
      tensor product of the 2 state spinor space with itself:
    
$$ {\bf 2} \otimes {\bf 2} = {\bf 3} \oplus {\bf 1} $$
     
being the tensor sum of a vector-space and a scalar space.

But: this is a bit much for the dear user to worry about, hence,
      ALL operations are overloaded via function emulation ("__call__");
      thus, all transformation representations act like functions that
      can transform whatever kind of object you give them.
      
For example, for quaternions, one might code"


	v' = M * v
	v' = (q * v * (~q)).vector  # quat. mul does v --> (0, v).

and get the right answer, but would you rather use the same form
for each transformation:

	v' = M(v)
	v' = q(v)

Likewise for vectors and tensors (and anything else). You transform
      your angular momentum vector, l, and your inertia tensor, I,  via:

	l' = M * l
	I' = M * I * M.T
	I' = M * I * ~M

but isn't it sweeter to just use:

	l' = M(l)
	I' = M(I)

      
and let the object figure out the right answer?
      
#### Euler Angle Classes
     
There is a base class, EulerAngleBase, for transformations represented as
      Euler angles. Like matrices and versors, you can compose and call them.
      Nevertheless, there are 2 points of confusion:
      
(1) What are the axes
      
(2) What are the units.
      
The answer:

I don't know. No, really, the EulerAngle class doesn't know. It uses its
[static class attributes](https://en.wikipedia.org/wiki/Class_variable),
[in the abstract base class](https://en.wikipedia.org/wiki/Method_(computer_programming)#Abstract_methods"):
geo.metric.euler.charts._ElementalRotationTriplet.AXES
( a length 3 tuple of vectors representing ordered intrinsic
      	      rotations), to figure it out. So, to support common radar problems like platform
      motion, there is a subclass:
      

		geo.metric.euler.charts.YPR
      
which has AXES set to (z, y, x)-- so that you get a yaw rotation followed
      by a pitch rotation followed by a roll rotation; Circumference=360, so
      that if you have an ASCII motion file from a gyroscope (in degrees),
      for instance, you can do:
      
      >>>attitude = YPR(*numpy.loadtxt("gyro.txt")
  

#### Rotation Summary:
      
In the spirit of OO: you do need to know which object you have when
      performing transformation. If you have "T", then:

		T(v)

Transforms v

		 ~T

Inverts the transformation


 		T*T'

composes 2 transformations, so that

		T*(~T)

is the identity transformation (for any rotation object). Specific forms
of the rotation are as follows:

		T.dcm()

returns the equivalent direction cosine matrix (as a Tensor)

		T.versor()

returns the equivalent versor

		T.ypr()

returns the equivalent YPR triplet 

		T.rpy()

returns the equivalent RPY triplet .

### Quaternions
The full quaternion algebra, __H__, is represented in desic/hamilton.py.
      While it has no use for transformations in __E3__, it is included
      for completeness. While versors only support the Grassmann product,
      quaternions allow for addition, subtraction, and dilation operations
      as well as inner, outer, odd, and even products.

### Affine
The Affine transformation is the most general linear function of
      a vector:
    
$$ A(\vec x) = {\bf M}\vec x + \vec b $$
      
For proper euclidean transformation (det(M) = 1), than any rotating
      object will do. Dilations, skews, reflections, etc.. will require a
      fully general matrix (rep'd as a rank-2 tensor).
    
    
The magic methods support:


	A(v)                    # function emulation
	(A * A')(v) = A(A'(v))  # composition
	~A(A(v)) = v            # inversion
	A ** (m/n)              # interpolation
	
	
### MODULE FUNCTIONS 
      
First of, I believe all functions should be PURE functions: that is,
      they do not change state. Methods can change the state of an objects,
      but functions should not. This is an important practice for proper OO
      design.
      
You will notice that the euclid module is full of functions for doing the
      vector operations. They are there so that methods can call them (i.e.,
      private, but that's not enforced). You can use them if you want.
      (An example
      is doing X^Y. The super says that is X.wedge(Y), which is laborious
      summation
      over the Levi-Civita tensor, and that is not acceptable for real
      computations.
      
Hence, it is overridden to call vector._cross_product(X, Y)-- which does
      the usual antisymmetric combination over the non-zero elements of the
      generalized sum).
   
Some functions are there for public use:
      
      geo.metric.euler.charts.roll  (aliased to geo.roll)
      geo.metric.euler.charts.pitch (aliased to geo.pitch)
      geo.metric.euler.charts.yaw   (aliased to geo.yaw)
      
These function are constructor-helpers for defining rotation
      objects.
      They are bumped up to the top namespace, which is supposed to
      just contain things you need access too- but that may need a
      little code review. The return Versors (unit Quaternions).
      
(Note: the geo.metric.euclid.euclid.LinearMap.aschart() method allows
      the user to convert any SO(3) to any form via an name, instance, or
      Type).
      
### Ellipsoids
THE WGS84 elipsoid, and an almanac of other ellipsoids are available:

      geo.detic.newton.wgs84.wgs84.WGS84 
      geo.detic.newton.almanac.ALMANAC are 

## Coordinates
Coordinates are a lot like vectors, since they can be represented by 3
      frame-dependent components. Both Vector and the coordinate base classes
      inherit from the ABCNumpy interface. Hence, all coordinate
      instances can be singletons or numpy arrays, and they all have
      broadcasting, mean, "__getitem__", "__getslice__",
      "__iter__", tolist(), and
      iter() and next() behavior at hand.
      
Moreover, they're all instantiated with 3 arguments and a possible peg
      point keyword or positional argument, but I'm getting ahead of myself.
     
At this time, there are 4 (four) coordinate systems:
     


      ECEF             Earth Centered Earth Fixed Cartesian (x, y, z)

      LLH              Geodetic   (lat, lon, hgt)
      
      LTP              Tangent Plane Cartesian (x, y, z, peg_point=peg_point)
      
      SCH              offset Tangent Sphere   (s, c, h, peg_point=peg_point)
      
      
A PegPoint is just a namedtuple consisting of (lat, lon, hdg). Of course,
      all arguments are in degrees-- radians are simply not for geolocation.
      Presumably, there are a bunch of ships off the coast of west Africa, full
      of people using radians for navigation.
      
#### Transformations:
     
The idea here is to be as polymorphic as possible, so for instance, if
      you have a coordinate instance, "r":

	r.ecef()	
	r.llh()	
	r.ltp()  or  r.ltp(peg)
	r.sch()  or  r.sch(peg)

will convert r to ECEF, LLH, LTP, SCH (or LTP, SCH with a different peg).
      Moreover, non-Cartesian coordinates can be readily converted to Cartesian
      coordinates, via their "cartesian_counter_part" method, so that:

	llh.cartesian_counter_part() --> llh.ecef()
	sch.cartesian_counter_part() --> llh.ltp()

Moreover, all classes have a "vector()" method that converts the points
      into a euclid.Vector instance relative to the origin of the coordinate
      system. Thus:
      
      	      r.vector()

always works, yielding a Vector relative to the center of a coordinate
system or relative to a peg point.
     
The point is that all instances have the same methods, so you don't need
      to know what you're dealing with to get the result you want. If I have a
      motion history, r,  in LTP coordinates and a target point on the ground,
      p, I get the vector in the peg pointing from the platform to the target
      via:

		>>>d = p.vector() - r.vector()

Now if the whole thing is in SCH coordinates, the computation is:

		>>>d = p.vector() - r.vector()

(because non-Cartesian coordinates call their Cartesian_counter_part()
      transformation via vector()).

Now if you have the motion in SCH and the target in LTP, you DO NOT have
      to type-check, then the formula is 
      
	        >>>d = p.vector() - r.vector()

Clearly, a pattern is forming-- but all good thing must come to an end.
      What if you motion is in LLH (or ECEF, or a different SCH or LTP system)?
      Then:
      
	       >>>d = p.vector() - r.ltp(p.peg_point).vector()

The strategy to overload the "__sub__" operator as follows:

   	       >>>d = p - r

will make a euclid.Vector instance relative to p's origin and Cartesian
      orientation, even p is not Cartesian. Meanwhile, r can be in any
      coordinate system.

Further Note: numpy will handle all the broadcasting, you can have 1
      coordinate be a bunch of arrays, and the other be a like shaped array,
      or a singleton. It's going to work.

Of course, "\_\_add\_\_", inverts the whole thing: you can add a vector to a
      point:

		>>>p' = p + v

and p' will be in p's coordinates, with v interpreted relative to p's
      Cartesian coordinates.
      
### The Ellipsoid Problem
The ellipsoid is defined by a semi-major axis, an inverse
      flattening, and and optional model name (ellipsoid.py)
     
Its method/properties convert to many of the various metrics
      associated with an ellipsoid of revolution. It also converts
      between all the various latitude definitions.
      
COORDINATES: Here is the key fact: coordinates only work on ellipsoids.
      So, the ellipsoid has methods that construct coordinate objects.
      
Suppose you have the WGS-84 ellipsoid:
  
  	>>>wgs84 = ellipsoid.Ellipsoid(6378137.0,c298.25722293286969,cmodel="WGS84")
						 
has four methods:
      
      	  wgs84.ECEF
      	  wgs84.LLH
     	  wgs84.LTP
      	  wgs84.SCH
	  
Transformation are then done using the coordinate methods, so to get
      JPL's coordinates in ECEF:
     	    
	jpl = wgs84.LLH(34.197, -118.175, 345.).ecef()

So, this makes an LLH instance and then converts it to an ECEF instance.
      Moreover, the ellipsoid has functions for conversion:

	x, y, z = wgs84.llh2ecef(34.197, -118.175, 345.)

That is, a triplet goes in, a triplet comes out.
      
#### WGS84:
      
The wgs84.py module has the classes and functions for WGS84 as module
      constants and functions, thereby allowing an entire procedure usage.

## Platform Motion 
After all that, you still don't have platform motion. Enter the motion.py
      It requires numpy, since you will have array_like attributes.
      The SpaceCurve class is basically Vector which takes that into
      consideration. 
      (Recall the fundamental theorem of Space Curves? That curvature
      and torsion define a unique space curve? Yes, that one-- well space
      curves define all that: velocity, normal, acceleration, angular velocity,
      yadda yadda. They key property is you can define a local tangent frame,
      with:
     
      x   parallel to the curve's velocity
      y   = z ^ x
      z   is some "z" orthogonal to x. The default "z" is DOWN, but you can
      	     	      		        make it UP, or something else.
					
Hence, given a 3-DoF motion history, you get the transformation from
      level Cartesian space to the tangent frame. Now if you throw in attitude,
      represented by any kind of rotation, boom, you have an affine
      transformation to body coordinates.
     
But wait: these things all have array_like attributes, that means in
      one object, you have the transformation from a local pegged coordinate
      system the body frame AT EVERY POINT IN THE MOTION HISTORY.
      	     
Now stop and think about that for a minute. IF you were still suffering
      primitive obsession, using arrays for vectors, and using stand-alone
      functions for coordinate transformations--you be in a real pickle.
      All these arrays, which are JUST NUMBERS and have NO intrinsic meaning-
      no you the developer has to keep it straight. Then, you have to pass
      them to functions, and then to other functions-- how you deal with the
      fact that the functions are not constant--- I don't know- but you do.

None of that. You got a GPS history and an attitude history:

	  f = SpaceCurve(*gps).tld2body(imu)
	

	 f(look)

does the whole affine transformation at every point.
     
### Connections
LLH and SCH are local coordinates: their (co)tangent basis vectors depend
      on position, and we may need to know them in terms of their canonical
      global bases: ECEF and LTP, respectively. Hence: christoffel.py.
      It adds mix-ins that computer the jacobians, covariant and
      	contravariant bases
       and connection coefficients between the cannonical frames.

      
## Utilities

### Trigonometry
The proliferation of "degrees" necessitates basic trig functions
      that take and return degree arguments. They are implemented in
      trig.py. Moreover, it is in this module where you try to use numpy,
      but get python's built in math library should numpy be unavailable.

### Exceptions
Various errors specific to the package are defined in exceptions.py.
      They're basically wonky ways to catch nonsense operations.
      

## Experimental Sub Packages: __Geodesic__

### Clifford Algebra
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
      $\mathbb R^3$. It is implemented in blades.py
      
### Quaternions
Versors are limited, in that they must have unit norm, thus only the
      Grassmann product is allowed. The space is not closed under the inner,
      outer, odd, and even products as well as addition and subtraction.
      Hence, for fun, they are implemented in
      geo/desic/hamilton.py.

### Cayley-Dickson Construction
Once you have full blown quaternions, there really is no reason to
      not have octonions and sedinions, tesserines, bi-quaternions,
      co-quaternions, and their "split" varieties (including the split
      complex numbers) They can all be implemented via the Cayley Dickson
      construction in extension.py

### Minkowski Space
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



### Higher Order Tensors
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
      tensor has only 2 obvious symmetries: $ {\bf T} \pm {\bf T}^T $,
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


## Index Gymnastics
    
The 'e' method allows you to select or run on indices, as follows
      for a rank 5 tensor $ T_{ijkmn} $:

      	>>>M = T.e(0, Ellipsis, 2, 1, Ellipsis)\endverbatim

$$ M_{jn} = T_{0j21n} $$

While a contraction to a rank 3 tensor:
     
 	>>>A = T.contraction(0, 3)
   
     
$$ A_{jkm} = T_{ijkim} $$

The inner and outer products, respectively are, for example:

	>>>C = A.inner(B)
	        
means:

$$ C_{ijm} = A_{ijk}B_{km} $$

and

	>>>C = A.outer(B)

$$ C_{ijkmn} = A_{ijk}B_{mn} $$

In the later case, the rank 5 tensor is created at run-time and added
      to the Zoo. An antisymmetric product is also possible:


	>>>C = A ^ B

which becomes:
      
$$ C_{ijkln} = A_{ijk}\epsilon_{klm}B_{mn} $$
      
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

 
      
## Developer's Tutorial
This is [SOLID](http://en.wikipedia.org/wiki/SOLID_(object-oriented_design))
      OO, pythonic code. It attempts to maximize its maintainability
      score as determined by [pylint](http://www.pylint.org), the current standard for
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

## Abstract Base Classes.
That so many mathematical object share behavior has lead to base
      classes (many base classes), where behavior is implemented once,
      abstractly, and the made concrete in the leaf classes. The python
      [abstract base class](https://docs.python.org/2/library/abc.html)
      library is used extensively to signal to the developer:
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

### Polymorphism
What is it? Well having to use "print", "printf", or "sprintf"
      is NOT IT.
      Polymorphism means "do the right thing without being told."
      Hence, if you code:

	   >>>c = a * b
			   
"c" better be the product of a and b. But what are a and b?
      	 
 $$ c_i = ab_i $$
      
differs from

$$ c = a_ib_i $$

differs from

$$ c_i = a_{ij}b_j $$

differs from

$$ c_{ij} = a_{ijk}b_k $$

and so on. In these cases, internal polymoprhism--the
      highest--form works:
      The class dispatches the correct "a\.\__mul\_\_" method to
      deal with the "*" operand.

What about:

$$c_i = a_i b $$

differs from

$$ c = a_i b_j $$

differs from

$$ c_j = a_i b_{ij} $$

differs from

$$ c_{ij} = a_i b_{ijk} $$

Here, each case calls

      geo.metric.euclid.vector.Vector.__mul__, which there is none,
      it's kicked up to the super:
      geo.metric.euclid.euclid.Tensor_.__mul__. Now what?
      
So one solution is to use case-switch like _if_-_then_ blocks to
      call the right function. That's too primitive. Instead, a hash
      table bound to the left operand looks up the rank of the
      right operand to get the functions (e.g.,
      geo.metric.euclid.vector.Vector._dispatch_mul). (Aside:
      That shows that
      things are further complicated by the quaternion and the
      geo.metric.frenet_serret.SpaceCurve). That is  external
      polymorphism, where explicit code decides the function to call
      based on the argument-- this can only be solved elegantly
      with a language that supports
      [multiple dispatch](http://en.wikipedia.org/wiki/Multiple_dispatch)
      aka: multimethods.

Well, we don't have that in python, so this is what you get.
      The solution presented is an attempt to avoid branchy
      dynamic code.
      The complexity is represented in a complex __static__ data
      structure--with one point of evaluation. It takes some
      getting used to, but it the right thing to do. If you disagree-
      read the next section. 


## Tensor Symmetries Revisited: Young Tableaux
We all know the symmetries of a rank 2 tensors:
      
$$ S_{ij} = \frac 1 2 [T_{ij} + T_{ji}] $$
     
$$ S_{ij} = \frac 1 2 [T_{ij} - T_{ji}] $$      
      
are the symmetric and antisymmetric parts of the tensor, $T_{ij}$.
      They are not just (anti)symmetric under interchange of the indicies,
      they are closed under rotations. So the question is: how is this
      extended
      to higher rank tensors? What is the math?

Well, the math is a deep dive. The invariant subspaces,
      aka irreducible representations (heretofore: irreps), of a rank $ N $
      tensor are related to the symmetric group $ S_N$, which is the
      group of permutation on $N$-letter, via Schur-Weyl Duality (SWD). These
      are computed from Young Tableuax with $N$ boxes via the
      Robinson-Schensted Correspondence (RSC) which in turn are related
      to the integer partitions of $N$. The following seciotn works through
      the lowest nontrivial example, rank 3.

What are the integer paritions of $ N=3 $? While you can work these out
      by hand, the pascal.py package will do it for you:
    
       >>>p = pascal.P(3)

       >>>print p
       P(3, k)

       >>>print type(p).mro()
       [<class 'geo.metric.schur_weyl.pascal.P'>,
       <class 'geo.metric.schur_weyl.pascal._Partition'>,
       <class 'geo.metric.schur_weyl.pascal.Triangle'>,
       <type 'long'>,
       <type 'object'>]
       
Note that the
[Parition Function](http://mathworld.wolfram.com/PartitionFunctionP.html)
       is fact just an extended long
       integer.

The geo.metric.schur_weyl.pascal.P.partitions() method generates
       the paritions:
       
     
       >>>list(p.partitions())
       [[3], [2, 1], [1, 1, 1]]
       
of which there are three.

Each parition can be represented by a
geo.metric.schur_weyl.young.Diagram, respectively as follows:
     
       In [169]: for d in p.young():
     ...:     print unicode(d)
     ...:     print
     ...:     
	?????????
	??????
	???

	???
	???
	???
	
Let's look at the seconds one 
	
	d = young.Diagram(2, 1)
	
To get to the permutation group, one get a Standard Tableau
	corresponding to the diagram. A standard tableau has numbers from
	$ 1, \ldots, N$ in the boxes, with each row and each column strictly
	increasing. There are 2 ways to fill the diagram that meet those
	criteria:
	
	>>>for t in d.standard_tableaux():
	     ...:     print t
	     ...:     print
	     ...:     
	[i][j]
	[k]

	[i][k]
	[j]
	
Note that I have filled the boxes not with number, but with ordered tensor
indices _i_, _j_, and _k_.
Now we have to pick one of those tableaux:
	
	>>>t = d.fill(0, 1 ,2)

	>>>t
	Tableau_[[0, 1], [2]]

	>>>print t
	[0][1]
	[2]
	
(Note that the numeric representation starts at $ 0$, not $ f$. This
	just makes it easier to work in a computer language that starts
	its indexing at zero.)

At this point, we move into Schur-Weyl duality: what does this
	tableau have to do with the permutation group? Let start with the
	symmetric group on 3 letters:
	
	>>>S3 = monte.Sym(3)
	
It has 6 elements:
	
	>>>:for count, perm in enumerate(S3):
	     ...:     print count, repr(perm)
	     ...:     
	0 
	1 (12)
	2 (01)
	3 (012)
	4 (021)
	5 (02)
		
Each geo.metric.schur_weyl.monte.Perm permutation is represented
	by its cycle structure, which shows the orbit of an element. (Note:
	the indentity, aka the nuetral elememt, is an empty cycle) They
	can also be show in two-line-format. For the last permutation, that
	is:
	
	>>>perm
	(02)

	>>>print perm
	 | (0, 1, 2) | 
 	 | (2, 1, 0) | 
       
where the latter format shows the trajectory of each element.
       
Note that the Young tableau represents a permutation:
   
	>>>print monte.Perm.fromcycles(*t.cycles())
	 | (0, 1, 2) | 
	 | (1, 0, 2) | 

But there is more. One can also consider the set of permutations
 that don't mix the letters in different rows:

	>>>for count, perm in enumerate(t.Row()):
	...:     print count, repr(perm)
	...:     
	...:     
	0 
	1 (01)
	
The Row method yields permutations.
Likewise, one can consider the set of permutations that don't mix the
colums; however, there is a twist: when collecting these, we must keep
track of the parity of the permutation:

	>>>for count, perm in enumerate(t.Col()):
	...:     print count, repr(perm)
	...:     
	...:     
	0 (1, )
	1 (-1, (02))

The Col method yields 2-ples of (parity, permutation).
	The parity is either +1 or -1, and depends on the parity of the number of
transpositions in the permutation.

These two sets are called the Row and Column symmetrizers, respectively.
The Young Symmetrizer for the tableeua is constructed my taking the
product of these 2 sets (along with the parity):

	for count, (parity, perm) in enumerate(t.S()):
	     ...:     print count, ":",  parity, repr(perm)
	     ...:     
	0: 1 
	1: 1 (01)
	2: -1 (02)
	3: -1 (021)
	
Finally: if you apply those permutations to the tensor indices and
add the results up, you get an irreducible subspace of the tensor, which
can be show lexigraphically as follows:


	>>>t.lexigraphicS()
	'(+T.ijk +T.jik -T.kji -T.kij) / 3'


This represents a mixed symmetry iredducible subspace.

There are several way to apply that symetrizer to a tensor. We can use
the tableuax's "symmetrize" method as follows:

Start with a rank three tensor:

	>>>T = Three(*(3*arange(27.)**2))
	>>>print T.broadcast(lambda x: str(int(x)).zfill(4)) # nice up for display
	[0000, 0003, 0012] [0243, 0300, 0363] [0972, 1083, 1200]
	[0027, 0048, 0075] [0432, 0507, 0588] [1323, 1452, 1587]
	[0108, 0147, 0192] [0675, 0768, 0867] [1728, 1875, 2028]

and use the method:

	q = t.symmetrize(T)
	>>>print q.broadcast(lambda x: str(int(x)).zfill(3))  # again, nicen up
	[000, -88, -352] [088, 000, -264] [352, 264, 000]
	[000, -172, -520] [172, 000, -348] [520, 348, 000]
	[000, -256, -688] [256, 000, -432] [688, 432, 000]

OK. So how do we know that that is rotationally closed? Let's rotate it
by an arbitrary axis:

	R = roll(34)*pitch(12)*yaw(-56)
	qprime = R(q)

Note that the versor $R$ knows exactly how to rotate a rank 3 tensor--you,
the user do not have be concerned with all the index gymnastics.

	>>>print qprime.broadcast(lambda x: round(x, 3))
	[0.0, 222.464, -877.141] [-222.464, 0.0, 281.535] [877.141, -281.535, 0.0]
	[-0.0, -86.987, 352.2] [86.987, 0.0, -113.674] [-352.2, 113.674, -0.0]
	[0.0, 70.665, -516.829] [-70.665, -0.0, 182.12] [516.829, -182.12, -0.0]

So I guess that's closed? The zeros are in the same spot, and, e.g.
$q'_{yzx} = -q'_{xzy}$, but I am not convinced. Let's take the tensor
$q'$ and ask it to show us its non-zero modules:

	>>>qprime.show_modules()
	[i][j]
	[k] 
	w=1641.91108164 

which matches the orginal:

	>>>q.show_modules()
	[i][j]
	[k] 
	w=1641.91108164 


Note that that procedure was deepy complicated. Neverthelss, it reduces
to the simple stuff that we understand from lower rank tensors:

Scalars are less than trivial, they have ther permutation group on nothing:
there
is simply no index to permute. Meanwhile, vectors have one index, and
the permutation group on 1 element is trivial.

$N=2$ has two paritiotns ( $2=2$ and $2=1+1$),
corresponding to the following diagrams:

	>>>p = pascal.P(2)
	>>>for d in p.young():
	     ...:     print unicode(d)
	     ...:     print
	??????
	???
	???

It should be easy to see that the first (last) diagram's Young symmetrizer
is total (anti)symmetric, corresponding to the familiar symmetric and
antisymetric tensors.

## Dimensionality

But wait, there is more. The Young diagrams have a remarkable formula,
called the "Hook Length Formula". The hook length is computed from
the arm and leg lengths. See:

	geo.metric.schur_weyl.young.Diagram.arm()
	geo.metric.schur_weyl.young.Diagram.leg()
	geo.metric.schur_weyl.young.Diagram.hook()

for more. These formula can calculate the dimensions of the closed
subspaces, for any rank of tensor, over _any_ field.

Let's start with the basis diagram:

	>>>d = Diagram(1)
	>>>print unicode(d)
	???


It is just a box. We can assign it any dimension and meaning we want:

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

	??????
	3	6	10

	???
	???
	1	3	6

So for combinig 2 quantum spinors ( $D=2$ ), we get a symmetric triplet and
an antisymmetric singlet, ala:

$$ |1, 1\rangle = \uparrow\uparrow $$

Meanwhile $D=3$ tells us:

$$ {\bf 3} \otimes {\bf 3} = {\bf 6}_S \oplus {\bf 3}_A $$

which says that (anti)symmetric rank 2 tensors have 6 (3) components in three
dimensions. In special relativity:

$$ {\bf 4} \otimes {\bf 4} = {\bf 10}_S \oplus {\bf 6}_A $$

so that the antisymmetric electromagnetic field strength tensor,
$ F_{\mu\nu} = \partial_{\mu}A_{\nu} - \partial_{nu}A_{\mu} $
has six componets (3 for the Electric field, and 3 for the magnetic field).
The symmetric stress-energy tensor has 10 components (1 for mass/energy
density, 3 for energy flux or momentum time derivative, and 6 for
spatial stress).

Now we can take it to 3rd rank:

	In [378]: for item in d ^ d ^ d:
	     ...:     print unicode(item)
	     ...:     print "{}\t{}\t{}".format(*map(item.dimW, [2, 3]))
	     ...:     print

	?????????
	4	10

	??????
	???
	2	8

	??????
	???
	2	8

	???
	???
	???
	0	1

From the $D=2$ columns we learn: there is no antisymmetric combination
of 3 electrons. There is one symmetric combination of 4 dimensions
$J = 3/2 $ ), and 2 mixed symmetry $J=1/2$ combinations.

For $D=3$ , there is one antisymmetric tensor, the Levi-Civita tensor.
While in the quark model, this irrep is the $\Omega^-$ baryon,
famoulsy predicted by Gell-Mann.

The application of Young tableau and symmetry spans a phenomenal amount
of physics.

### Quantum Spins

#S pin 1/2:

The addition of quantum spins is a special case of __SU__(2) irreducible
representations. For instance, a spin 1/2 particle has 2 eigenstates of
"alignment": the spin along and arbitary axis (taken to $J_z$)
has 2 eigenvalues:
$m = \pm frac 1 2$ . When the spins of 2 indentical particales are combined,
the eigenstates of total $J_z$ are not eigenstates the individual particles'
$J_z$ .

	>>>u = racah.Ket(0.5, 0.5)
	>>>print unicode(u)
	|??, ?????
 
	>>>print unicode(d)
	|??, -?????

	>>>d = racah.Ket(0.5, -0.5)

	>>>print unicode(u*u)
	[ |1, 1???]

	>>>print unicode(u*d)
	[????????) |0, 0??? + ????????) |1, 0???]

	>>>print unicode(d*u)
	[-????????) |0, 0??? + ????????) |1, 0???]

	>>>print unicode(d*d)
	[ |1, -1???]

So that now, by inspection, we can construct linear combinations that a
(normalized) eigenstates of $J^2$ and $J_z=0$ :

	>>>print (u*d + d*u)/sqrt(2)
	[ |1, 0>]

	>>>print (u*d - d*u)/sqrt(2)
	[ |0, 0>]


### Spin 1

The racah.py has the basis kets for spin 1:

	>>>Z, P, M = [Ket(1, m) for m in (0, 1, -1)]

so that:

	>>>print racah.Z, racah.P, racah.M
	 |1, 0>  |1, 1>  |1, -1>

Then we can take linear combination:

### Spherical Vectors:

With the above tools we can now delve in to spherical vectors.
Mathematically, the basis vectors transform as the fundemental representation.
[To be written]

