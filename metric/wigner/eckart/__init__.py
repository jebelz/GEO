"""The wigner/eckart package builds full-blown 3**N DoF rank-N tensors
in the Spherical Basis: the complex basis of rotation eigenstates.

The break down differs significantly from the Cartesian method of
building on the 3**N combinations of N-length outer products of unit vectors.

Rather, the basis states are identified by the Casimir invariant, or weight:
j, which is the usual angular momentum j meaning:

j  meaning
-----------
0  Scalar
1  Vector
2  Symmetric trace-free tensor
3  Symmetric trace-free rank-3 tensor
4  Symmetric trace-free rank-4 tensor

Of course, each j admits 2j+1 orders, which are rotational eigenstates (about
the z-axis with eigenvalue exp(im phi), for -j <= m <= j.

Beyond rank 3, there are also degeneracies, that is, j may not be unique.
In this case, there is a seniority index, q, where, there appears to be
no explicit convention for ordering, geo will order the smallest q with
state of highest symmetry (in the sense of Young tableaux-- see schur_weyl/).

On notation, a "j" weight is labeled by its multiplicity, and the expansion
from Cartesian (e.g. 3 x 3 for rank-2 tensors) into 5 + 3 +1 is standard
representation theory (group theory) where a reducible product representation
is expressed as a sum of irreducuble representations.

This is 100% entirely equivlent to expressing the nice neat Cartesian
polynomial (on the unit sphere):

(x**a)*(y**b)**(z**c)

as some linear combination of spherical harmonics:

Y_lm,  Y_l'm', ...,  Y_00

with the maximum l = a + b + c
(Note, l is j, and in general l is a non-negative integer-- same for j, but it
can also be 1/2 integer. Since we're dealing with tensors, and not spinors,
j is always integer, there are spinoral-harmonics for 1/2 integer j, but
we're not going there).

Word to the wise: it is not easy for large j. Here's the usual break down:

j  product           sum                     comments
-------------------------------------------------------------------------
0  1                  1               Scalar is the Trivial Rep.
1  3                  3               Fundmental Rep vs. Adjoint Rep.
2  3 x 3              5 + 3 + 1       Sym, Trace Free; Antisymmetrc; Trace.
3  3 x 3 x 3          7 + 5 + 5 + 3 + 3 + 3 + 3 + 3 + 1 + 1 + 1
4  3 x 3 x 3 x 3      9 + 7 + 7 + 7 + 5 + 5 + 5 + 5 + 5 + 5 + 3 + 3 + 3...
.                     + 3 + 3 + 3 + 1 + 1 + 1

Note that for rank-N, the highest weight is j=N. It has 2N+1 components are
is the symmetric trace-free tensor of rank-N. It's called the Natural
Form of the rank-N tensor, and can be contructed from the meta-class:

natural_form(0, 1, ..., j, -j, -j+1, ..., -1)

See the detracer() tensor method.

Now for implementations:

The classes:

Scalar
Vector
Tensor
Three
Four

are indeed the spherical tensors of obvious rank. Their construtor signatures
are not straight-forward: this is essential complexity.

They are built from component classes in component modules-- each module
has a different responsibility for its Scalar, Vector, ..., Four classes.
The break down is as follows:

ec.py   Python container behavior, keeping in mind degerneracies.
ka.py   Do algebraic operations, keeping in mind degenracies
rt.py   Handle conversion code to/from Cartesian representations[1].

The problem is complex enough that an architecture with such apparent non-
flatness is warranted. (Do not confuse apparent flatness with true flatness.)


[1] The actual formulae for conversion to/from natural form Cartesian
to the heighest weight irrep are tabulated in modules:

Module         Thingy      Comments
-------------------------------------------------------------------------
zero.py        Scalar     Trivial
one.py         Vector     This the basic Cartesian -> Spherical
                                or       Adjoint -> Fundamental Rep.
two.py         Tensor     Extract |j=2, m> for symmetric trace-free tensor
three.py       Rank-3     Rank 3 is where things get complicated
four.py        Rank-4     NotImplemented.

Lower weight Irreps should reuse code in lesser modules.
"""
## \namespace geo.metric.wigner.eckart All together: ec.py + ka.py + rt.py
from . import ec
from . import ka
from . import rt

from . import zero, one, two, three, four

## Class-decorator adds a degeneracy (class)method from pascal.Degen; since
# the [n, j=0, ... n] for a Number Triangle, and each has, of course a
# fixed n.
# \param cls A spherical tensor class (with a maximum weight cls.L)
# \returns cls with a degeneracy classmethod added.
# \sideeffect Adds "degeneracy" classmethod to cls.
def weights(cls):
    """@weights
    class Class(....:

    adds the degeneracy method to a class, which is a function of cls.L
    """
    from ...schur_weyl.pascal import Degen
    # This becomes a class method for computing degeneracy of weight j.
    cls.degeneracy = Degen(cls.L)
    # Add the correct natural form projection hash table
    cls.projectors = getattr(
        {0: zero, 1: one, 2: two, 3: three, 4: four}[cls.L],
        'PROJECTIONS')
    return cls


## \f$ {\bf 1} = {\bf 1} \f$
@weights
class Scalar(rt.Scalar, ka.Scalar, ec.Scalar):
    """s = Scalar(|0, 0>)"""


## \f$ {\bf 3} = {\bf 3} \f$
@weights
class Vector(rt.Vector, ka.Vector, ec.Vector):
    """v = Vector(|1, 0>, |1, +1>, |1, -1>)"""


## \f$ {\bf 3} \otimes {\bf 3} = {\bf 5} \oplus {\bf 3} \oplus {\bf 1} \f$
@weights
class Tensor(rt.Tensor, ka.Tensor, ec.Tensor):
    """t = Tensor(|0, 0>,
    ...           [|1, 0>, |1, +1>, |1, -1>],
    ...           [|2, 0>, |2, +1>, |2, +2>, |2, -2>, |2, -2>])
    """


## \f$ {\bf 3} \otimes {\bf 3} \otimes {\bf 3} = {\bf 7} \oplus 2\cdot{\bf 5}
# \oplus 3\cdot{\bf 3} \oplus {\bf 1} \f$
@weights
class Three(rt.Three, ka.Three, ec.Three):
    """
    t = Three([(|0, 0>,)],
    ...       [(|1, 0>, |1, 1>, |1, -1>)_0,
    ...        (|1, 0>, |1, 1>, |1, -1>)_1,
    ...        (|1, 0>, |1, 1>, |1, -1>)_2],
    ...       [(||2, 0>, |2, 1>, |2, 2>, |2, -2>, |2, -1>)_0,
    ...         (||2, 0>, |2, 1>, |2, 2>, |2, -2>, |2, -1>)_1],
    ...       [(||3, 0>, |3, 1>, |3, 2>, |3, 3>, |3, -3>, |3, -2>, |3, -1>)])
    """

    

## \f$ {\bf 3} \otimes {\bf 3} \otimes {\bf 3} \otimes {\bf 3} =
# {\bf 9} \oplus 3\cdot{\bf 7} \oplus 6\cdot{\bf 5}
# \oplus 6\cdot{\bf 3} \oplus 3\cdot{\bf 1} \f$
@weights
class Four(rt.Four, ka.Four, ec.Four):
    """f = Four(
    [(0,), (0,), (0,)],
    [(0, 1, -1), (0, 1, -1), (0, 1, -1), (0, 1, -1), (0, 1, -1), (0, 1, -1)],
    [(0, 1, 2, -2, -1),
     (0, 1, 2, -2, -1),
     (0, 1, 2, -2, -1),
     (0, 1, 2, -2, -1),
     (0, 1, 2, -2, -1),
     (0, 1, 2, -2, -1)],
    [(0, 1, 2, 3, -3, -2, -1),
     (0, 1, 2, 3, -3, -2, -1),
     (0, 1, 2, 3, -3, -2, -1)],
    [(0, 1, 2, 3, 4, -4, -3, -2, -1)])

    That is:
    position   list of...
    ---------------------
    0      [[s0], [s1], [s2]]        three length 1 lists of scalars
    1      [v0, v1, v2, v3, v4, v5]  five length 3 vectors
    2      [t0, t1, t2, t3, t4, t5]  five length 5 natural form rank-2 tensors
    3      [T0, T1, T2]              three length 7 natural form rank-3 tensors
    4      [[F]]                     One length 9 natural form rank-4 tensor
    """

    

## Natural Form Rank N Tensor Factory from 2N+1 arguments.
class natural_form(type):
    """spherical_tesnor = natural_form(*args)

    For len(args) == 2N+1, you get a natural form (symmetric in all indices
    and trace-free in all indices), stuffed into a full rank-N spherical
    tensor argument from the ZOO.
    """

    ## Construct a Natural Form Tensor from arguments.
    # \param args 2N+1 values
    # \returns tensor (rank-N, spherical, with only |j=N> non-zero
    def __new__(mcs, *args):
        j, zero = divmod(len(args)-1, 2)
        if zero:
            from ....utils.exceptions import RepresentationError as err
            raise err("Got {} arguments.".format(len(args)))
        try:
            cls = ZOO[j]
        except KeyError:
            raise ValueError(
                "Can't pack: \n {}\n into a Natural Tensor".format(args))
        if j <= 1:  # special case of singular j
            self = cls(*args)
        elif j == 2:  # special case of singular q
            self = cls(0, [0]*3, list(args))
        else:  # general case of multiple j, q.
            data = [[tuple(0.*m for m in q) for q in j] for j in cls.jqm()]
            data[-1][-1] = args
            self = cls(*data)
        return self

## The Spherical Zoo
ZOO = {cls.L: cls for cls in (Scalar, Vector, Tensor, Three, Four)}


## spherical tensor factory
class spherical(type):
    """T = spherical(*args)

    will pack 3**N args is a rank N tensor, or 2J+1 args in a
    natural form tensor.
    """
    def __new__(mcs, *args):
        for l, cls in ZOO.iteritems():
            if 3 ** l == len(args):
                return cls._pack(*args)
        return natural_form(*args)


## Spherical Vector Basis: VECTOR[m] --> \f$|1, m\rangle \equiv
# {\bf \hat{e}}^m\f$
VECTOR = Vector.basis()[1][0]

## Rank 2 Basis:  TENSOR[j][m] --> \f$|j, m\rangle \equiv \hat{Y}_j^m\f$
TENSOR = [[m for m in item[0]] for item in Tensor.basis()]

## Rank 3 Basis THREE[j][q][m]
THREE = Three.basis()


## Rank 4 Basis FOUR[j][q][m]
FOUR = Four.basis()
