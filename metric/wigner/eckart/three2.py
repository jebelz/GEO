"""This module converts rank-3 Spherical Tensors into Cartesian
Tensors. While rank-1 and 2 are simple enough to allow inversion of
the "from" algebraic relations, inverting 27 equations is a bit much.

Instead, the Cartesian forms of the Eigenstates are computed directly
(i.e, by fiat) from the spherical basis vectors and Racah algebra."""
## \namespace geo.metric.wigner.eckart.three2 Rank-3 Cartesian forms of
# spherical eigenstates.
from . import VECTOR
from geo import EPSILON, DELTA

# \f$ {\bf \hat{e}}^q} \ \forall \ q \in (0, 1, -1) \f$ in Cartesian
# coordinates.
Z, P, M = [item.tocart() for item in VECTOR]



T33 = P & P & P
T3m3= M & M & M

T32 = (P & P & Z).natural_form() / 3
T3m2= (P & P & Z).natural_form() / 3

T31a = ((P & P & M)  + 0*(P & Z & Z)).natural_form() 
T3m1a= ((M & M & P)  + 0*(M & Z & Z)).natural_form()

T31b = (0*(P & P & M)  + (P & Z & Z)).natural_form() 
T3m1b= (0*(M & M & P)  + (M & Z & Z)).natural_form()


T31 = T31a + T31b  # TODO: weight
T3m1 = T3m1a + T3m1b  # TODO: weight

T30 = -(Z & Z & Z).natural_form() * (251. /81.) ** 0.5


T11_0 = (DELTA & VECTOR[1].tocart()).natural_form()
T1m1_0 = (DELTA & VECTOR[-1].tocart()).natural_form()
T10_0 = (DELTA & VECTOR[0].tocart()).natural_form()


T11_1 = NotImplemented
T1m1_1 = NotImplemented
T10_1 = NotImplemented


T11_2 = NotImplemented
T1m1_2 = NotImplemented
T10_2 = NotImplemented


T22_0 = NotImplemented
T2m2_0 = NotImplemented
T21_0 = NotImplemented
T2m1_0 = NotImplemented
T20_0 = NotImplemented

T22_1 = NotImplemented
T2m2_1 = NotImplemented
T21_1 = NotImplemented
T2m1_1 = NotImplemented
T20_1 = NotImplemented


## Dictionary of Normalized Cartesian spherical
# eigenstates: \f$ T_{j}^{(m; q)} \f$
THREE = {(3, 0, 3): T33,
         (3, 0, 2): T32,
         (3, 0, 1): T32,
         (3, 0, 0): T30,
         (3, 0, -1): T3m1,
         (3, 0, -2): T3m2,
         (3, 0, -3): T3m3,
         (2, 0, 2): T22_0,
         (2, 0, 1): T22_0,
         (2, 0, 0): T20_0,
         (2, 0, -1): T2m1_0,
         (2, 0, -2): T2m2_0,
         (2, 1, 2): T22_1,
         (2, 1, 1): T22_1,
         (2, 1, 0): T20_1,
         (2, 1, -1): T2m1_1,
         (2, 1, -2): T2m2_1,
         (1, 0, 1): T21_0,
         (1, 0, 0): T21_0,
         (1, 0, -1): T1m1_0,
         (1, 1, 1): T11_1,
         (1, 1, 0): T10_1,
         (1, 1, -1): T1m1_1,
         (1, 2, 1): T11_2,
         (1, 2, 0): T10_2,
         (1, 2, -1): T1m1_2,
         (0, 0, 0): EPSILON}
