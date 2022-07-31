"""Rank-0 (Trival) Rep Conversion

Spherical vs. Cartesian scalars.

Weyl Module  dimW Irreps (j; q)         
Tableau                                 
=====================================================================
[]            1     1    (0; 0)                 s = s"""
## \namespace geo.metric.wigner.eckart.zero Trival Reps of SO(3).


## S0
# \param s A Cartesian Scalar
# \returns \f$ s^0_0 = s \f$
def S0(s):
    return s.w

## The one projection onto |0, 0>
PROJECTIONS = {0: S0}


