"""Tucker decomposition of higher-order objects. It's clifford related,
but beyond scope.

Doing this right would require a more general clifford algebra implementation
followed by Spin(n) and Pin(n) groups.
"""
## \namespace geo.desic.tucker HOSVD Place Holder
# @image html tucker.tiff
# @image latex tucker.tiff


## A popular thing in HOSVD.
def fourth_order_cumulant(vector_):
    """Fourth order Cummulant of a random vector (array)."""
    a = vector_.centralized()
    b = a.C
    ee = (a & b).mean() & (b & a).mean()
    return (
        (a & b & b & a).mean() -
        ee -
        ee.transpose(0, 2, 1, 3) -
        ee.transpose(0, 3, 1, 2)
        )  # maybe.


class EMP(object):
    """Elementary Multilinear Projection."""


class TTP(EMP):
    """Tensor-Tensor Projection."""


class TVP(EMP):
    """Tensor-Vector Projection."""
