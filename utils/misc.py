"""Linear Algebra Accelerators"""
## \namespace geo.utils.misc Linear Algebra Helpers
import numpy as np


## Use lapck_lite to faster matrix inverse for fixed dimensions
class SquareMatrixInverter(object):
    """n -->matrix size
    callable with n x n numpy matrix argument
    """

    ## n x n square matrix size n
    #  \param n Square matrix (root) size.
    def __init__(self, n):
        ## n x n matrix size n
        self.n = n
        ## <a href="http://mathworld.wolfram.com/Pivoting.html">pivots</a>
        self.pivots = np.zeros(n, np.intc)
        ## <a href="http://en.wikipedia.org/wiki/Identity_matrix">
        #  \f$ I_n \f$ </a>
        self.eye = np.eye(n, dtype=np.float64)

    ## Call takes a memory contiguous array
    #  \param a a specialized contiguous memory array
    #  \retval result LAPACK.dgesv
    def fast__call__(self, a):
        """smi(a) -->a**(-1)"""
        #pylint: disable=W0612
        from numpy.linalg import lapack_lite as LAPACK
        b = np.copy(self.eye)
        try:
            func = LAPACK.dgesv
        except AttributeError:
            raise RuntimeError("Upgrayyedd numpy version.")
        result = func(self.n,
                      self.n,
                      a.astype(np.float64),
                      self.n,
                      self.pivots,
                      b,
                      self.n,
                      0)

        ## CHECK FOR LINALG ERROR  HERE
        return b

    #pylint: disable=R0201
    ##  Nominal matrix inversion
    # \param a An array
    # \returns  \f$ a^{-1} \f$
    def nominal__call__(self, a):
        """Executes standard linalg interface"""
        from numpy import linalg
        return linalg.inv(a.astype(np.float64))

    ## set call method
    __call__ = fast__call__


## Invert a 6 x 6 matrix
INV_6X6_FAST = SquareMatrixInverter(6)

## Invert a 2 x 2 matrix - this is not as fast as an algebraic solution
INV_2X2_FAST = SquareMatrixInverter(2)


#pylint: disable=R0914,C0103
## Compute 2 x 2 eigenvalues of an array of 2 x 2's-- \n
# straight old school procedural.
# \param A (n, 2, 2,) array of (2, 2) to diagonalize
# \retval vals  array of (n, 2) eigenvalues
# \retval vecs  array of (n, 2, 2) eigenvectors
# \throws RuntimeWarning invalid value encountered in divide
def EIG2X2(A):
    """A direct computation of:

    return np.array(map(linalg.eigvals, A))

    for a speed bump of 60x
    """
    a = A[Ellipsis, 0, 0]
    b = A[Ellipsis, 0, 1]
    c = A[Ellipsis, 1, 0]
    d = A[Ellipsis, 1, 1]

    D = a*d - b*c
    T = a + d

    R2 = (T/2)**2 - D

    R = np.sqrt(R2)

    lam1 = T/2 + R
    lam2 = T/2 - R

    vals = np.array((lam1, lam2)).T

    v1 = np.array([lam1-d, c])
    v2 = np.array([lam2-d, c])

    v1 /= (v1**2).sum(axis=0)**0.5
    v2 /= (v2**2).sum(axis=0)**0.5

    vecs = np.array([[v1[0], v1[1]],
                     [v2[0], v2[1]]]).T

    ## Find c=0 bum solutions
    if len(c) != len(c.nonzero()[0]):  # algorithm degeneracy
        idx = np.where(c == 0)[0]

        u1 = np.array([b[idx], lam1[idx] - a[idx]])
        u2 = np.array([b[idx], lam2[idx] - a[idx]])

        u1 /= (u1**2).sum(axis=0)**0.5
        u2 /= (u2**2).sum(axis=0)**0.5

        vecs[idx, Ellipsis, :] = np.array([[u1[0], u1[1]],
                                           [u2[0], u2[1]]]).T

        ## Find c=0, b=0 bum solutions
        if not all(b[idx]):
            idx = np.where((c == 0) * (b == 0))
            vecs[idx, Ellipsis, :] = np.array([[1, 0],
                                               [0, 1]])

    return vals, vecs


__all__ = ('INV_6X6_FAST', 'INV_2X2_FAST', 'EIG2X2')
