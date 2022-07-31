"""The cauchy module could very-well be called hilbert.py, but that's for
use only with quantum mechanic a relativity- for me. So here we have functions
that are relevant for Cauchy spaces, by which I mean an inner product space.

The implementation is such that all objects bring their own dot product to the
show, and these module functions take it from there, e.g:

norm
hat
scalar_projection
projection
rejection
reflection

We all know what that means for vectors; the idea is that these functions
also work for rank N tensors, quaternions (not versors, though, they only
rotate), 4-vectors, and we will see what else.

For Vectors:
These functions work for Real Euclidean (vector.Vector) and Complex
(gibbs.Gibbs) vectors.


"""
## \namespace geo.utils.cauchy Generalized
# <a href="http://en.wikipedia.org/wiki/Inner_product_space">Cauchy inner
# product space</a> functions.


## Dot is responsible invoking the correct inner product (Vector vs. Gibbs),
# for higher rank, dot contracts on all indices.
# \returns Scalar \f$ \sqrt{\langle a, a \rangle} \f$
def norm(a):
    """scalar = norm(a)"""
    return a.dot(a) ** 0.5


## Normalizer
# \returns \f$ \frac{\bf a}{||a||} \f$
def hat(a):
    """z --> z / ||z||"""
    #    print "DD: cauchy.hat:", type(a)
    return a / norm(a)


## Scalar Projection <a, b> / ||b||**2
# \return \f$ \langle {\bf a}, \hat{\bf b} \rangle \f$
def scalar_projection(a, b):
    """a * b / ||b||"""
    # reversed order is b/c complex space defines <a|b> = a.C * b
    return hat(b).dot(a)


## Projection of a >> b.
# \returns \f$ {\rm proj}_a(b) = \langle a, \hat{b} \rangle  \hat{b}\f$
def projection(a, b):
    """a * b / ||b||**2"""
    # (bb*)a
    return scalar_projection(a, b) * hat(b)


## Rejection of a << b.
# \returns \f$ {\rm rej}_a(b) \equiv {\bf a} - {\rm proj}_a(b) \f$
def rejection(a, b):
    """a - [a * b / ||b||**2]"""
    return a - projection(a, b)


def rejection_vector_only(a, b):
    c = b.hat()
    return -(a ^ c.C) ^ c


## Reflection of a | b.
# \returns \f$ {\rm ref}_a(b) \equiv {\rm rej}_a(b) - {\rm proj}_a(b) \f$
def reflection(a, b):
    """a - 2[a * b / ||b||**2]"""
    return projection(a, b) - rejection(a, b)
