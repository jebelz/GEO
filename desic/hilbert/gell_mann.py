"""This is an extemporaneous module with no unifying idea"""

#pylint: disable=C0301,C0103

## \namespace geo.desic.hilbert.gell_mann Eightfold Way and/or
# <a href="http://en.wikipedia.org/wiki/Gell-Mann_matrices">
# __SU(3)__ </a>
import numpy as np
from ...metric.euclid import tensor


# prolly need to attach this to a metaclass to get it to work.
axes = ('r', 'g', 'b')


## Lie Bracket \n
# [A, B] = \f$ \frac{1}{2i}(AB-BA) \f$
def Lie(a, b):
    """Lie Bracket."""
    return -(1j) * (a*b - b*a) / 2


## \f$ \frac{\langle A|B \rangle}
# {\sqrt{\langle A|A \rangle \langle B|B \rangle}} \f$
def proj(a, b):
    """projection."""
    return a.dot(b.C).w / (abs(a).w*abs(b).w or 1)


## An SU(3) Matrix as a complexified tensor (in development).
class SU3(tensor.Tensor):
    """SU3 is a complex 3 x 3 with a Lie product bound to xor."""

    def __xor__(self, other):
        return Lie(self, other)


## \f$ i = \sqrt{-1} \f$
i = 1j


##\f$\lambda_1=\left(\begin{array}{ccc}0&1&0\\1&0&0\\0&0&0\end{array}\right)\f$
Lambda1 = SU3(0, 1, 0,
              1, 0, 0,
              0, 0, 0)


##\f$\lambda_2=\left(\begin{array}{ccc}0&-i&0\\i&0&0\\0&0&0\end{array}\right)\f$
Lambda2 = SU3(0, -i, 0,
              i, 0, 0,
              0, 0, 0)


##\f$\lambda_3=\left(\begin{array}{ccc}1&0&0\\0&-1&0\\0&0&0\end{array}\right)\f$
Lambda3 = SU3(1, 0, 0,
              0, -1, 0,
              0, 0, 0)


##\f$\lambda_4=\left(\begin{array}{ccc}0&0&1\\0&0&0\\1&0&0\end{array}\right)\f$
Lambda4 = SU3(0, 0, 1,
              0, 0, 0,
              1, 0, 0)


##\f$\lambda_5=\left(\begin{array}{ccc}0&0&-i\\0&0&0\\i&0&0\end{array}\right)\f$
Lambda5 = SU3(0, 0, -i,
              0, 0, 0,
              i, 0, 0)


##\f$\lambda_6=\left(\begin{array}{ccc}0&0&0\\0&0&1\\0&1&0\end{array}\right)\f$
Lambda6 = SU3(0, 0, 0,
              0, 0, 1,
              0, 1, 0)


##\f$\lambda_7=\left(\begin{array}{ccc}0&0&0\\0&0&i\\0&i&0\end{array}\right)\f$
Lambda7 = SU3(0, 0, 0,
              0, 0, -i,
              0, i, 0)


##\f$\lambda_8=\frac{1}{\sqrt{3}}\left(\begin{array}{ccc}1&0&0\\0&0&0\\0&0&-2
# \end{array}\right)\f$
Lambda8 = SU3(1, 0, 0,
              0, 0, 0,
              0, 0, -2) / 3**0.5


## \f$ \lambda_i \f$
LAMBDA = (tensor.DELTA,
          Lambda1,
          Lambda2,
          Lambda3,
          Lambda4,
          Lambda5,
          Lambda6,
          Lambda7,
          Lambda8)

## \f$ T^k_{ij} = \langle [\lambda_i, \lambda_j] | \lambda_k \rangle \f$ \n
# The 8 Dimensional Adjoint Representation - in 1 line: viva python.
EIGHTFOLD_WAY = np.array([[[proj(l1 ^ l2, l3).real for l1 in LAMBDA[1:]] for l2 in LAMBDA[1:]] for l3 in LAMBDA[1:]])
