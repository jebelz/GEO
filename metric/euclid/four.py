"""Rank 4 Tensors.
"""
#pylint: disable=C0301, R0902, R0913, R0904
## \namespace geo::metric::euclid::four Rank Four cartesian
# tensors.

from . import euclid
from .three import HigherOrder

__all__ = ('Four',)


## Rank 4 tensor, Transforms as
# \f$ T_{ijkl}' = M_{ii'}M_{jj'}M_{kk'}M_{ll'}
# T_{i'j'k'l'} \f$
class Four(HigherOrder):
    """Four(xxxx, xxxy, ..., zzzz)"""

    __metaclass__ = euclid.ranked

    rank = 4

    ## Transpose and descend rank to make a nice string, w/ borders.
    def __str__(self):
        """somethig like this:
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
        """
        # value
        v = "\n%s\n".join(map(str, map(self.e, xrange(3))))
        l = max(map(len, v.split('\n')))   # max length
        b = "=" * l                        # border
        s = "-" * l                        # separator
        return "%s\n%s\n%s" % (b, v % (s, s), b)

    ## Invert
    # \returns \f$ T_{ijmn}:{\bar{T}}_{mnkl} = \delta_{ij}\delta_{kl} \f$
    def __invert__(self):
        #
        return type(self)(*inverse(self))

    ## Right minor symmetry:
    # \return \f$ C_{ij\{kl\}} \equiv \frac{1}{2}(C_{ijkl}+C_{ijlk}) \f$
    def Sright(self):
        """T.ij{kl}"""
        return 0.5 * (self + self.transpose(0, 1, 3, 2))

    ## Left minor symmetry:
    # \return \f$ C_{\{ij\}kl} \equiv \frac{1}{2}(C_{ijkl}+C_{jikl}) \f$
    def Sleft(self):
        """T.{ij}kl"""
        return 0.5 * (self + self.transpose(1, 0, 2, 3))

    ## Major symmetry:
    # \return \f$ \frac{1}{2}(C_{ijkl}+C_{klij}) \f$
    def Smajor(self):
        """T.{(ij)(kl)}"""
        return 0.5 * (self + self.transpose(2, 3, 0, 1))

    ## Right minor (anti)symmetry:
    # \return \f$ C_{ij[kl]} \equiv \frac{1}{2}(C_{ijkl}-C_{ijlk}) \f$
    def Aright(self):
        """T.ij[kl]"""
        return 0.5 * (self - self.transpose(0, 1, 3, 2))

    ## Left minor (anti)symmetry:
    # \return \f$ C_{[ij]kl} \equiv \frac{1}{2}(C_{ijkl}-C_{jikl}) \f$
    def Aleft(self):
        """T.[ij]kl"""
        return 0.5 * (self - self.transpose(1, 0, 2, 3))

    ## Major (anti)symmetry:
    # \return \f$ \frac{1}{2}(C_{ijkl}-C_{klij}) \f$
    def Amajor(self):
        """T.[(ij)(kl)]"""
        return 0.5 * (self - self.transpose(2, 3, 0, 1))

    ## Extract 21 components that are left/right-minor and major
    # symmetric.
    # \returns \f$ T_{\{\{ij\}\{kl\}\}} \f&
    def elastic(self):
        """Get components that are valid elasticity tensor:
        T.{(ij)(kl)}.{ij}kl.ij{kl}"""
        return self.Smajor().Sleft().Sright()

    ## The Natural Form (Symmetric and Trace Free in All indices)
    # \returns \f$
    # T^{(4,1)}_{ijkl} = {\bf T^{(S)}}_{ijkl} -
    # \frac{1}{84}[
    # \delta_{ij}(T_{mmkl}+T_{mmlk}+T_{mkml}+T_{mlmk}+T_{mklm}+T_{mlkm}+
    # T_{kmml}+T_{lmmk}+T_{kmlm}+T_{lmkm}+T_{klmm}+T_{lkmm}) +
    # \delta_{ik}(T_{mmjl}+T_{mmlj}+T_{mjml}+T_{mlmj}+T_{mjlm}+T_{mljm}+
    # T_{jmml}+T_{lmmj}+T_{jmlm}+T_{lmjm}+T_{jlmm}+T_{ljmm}) +
    # \delta_{il}(T_{mmkj}+T_{mmjk}+T_{mkmj}+T_{mjmk}+T_{mkjm}+T_{mjkm}+
    # T_{kmmj}+T_{jmmk}+T_{kmjm}+T_{jmkm}+T_{kjmm}+T_{jkmm}) +
    # \delta_{jk}(T_{mmil}+T_{mmli}+T_{miml}+T_{mlmi}+T_{milm}+T_{mlim}+
    # T_{imml}+T_{lmmi}+T_{imlm}+T_{lmim}+T_{ilmm}+T_{limm}) +
    # \delta_{jl}(T_{mmki}+T_{mmik}+T_{mkmi}+T_{mimk}+T_{mkim}+T_{mikm}+
    # T_{kmmi}+T_{immk}+T_{kmim}+T_{imkm}+T_{kimm}+T_{ikmm}) +
    # \delta_{kl}(T_{mmij}+T_{mmji}+T_{mimj}+T_{mjmi}+T_{mijm}+T_{mjim}+
    # T_{immj}+T_{jmmi}+T_{imjm}+T_{jmim}+T_{ijmm}+T_{jimm})]+
    # \frac{1}{105}[ (\delta_{ij}\delta_{kl}+ \delta_{ik}\delta_{jl} +
    # \delta_{il}\delta_{jk})(T_{mmnn} + T_{mnmn} + T_{mnnm})] \f$
    def natural_form(self):
        """Symmetric, trace-free part:

        (1) Symmetrize with Young Diagram: [][][][].
        (2) Subtract off all combos of delta_ij and T.contract(j,k) / 84
        (3) Subtract add sum(ISOtropic rank-4 tensors) *
                          sum(all double contractions) /105.
        Easy peasy.
        """
        from .syzygy import ISO
        D, = ISO[2]
        f1, f2, f3 = ISO[4]

        S = self.S()

        s2 = 6*(S.iijk & D).S()

        s4 = (3*S.iijj * (f1 + f2 + f3))
        return S - s2/7. +  s4/105.

    ## Cast in poly basis as a 9 x 9 so that
    # \f$ \delta_{ij}{kl} \f$ is diagonal
    # \return matrix 9 x 9.
    def nineXnine(self):
        """numpy matrix as:
        [[ 0  1  2  9 10 11 18 19 20]
         [ 3  4  5 12 13 14 21 22 23]
         [ 6  7  8 15 16 17 24 25 26]
         [27 28 29 36 37 38 45 46 47]
         [30 31 32 39 40 41 48 49 50]
         [33 34 35 42 43 44 51 52 53]
         [54 55 56 63 64 65 72 73 74]
         [57 58 59 66 67 68 75 76 77]
         [60 61 62 69 70 71 78 79 80]]

        or

        for item in array(map("".join, product(*repeat('xyz', 4)))
         ).reshape(3,3,3,3).transpose(0,2,1,3).ravel():
            print item

        """
        from numpy import matrix
        return matrix(
            self.asarray().reshape(3, 3, 3, 3).transpose(0, 2, 1, 3).reshape(9, 9)
        )


_dummy = lambda _: NotImplemented

## Dictionary of (j, q) projecting functions.
IRREPS = {(0, 0): _dummy,
          (0, 1): _dummy,
          (0, 2): _dummy,
          (1, 0): _dummy,
          (1, 1): _dummy,
          (1, 2): _dummy,
          (1, 3): _dummy,
          (1, 4): _dummy,
          (1, 5): _dummy,
          (2, 0): _dummy,
          (2, 1): _dummy,
          (2, 2): _dummy,
          (2, 3): _dummy,
          (3, 0): _dummy,
          (3, 1): _dummy,
          (3, 3): _dummy,
          (4, 0): Four.natural_form}


# doesn't work, solve d.ikjl
def inverse(t):
    """Non functional"""
    from numpy import array, matrix
    from numpy import linalg
    a = array(t.tolist())
    m = matrix(a.reshape(9, 9))
    mbar = linalg.inv(m)
    return array(mbar).ravel().tolist()
