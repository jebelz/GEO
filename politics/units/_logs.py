## \namespace geo.politics.units._logs
# <a href="http://en.wikipedia.org/wiki/Logarithmic_scale#Logarithmic_unit">
# Logrithmic units</a>

from ._bases import *

__all__ = ('bel', 'dB', 'fromdB', 'Np', 'Np2dB', 'dB2Np', 'dBi2dBd', 'dBd2dBi', 'dBi2dBq', 'dBq2dBi')


class _Log(object):

    def __init__(self, base, prefix):
        self.base = base
        self.prefix = prefix
        # work in __init__, a sin.
        self.func = functools.partial(pow, base)
        self._logbase = scipy.log(base)

    def __rmul__(self, linear):
        return self.func(linear * self.prefix)

    def __call__(self, other):
        return scipy.log(other)/self._logbase/self.prefix

bel = _Log(base=10, prefix=1)
dB = _Log(base=10, prefix=constants.deci)
fromdB = lambda x: x * dB

Np = _Log(base=scipy.e, prefix=1)

## <a href="http://en.wikipedia.org/wiki/Neper">Nepers</a> vs
# <a href="http://en.wikipedia.org/wiki/Decibel">Decibels</a>.
Np2dB = SI(10 * math.log10(math.exp(1)))

## <a href="http://en.wikipedia.org/wiki/Neper">Nepers</a> vs
# <a href="http://en.wikipedia.org/wiki/Decibel">Decibels</a>.
dB2Np = ~Np2dB


## dB(dipole) http://en.wikipedia.org/wiki/DBi#Antenna_measurements  .
dBi2dBd = scipy.poly1d([1, -2.15])
dBd2dBi = scipy.poly1d([1, 2.15])

## dB(quaterwave) http://en.wikipedia.org/wiki/DBi#Antenna_measurements .
dBi2dBq = scipy.poly1d([1, 0.15])
dBq2dBi = scipy.poly1d([1, 0.15])
