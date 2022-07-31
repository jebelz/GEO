"""Units Conversion Functions, which live in their dimensionally named
modules.

The sub-package is designed from the user using tab-completion, e.g:

>>>from geo.politics import units as u
>>> u.<TAB>
u.acceleration  u.force         u.logs          u.radiation
u.area          u.frequency     u.mass          u.speed
u.energy        u.intensity     u.power         u.times
u.entropy       u.length        u.pressure      u.volume

Then, you choose your dimensions and do it again:

In [2]: u.energy.
u.energy.BTU2J              u.energy.eV2J
u.energy.J2BTU              u.energy.foe
u.energy.J2eV               u.energy.foe2J
u.energy.J2kWhours          u.energy.kWhours2J
u.energy.Ry                 u.energy.photons_per_joule
u.energy.eV                 u.energy.planck_E

Then you can see how many Rydbergs are in a FOE:

>>> (u.energy.Ry/u.energy.foe)(1)
4.5874249581853503e+61

Wow, that's a lot of Rydbergs.
"""
## \namespace geo.politics.units Unit Converting Functions.

__version__ = "1.0"
print "importing %s version::%s"%(__name__, __version__)

from ._bases import SI

from . import times
from . import frequency


from . import length
from . import area
from . import volume


from . import speed
from . import acceleration

from . import mass

from . import force
from . import pressure
from . import energy

from . import power

from . import entropy

from . import intensity

from . import logs

from . import radiation

