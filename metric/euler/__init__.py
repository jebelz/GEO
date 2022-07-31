"""The Euler subpackage is about reppin' rotations:

CORE:

charts/    Are various angle-base rotation representations (ACTIVE/INTRINSIC)
hamilton/  The unit-quaternion (versor) rep: (ACTIVE/INTRINISC when combined)


Bonus:

van_elfrinkhof.py  Quaternions as real 4x4 matrices.
"""

## \namespace geo.metric.euler
# <a href="http://en.wikipedia.org/wiki/Euler_angles">Euler</a> angles
# and other rotations.
from .charts import *
from .hamilton import *
