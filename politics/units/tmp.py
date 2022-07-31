I = ("acceleration_.py",
     "area_.py",
     "energy_.py",
     "entropy_.py",
     "force_.py",
     "frequency_.py",
     "intensity_.py",
     "logs_.py",
     "mass_.py",
     "power_.py",
     "pressure_.py",
     "radiation_.py",
     "speed_.py",
     "times_.py")


O = ("_acceleration.py",
     "_area.py",
     "_energy.py",
     "_entropy.py",
     "_force.py",
     "_frequency.py",
     "_intensity.py",
     "_logs.py",
     "_mass.py",
     "_power.py",
     "_pressure.py",
     "_radiation.py",
     "_speed.py",
     "_times.py")

import os
if __name__ == '__main__':
    for i, o in zip(I, O):
        os.rename(i, o)
