#------------------------------------------------
#
# This code generates different model parameters
#   value0a: initial crss
#   value0b: saturated crss
#   value0c: hardening rate
#   value0d: saturated back stress
#
#------------------------------------------------

import os
import math
import sys 
import random

caseNM0="cycDefModelling0"

f=open("err.txt","w+")
f.close()
f=open("parABCD.dat","w+")
f.close()



for iloop in range(1, 100):
    value0a = 20.
    value0b = 180.
    value0c = 230.
    value0d = 50.

    Ia = random.randint(1, 10)
    Ib = random.randint(1, 10)
    Ic = random.randint(1, 10)
    Id = random.randint(1, 10)

    valueXa = value0a / 2.0 + value0a / 10.0 * Ia
    valueXb = value0b / 2.0 + value0b / 10.0 * Ib
    valueXc = value0c / 2.0 + value0c / 10.0 * Ic
    valueXd = value0d / 2.0 + value0d / 10.0 * Id

    with open("parABCD.dat", "a") as f:
        f.write("{:.{}f}".format(valueXa, 4) + "  ")
        f.write("{:.{}f}".format(valueXb, 4) + "  ")
        f.write("{:.{}f}".format(valueXc, 4) + "  ")
        f.write("{:.{}f}".format(valueXd, 4) + "\n")
        f.close()



