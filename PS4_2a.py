#Problem 2 - part a
#should first do "pip install scipy" if scipy is not installed
import numpy
import math
from scipy.optimize import *

Vmax1 = 5.0;
Vmax2 = 5.0;
Vmax3 = 1.0;
Vmax4 = 1.0;
Ks1 = 5.0;
Ks2 = 5.0;
Ks3 = 5.0;
Ks4 = 5.0;
Ki1 = 1.0;
Ki2 = 1.0;
Stot = 100;

def system(E):
    A = E[0]
    B = E[1]
    C = E[2]
    Z = numpy.empty((3))
    Z[0] = (Vmax1*A)/(Ks1+A)-(Vmax3*B)/(Ks3+B)
    Z[1] = (Vmax2*A)/(Ks2+A)-(Vmax4*C)/(Ks4+C)
    Z[2] = Stot-A-B-C
    return Z

P = numpy.array([1,1,1])
A,B,C = fsolve(system,[1,1,1])
print(A,B,C)
