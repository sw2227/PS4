#Problem 2 - part d
import numpy
from scipy.optimize import *
from matplotlib import *
from mpl_toolkits import *

Vmax1 = 5.0;
Vmax2 = 5.0;
Vmax3 = 1.0;
Vmax4 = 1.0;
Ks1 = 35.0;
Ks2 = 35.0;
Ks3 = 35.0;
Ks4 = 35.0;
Ki1 = 1.0;
Ki2 = 1.0;
Stot = 100;

mesh = 190 #could not achieve a plot above this mesh number with my PC
i1 = numpy.logspace(-2,3,num = n)
i2 = numpy.logspace(-2,3,num = n)
A = numpy.zeros((n,n))
B = numpy.zeros((n,n))
C = numpy.zeros((n,n))
P = numpy.array([1,1,1])
for x1 in range(n):
    P = numpy.array([1,1,1])
    for x2 in range(n):
        def with_inhibitor(E):
            A = E[0]
            B = E[1]
            C = E[2]
            d1 = i1[x1]
            d2 = i2[x2]
            Z = numpy.empty((3))
            Z[0] = (Vmax1*A)/((Ks1+A)*(1+(d1/Ki1)))-(Vmax3*B)/(Ks3+B)
            Z[1] = (Vmax2*A)/((Ks2+A)*(1+(d2/Ki2)))-(Vmax4*C)/(Ks4+C)
            Z[2] = Stot-A-B-C
            return Z
        A[x1,x2],B[x1,x2],C[x1,x2] = fsolve(with_inhibitor, P)
        P[0] = A[x1,x2]
        P[1] = B[x1,x2]
        P[2] = C[x1,x2]


plot1 = matplotlib.pyplot.figure(figsize=(20,10))
pl = matplotlib.pyplot.axes(projection="3d")
x, y = numpy.meshgrid(i1,i2)
pl.plot_surface(numpy.log10(x), numpy.log10(y), A, cmap='cool')
pl.set_title('Low Ks values', fontsize = 30)
pl.set_xlabel('[I1] (in log10)', fontsize = 15)
pl.set_ylabel('[I2] (in log10)', fontsize = 15)
pl.set_zlabel('[A]', fontsize = 15)

matplotlib.pyplot.show()
