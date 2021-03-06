#Problem 2 - part b
import numpy
from scipy.optimize import *
from matplotlib import *
from mpl_toolkits import *

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

m = 190 #could not achieve a plot above this mesh number with my PC
i1 = numpy.logspace(-2,3,num = m)
i2 = numpy.logspace(-2,3,num = m)
A = numpy.zeros((m,m))
B = numpy.zeros((m,m))
C = numpy.zeros((m,m))
P = numpy.array([1,1,1])
for x1 in range(m):
    P = numpy.array([1,1,1])
    for x2 in range(m):
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
pl.plot_surface(numpy.log10(x), numpy.log10(y), A, cmap='plasma')
pl.set_title('Low Ks values', fontsize = 30)
pl.set_xlabel('[I1] (in log10)', fontsize = 15)
pl.set_ylabel('[I2] (in log10)', fontsize = 15)
pl.set_zlabel('[A]', fontsize = 15)

matplotlib.pyplot.show()
