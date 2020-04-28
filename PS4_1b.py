#Problem 1 - part b
import numpy
from scipy.optimize import *
from matplotlib import *

K1 = 0.1
K2 = 10
mesh = 1000
KDi = numpy.linspace(0.1,100,mesh)
sitaB = KDi/(1+KDi)

# For x*
xnorm_K1 = numpy.zeros(mesh)
xnorm_K2 = numpy.zeros(mesh)
for a in range(mesh):
        def for_x(E):
          x1 = E[0]
          x2 = E[1]
          Z = numpy.empty((2))
          Z[0] = 5*sitaB[a]*(1-x1)/(K1+1-x1)-x1/(K1+x1)
          Z[1] = 5*sitaB[a]*(1-x2)/(K2+1-x2)-x2/(K2+x2)
          return Z
        xnorm_K1[a],xnorm_K2[a] = fsolve(for_x, numpy.array([1,1]))

# For y*
ynorm_K1 = numpy.zeros(mesh)
ynorm_K2 = numpy.zeros(mesh)
for b in range(mesh):
    def for_y(F):
        y1 = F[0]
        y2 = F[1]
        Z = numpy.empty((2))
        Z[0] = 10*xnorm_K1[b]*(1-y1)/(K1+1-y1)-y1/(K1+y1)
        Z[1] = 10*xnorm_K2[b]*(1-y2)/(K2+1-y2)-y2/(K2+y2)
        return Z
    ynorm_K1[b],ynorm_K2[b] = fsolve(for_y, numpy.array([1,1]))

#Plot
matplotlib.pyplot.figure(figsize = (10,8))
matplotlib.pyplot.title(
    'Responses vs. Non-dimensional Input',fontsize = 30)
matplotlib.pyplot.xlabel('$1 /'+'\kappa_D$',fontsize = 15)
matplotlib.pyplot.ylabel('',fontsize = 15)
matplotlib.pyplot.xlim(0,100)
matplotlib.pyplot.ylim(0,1.1)
matplotlib.pyplot.plot(KDi,sitaB,color='black',label='$\Theta _b$ ($\kappa$ = 0.1)');
matplotlib.pyplot.plot(KDi,xnorm_K1,color='red',label='x* ($\kappa$ = 0.1)');
matplotlib.pyplot.plot(KDi,xnorm_K2,color='red',linestyle='dashed',label='x* ($\kappa$ = 10)');
matplotlib.pyplot.plot(KDi,ynorm_K1,color='blue',label='y* ($\kappa$ = 0.1)');
matplotlib.pyplot.plot(KDi,ynorm_K2,color='blue',linestyle='dashed',label='y* ($\kappa$ = 10)');
matplotlib.pyplot.legend(loc='lower right',fontsize = 20);
