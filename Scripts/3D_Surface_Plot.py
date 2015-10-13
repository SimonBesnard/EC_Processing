from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import *

a0=0.387217066
a1=0.328086337
a2=0.004177033
p0=1.117997
p1=5.426564
p2=-220.745499

df = np.loadtxt('/home/simonbesnard/Desktop/3plot_Test.csv', delimiter=',',skiprows=1)
obs= np.loadtxt('/media/simonbesnard/External_SB/SimonBesnard/PhD_MPI/RO1/Rproject/EC_Processing/obs.csv', delimiter=',',skiprows=1)
X = df[:,0]
Y = df[:,1]
X, Y = np.meshgrid(X, Y)
Z1=a0*(X**a1)*exp(-a2*X)*p0*(1-exp(-exp(-p1)*(Y-p2)))
#Z1 = dat[:,3]

fig = plt.figure()
ax = fig.gca(projection='3d',elev=30)

surf1 = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
ax.plot3D(obs[:,1], obs[:,2], obs[:,0], c='k', marker='o',ms=5,lw=0,mfc='None')

fig.colorbar(surf1, shrink=0.75, aspect=30)
plt.xlabel('Age',fontsize=8)
plt.ylabel('Precipitation',fontsize=8)

plt.savefig('prcp.png',dpi=250)
plt.show()
'''
#plt.zlabel('f(Prcp)',fontsize=9)
#ax.set_zlim(-1.01, 1.01)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))Z2 = dat[:,4]
fig = plt.figure()
ax = fig.gca(projection='3d')

surf2 = ax.plot_surface(X, Y, Z2, rstride=1, cstride=1, cmap=cm.jet,
        linewidth=0, antialiased=False)
#ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))


fig.colorbar(surf2, shrink=0.5, aspect=5)

plt.xlabel('Prec',fontsize=9)
plt.ylabel('Age',fontsize=9)
#plt.zlabel('f(Age)',fontsize=9)
plt.savefig('age.png',dpi=250)

'''
