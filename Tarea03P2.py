#! /home/jynd/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#Fuente: http://www.codeproject.com/Tips/792927/Fourth-Order-Runge-Kutta-Method-in-Python#
def rKN(x, fx, n, hs):
    k1 = []
    k2 = []
    k3 = []
    k4 = []
    xk = []
    for i in range(n):
        k1.append(fx[i](x)*hs)
    for i in range(n):
        xk.append(x[i] + k1[i]*0.5)
    for i in range(n):
        k2.append(fx[i](xk)*hs)
    for i in range(n):
        xk[i] = x[i] + k2[i]*0.5
    for i in range(n):
        k3.append(fx[i](xk)*hs)
    for i in range(n):
        xk[i] = x[i] + k3[i]
    for i in range(n):
        k4.append(fx[i](xk)*hs)
    for i in range(n):
        x[i] = x[i] + (k1[i] + 2*(k2[i] + k3[i]) + k4[i])/6
    return x

#------------------------Se definen las funiones a resolver--------------------#

def dxdt(k):return sigma*(k[1]-k[0])
def dydt(k):return k[0]*(rho-k[2])-k[1]
def dzdt(k):return k[0]*k[1]-beta*k[2]

sigma=10
beta=8/3
rho=28

N_steps = 10000
x=np.zeros(N_steps)
y=np.zeros(N_steps)
z=np.zeros(N_steps)
h=(100) / N_steps
x[0]=1
y[0]=1
z[0]=1
f=[dxdt,dydt,dzdt]

for i in range(1, N_steps):
    xyz=[x[i-1],y[i-1],z[i-1]]
    x[i], y[i],z[i]=rKN(xyz,f,3,h)


#------------------------ se crea el grafico----------------------------------#


fig = plt.figure(1)
fig.clf()

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

t = [h * i for i in range(N_steps)]


ax.plot(x, y, z)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
