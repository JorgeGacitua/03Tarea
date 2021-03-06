#! /home/jynd/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import ode

from mpl_toolkits.mplot3d import Axes3D

def Lorenz_system(t,k):
    return [sigma*(k[1]-k[0]), k[0]*(rho-k[2])-k[1], k[0]*k[1]-beta*k[2]]

sigma=10
beta=8/3
rho=28

N_steps = 10000
x=np.zeros(N_steps)
y=np.zeros(N_steps)
z=np.zeros(N_steps)
xyzi=[15,3,10]
t_values = np.linspace(1e-3, 30, N_steps)

r = ode(Lorenz_system)
r.set_integrator('dopri5')
r.set_initial_value(xyzi)


for i in range(N_steps):
    r.integrate(t_values[i])
    x[i], y[i], z[i] = r.y

#------------------------ se crea el grafico----------------------------------#

plt.figure(1)
plt.clf()

plt.plot(x,y)

plt.xlabel('$x(t)$', fontsize=18)
plt.ylabel('$y(y)$', fontsize=18)
plt.show()
plt.draw()

plt.figure(2)
plt.clf()

plt.plot(y,z)

plt.xlabel('$y(t)$',fontsize=18)
plt.ylabel('$z(t)$', fontsize=18)
plt.show()
plt.draw()

plt.figure(3)
plt.clf()

plt.plot(x,z)

plt.xlabel('$x(t)$', fontsize=18)
plt.ylabel('$z(t)$', fontsize=18)
plt.show()
plt.draw()





fig = plt.figure(4)
fig.clf()

ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')


ax.plot(x, y, z)

ax.set_xlabel('$x(t)$', fontsize=18)
ax.set_ylabel('$y(t)$', fontsize=18)
ax.set_zlabel('$z(t)$', fontsize=18)

plt.show()
plt.draw()
