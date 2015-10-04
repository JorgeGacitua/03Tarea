#! /home/jynd/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt



def get_k1(f,h,g,g_prima,xi):
    '''
    Recive la funcion vectorial 'f', el paso 'h', la funcion original g y su derivada y
    entrega el verctor k1 para utilizar en el metodo Runge Kutta 3
    '''
    f_eval=f(xi,g,g_prima)
    return f_eval[0],f_eval[1]


def get_k2(f,h,g,g_prima,xi):
    '''
    Recive la funcion vectorial 'f', el paso 'h', la funcion original g y su derivada y
    entrega el verctor k1 para utilizar en el metodo Runge Kutta 3
    '''
    k1 = get_k1(f,h,g,g_prima,xi)
    f_eval = f(xi+h/2,g + k1[0]/2, g_prima + k1[1]/2)
    return h * f_eval[0], h * f_eval[1]

def get_k3(f,h,g,g_prima,xi):
    '''
    Recive la funcion vectorial 'f', el paso 'h', la funcion original g y su derivada y
    entrega el verctor k2 para utilizar en el metodo Runge Kutta 3
    '''
    k1 = get_k1(f,h,g,g_prima,xi)
    k2 = get_k2(f,h,g,g_prima,xi)
    f_eval = f(xi+h,g - k1[0]+2*k2[0], g_prima - k1[0]+2*k2[0])
    return h * f_eval[0], h * f_eval[1]

def RK3(f,h,g,g_prima,xi):
    '''
    Recive la funcion vectorial 'f', el paso 'h', la funcion original g y su derivada y
    el resultado usando el metodo Runge Kutta 3
    '''
    k1 = get_k1(f,h,g,g_prima,xi)
    k2 = get_k2(f,h,g,g_prima,xi)
    k3 = get_k1(f,h,g,g_prima,xi)

    Y1n1 = g + (k1[0]+4*k2[0]+k3[0])/6.0
    Y2n1 = g_prima + (k1[1]+4*k2[1]+k3[1])/6.0
    return Y1n1,Y2n1

#----------------------------------Problema 1----------------------------------#
mu=1.560 #rut=17471560-2
N_steps = 10000
dyds=np.zeros(N_steps)
y=np.zeros(N_steps)
h=(20*np.pi)/N_steps

def f(xi,g,g_prima):
    f1=g_prima
    f2=-g-mu*((g**2)-1)*g_prima
    return f1,f2

dyds[0]=0
y[0]=4.0

for i in range(1, N_steps):
    y[i], dyds[i]=RK3(f,h,y[i-1],dyds[i-1],0)


s=[h * i for i in range(N_steps)]


plt.figure(1)
plt.clf()

plt.plot(s,y)

plt.xlabel('s ()')
plt.ylabel('$y(s)$', fontsize=18)
plt.show()
plt.draw()

plt.figure(2)
plt.clf()

plt.plot(y,dyds)

plt.xlabel('y(s)',fontsize=18)
plt.ylabel('$dy/ds$', fontsize=18)
plt.show()
plt.draw()
