#! /home/jynd/anaconda/bin/python

import numpy as np
import matplotlib.pyplot as plt

mu=1.560 #rut=17471560-2

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
    f_eval = f(xi+h/2,phi_n + k1[0]/2, w_n + k1[1]/2)
    return h * f_eval[0], h * f_eval[1]

def get_k3(f,h,g,g_prima,xi):
    '''
    Recive la funcion vectorial 'f', el paso 'h', la funcion original g y su derivada y
    entrega el verctor k2 para utilizar en el metodo Runge Kutta 3
    '''
    k1 = get_k1(f,h,g,g_prima,xi)
    k2 = get_k2(f,h,g,g_prima,xi)
    f_eval = f(xi+h,phi_n + k1[0]+2*k2[0], w_n + k1[0]+2*k2[0])
    return h * f_eval[0], h * f_eval[1]

def RK3(f,h,g,g_prima,xi):
    '''
    Recive la funcion vectorial 'f', el paso 'h', la funcion original g y su derivada y
    el resultado usando el metodo Runge Kutta 3
    '''
    k1 = get_k1(f,h,g,g_prima,xi)
    k2 = get_k2(f,h,g,g_prima,xi)
    k3 = get_k1(f,h,g,g_prima,xi)

    Y1n1 = g + (k1[0]+4*k2[0]+k3[0])/6
    Y2n1 = g_prima + (k1[1]+4*k2[1]+k3[1])/6
    return phi_n1, w_n1
