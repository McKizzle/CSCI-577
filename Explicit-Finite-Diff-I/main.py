#!/usr/bin/env python
import scipy.signal as sig
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import time as time

# Renalds: The characteristic length / ch

def main():
    a = 0.0
    b = 1.0
    u = np.array(np.zeros(21)).astype(float)
    u[0] = a
    u[-1] = b
    x = np.arange(0.0, float(len(u)), 1.0) / float(len(u) - 1)
    u = relax1D(u)

    animate1Dframes(x, u, 12.0)
    #plt.plot(x, u)
    #plt.show()

def relax1D(u, a=0.0, b=1.0,
            num=np.array([1.0, 0.0, 1.0]),
            den=2.0, tol=1E-4):
    """ Relaxes a 1D array of data given. 
        :param u: the initial array of values to relax.  
        :param a: constant value to set the initial index of u
        :param b: constant value to set the terminal index of u
        :param num: the array to pass into the convolver. 
        :param den: the default value to divide u_new by at each iteration. 
    """
    
    data = []

    d_sol = 1
    u[0] = a
    u[-1] = b
    while(d_sol > tol):
        data.append(u)
        u_new = sig.convolve(u, num, 'same') / den
        u_new[0] = a
        u_new[-1] = b
        d_sol = Delta_solution(u, u_new)
        u = u_new
    return data

def animate1Dframes(x, data, fps):
    """ Animates a 2D array of data using pyplot. 
        :param x: the x values to plot over
        :param y: the y value 'frames' to plot at each iteration. 
        :param fps: the number of frames per second. 
    """
    
    #print data

    for u in data:
        plt.plot(x, u)
        plt.draw()
        time.sleep(1.0 / fps)


def Delta_solution(u_n, u_nm1):
    """ Expects np.arrays to calculate the Delta solution """
    return np.sqrt(np.sum((u_n - u_nm1)**2))

if __name__ == "__main__":
    main()
