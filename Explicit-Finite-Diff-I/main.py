#!/usr/bin/env python
import scipy.signal as sig
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import time as time

# Renalds: The characteristic length / ch
def main():
    # Lets relax a 1D array of values. 
    N = 20
    a = 0.0
    b = 1.0
    u = np.zeros(20)
    u[0] = a
    u[-1] = b
    x = np.arange(0.0, float(len(u)), 1.0) / float(len(u) - 1)
    u = relax1D(u)
    #animate1Dframes(x, u)
    plt.plot(x, u[-1], "-sk")
    plt.savefig("Relaxed1D.png")
    plt.close()

    # Now lets relax a 2D array of values. Such that the top and bottom edges are set to zero
    # and the left and right edges are set to one. 
    u2D = np.zeros((N, N))
    u2D[0,:]  = a
    u2D[-1,:] = a
    u2D[:,0]  = b
    u2D[:,-1] = b
    u2D = relax2D(u2D)
    u2D_cf = plt.contour(u2D[-1], levels=np.arange(0, 1, 0.1))
    plt.clabel(u2D_cf, colors='k')
    plt.savefig("Relaxed2D.png")
    plt.close()

def relax2D(u, a=0.0, b=1.0, 
            num=np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]]), 
            den=4.0, tol=1E-4):
    """ 
        :param u: the initial array of values to relax.  
        :param a: constant value to set the initial index of u
        :param b: constant value to set the terminal index of u
        :param num: the stencil array. (needs to be 2D and of equal dimensions)
        :param den: the default value to divide u_new by at each iteration. 
    """
    data = []
    
    d_sol = 1
    u[0] = a
    u[-1] = b
    while(d_sol > tol):
        data.append(u)
        u_new = sig.convolve2d(u, num, 'same') / den
        u_new[0, :] = a
        u_new[-1,:] = a
        u_new[:, 0] = b
        u_new[:,-1] = b
        d_sol = Delta_solution(u, u_new)
        u = u_new
        #print d_sol
    return data

def animate2Dframes(data):
    """ Animates a 2D array of data using pyplot. 
        :param x: the x values to plot over
        :param y: the y value 'frames' to plot at each iteration. 

        Follows example at http://www.lebsanft.org/?p=48 but uses imshow
    """
    plt.ion()
    img = plt.imshow(data[0])

    for u in data:
        img.set_data(u)
        plt.draw()

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

def animate1Dframes(x, data):
    """ Animates a 2D array of data using pyplot. 
        :param x: the x values to plot over
        :param y: the y value 'frames' to plot at each iteration. 

        Follows example at http://www.lebsanft.org/?p=48
    """
    plt.ion() # Set the plot to animated.    
    ax1 = plt.axes()
    line, = plt.plot(x, data[0])

    for u in data:
        line.set_ydata(u)
        plt.draw()

def Delta_solution(u_n, u_nm1):
    """ Expects np.arrays to calculate the Delta solution """
    return np.sqrt(np.sum((u_n - u_nm1)**2))

if __name__ == "__main__":
    main()
