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
    x = np.linspace(a, b, num=N)  
    y = np.linspace(a, b, num=N) 
    u = np.zeros(N)
    u[0] = a
    u[-1] = b
    u = relax1Dpartial(u, x, dt=0.0005, tol=1e-4, a=a, b=b)
    plt.plot(x, u[-1], "-sk")
    plt.savefig("partial1D.png")
    plt.close()
    #animate1Dframes(x, u)

    u2D = np.zeros((N, N))
    u2D[:,-1]= b
    u2D[:,0] = b
    u2D[-1,:] = b
    u2D[0,:]  = a
    u2D = relax2Dpartial(u2D, x, y, dt=0.0005, tol=1e-4, a=a, b=b)
    u2D_cf = plt.contour(x, y, u2D[-1], levels=np.arange(0, 1, 0.1))
    plt.clabel(u2D_cf, colors='k')
    plt.savefig("partial2D.png")
    plt.close()


def d2x(u, dx, t):
    return sig.convolve(u, [1.0, -2.0, 1.0], 'same') / dx**2.0

def euler_partial(u, x, dt, t):
    dx = x[1] - x[0]
    return u  + d2x(u, dx, t) * dt

def relax1Dpartial(u, x, dt=0.005, tol=1e-4, a=0.0, b=1.0):
    data = []
    d_sol = 1.0
    while(d_sol > tol):
        data.append(u)
        u_new = euler_partial(u, x, dt, 0.0) 
        u_new[0] = a
        u_new[-1] = b
        d_sol = Delta_solution(u, u_new)
        u = u_new
    return data

def d2xv2(u, dx, dy, t):
    return sig.convolve2d(u, [[0.0, 1.0, 0.0], [1.0, -2.0, 1.0], [0.0, 1.0, 0.0]], 'same') / (4 *(dy * dx)**2.0)

def partial(u, dx, t, num=[1.0, -2.0, 1.0]):
    dudx = []

    for i in range(0, u.shape[0]):
        dudx.append(sig.convolve(u[i,:], num, 'same') / dx**2.0)

    return np.array(dudx)

def euler_partial_v2(u, x, y, dt, t):
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    dudx = partial(u, dx, t)
    dudy = np.transpose(partial(np.transpose(u), dy, t))
    return u + (dudx + dudy) * dt

def relax2Dpartial(u, x, y, dt=0.00005, tol=1e-4, a=0.0, b=1.0):
    """ testing the relax2D """ 
    data = []
    d_sol = 1.0
    while(d_sol > tol):
        data.append(u)
        u_new = euler_partial_v2(u, x, y, dt, 0.0) 
        u_new[:,-1]= b
        u_new[:,0] = b
        u_new[-1,:] = b
        u_new[0,:]  = a
        d_sol = Delta_solution(u, u_new)
        u = u_new
    return data


#def relax2D(u, a=0.0, b=1.0, 
#            num=np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]]), 
#            den=4.0, tol=1E-4):
#    """ 
#        :param u: the initial array of values to relax.  
#        :param a: constant value to set the initial index of u
#        :param b: constant value to set the terminal index of u
#        :param num: the stencil array. (needs to be 2D and of equal dimensions)
#        :param den: the default value to divide u_new by at each iteration. 
#    """
#    data = []
#    
#    d_sol = 1
#    u[0] = a
#    u[-1] = b
#    while(d_sol > tol):
#        data.append(u)
#        u_new = sig.convolve2d(u, num, 'same') / den
#        u_new[0, :] = a
#        u_new[-1,:] = a
#        u_new[:, 0] = b
#        u_new[:,-1] = b
#        d_sol = Delta_solution(u, u_new)
#        u = u_new
#    return data

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
    line, = plt.plot(x, data[0], '-ok')

    for u in data:
        line.set_ydata(u)
        plt.draw()

def Delta_solution(u_n, u_nm1):
    """ Expects np.arrays to calculate the Delta solution """
    return np.sqrt(np.sum((u_n - u_nm1)**2.0))

if __name__ == "__main__":
    main()
