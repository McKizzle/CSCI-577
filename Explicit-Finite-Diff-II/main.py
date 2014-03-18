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
    u = relax1D(u, x, dt=0.0001, tol=1e-4, a=a, b=b)
    plt.plot(x, u[-1], "-sk")
    plt.savefig("partial1D.png")
    plt.close()
    #animate1Dframes(x, u)

    u2D = np.zeros((N, N))
    u2D[:,-1]= b
    u2D[:,0] = b
    u2D[-1,:] = b
    u2D[0,:]  = a
    u2D = relax2D(u2D, x, y, dt=0.0001, tol=1e-4, a=a, b=b)
    u2D_cf = plt.contour(x, y, u2D[-1], levels=np.arange(0, 1, 0.1))
    plt.clabel(u2D_cf, colors='k')
    plt.savefig("partial2D.png")
    plt.close()
    #animate2Dframes(u2D)


def d2dx2(u, dx, t):
    return sig.convolve(u, [1.0, -2.0, 1.0], 'same') / dx**2.0

def nabla_x(u, x, dt, t):
    dx = x[1] - x[0]
    return u  + d2dx2(u, dx, t) * dt

def relax1D(u, x, dt=0.005, tol=1e-4, a=0.0, b=1.0):
    data = []
    d_sol = 1.0
    while(d_sol > tol):
        data.append(u)
        u_new = nabla_x(u, x, dt, 0.0) 
        u_new[0] = a
        u_new[-1] = b
        d_sol = Delta_solution(u, u_new)
        u = u_new
    return data

def d2dx2_d2y2(u, dx, dy, t):
    return sig.convolve(u, [[0.0, 1.0, 0.0], [1.0, -4.0, 1.0], [0.0, 1.0, 0.0]], 'same') / (2.0 * (dy * dx))

def nabla_xy(u, x, y, dt, t):
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    dudx = d2dx2_d2y2(u, dx, dy, t)
    return u + dudx * dt

def relax2D(u, x, y, dt=0.00005, tol=1e-4, a=0.0, b=1.0):
    """ testing the relax2D """ 
    data = []
    d_sol = 1.0
    while(d_sol > tol):
        data.append(u)
        u_new = nabla_xy(u, x, y, dt, 0.0) 
        u_new[:,-1]= b
        u_new[:,0] = b
        u_new[-1,:] = b
        u_new[0,:]  = a
        d_sol = Delta_solution(u, u_new)
        u = u_new
    return data

def animate2Dframes(data):
    """ Animates a 2D array of data using pyplot. 
        :param x: the x values to plot over
        :param y: the y value 'frames' to plot at each iteration. 

        Follows example at http://www.lebsanft.org/?p=48 but uses imshow
    """
    plt.ion()
    img = plt.imshow(data[0])
    plt.colorbar()

    for u in data:
        img.set_data(u)
        plt.draw()

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
