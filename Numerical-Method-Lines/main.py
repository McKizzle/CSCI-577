#!/usr/bin/env python
import scipy.signal as sig
import scipy as scpy
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import time as tm
import scipy.integrate as scint

def main():
    N = 100
    dt = 0.1
    t_max = 15
    time_steps = t_max / dt
    v = 5.0
    k = 1.0
    a = 0
    b = 1
    x_0 = 0.0
    x_N = 10.0
    x_stp_size = x_N / N 
    X = np.linspace(x_0, x_N, N)
    x_sw0 = 20 
    x_swk = x_sw0 + int(0.05 * float(N))
    dx = X[1] - X[0]
    u = np.zeros(N)
    u[x_sw0:x_swk] = b
    u[0] = a
    alpha = dt / dx
    u[-1] = a
    U = [u]

    A = convection_diffusion(dx, N, v=v, k=k);
    
    # set up the integrtn env. 
    intgr = scint.ode(rhs)
    intgr.set_integrator('vode', method='bdf')  
    intgr.set_f_params(A)
    intgr.set_initial_value(u, 0.0)
    # Outer integration loop over times
    while intgr.t < 10:
        U.append(intgr.integrate(intgr.t+dt))

    animate1Dframes(X, U)
    #plt.plot(X, U[0])
    #plt.show()

    return 0

def rhs(t, y, A):
    return A * y

def wiki_upwinded_derivative(dx, N, v):
    A = scpy.sparse.lil_matrix((N, N))

    """ Upwind formula from the wiki """
    if v < 0:
        A.setdiag(np.zeros(N), k=-1)
        A.setdiag(np.ones(N) * -1.0)
        A.setdiag(np.ones(N), k=1)
    elif v == 0:
        A = A * 0
    elif v > 0:
        A.setdiag(np.ones(N) * -1.0, k=-1)
        A.setdiag(np.ones(N))
        A.setdiag(np.zeros(N), k=1)

    A = A * v / dx
 
    A[0,:] = np.zeros(N)
    A[-1,:] = np.zeros(N) 
    A[0, 0] = 1
    A[-1, -1] = 1

    print "A: ", A.todense()

    return A.tolil()

def convection_diffusion(dx, N, v=1.0, k=1.0, firstDFormula=wiki_upwinded_derivative):
    """ Convection diffusion formula
        :param dx: the step size.
        :param v: the velocity
        :param k: some constant
        :param firstDFormula: a function that returns a matrix.

        :return: the next state of the system. 
    """

    B = scpy.sparse.lil_matrix((N, N)) 
    B.setdiag(np.ones(N), k=1)
    B.setdiag(np.ones(N) * -2.0)
    B.setdiag(np.ones(N), k=-1)
    B = B * k / dx**2
    #print "B: ", B.todense()
    B = B - firstDFormula(dx, N, v)
    B[0,:] = np.zeros(N)
    B[-1,:] = np.zeros(N)
    B[0, 0] = 1
    B[-1, -1] = 1
    return B.tolil()

def animate1Dframes(x, data):
    """ Animates a 2D array of data using pyplot. 
        :param x: the x values to plot over
        :param y: the y value 'frames' to plot at each iteration. 

        Follows example at http://www.lebsanft.org/?p=48
    """
    plt.ion() # Set the plot to animated.    
    ax1 = plt.axes()
    line, = plt.plot(x, data[0], '-*k')

    for u in data:
        line.set_ydata(u)
        plt.draw()
        #tm.sleep(0.25)


if __name__ == "__main__":
    main()
