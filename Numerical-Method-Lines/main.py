#!/usr/bin/env python
import scipy.signal as sig
import scipy as scpy
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import time as tm

def main():
    N = 100
    dt = 0.001
    t_max = 4.0
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

    #plt.plot(X, u)
    #plt.show()
    
    # 
    # kBu^u + v2u/2x
    #

    M = scpy.sparse.lil_matrix((N, N)) 
    M.setdiag(np.ones(N), k=1)
    M.setdiag(np.ones(N) * -2.0)
    M.setdiag(np.ones(N), k=-1)
    
    # Construct the Matricies to reflexct the lax formula 
    J = scpy.sparse.lil_matrix((N, N))
    J.setdiag(np.zeros(N) - 1, k=-1)
    J.setdiag(np.zeros(N) + 1, k=1)
    J[0,:] = np.zeros(N)
    J[-1,:] = np.zeros(N)
    J[0, 0] = 1
    J[-1, -1] = 1
    return 0

def rhs(t, y, A):
    return A * y

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

    pass

def wiki_upwinded_derivative(dx, N, v):
    A = scpy.sparse.lil_matrix((N, N))

    """ Upwind formula from the wiki """
    if v < 0:
        A.setdiag(np.zeros(N), k=-1)
        A.setdiag(np.ones(N) * -1.0)
        A.setdiag(np.ones(N), k=1)
    elif v = 0:
        A = A * 0
    elif v > 0:
        A.setdiag(np.onos(N) * -1.0, k=-1)
        A.setdiag(np.ones(N))
        A.setdiag(np.zeros(N), k=1)
 
    A[0,:] = np.zeros(N)
    A[-1,:] = np.zeros(N)
    A[0, 0] = 1
    A[-1, -1] = 1

    return A * v / dx


def integral(x, u):
    """ Calculate the area underneath a curve.
        :param x: the x values of the curve. 
        :param u: the y values of the curve.
    """
    area = 0
    h = x[1] - x[0]
    for i in range(1, len(u)):
        p0 = np.array([x[i - 1], 0])
        p1 = np.array([x[i - 1], u[i - 1]])
        p2 = np.array([x[i], u[i]])
        p3 = np.array([x[i], 0])
        a = np.linalg.norm(p0 - p1)
        b = np.linalg.norm(p2 - p3)
        area = area + 0.5 * (a + b) * h
    
    return area

def simulate(u_0, A, boundries = 0.0, n=2000):
    """ Run the simulation for a defined number of steps 
        :param 
        Returns the u^{n}
    """
    U = [] 
    u_n = u_0
    for i in range(0, int(n)):
        U.append(u_n)
        u_np1 = scpy.sparse.linalg.spsolve(A,u_n) 
        u_np1[0] = boundries
        u_np1[-1] = boundries
        u_n = u_np1
        
    return U

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
        #tm.sleep(0.25)


if __name__ == "__main__":
    main()
