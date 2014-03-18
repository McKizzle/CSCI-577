#!/usr/bin/env python
import scipy.signal as sig
import scipy as scpy
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import time as tm

def main():
    N = 51
    a = 0
    b = 1
    x_0 = 0.0
    x_N = 1e3
    X = np.linspace(x_0, x_N, N)
    dx = X[1] - X[0]
    dt = 1.0
    u = np.zeros(N)
    u[int(N/2)] = b 
    alpha = dt / dx**2.0
    A = scpy.sparse.lil_matrix((N, N))
    A.setdiag(np.ones(N) + 2.0 * alpha)
    A.setdiag(np.zeros(N) - alpha, k=1)
    A.setdiag(np.zeros(N) - alpha, k=-1)
    A = A.tocsr()
    U = simulate(u, A, dx, dt=0.001)
    plt.plot(X, U[-1], '-ok')
    plt.show()
    return 0

def simulate(u_0, A, dx, boundries = 0.0, dt=0.001, n=2000):
    """ Run the simulation for a defined number of steps 
        :param 
        Returns the u^{n}
    """
    U = [] 
    u_n = u_0
    for i in range(0, n):
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


def wiki_code():
    N = 51
    # Create a sparse matrix in the "list of lists" format. 
    # This provides a flexible syntax for creating sparse matrices.
    A = scpy.sparse.lil_matrix((N, N))
    # Finite difference matrices will always have some 
    # regular structure involving diagonals and off diagonals.
    # The lil_matrix object has a function to help you with this.
    A.setdiag(np.ones(N))              # The diagonal
    A.setdiag(-2*np.ones(N-1),k=1)       # The fist upward off-diagonal.
    A.setdiag(2*np.ones(N-1), k=-1)
    # etc. observe that I leave it to you to insert correct values along the diagonals...

    # To fix the boundaries, you'll have to reset some rows.
    # The numpy indexing you expect to work will work.
    #A[0,:] = np.zeros(N) # zero out the first row
    #A[0,0] = 1.0       # set diagonal in that row to 1.0


    # For performance, other sparse matrix formats are preferred.
    # Convert to "Compressed Row Format" CR
    A = A.tocsr()

    # Helpful for diagnostics

    print A.todense()

    # and 

    plt.spy(A.todense())
    plt.show()

if __name__ == "__main__":
    main()
