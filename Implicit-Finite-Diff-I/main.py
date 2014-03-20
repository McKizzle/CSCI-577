#!/usr/bin/env python
import scipy.signal as sig
import scipy as scpy
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import time as tm

def main():
    N = 51
    time_steps = 2000
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

    I = scpy.identity(N) 
    
    # Build the A matrix
    A = scpy.sparse.lil_matrix((N, N))
    A.setdiag(np.zeros(N) - 2.0 * alpha)
    A.setdiag(np.zeros(N) + alpha, k=1)
    A.setdiag(np.zeros(N) + alpha, k=-1)

    S = (I - A)
    S = scpy.sparse.csr_matrix(S)
    
    U = simulate(u, S, dx, dt=dt, n=time_steps)
    
    #animate1Dframes(X, U)
    # Plot the final curve
    plt.plot(X, U[-1], '-ok')
    plt.title("Implicit Method")
    plt.ylabel("x")
    plt.xlabel("u")
    plt.savefig("implicitDiff.png")
    plt.close()
    
    # Plot the area under the curve at each time step
    intU = np.array([np.array([i, integral(X, U[i])]) for i in range(0, time_steps)])
    plt.plot(intU[:,0], intU[:,1])
    plt.title("Implicit Method Area")
    plt.ylabel("time")
    plt.xlabel("area")
    plt.savefig("implicitDiff_Area.png")
    plt.close()

    ############### Crank Nicolson Method. ######
    # Build the identity (I) matrix
    alpha = dt / (2.0 * dx**2.0)
    I = scpy.identity(N)
    
    # Build the A matrix
    A = scpy.sparse.lil_matrix((N, N))
    A.setdiag(np.zeros(N) - 2.0 * alpha)
    A.setdiag(np.zeros(N) + alpha, k=1)
    A.setdiag(np.zeros(N) + alpha, k=-1)

    Q = (I - A)
    R = (I + A)
    
    S = np.mat(scpy.linalg.inv(R)) * Q
    S = scpy.sparse.csr_matrix(S)
    
    U = simulate(u, S, dx=dx, dt=dt, n=time_steps)


    # Plot the final curve
    #animate1Dframes(X, U)
    plt.plot(X, U[-1], '-ok')
    plt.title("Crank Nicolson Method")
    plt.ylabel("x")
    plt.xlabel("u")
    plt.savefig("CrankNicolsonMethod.png")
    plt.close()
    
    # Plot the area under the curve at each time step
    intU = np.array([np.array([i, integral(X, U[i])]) for i in range(0, time_steps)])
    plt.plot(intU[:,0], intU[:,1])
    plt.title("Implicit Method Area")
    plt.ylabel("time")
    plt.xlabel("area")
    plt.savefig("CrankNicolsonMethod_Area.png")
    plt.close()

    return 0

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
