#!/usr/bin/env python
import scipy.signal as sig
import scipy as scpy
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import time as tm

def main():
    N = 100
    dt = 0.02
    t_max = 10.0
    time_steps = t_max / dt
    v = 1.0
    a = 0
    b = 1
    x_0 = 0.0
    x_N = 10.0
    x_stp_size = x_N / N 
    X = np.linspace(x_0, x_N, N)
    x_sw0 = int(4.5 / x_stp_size) 
    x_swk = int(5.5 / x_stp_size)
    dx = X[1] - X[0]
    u = np.zeros(N)
    u[x_sw0:x_swk] = b
    u[0] = a
    alpha = dt / dx
    print "|v|dt/dx: ", v * dt / dx
    u[-1] = a

    # Construct the Matricies to reflexct the lax formula 
    J = scpy.sparse.lil_matrix((N, N))
    J.setdiag(np.zeros(N) - 1, k=-1)
    J.setdiag(np.zeros(N) + 1, k=1)
    J[0,:] = np.zeros(N)
    J[-1,:] = np.zeros(N)
    J[0, 0] = 1
    J[-1, -1] = 1

    K = scpy.sparse.lil_matrix((N, N))
    K.setdiag(np.zeros(N) + 1, k=-1)
    K.setdiag(np.zeros(N) + 1, k=1)
    K[0,:] = np.zeros(N)
    K[-1, :] = np.zeros(N)
    K[0, 0] = 1
    K[-1, -1] = 1
    
    # Simulate where CFL is satisfied
    A = -v * 0.5 * alpha * J + 0.5 * K
    A[0,:] = np.zeros(N)
    A[-1, :] = np.zeros(N)
    A[0, 0] = 1
    A[-1, -1] = 1
    A = scpy.sparse.csr_matrix(scpy.linalg.inv(A.todense()))
    U = simulate(u, A, n=time_steps)
    animate1Dframes(X, U)
    
    return 0
    
    plt.title("Lax Method: dt = %0.2f, dx = %0.2f, and v = %0.2f" % (dt, dx, v))
    plt_0,  = plt.plot(X, U[0], "-ks")
    plt_32, = plt.plot(X, U[int(len(U)/32)], "-bs")
    plt_16, = plt.plot(X, U[int(len(U)/16)], "-gs")
    plt_8,  = plt.plot(X, U[int(len(U)/8)], "-rs")
    plt.legend([plt_0, plt_32, plt_16, plt_8], 
            [r"t = 0.0", r"t = %0.2f" % (int(len(U)/32) * dt), 
             r"t = %0.2f" % (int(len(U)/16) * dt),
             r"t = %0.2f" % (int(len(U)/8) * dt)]
            )

    plt.savefig("laxMethod_stable.png")
    plt.show()
    plt.close()
    #animate1Dframes(X, U)

    # Simulate without CFL being satisfied. 
    dt = 0.5
    time_steps = t_max / dt
    alpha = dt / dx
    print "|v|dt/dx: ", v * dt / dx

    A = -v * 0.5 * alpha * J + 0.5 * K
    A = scpy.sparse.csr_matrix(scpy.linalg.inv(A.todense()))
    U = simulate(u, A, n=time_steps)
    #animate1Dframes(X, U)
   
    plt.title("Lax Method: dt = %0.2f, dx = %0.2f, and v = %0.2f" % (dt, dx, v))
    plt_0,  = plt.plot(X, U[0], "-ks")
    plt_32, = plt.plot(X, U[1], "-bs")
    plt_16, = plt.plot(X, U[2], "-gs")
    plt_8,  = plt.plot(X, U[4], "-rs")
    plt.legend([plt_0, plt_32, plt_16, plt_8], 
            [r"t = 0.0", r"t = %0.2f" % (1 * dt), 
             r"t = %0.2f" % (2 * dt),
             r"t = %0.2f" % (3 * dt)]
            )
    plt.savefig("laxMethod_unstable.png")
    plt.show() 
    plt.close()

    # Leapfrog method
    dt = 0.01
    time_steps = t_max / dt
    alpha = dt / dx
    
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

def leap_frog(u_n, u_nm1=None, A=None):
    """ Leap frog method. Takes in u^n and u^n-1 and returns u^{n + 1}

        :param u_n: The current state of the system. 
        :param u_nm1: The next state of the system. If left as None
            then uses the lax method to calculate the next state. 
        :param A: Defaults to the identity matrix.

        :return: n^{n + 1} the next state of the system.
    """
    if not A:
        A=scpy.identity(len(u_n));
    


    return 0

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
