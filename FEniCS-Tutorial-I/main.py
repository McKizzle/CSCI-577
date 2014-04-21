#!/usr/bin/env python
from dolfin import *
import numpy

def main():
    test2()

def test1():
    """
    FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
    Simplest example of computation and visualization with FEniCS.

    -Laplace(u) = f on the unit square.
    u = u0 on the boundary.
    u0 = u = 1 + x^2 + 2y^2, f = -6.
    """

    # Create mesh and define function space
    mesh = UnitSquare(6, 4)
    #mesh = UnitCube(6, 4, 5)
    V = FunctionSpace(mesh, 'Lagrange', 1)

    # Define boundary conditions
    u0 = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]')
   
    bc = DirichletBC(V, u0, u0_boundary)

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = inner(nabla_grad(u), nabla_grad(v))*dx
    L = f*v*dx

    # Compute solution
    u = Function(V)
    solve(a == L, u, bc)

    # Plot solution and mesh
    plot(u)
    plot(mesh)

    # Dump solution to file in VTK format
    file = File('poisson.pvd')
    file << u

    # Hold plot
    interactive()
    return 0

def test2():
    T = 10.0 # tension
    A = 1.0 # pressure amplitude
    R = 0.3 # radius of domain
    theta = 0.2
    x0 = 0.6 * R * cos(theta)
    y0 = 0.6 * R * sin(theta)
    sigma = 50.0
    n = 40 
    mesh = UnitCircle(n)
    V = FunctionSpace(mesh, 'Lagrange', 1)

    bc = DirichletBC(V, Constant(0.0), boundary)
    
    # Define variational Problem
    w = TrialFunction(V)
    v = TestFunction(V)
    a = inner(nabla_grad(w), nabla_grad(v)) * dx
    f = Expression('4*exp(-0.5*(pow((R*x[0] - x0)/sigma, 2)) '
                   '      -0.5*(pow((R*x[1] - y0)/sigma, 2)))',
                   R=R, x0=x0, y0=y0, sigma=sigma)
    L = f * v * dx

    # Compute Solution
    w = Function(V)
    problem = LinearVariationalProblem(a, L, w, bc)
    solver = LinearVariationalSolver(problem)
    solver.parameters['linear_solver'] = 'cg'
    solver.parameters['preconditioner'] = 'ilu'
    solver.solve()

    # Plot scaled solution, mesh and pressure
    #plot(mesh, title='Mesh over scaled domain')
    #plot(w, title='Scaled deflection')
    #f = interpolate(f, V)
    #plot(f, title='Scaled pressure')

    # Find maximum real deflection
    max_w = w.vector().array().max()
    max_D = A*max_w/(8*pi*sigma*T)
    print 'Maximum real deflection is', max_D

    # Verification for "flat" pressure (large sigma)
    if sigma >= 50:
        w_e = Expression("1 - x[0]*x[0] - x[1]*x[1]")
        w_e = interpolate(w_e, V)
        dev = numpy.abs(w_e.vector().array() - w.vector().array()).max()
        print 'sigma=%g: max deviation=%e' % (sigma, dev)

    # Should be at the end
    #interactive()

    V_g = VectorFunctionSpace(mesh, 'Lagrange', 1)
    w = TrialFunction(V_g)
    v = TestFunction(V_g)
    u = 

    a = inner(w, v) * dx
    L = inner(grad(v), v) * dx
    grad_v = Function(V_g)
    solve(a == L, grad_u)

    plot(grad_u, title='grad(u)')
    interactive()

    return 0


def boundary(x, on_boundary):
   return on_boundary

def u0_boundary(x, on_boundary):
   return on_boundary


def animate1Dframes(x, data):
    import time as tm
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
