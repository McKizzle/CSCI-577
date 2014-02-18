import matplotlib.pyplot as plt
import numpy as np
import ode as ode
import system as sys
import math as math

def run():
    sim_e = sys.simulate([np.array([0])], 0.1, 0, 1.5, erf_ode, ode.euler)
    sim_er = sys.simulate([np.array([0])], 0.1, 0, 1.5, erf_ode, ode.euler_richardson)
    sim_rk4 = sys.simulate([np.array([0])], 0.1, 0, 1.5, erf_ode, ode.runge_kutta) 
    anl_data = sys.simulate([np.array([0])], 0.1, 0, 1.5, erf, ode.analytical)

    # TODO Compute the sum squared error of your erf and pythons erf.

    plt.plot(sim_e[:, 0], sim_e[:, 1])
    plt.plot(sim_er[:, 0], sim_er[:, 1])
    plt.plot(sim_rk4[:, 0], sim_rk4[:, 1])
    plt.plot(anl_data[:, 0], anl_data[:, 1])
    plt.show()

def erf_ode(x, t):
    """ ODE system for erf(t) """
    return np.array([(2 / np.sqrt(np.pi)) * np.exp(-t**2)])

def erf(x, dt, t):
    """ Simply a wrapper function for erf(t) """
    return np.array([math.erf(t)])

