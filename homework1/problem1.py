import matplotlib.pyplot as plt
import numpy as np
import ode as ode
import system as sys
import math as math
import utils as ut

def run(): 
    dt = 0.1
    t_0 = 0.0
    t_max = 1.5
    odes = [ode.euler, ode.euler_richardson, ode.runge_kutta] 
    ode_ldgns = ["Euler", "Euler Richardson", "Runge Kutta"]
    plt_files = ["erf_euler", "erf_euler-richardson", "erf_rk4"]

    anl_data = sys.simulate([np.array([0])], 0.1, 0, 1.5, erf, ode.analytical)
    for i in range(0, len(odes)):
        sim_d = sys.simulate([np.array([0])], 0.1, 0, 1.5, erf_ode, odes[i])
        
        # compute the sum squared error
        sum_sqrs = np.sum((sim_d[:, 1] - anl_data[:, 1])**2)
        
        ut.plt2dcmpr(
                sim_d[:,0], sim_d[:,1], anl_data[:,0], anl_data[:,1],
                ["bs-", "ro:"], [ode_ldgns[i], "Python erf(x)"],
                "t", "x", "Pyrthon erf(x) vs ODE erf(x) %s Integration" % ode_ldgns[i]
                )
        plt.text(0.8, 0.4, "Sum Squares Error = %0.4f" % sum_sqrs)
        plt.savefig("%s.png" % plt_files[i], bbox_inches='tight')
        plt.close()

def erf_ode(x, t):
    """ ODE system for erf(t) """
    return np.array([(2 / np.sqrt(np.pi)) * np.exp(-t**2)])

def erf(x, dt, t):
    """ Simply a wrapper function for erf(t) """
    return np.array([math.erf(t)])

