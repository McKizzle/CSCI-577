#!/usr/bin/env python
import scipy.integrate as scint
import numpy as np
import system as sys
import math as mt
import ode as ode
import utils as ut
import matplotlib.pyplot as plt

# Constants
GM = 4 * mt.pi ** 2
M_EAR = 0.001
M_JUP= 0.04

STAR1_POS = np.array([0.0, 0.0])
STAR2_POS = np.array([2.0, 0.0])

def main(): 
    dt = 0.001
    t_0 = 0.0
    t_max = 100

    # Simulate the two planets. 
    sim_d = simulate(binary_initial_cond1(), dt, t_0, t_max) 
    plt.subplots(nrows=2, ncols=1)
    plt.tight_layout()
    plt.subplot(211)
    plt.axis("equal")
    plt.title("Binary Star System with Single Planet (dt = %0.3f)" % dt)
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    plt.plot(sim_d[:, 0], sim_d[:, 1])

    sim_d = simulate(binary_initial_cond2(), dt, t_0, t_max) 
    plt.subplot(212)
    plt.axis("equal")
    plt.plot(sim_d[:, 0], sim_d[:, 1])
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    plt.savefig("kepler_orbits.png", bbox_inches="tight")

def simulate(init_cond, dt, t_0, t_max):
    # Simulate the two planets. 
    intgr = scint.ode(binary_system)
    intgr.set_integrator('vode', method='adams')
    intgr.set_initial_value(init_cond, t_0)
    times = np.arange(t_0, t_max, dt)
    sim_d = [init_cond]
    for t in times[1:]:
        intgr.integrate(intgr.t + dt)
        sim_d.append(intgr.y)
    sim_d = np.array(sim_d)

    return sim_d

def binary_system(t, x, s1_p=STAR1_POS, s2_p=STAR2_POS, GM=GM): 
    dxdt = np.zeros(np.size(x))
    
    vr_1 = x[0:2] - s1_p
    vr_23 = s2_p - x[0:2]

    r_1 = np.linalg.norm(vr_1)
    r_23 = np.linalg.norm(vr_23)

    dvdt = -GM / r_1**3 * vr_1 + GM / r_23**3 * vr_23 

    dxdt[0] = x[2]
    dxdt[1] = x[3]
    dxdt[2] = dvdt[0]
    dxdt[3] = dvdt[1]

    return dxdt

def binary_initial_cond1(GM=GM):
    # planet initial conditions
    x = 1.1
    y = 1.0
    v_x = -mt.sqrt(GM / y)
    v_y = mt.sqrt(GM / x)

    return np.array([x, y, v_x, v_y]) 


def binary_initial_cond2(GM=GM):
    # planet initial conditions
    x = 1.1
    y = 1.0
    v_x = mt.sqrt(GM / y)
    v_y = 1.05 * mt.sqrt(GM / x)

    return np.array([x, y, v_x, v_y]) 

if __name__ == '__main__':
    main()


