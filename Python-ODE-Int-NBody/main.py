#!/usr/bin/env python
import numpy as np
import system as sys
import math as mt
import ode as ode
import utils as ut
import matplotlib.pyplot as plt

XE = 0
YE = 1
VXE= 2
VYE= 3
XJ = 4
YJ = 5
VXJ= 6
VYJ= 7

# Constants
GM = 4 * mt.pi ** 2
M_EAR = 0.001 * GM
M_JUP= 0.04 * GM

def main():
    # Earth initial Cond
    x_e = 2.52
    y_e = 0.0
    v_xe = 0.0
    v_ye = mt.sqrt(GM / x_e)

    # Jupiter Initial Conditions
    x_j = 5.24
    y_j = 0.0
    v_xj = 0.0
    v_yj = mt.sqrt(GM / x_j) 


    init_cond = np.array([x_e, y_e, v_xe, v_ye, x_j, y_j, v_xj, v_yj]) 
    
    # System Settings
    dt = 0.01
    t_0 = 0.0
    t_max = 100
    ode_sys = planets
    integrator = ode.runge_kutta

    print "Simulating NBody..."
    sim_d = sys.simulate([init_cond], dt, t_0, t_max, ode_sys, integrator)
    print "Done"
    
    ut.plt2dcmpr(sim_d[:, XE + 1], sim_d[:, YE + 1], sim_d[:, XJ + 1], sim_d[:, YJ + 1],
        ["r", "b"], ["Earth", "Jupiter"], "X", "Y", "Jupiter and Earth Orbit Simulation Results")
    plt.show()


def planets(x, t, GM=GM, m_e=M_EAR, m_j=M_JUP):
    """ Planet ODE system """
    dxdt = np.zeros(np.size(x))
     
    r_se = np.sqrt(np.sum(x[XE:(YE + 1)]**2))                           # r_se between sun and earth.
    r_ej = np.sqrt(np.sum(x[XJ:(YJ + 1)]**2 + x[XE:(YE + 1)]**2))       # r_ej between earth and jupiter.
    r_sj = np.sqrt(np.sum(x[XJ:(YJ + 1)]**2))                           # r_sj between sun and jupiter.
     
    # Earth
    dxdt[XE]    = x[VXE]
    dxdt[YE]    = x[VYE]
    dxdt[VXE]   = -GM / r_se**3 * x[XE] + GM * m_j / r_ej**3 * (x[XJ] - x[XE])
    dxdt[VYE]   = -GM / r_se**3 * x[YE] + GM * m_j / r_ej**3 * (x[YJ] - x[YE]) 
    # Jupiter
    dxdt[XJ]    = x[VXJ]
    dxdt[YJ]    = x[VYJ]
    dxdt[VXJ]   = -GM / r_sj**3 * x[XJ] - GM * m_e / r_ej**3 * (x[XE] - x[XJ])
    dxdt[VYJ]   = -GM / r_sj**3 * x[YJ] - GM * m_e / r_ej**3 * (x[YE] - x[YJ])

    #print dxdt

    return dxdt

if __name__ == '__main__':
    main()

