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
    x_1 = 2.52
    y_1 = 0.0
    v_x1 = 0.0
    v_y1 = mt.sqrt(GM / x_1)

    # Jupiter Initial Conditions
    x_2 = 5.24
    y_2 = 0.0
    v_x2 = 0.0
    v_y2 = mt.sqrt(GM / x_2) 
    init_cond = np.array([x_1, y_1, v_x1, v_x1, x_2, y_2, v_x2, v_y2]) 
    
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
    
    r_se = np.sqrt(np.sum(x[XE:(YE + 1)]**2))                      # r_1 between sun and earth.
    r_ej = np.sqrt(np.sum((x[XJ:(YJ + 1)] - x[XE:(YE + 1)])**2))   # r_21 between earth and jupiter.
    r_sj = np.sqrt(np.sum(x[XJ:(YJ + 1)]**2))                      # r_2 between sun and jupiter.

    #print r_se, r_ej, r_sj
    #exit()
    
    # Earth
    dxdt[XE]    = x[VXE]
    dxdt[YE]    = x[VYE]
    dxdt[VXE]   = (GM / r_se**3 * x[XE]) + (GM * m_e / r_ej**3 * (x[XJ] - x[XE]))
    dxdt[VYE]   = (-1 * GM / r_se**3 * x[YE]) + (GM * m_j / r_ej**3 * (x[YJ] - x[YE]))
    # Jupiter
    dxdt[XJ]    = x[VXJ]
    dxdt[YJ]    = x[VYJ]
    dxdt[VXJ]   = (GM / r_sj**3 * x[XJ]) + (GM * m_j / r_ej**3 * (x[XE] - x[XJ]))
    dxdt[VYJ]   = (-1 * GM / r_sj**3 * x[YJ]) + (GM * m_j / r_ej**3 * (x[YE] - x[YJ]))

    #print dxdt

    return dxdt

if __name__ == '__main__':
    main()

