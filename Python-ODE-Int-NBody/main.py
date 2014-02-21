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
G = 6.67384 * 10 ** -11
GM = 4 * mt.pi ** 2
M_EAR = 0.001
M_JUP= 0.04

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
    integrator = ode.euler_richardson

    print "Simulating nbody..."
    sim_d = sys.simulate([init_cond], dt, t_0, t_max, ode_sys, integrator) # System simulation
    print "Done"
    print "Calculating the angular momentum and the energy..."
    sim_am = [angular_momentum(sim_d[0, 1:sim_d.shape[1]])]
    sim_en = [energy(sim_d[0, 1:sim_d.shape[1]])]
    for i in range(1, sim_d.shape[0]):
        sim_am.append(angular_momentum(sim_d[i, 1:sim_d.shape[1]]))
        sim_en.append(energy(sim_d[i, 1:sim_d.shape[1]]))
    print "Done"
 
    # Plot the orbits, angular momentum, and energy.
    plt.subplot(311)
    plt.axis('equal')
    ut.plt2dcmpr(sim_d[:, XE + 1], sim_d[:, YE + 1], sim_d[:, XJ + 1], sim_d[:, YJ + 1],
        ["r", "b"], ["Earth", "Jupiter"], "X", "Y", "Jupiter and Earth Orbit Simulation Results")
    plt.subplot(312)
    plt.plot(sim_d[:,0], sim_am, "-")
    plt.xlabel("time (years)")
    plt.ylabel(r"% $\Delta$ E / M")
    plt.subplot(313)
    plt.plot(sim_d[:,0], sim_en, "-")
    plt.xlabel("time (years)")
    plt.ylabel(r"% $\Delta$ L / M")
    plt.show()

def energy(x, GM=GM, G=G, m_e=M_EAR, m_j=M_JUP):
    r_se = np.sqrt(np.sum(x[XE:(YE + 1)]**2))                           # r_se between sun and earth.
    r_ej = np.sqrt(np.sum((x[XJ:(YJ + 1)] - x[XE:(YE + 1)])**2))        # r_ej between earth and jupiter.
    r_sj = np.sqrt(np.sum(x[XJ:(YJ + 1)]**2))                           # r_sj between sun and jupiter.

    E_x =  1.0 / 2.0 * m_e * x[VXE] + 1.0/2.0 * m_j * x[VXJ]
    E_x = E_x - GM * m_e / r_se - GM * m_j / r_sj - G * m_j * m_e / r_ej

    E_y =  1.0 / 2.0 * m_e * x[VYE] + 1.0/2.0 * m_j * x[VYJ]
    E_y = E_y - GM * m_e / r_se - GM * m_j / r_sj - G * m_j * m_e / r_ej

    return np.array([E_x + E_y])

def angular_momentum(x, GM=GM, m_e=M_EAR, m_j=M_JUP):
    L_e = m_e*(x[XE] * x[VXE] - x[YE] * x[VYE])
    
    return np.array([L_e + m_e * (x[XJ] * x[VXJ] - x[YJ] * x[VYJ])])
    


def planets(x, t, GM=GM, m_e=M_EAR, m_j=M_JUP):
    """ Planet ODE system """
    dxdt = np.zeros(np.size(x))

    #print "x:\t", x 
    r_se = np.sqrt(np.sum(x[XE:(YE + 1)]**2))                           # r_se between sun and earth.
    r_ej = np.sqrt(np.sum((x[XJ:(YJ + 1)] - x[XE:(YE + 1)])**2))        # r_ej between earth and jupiter.
    r_sj = np.sqrt(np.sum(x[XJ:(YJ + 1)]**2))                           # r_sj between sun and jupiter.

    #print "r_se:\t%0.4f"% r_se
    #print "r_ej:\t%0.4f"% r_ej
    #print "r_sj:\t%0.4f"% r_sj
     
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

    #print "DxDt:\t", dxdt

    #exit()

    return dxdt

if __name__ == '__main__':
    main()

