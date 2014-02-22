#!/usr/bin/env python
#import scipy.integrate as spint
import scipy.integrate as scint
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
G = 6.67384e-11
GM = 4 * mt.pi ** 2
#M = 1.989e30
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
    dt = 0.001
    t_0 = 0.0
    t_max = 100
    ode_sys = planets
    integrator = ode.euler_richardson

    print "Simulating nbody..."
    intgr = scint.ode(planets_scipy)
    intgr.set_integrator('vode', atol=0.00001, rtol=0.00001)
    intgr.set_initial_value(init_cond, t_0)
    times = np.arange(t_0, t_max, dt)

    sim_d = [init_cond]
    for t in times[1:]:
        intgr.integrate(intgr.t+dt)
        sim_d.append(intgr.y)
    sim_d = np.array(sim_d)
    #times = np.arange(t_0, t_max, dt)
    #sim_d = scint.odeint(planets_scipy, init_cond, times)
    #sim_d = sys.simulate([init_cond], dt, t_0, t_max, ode_sys, integrator) # System simulation
    print "Done"

    print "Calculating the angular momentum and the energy..."
    print "TODO"
    print "Done"
 
    # Plot the orbits, angular momentum, and energy.
    #plt.subplot(311)
    plt.axis('equal')
    ut.plt2dcmpr(sim_d[:, XE + 1], sim_d[:, YE + 1], sim_d[:, XJ + 1], sim_d[:, YJ + 1],
        ["r", "b"], ["Earth", "Jupiter"], "X", "Y", "Jupiter and Earth Orbit Simulation Results")
    #plt.subplot(312)
    #plt.plot(sim_d[:,0], sim_am, "-")
    #plt.xlabel("time (years)")
    #plt.ylabel(r"% $\Delta$ E / M")
    #plt.subplot(313)
    #plt.plot(sim_d[:,0], sim_en, "-")
    #plt.xlabel("time (years)")
    #plt.ylabel(r"% $\Delta$ L / M")
    plt.show() 
 
def planets_scipy(t, y):
    return planets(y, t)

def planets2(x, t, GM=GM, m_eM=M_EAR, m_jM=M_JUP):
    dxdt = np.zeros(np.size(x))

    r_1     = mt.sqrt(x[0]**2 + x[1]**2)
    r_21    = mt.sqrt((x[4] - x[0])**2 + (x[5] - x[1])**2)
    r_2     = mt.sqrt(x[4]**2 + x[5]**2)

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

#def percent_delta(a_0, a_n):
#    if a_0:
#        return np.array((a_n - a_0) / a_0) * 100
#    else:
#        return 0
#
#def energy(x, GM=GM, G=G, m_e=M_EAR, m_j=M_JUP):
#    r_se = np.sqrt(np.sum(x[XE:(YE + 1)]**2))                           # r_se between sun and earth.
#    r_ej = np.sqrt(np.sum((x[XJ:(YJ + 1)] - x[XE:(YE + 1)])**2))        # r_ej between earth and jupiter.
#    r_sj = np.sqrt(np.sum(x[XJ:(YJ + 1)]**2))                           # r_sj between sun and jupiter.
#
#    v_e = np.sqrt(np.sum(x[VXE:VYE + 1]**2))**2
#    v_j = np.sqrt(np.sum(x[VXJ:VYJ + 1]**2))**2
#    
#    E = 1.0 / 2.0 * m_e * v_e**2
#    E = E + 1.0 / 2.0 * m_j * v_j**2
#    E = E - GM * m_e / r_se - GM * m_j / r_sj - GM * m_j * m_e / r_ej
#    return np.array(E)
#
#
#def angular_momentum(x, GM=GM, m_e=M_EAR, m_j=M_JUP):
#    L_e = m_e*(x[XE] * x[VXE] - x[YE] * x[VYE]) 
#    L_j = m_j * (x[XJ] * x[VXJ] - x[YJ] * x[VYJ])
#    return np.array([L_e + L_j])
