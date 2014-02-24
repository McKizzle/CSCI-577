#!/usr/bin/env python
import scipy.integrate as scint
import numpy as np
import system as sys
import math as mt
import ode as ode
import utils as ut
import matplotlib.pyplot as plt

# Constants
G = 6.67384e-11
GM = 4 * mt.pi ** 2
M_EAR = 0.001
M_JUP= 0.04

def main(): 
    dt = 0.001
    t_0 = 0.0
    t_max = 100

    # Simulate the two planets using the scipy integrator
    init_cond = two_planets_initial_cond()
    intgr = scint.ode(two_planets)
    intgr.set_integrator('dopri5', method = 'adams', atol=0.000001, rtol=0.000001)
    intgr.set_initial_value(init_cond, t_0)
    times = np.arange(t_0, t_max, dt)
    sim_d = [init_cond]
    for t in times[1:]:
        intgr.integrate(intgr.t + dt)
        sim_d.append(intgr.y)
    sim_d = np.array(sim_d)
    
    # Plot the orbits
    plt.subplots(nrows=3, ncols=1)
    plt.tight_layout()
    plt.subplot(311)
    plt.axis("equal")
    ut.plt2dcmpr(sim_d[:, 0], sim_d[:, 1], sim_d[:, 4], sim_d[:, 5],
        ["r", "b"], ["Earth", "Jupiter"], "X (AU)", "Y (AU)", "Jupiter and Earth Orbit Simulation Results (dt = %0.3f)" % dt)

    # Now simulate the system using the class built integrator.
    sim_c = sys.simulate([two_planets_initial_cond()], dt, t_0, t_max, two_planets_wrapper, ode.runge_kutta)
    #ut.plt2dcmpr(sim_c[:, 0], sim_c[:, 1], sim_c[:, 4], sim_c[:, 5],
    #    ["r", "b"], ["Earth", "Jupiter"], "X (AU)", "Y (AU)", "Jupiter and Earth Orbit Simulation Results")

    # Calcualte the percent delta change for the energy and angular momentum.
    sim_perde = []
    sim_perdam = []
    for i in range(0, sim_d.shape[0]):
        sim_perde.append(per_diff(two_planets_energy(sim_d[i]), two_planets_energy(sim_d[0])))
        sim_perdam.append(per_diff(two_planets_angular_momentum(sim_d[i]), two_planets_angular_momentum(sim_d[0])))
    
    # Plot 
    sim_perde = np.array(sim_perde)
    sim_pardam = np.array(sim_perdam)
    plt.subplot(312)
    plt.plot(times, sim_perde)
    plt.xlabel("time (years)")
    plt.ylabel(r"% $\Delta$ E / M")
    plt.subplot(313)
    plt.plot(times, sim_perdam)
    plt.xlabel("time (years)")
    plt.ylabel(r"% $\Delta$ L / M")
    plt.savefig("two_planets.png", bbox_inches="tight", dpi=120)
    plt.close()


def per_diff(x_1, x_2):
    return (x_2 - x_1) / x_1 * 100

def two_planets_energy(x, GM=GM, m_1overM=M_EAR, m_2overM=M_JUP):
    # Build the vectors. 
    vr_1  = np.array([x[0], x[1]])
    vr_21 = np.array([x[4] - x[0], x[5] - x[1]])
    vr_2  = np.array([x[4], x[5]])
    # Magnitudes
    r_1     = np.linalg.norm(vr_1) 
    r_21    = np.linalg.norm(vr_21)
    r_2     = np.linalg.norm(vr_2)
    # Construct the Velocity magnitudes.
    v_1 = np.linalg.norm(x[2:4])
    v_2 = np.linalg.norm(x[6:8])

    EoverM = 1.0/2.0*m_1overM*v_1**2 + 1.0/2.0*m_2overM*v_2**2 \
            - GM * m_1overM / r_1 - GM * m_2overM / r_2 \
            - GM * m_1overM * m_2overM / r_21
    return EoverM

def two_planets_angular_momentum(x, GM=GM, m_1overM=M_EAR, m_2overM=M_JUP):
    # Build the vectors. 
    L_1overM = m_1overM*(x[0]*x[3] - x[1]*x[2])
    L_2overM = m_1overM*(x[4]*x[7] - x[5]*x[6])
    return L_1overM + L_2overM

def two_planets_wrapper(x, t):
    return two_planets(t, x)
 
 
def two_planets(t, x, GM=GM, m_1overM=M_EAR, m_2overM=M_JUP):
    dxdt = np.zeros(np.size(x))

    # Build the vectors. 
    vr_1  = np.array([x[0], x[1]])
    vr_21 = np.array([x[4] - x[0], x[5] - x[1]])
    vr_2  = np.array([x[4], x[5]])
    # Magnitudes
    r_1     = np.linalg.norm(vr_1) 
    r_21    = np.linalg.norm(vr_21)
    r_2     = np.linalg.norm(vr_2)

    dvdt_1 = -GM / r_1**3 * vr_1 + GM * m_2overM / r_21**3 * vr_21
    dvdt_2 = -GM / r_2**3 * vr_2 - GM * m_1overM / r_21**3 * vr_21
    
    dxdt[0] = x[2]
    dxdt[1] = x[3]
    dxdt[2] = dvdt_1[0]
    dxdt[3] = dvdt_1[1]
    dxdt[4] = x[6]
    dxdt[5] = x[7]
    dxdt[6] = dvdt_2[0]
    dxdt[7] = dvdt_2[1]

    return dxdt

def two_planets_initial_cond(GM=GM):
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

    return np.array([x_1, y_1, v_x1, v_y1, x_2, y_2, v_x2, v_y2]) 

if __name__ == '__main__':
    main()


