#!/usr/bin/env python
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
import ode as ode



""" Runs a harmonic simulator
   
   Uses the following formulas to simulate a harmonic simulator. 

   F = -kx
   x(t) = x_0 cos(\sqrt{\frac{k}{m}} t)

   Energy conservation is then checked with 

   E= \frac{1}{2} k x^2 + \frac{1}{2} m v^2
"""
def main():
    x_0 = 1.0 # starting position
    v_0 = 0.0 # starting velocity
    t_0 = 0.0

    use_predictor_corrector = False
    intgtr = ode.euler_richardson # The default funtion to use. (ode.euler, ode.euler_richardson, runge_kutta)

    sim_d = [np.array([None, None, None])]
    sim_d.append(np.array([t_0, v_0, x_0]))
    
    sim_d_e = [None]
    sim_d_e.append(sho_anal_ener(x_0, v_0))
    
    sim_a = [np.array([None, None, None])]
    sim_a.append(np.array([t_0, sho_anal_vel(x_0, t_0), sho_anal_pos(x_0, t_0)]))
    
    sim_a_e = [None]
    sim_a_e.append(sho_anal_ener(x_0, v_0))

    dt = 0.01
    t_max = int( 5 * 2 * np.pi + 1) # only go out to five cycles. 
    dt_steps = int(t_max / dt) + 1
    for step in range(2, dt_steps):
        t = step * dt + t_0
        tp1 = np.array([t])
        
        # Integrate and calcuate the energy in the integrated system.
        if use_predictor_corrector:
            sim_d.append(np.append(tp1, ode.predictor_corrector(sho, sim_d[step - 1][1:3], dt, t, intgtr, sim_d[step - 2][1:3])))
        else:
            sim_d.append(np.append(tp1, intgtr(sho, sim_d[step - 1][1:3], dt, t)))  
        sim_d_e.append(sho_anal_ener(sim_d[step][2], sim_d[step][1]))
    
        # Calculate analytic values. 
        sim_a.append(np.append(tp1, [sho_anal_vel(x_0, t), sho_anal_pos(x_0, t)]))
        sim_a_e.append(sho_anal_ener(sim_a[step][2], sim_a[step][1])) 

    # Calclate the error at the fifth cycle. 
    t_fcyl = 5 * 2 * np.pi
    fcyl_x_d = sim_d[dt_steps - 1][2]
    fcyl_x_a = sim_a[dt_steps - 1][2]
    print "The error is: %0.2f%% percent" % (error_calc([fcyl_x_d, fcyl_x_a]) * 100)

    # Massage the data
    sim_d.pop(0)
    sim_d_e.pop(0)
    sim_a.pop(0)
    sim_a_e.pop(0)
    sim_d = np.array(sim_d)
    sim_a = np.array(sim_a)

    # Plot stuff
    plt2dcmpr(sim_d[:,0], sim_d[:,2], sim_a[:,0], sim_a[:,2], 
            ['Body Position', 'Body Position Analytical'], "Time (s)", 
            "Position (m)", "Position vs Time of Simple Harmonic Oscillator")
    plt.show()

    # Plot the energies of the two systems. 
    sim_d_e_plt, = plt.plot(sim_d[:,0], sim_d_e, "b")
#    sim_a_e_plt, = plt.plot(sim_a[:,0], sim_a_e, "b")
    #plt.legend([sim_d_e_plt, sim_a_e_plt], ['Energy Euler', 'Energy Analytical'])
    #plt.legend([sim_d_e_plt], ['Energy'])
    plt.show()

    # Calculate the diffence between the two systems.
    e_diff = abs(sim_d_e[len(sim_d_e) - 1] - sho_anal_ener(sho_anal_pos(x_0, t_max), sho_anal_vel(x_0, t_max)))
    print "Energy difference is %0.2fJ" % e_diff

# Plot two systems together for comparing.  
def plt2dcmpr(X1, Y1, X2, Y2, legend_labels, xlab, ylab, title):
    plt1, = plt.plot(X1, Y1, "b")
    plt2, = plt.plot(X2, Y2, "r")
    plt.legend([plt1, plt2], legend_labels)
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.title(title)

# Calculate the error between two points. 
def error_calc(x): 
    x_max = max(x)
    x_min = min(x)
    return (max(x) - min(x))/max(x)

# Simple harmonic oscilator slope function.
def sho(x, t, k=1.0, m=1.0):
    return np.array([ -k * x[1] / m, x[0]])

# Simple harmonic oscilator analytic position function.  
def sho_anal_pos(x_0, t, k=1.0, m=1.0):
    return x_0 * np.cos((k/m)**0.5 * t)

# Simple harmonic oscilator analytic velocity function. 
def sho_anal_vel(x_0, t, k=1.0, m=1.0):
    return -0.5 * x_0 * np.sin((k/m)**0.5 * t)

# Simple harmonic oscilator energy function. 
def sho_anal_ener(x, v, k=1.0, m=1.0):
    return 0.5 * (k * x**2 + m * v**2)

if __name__ == '__main__':
    main()



