#!/usr/bin/env python
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as ax3d
import ode as ode
import finite_diff as fd
import utils as ut

v_t = 0.0
VEL = 1
POS = 0

def main():
    data = get_data()

    data_veloc = np.array([
            [data[0, index] for index in range(1, data.shape[1] - 1)],
            [fd.velocity(index, data) for index in range(1, data.shape[1] - 1)]
        ])
    data_accel = np.array([ 
            [data[0, index] for index in range(1, data.shape[1] - 1)],
            [fd.acceleration(index, data) for index in range(1, data.shape[1] - 1)]
        ])
     
    # Apply the guassian blend to the acceleration data. 
    window = gaussian_window(7, 1)
    data_agss = gaussian_blend_row(data_accel, 1, window)
    data_vgss = gaussian_blend_row(data_veloc, 1, window)
   
    # Find when the acceleration first reaches zero then find the terminal
    # velocity from the smoothed velocity data.
    term_acc = [first_x_intercept(data_agss), 0]
    term_vel = [term_acc[0], y_for_x(data_vgss, term_acc[0])]

    # Render a plot that compares the velocity to the acceleration.
    plt.plot(data_veloc[1, :], data_accel[1, :], "sb-") 
    plt.xlabel("Velocity (m/s)")
    plt.ylabel("Acceleration (m/s^2)")
    plt.title("Acceleration as a Function of Velocity")
    plt.show()

    global v_t
    v_t = term_vel[1]
    print "Estimated terminal velocity is %0.2f m/s" % v_t

    # Simulate the falling object. 
    sim_lin = [np.array([data[0, 0], data[1, 0], 0.0])] # linear drag forces simulation
    sim_qua = [np.array([data[0, 0], data[1, 0], 0.0])] # quadratic drag forces simulation

    #print sim_lin
    #print data
    
    dt = 0.001
    t_0 = data[0, 0]
    t_max = np.max(data[0,:]) # max seconds
    dt_steps = int(((t_max - t_0) + dt) / dt)
    for step in range(1, dt_steps + 1):
        t = step * dt + t_0
        tp1 = [t]
        sim_lin.append(np.append(tp1, ode.runge_kutta(df_linear, sim_lin[step - 1][1:3], dt, t)))
        sim_qua.append(np.append(tp1, ode.runge_kutta(df_quadratic, sim_qua[step - 1][1:3], dt, t)))

    sim_lin = np.array(sim_lin)
    sim_qua = np.array(sim_qua)

    # Compare real pos vs sim pos
    plt.subplot(211)
    ut.plt2dcmpr(sim_lin[:,0], sim_lin[:, 1], data[0, :], data[1, :], 
            ["r", "sb"],
            ["Simulated Position", "Real Position"], "Time (s)", 
            "Position (m)", "Simulated vs Real Position of Coffee Filter Using the Linear Drag Model"
            )
    plt.subplot(212)
    ut.plt2dcmpr(sim_qua[:,0], sim_qua[:, 1], data[0, :], data[1, :], 
            ["r", "sb"],
            ["Simulated Position", "Real Position"], "Time (s)", 
            "Position (m)", "Simulated vs Real Position of Coffee Filter Using the Quadratic Drag Model"
            )
    plt.show()

    # Compare real vel vs sim vel
    plt.subplot(211)
    ut.plt2dcmpr(sim_lin[:,0], sim_lin[:, 2], data_veloc[0, :], data_veloc[1, :], 
            ["r", "sb"],
            ["Simulated Velocity", "Real Velocity"], "Time (s)", 
            "Velocity (m/s)", "Simulated vs Real Velocity of Coffee Filter Using the Linear Drag Model"
            )
    plt.subplot(212)
    ut.plt2dcmpr(sim_qua[:,0], sim_qua[:, 2], data_veloc[0, :], data_veloc[1, :], 
            ["r", "sb"],
            ["Simulated Velocity", "Real Velocity"], "Time (s)", 
            "Velocity (m/s)", "Simulated vs Real Velocity of Coffee Filter Using the Quadratic Drag Model"
            )
    plt.show()


    # Compare real vel vs sim vel
    
    # Plot results 
    #plt.subplot(211) 
    #plt.plot(data_veloc[0,:], data_veloc[1,:], "go-")
    #plt.plot(data_vgss[0,:], data_vgss[1,:], "ro-")
    #plt.plot(term_vel[0], term_vel[1], "bs")
    #plt.ylabel(r"Velocity ($\frac{m}{s}$)")
    #plt.xlabel("Time (sec)")
    #plt.title("Velocity vs Time")
    #plt.subplot(212) 
    #plt.plot(data_accel[0,:], data_accel[1,:], "go-")
    #plt.plot(data_agss[0,:], data_agss[1,:], "ro-")
    #plt.plot(term_acc[0], term_acc[1], "bs")
    #plt.ylabel(r"Acceleration ($\frac{m}{s^2}$)")
    #plt.xlabel("Time (sec)")
    #plt.title("Acceleration vs Time")
    #plt.show()

# Drag forces linear ODE.
def df_linear(x, t, g=9.8):
    return np.array([x[VEL], -g * (1 - x[VEL] / v_t)])

# Drag forces quadratic ODE
def df_quadratic(x, t, g=9.8):
    return np.array([x[VEL], -g * (1 - (x[VEL] / v_t)**2)])

# Finds the first x intercept.
# Takes in a 2d numpy array assumes that the first row 
# is the x-values and that the second row is the y-vals.
# 
# Using linear interpolation returns the x value for which 
#       the y values first reach zero.
def first_x_intercept(data):
    for i in range(1, data.shape[1]):
        xc, yc = data[:, i]
        xp, yp = data[:, i - 1]
         
        # Check for intersection and perform interpolation.
        hit = False
        if yc >= 0 and yp < 0:
            hit = True
        elif yc <= 0 and yp > 0:
            hit = True

        if hit:
            m = (yc - yp) / (xc - xp) 
            b = yc - m * xc
            x = -b / m
            return x

    return None

# Finds a corresponding y value given an x value in a dataset. 
# assumes that the argument data is a 2d numpy array where the
# first row is the x values and the second row is the y values.
#
# Returns the associated y value using linear interpolation.
def y_for_x(data, x):
    for i in range(1, data.shape[1]):
        xc, yc = data[:, i]
        xp, yp = data[:, i - 1]

        if xp <= x and xc > x: 
            m = (yc - yp) / (xc - xp) 
            b = yc - m * xc
            return m * x + b

    return None

# First column is time. Second column is position.
def get_data():
    return np.array([ [.2055,.2302,.2550,.2797,.3045,.3292,.3539,.3786,.4033,.4280,
                       .4526,.4773,.5020,.5266,.5513,.5759,.6005,.6252,.6498,.6744,
                       .6990,.7236,.7482,.7728,.7974,.8220,.8466],
                      [.4188,.4164,.4128,.4082,.4026,.3958,.3878,.3802,.3708,.3609,
                       .3505,.3400,.3297,.3181,.3051,.2913,.2788,.2667,.2497,.2337,
                       .2175,.2008,.1846,.1696,.1566,.1393,.1263]])

# Used to calculate the gauss weights for the sliding window. 
def gauss_weights(r, sigma):
    return np.exp(-(r**2) / (2 * sigma**2))

# Takes in a numpy array and applies a guassian blend to the specified cell
def gaussian_blend(data, index, window):
    deg = int(len(window)/2)
    start = index - deg
    stop = index + deg 
    subset = np.array([data[i] for i in range(start, stop)])
    gaussian_sum = np.sum([subset[i] * window[i] for i in range(0, len(subset))])
    
    return gaussian_sum 

# Takes in a 2d numpy array and applied a gaussian blend to the specified row.
def gaussian_blend_row(data, row, window): 
    data_gauss = np.copy(data)
    wbounds = int(len(window)/2)
    for index in range(wbounds, data.shape[1] - wbounds):
        data_gauss[row,index] = gaussian_blend(data[row,:], index, window)
    
    i = wbounds
    j = data.shape[1] - wbounds
    return data_gauss[:, i:j]

# Creates a window for a gaussian blender. Assumes that wlen is an odd number.
def gaussian_window(wlen, sigma):
    rad = int(wlen / 2)
    distances = np.abs(np.array(range(-rad, rad + 1)))
    window = np.array([gauss_weights(r, sigma) for r in distances])
    return window/np.sum(window)

if __name__ == "__main__":
    main()

