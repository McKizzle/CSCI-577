#!/usr/bin/env python
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as ax3d
import ode as ode
import finite_diff as fd

# First column is time. Second column is position.
def get_data():
    return np.array([ [.2055,.2302,.2550,.2797,.3045,.3292,.3539,.3786,.4033,.4280,
                       .4526,.4773,.5020,.5266,.5513,.5759,.6005,.6252,.6498,.6744,
                       .6990,.7236,.7482,.7728,.7974,.8220,.8466],
                      [.4188,.4164,.4128,.4082,.4026,.3958,.3878,.3802,.3708,.3609,
                       .3505,.3400,.3297,.3181,.3051,.2913,.2788,.2667,.2497,.2337,
                       .2175,.2008,.1846,.1696,.1566,.1393,.1263]])

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
    

    #plt.subplot(221)
    #plt.plot(data[0,:], data[1,:], "ro-")
    #plt.ylabel("Position (m)")
    #plt.xlabel("Time (sec)")
    #plt.title("Position vs Time")
    #plt.subplot(222)
    #plt.plot(data_veloc[0,:], data_veloc[1,:], "go-")
    #plt.ylabel(r"Velocity ($\frac{m}{s}$)")
    #plt.xlabel("Time (sec)")
    #plt.title("Velocity vs Time")
    #plt.subplot(223)
    #plt.plot(data_accel[0,:], data_accel[1,:], "ko-")
    #plt.ylabel(r"Acceleration ($\frac{m}{s^2}$)")
    #plt.xlabel("Time (sec)")
    #plt.title("Acceleration vs Time")
    #plt.subplot(224)
    #plt.plot(data_accel[1,:], data_veloc[1,:], "ko-")
    #plt.ylabel(r"Velocity ($\frac{m}{s}$)")
    #plt.ylabel(r"Acceleration ($\frac{m}{s^2}$)")
    #plt.title("Acceleration vs Velocity")
    #plt.show()
    
    # Calculate the window wieghts for the guassian blend. 
    wwghts = np.array([gauss_weights(r, 1) for r in [2, 1, 0, 1, 2]])
    wwghts = wwghts / np.sum(wwghts)
    bnds = int(len(wwghts)/2)
    
    # Apply guassian blend to the acceleration and velocity data. 
    data_vgss = np.array([
        [data_veloc[0, index] for index in range(bnds, data_veloc.shape[1] - bnds)],
        [gaussian_blend(index, data_veloc[1,:], wwghts) for index in range(bnds, data_veloc.shape[1] - bnds)]
    ])

    data_agss = np.array([
        [data_accel[0, index] for index in range(bnds, data_accel.shape[1] - bnds)],
        [gaussian_blend(index, data_accel[1,:], wwghts) for index in range(bnds, data_accel.shape[1] - bnds)]
    ])

    # Determine the time the acceleration first reaches zero. 
    tv_time = 0.0 # the point the object stops accelerating (terminal velocity)
    tv_index = 0.0  
    for i in range(0, data_agss.shape[1]):
        if data_agss[1, i] >= 0 and i > 0:
            tv_time = data_agss[0,i]
            break
        elif data_agss[1, i] >= 0:
            tv_time = data_agss[0,i]
            break

    print zerodex
    
    #plt.subplot(221) 
    #plt.plot(data_veloc[0,:], data_veloc[1,:], "go-")
    #plt.plot(data_vgss[0,:], data_vgss[1,:], "ro-")
    #plt.ylabel(r"Velocity ($\frac{m}{s}$)")
    #plt.xlabel("Time (sec)")
    #plt.title("Velocity vs Time")
    #plt.subplot(222) 
    #plt.plot(data_accel[0,:], data_accel[1,:], "go-")
    #plt.plot(data_agss[0,:], data_agss[1,:], "ro-")
    #plt.ylabel(r"Acceleration ($\frac{m}{s^2}$)")
    #plt.xlabel("Time (sec)")
    #plt.title("Acceleration vs Time")
    #plt.subplot(224)
    #plt.plot(data_vgss[1,:], data_agss[1,:], "ko-")
    #plt.xlabel(r"Velocity ($\frac{m}{s}$)")
    #plt.ylabel(r"Acceleration ($\frac{m}{s^2}$)")
    #plt.title("Acceleration vs Velocity")
    #plt.show()
    #plt.show()


    #gaussian_blend(2, data[1,:], window_weights)

    #ax = ax3d(plt.gcf())
    #ax.plot(data_veloc[0,:], data_accel[1,:], data_veloc[1,:])
    #plt.show()

# Used to calculate the gauss weights for the sliding window. 
def gauss_weights(r, sigma):
    return np.exp(-(r**2) / (2 * sigma**2))

# Assumes that the data is an numpy array.
def gaussian_blend(index, data, window):
    deg = int(len(window)/2)
    start = index - deg
    stop = index + deg 
    subset = [data[i] for i in range(start, stop)]
    gaussian_sum = np.sum([subset[i] * window[i] for i in range(0, len(subset))])

    return gaussian_sum
     


if __name__ == "__main__":
    main()

