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

""" 

"""
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

    plt.subplot(221)
    plt.plot(data[0,:], data[1,:], "b")
    plt.subplot(222)
    plt.plot(data_veloc[0,:], data_veloc[1,:], "g")
    plt.subplot(223)
    plt.plot(data_accel[0,:], data_accel[1,:], "g")
    plt.show()

    ax = ax3d(plt.gcf())
    ax.plot(data_veloc[0,:], data_accel[1,:], data_veloc[1,:])
    plt.show()

    


if __name__ == "__main__":
    main()

