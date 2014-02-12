import numpy as np

"""
    v(t) = \frac{y(t+\Delta t) - y(t-\Delta t)}{2\Delta t}
    a(t) = \frac{y(t+\Delta t )- 2 y(t) + y(t - \Delta t)}{\Delta t^2}
"""

""" Estimates the velocity

Takes in a estimates the velocity at a given index in the data. 
The data needs to be in the following form. 
The first row contains the time data and the second row contains the position data. 

"""
def velocity(index, data):
    y_tpdt = data[1, index + 1]
    y_tmdt = data[1, index - 1]
    pdt = data[0, index + 1] - data[0, index]
    mdt = data[0, index] - data[0, index - 1]
    return (y_tpdt - y_tmdt) / (pdt + mdt)

""" Estimates the acceleration

Takes in and estimates the acceleration values in a data set.

"""
def acceleration(index, data):
    y_tpdt = data[1, index + 1]
    y_tmdt = data[1, index - 1]
    y_t = data[1, index]
    pdt = data[0, index + 1] - data[0, index]
    mdt = data[0, index] - data[0, index - 1]

    return (y_tpdt - 2 * y_t + y_tmdt) / (pdt * mdt)


