#!/usr/bin/env python
import numpy as np
import system as sys
import math as mt

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
    print init_cond
    
    m1_M = 0.001
    m2_M = 0.04

    print "Simulating NBody"

def planets(x, t):
    """ Planet ODE system """
    dxdt = np.zeros(np.size(x))
    
    r_se = np.sqrt(np.sum(x[XE:(YE + 1)]**2))                      # r_1 between sun and earth.
    r_ej = np.sqrt(np.sum((x[XJ:(YJ + 1)] - x[XE:(YE + 1)])**2))   # r_21 between earth and jupiter.
    r_sj = np.sqrt(np.sum(x[XJ:(YJ + 1)]**2))                      # r_2 between sun and jupiter.

    dxdt[XE]    = x[VXE]
    dxdt[YE]    = x[VYE]
    dxdt[VXE]   = 1 #TODO
    dxdt[VJE]   = 1 #TODO


    return dxdt

if __name__ == '__main__':
    main()

