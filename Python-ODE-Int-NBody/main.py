#!/usr/bin/env python
import numpy as np
import system as sys
import math as mt

X1 = 0
Y1 = 1
VX1= 2
VY1= 3
X2 = 4
Y2 = 5
VX2= 6
VY2= 7

def main():
    GM = 4 * mt.pi ** 2
    x_1 = 2.52
    y_1 = 0.0
    v_x1 = 0.0
    v_y1 = mt.sqrt(GM / x_1)
    x_2 = 5.24
    y_2 = 0.0
    v_x2 = 0.0
    v_y2 = mt.sqrt(GM / x_2)
    
    m1_M = 0.001
    m2_M = 0.04

    init_cond = np.array([])
    sim_d = sys.simulate(
    print "Simulating NBody"

def planets(x, t):
    dxdt = np.zeros(np.size(x))
    dxdt[X1] = x[VX1]
    dxdt[Y1] = x[VY1]
    #dxdt[X2] = 
    #dxdt[Y2] = 

    return dxdt

if __name__ == '__main__':
    main()

