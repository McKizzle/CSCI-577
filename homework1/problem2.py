import matplotlib.pyplot as plt
import numpy as np
import ode as ode
import system as sys
import utils as ut

NEA = 0
EMM = 1

A = 1.0 / 12.0
B = 1.0 / 35.0
D = 1.0 / 1000.0
s = 0.96

def run():
    dt = 10
    
    sim_d = sys.simulate([np.array([50, 1])], dt, 0, 10000, pop, ode.euler_richardson)

    ut.plt2dcmpr(
            sim_d[:, 0], sim_d[:, 1], sim_d[:, 0], sim_d[:, 2],
            ["b-", "r-"], ["Neanderthal Population", "Human Population"],
            "time (years)", "Population", 
            "Early Modern Human vs Neanderthal Populations"
            )
    plt.text(5000, 30, "s = %0.3f" % s)
    plt.savefig("clinton_nea-emm-competition.png", bbox_inches="tight")
    plt.close()

def pop(x, t):
    """ Early Modern Man and Neanderthal ODE System """
    global NEA, EMM, A, B, D, s

    N_v = x[NEA] * (A - B - D * (x[NEA] + x[EMM]))
    E_v = x[EMM] * (A - (s * B) - D * (x[NEA] + x[EMM]))

    return np.array([N_v, E_v])

