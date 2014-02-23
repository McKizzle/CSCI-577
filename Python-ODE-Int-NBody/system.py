import numpy as np

def simulate(init_cond, dt, t_0, t_max, ode_sys, intgrtr, verbose=False):
    """ Runs a system for a predifined timerange. 
    :param init_cond: The initial conditions of the simulation. Takes in a list
    of np.arrays of length one.
    :param dt: time step size in seconds.
    :param t_0: starting time.
    :param t_max: stop time.
    :param ode_sys: the ode system to solve. 
    :param intgrtr: the ode solver to use. Do not use the predictor corrector.
    :param verbose: 

    :return: Returns a numpy array of times and the simulation results as
        a numpy array. 
    """
    
    if verbose:
        print "dt: %0.4f" % dt
        print "Starting time: %0.2f" % t_0
        print "Stop time: %0.2f" % t_max
        print "Initial Conditions", init_cond

    sim_data = init_cond
    time_data = [[t_0]]

    dt_steps = int((t_max - t_0) / dt)
    if verbose: 
        print "Steps to take: %d" % dt_steps
    for step in range(1, dt_steps):
        t = step * dt + t_0
        sim_data.append(intgrtr(ode_sys, init_cond[step - 1], dt, t))
        time_data.append([t])
        #print "\t%0.4f seconds at step %d" % (t, step) and verbose
        if verbose:
            print "\t%0.4f seconds at step %d" % (t, step)
                
    sim_data = np.array(sim_data)
    return sim_data

def falling(x, t, g = -9.8): 
    """ Falling ode for testing """
    return np.array([x[1], g])

