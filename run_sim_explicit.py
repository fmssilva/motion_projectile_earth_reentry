
import numpy as np

from sim_common_fun import *


def forward_euler_step(Sk, Mk, p: Params, slope_f): 
    ''' Forward Euler step for the capsule reentry simulation. 
        Given previous state (Sk) and metrics (Mk), calculates the current state (Sk + 1) and metrics (Mk +1).
        returns Sk1 as a new vector by value, not altering the original Sk.'''
    slopes, Mk1 = slope_f(Sk, Mk, p)
    steps = slopes * p.dt
    Sk1 = Sk + steps
        
    # Mk1[A] and Mk1[CHUTE_OPEN] were already calculated in get_acceleration function
    Mk1[V] = np.sqrt(Sk1[VX]**2 + Sk1[VY]**2)
    Mk1[EARTH_ANGLE] = Mk[EARTH_ANGLE] + (Sk1[X] - Sk[X]) / Sk1[Y]
    return Sk1, Mk1


if __name__ == '__main__':
    run_all_simulations(forward_euler_step)
