# -*- coding: utf-8 -*-
"""
This example solves a plug-flow reactor problem with an ideal gas constant pressure mole reactor. It then compares and reports the timing of this with a preconditioned version.
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import time


def preconditioner_pfr(fuel="CH4:1", air="O2:1, N2:3.76", mech="gri30.yaml", phi=1):
    # Preconditioned solver
    T0 = 1500.0  # inlet temperature [K]
    pressure = ct.one_atm # constant pressure [Pa]
    # first time
    t1 = time.time()
    gas1 = ct.Solution(mech)
    gas1.TP = T0, pressure
    gas1.set_equivalence_ratio(phi, fuel, air)
    # create a new reactor
    reactor1 = ct.IdealGasConstPressureMoleReactor(gas1)
    reactor1.volume = 1.0
    # create a reactor network for performing time integration
    net1 = ct.ReactorNet([reactor1,])
    # setup preconditioner
    precon1 = ct.AdaptivePreconditioner()
    precon1.set_threshold(1e-8)
    net1.preconditioner = precon1
    net1.problem_type = "GMRES"
    # approximate a time step to achieve a similar resolution as in
    # the next method
    tf = 1.0
    curr_time = 0
    while curr_time < tf:
        # perform time integration
        curr_time = net1.step()
    t2 = time.time()
    gas2 = ct.Solution(mech)
    gas2.TP = T0, pressure
    gas2.set_equivalence_ratio(phi, fuel, air)
    # create a new reactor
    reactor2 = ct.IdealGasConstPressureMoleReactor(gas2)
    reactor2.volume = 1.0
    # create a reactor network for performing time integration
    net2 = ct.ReactorNet([reactor2,])
    # approximate a time step to achieve a similar resolution as in
    # the next method
    curr_time = 0
    while curr_time < tf:
        # perform time integration
        curr_time = net1.step()
    t4 = time.time()
    # Report speedup
    print("Speedup: {:0.8f}".format((t4-t2)/(t2-t1)))


if __name__ == "__main__":
    preconditioner_pfr()
