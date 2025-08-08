"""
Vibrational species declaration check
==============
Compute EEDF with two term approximation solver at constant E/N.

The goal of this file is to check whether the vibrational species declared within the phase in the YAML can be seen.

Requires: cantera >= XX.

.. tags:: Python, plasma
"""


import matplotlib.pyplot as plt
import cantera as ct

gas = ct.Solution('data/example_data/air-plasma-test-vib.yaml')
gas.TPX = 300., 101325., 'N2:0.79, O2:0.21, N2+:1E-10, Electron:1E-10'
gas.reduced_electric_field = 200.0 * 1e-21 # Reduced electric field [V.m^2]
gas.update_EEDF()

print("The number of declared vibrational species withint the load plasma phase is: ", gas.nsp_evib)
