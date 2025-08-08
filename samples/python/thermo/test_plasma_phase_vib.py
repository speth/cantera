import matplotlib.pyplot as plt
import cantera as ct

gas = ct.Solution('./air-plasma-test-vib.yaml')
gas.TPX = 300., 101325., 'N2:0.79, O2:0.21, N2+:1E-10, Electron:1E-10'
gas.reduced_electric_field = 200.0 * 1e-21 # Reduced electric field [V.m^2]
gas.update_EEDF()

print(gas.nsp_evib)
