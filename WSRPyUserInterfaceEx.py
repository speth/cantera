#WSR Python User Interface Example 1

import sys
import cantera as ct
import numpy as np

class PrimaryZoneReactor(ct.DelegatedIdealGasConstPressureReactor):

   def __init__(self, *args, neighbor, **kwargs): #takes arguments the user specifies when creating PrimaryZone Reactor
      super().__init__(*args, **kwargs) #calls init of superclass (aka DelegatedReactor and Reactor)
 #add variables and values to be initialized upon object creation
      self.soot_comp = 'C' #pure carbon
      self.coag_coeff = 1 #a number 1-9 chosen by user
      self.PAH_list = 'C12H8, C12H10, C13H10, C14H10, C16H10' #PAH species (soot nucleation precursors)
      self.neighbor = neighbor
      self.M = 0 #initialize soot variables to zero
      self.N = 0

   def after_initialize(self,t0): #set initial time to t0
      self.n_vars += 2 #increases neq in reactor by 2 (doesnt this just have the get function but no set in the property definition in reactor.pyx?)
      self.i_M = self.nvars - 1 #index for soot mass concentration
      self.i_N = self.n_vars - 2 #index for soot number concentration
 
   def after_get_state(self, y):
 # This method is used to set the initial condition used by the ODE solver
      y[self.i_M] = self.M
      y[self.i_N] = self.N
    
   def after_update_state(self, y):
 #save solution to state vector y
      self.M = y[self.i_M]
      self.N = y[self.i_N]
 #add other variables to update here? (if want to plot/save)
    
   def before_eval(self,): #can you define variables here and they be used in c++?
      #define user coefficients for gov eq, all equal to 1 or 0 initially, user can edit
      userCoeff_mcpdTdt = 1 #multiply the LHS by a term ex. = ((1/userCoeff_mcpdTdt) + msoot*cpsoot)^-1
      userAdd_mcpdTdt = 0 #add a term to the RHS of mcpdTdt

   def after_eval(self, t, ydot):
 # Calculate the time derivative for the additional equation
 #this is where I will ask cantera for all of the fun properties needed 
 #ex. the formation rates of the species

#Application of user created class
gas = ct.Solution('h2o2.yaml')

# Initial condition
P = ct.one_atm
gas.TPY = 920, P, 'H2:1.0, O2:1.0, N2:3.76'

# Set up the reactor network
res = ct.Reservoir(gas)
r = PrimaryZoneReactor(gas, neighbor=res)
w = ct.Wall(r, res)
net = ct.ReactorNet([r])

# Integrate the equations, keeping T(t) and Y(k,t)
states = ct.SolutionArray(gas, 1, extra={'t': [0.0], 'V': [r.volume]})
while net.time < 0.5:
    net.advance(net.time + 0.005)
    states.append(TPY=r.thermo.TPY, V=r.volume, t=net.time)