# `lattice` Phase Model

A simple thermodynamic model for a bulk phase, assuming a lattice of solid
atoms

The bulk consists of a matrix of equivalent sites whose molar density does not vary with
temperature or pressure. The thermodynamics obeys the ideal solution laws. The phase and
the pure species phases which comprise the standard states of the species are assumed to
have zero volume expansivity and zero isothermal compressibility.

The density of matrix sites is given by the variable $C_o$, which has SI units of kmol
m-3.

## Specification of Species Standard State Properties

It is assumed that the reference state thermodynamics may be obtained by a pointer to a
populated species thermodynamic property manager class (see
{ct}`ThermoPhase::m_spthermo`). However, how to relate pressure changes to the reference
state thermodynamics is within this class.

Pressure is defined as an independent variable in this phase. However, it has no effect
on any quantities, as the molar concentration is a constant.

The standard state enthalpy function is given by the following relation, which has a
weak dependence on the system pressure, $P$.

$$
     h^o_k(T,P) = h^{ref}_k(T) +  \left( \frac{P - P_{ref}}{C_o} \right)
$$

For an incompressible substance, the molar internal energy is independent of pressure.
Since the thermodynamic properties are specified by giving the standard-state enthalpy,
the term $\frac{P_{ref}}{C_o}$ is subtracted from the specified reference molar enthalpy
to compute the standard state molar internal energy:

$$
     u^o_k(T,P) = h^{ref}_k(T) - \frac{P_{ref}}{C_o}
$$

The standard state heat capacity, internal energy, and entropy are independent of
pressure. The standard state Gibbs free energy is obtained from the enthalpy and entropy
functions.

The standard state molar volume is independent of temperature, pressure, and species
identity:

$$
     V^o_k(T,P) = \frac{1.0}{C_o}
$$

## Specification of Solution Thermodynamic Properties

The activity of species $k$ defined in the phase, $a_k$, is given by the ideal solution
law:

$$
     a_k = X_k ,
$$

where $X_k$ is the mole fraction of species $k$. The chemical potential for species
$k$ is equal to

$$
     \mu_k(T,P) = \mu^o_k(T, P) + R T \ln X_k
$$

The partial molar entropy for species *k* is given by the following relation,

$$
     \tilde{s}_k(T,P) = s^o_k(T,P) - R \ln X_k = s^{ref}_k(T) - R \ln X_k
$$

The partial molar enthalpy for species *k* is

$$
     \tilde{h}_k(T,P) = h^o_k(T,P) = h^{ref}_k(T) + \left( \frac{P - P_{ref}}{C_o} \right)
$$

The partial molar Internal Energy for species *k* is

$$
     \tilde{u}_k(T,P) = u^o_k(T,P) = u^{ref}_k(T)
$$

The partial molar Heat Capacity for species *k* is

$$
     \tilde{Cp}_k(T,P) = Cp^o_k(T,P) = Cp^{ref}_k(T)
$$

The partial molar volume is independent of temperature, pressure, and species
identity:

$$
     \tilde{V}_k(T,P) =  V^o_k(T,P) = \frac{1.0}{C_o}
$$

It is assumed that the reference state thermodynamics may be obtained by a
pointer to a populated species thermodynamic property manager class (see
{ct}`ThermoPhase::m_spthermo`). How to relate pressure changes to the reference
state thermodynamics is resolved at this level.

Pressure is defined as an independent variable in this phase. However, it
only has a weak dependence on the enthalpy, and doesn't effect the molar
concentration.

## Application within Kinetics Managers

$C^a_k$ are defined such that $C^a_k = a_k = X_k$. $C^s_k$, the standard concentration,
is defined to be equal to one. $a_k$ are activities used in the thermodynamic functions.
These activity (or generalized) concentrations are used by kinetics manager classes to
compute the forward and reverse rates of elementary reactions. The activity
concentration, $C^a_k$, is given by the following expression.

$$
     C^a_k = C^s_k  X_k  =  X_k
$$

The standard concentration for species *k* is identically one

$$
    C^s_k =  C^s = 1.0
$$

For example, a bulk-phase binary gas reaction between species $j$ and $k$,
producing a new species l would have the following equation for its rate of
progress variable, $ R^1 $, which has units of kmol m-3 s-1.

$$
   R^1 = k^1 C_j^a C_k^a =  k^1  X_j X_k
$$

The reverse rate constant can then be obtained from the law of microscopic
reversibility and the equilibrium expression for the system.

$$
    \frac{X_j X_k}{ X_l} = K_a^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
$$

$  K_a^{o,1} $ is the dimensionless form of the equilibrium constant,
associated with the pressure dependent standard states $ \mu^o_l(T,P) $
and their associated activities,
$ a_l $, repeated here:

$$
     \mu_l(T,P) = \mu^o_l(T, P) + R T \ln a_l
$$

The concentration equilibrium constant, $ K_c $, may be obtained by
changing over to activity concentrations. When this is done:

$$
    \frac{C^a_j C^a_k}{ C^a_l} = C^o K_a^{o,1} = K_c^1 =
        \exp(\frac{\mu^{o}_l - \mu^{o}_j - \mu^{o}_k}{R T} )
$$

Kinetics managers will calculate the concentration equilibrium constant, $K_c$, using
the second and third part of the above expression as a definition for the concentration
equilibrium constant.
