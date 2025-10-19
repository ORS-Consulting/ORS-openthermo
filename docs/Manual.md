---
title: openthermo
subtitle: User guide and technical reference
author: Anders Andreasen
titlepage: true
toc-own-page: true
book: true
reference-section-title: References
bibliography: references.bib
listings: True
---

# Introduction
*openthermo* is an open source Python3 tool for calculation of vessel depressurization / blowdown. The main phenomena modelled are visualized in [@Fig:logo].
The thermodynamic state inside the vessel changes over time as seen from immediately observable variables temperature (T) and pressure (P).
This is caused by change in fluid inventory (density) due to flow of gas and/or liquid out of the vessel.
Further, heat is transferred from or to the surroundings via convective heat transfer on the in- and outside of the vessel with heat being conducted through the vessel wall. Due to differences in thermal resistance the vessel wall will obtain a temperature different from the fluid. Depending on the assumptions regarding the description of the fluid inside the vessel, the gas and liquid may have the same temperature (equilibrium assumption) or the two-phases may have different temperature (partial equilibrium assumption).

![openthermo main sketch](docs/img/vessel_sketch.png){#fig:logo}

## Citing *openthermo*
If you use *openthermo* please cite the following reference:

Andreasen, A., Stegelmann, C. (2025). Open source pressure vessel blowdown modelling under partial phase equilibrium. Process Safety Progress (Under review)

    @article{AndreasenStegelmann,
      year = {2025},
      publisher = {Wiley},
      author = {Anders Andreasen and Carsten Stegelmann},
      title = {Open source pressure vessel blowdown modelling under partial phase equilibrium (Under review)},
      journal = {Process Safety Progress}
    }

A preprint of the paper is available on ChemRxiv: https://doi.org/10.26434/chemrxiv-2025-00xzc-v2. 

## Background
### Early works
The foundation of *openthermo* was laid by Carsten Stegelmann, more than a decade ago, who developed code for running blowdown calculations in a spreadsheet relying heavily on VBA and a legacy flash calculation routine (DLL) coded in FORTRAN by late Prof. Michael Michelsen [@michelsen_isothermal_1982;@michelsen_isothermal_1982-1]. This worked surprisingly well and executed very efficiently. The short-comings where lacking heat transfer modelling as well as an *equilibrium* only approach i.e. two-phase fluids were in full equilibrium at all times. 

### Challenges
The VBA code provided tricky to maintain in a version control system, and in the meantime the availability of high-quality thermodynamic packages for Python increased significantly. However, reimplementing an entire codebase is a time-consuming task, and this proved difficult to manage working full time as engineering consultants, and years went by without being able to fully to the long haul required. 

### Proof of concept
Having worked together in the same company, Carsten and I split ways 5 years ago. Then the Covid-19 hit and that freed up some spare time for me, staying at home without a lot of activities being possible and that lead to the development of [HydDown](https://github.com/andr1976/HydDown) [@Andreasen2021]. It started as a small spare-time project for calculation of vessel filling and depressurization behaviour. At that time the expectation was that a lot of engineering work was expected in high-pressure storage and filling stations. The work on HydDown served as a proof of concept for an efficient implementation in Python mainly provided by the [Coolprop](http://www.coolprop.org/) back-end [@doi:10.1021/ie4033999]. Eventually HydDown matured and is now in a state where it can model heat transfer in both steel and dual-layer low thermal conductivity composites during depressurisation/pressurisation. However, it cannot manage two-phase (gas/liquid) behaviour due to limitations in the flash calculation. Thus, a change in thermodynamic back-end was inevitable.

### *openthermo* development
Recently, I joined ORS Consulting with Carsten, and we revived our plans for a rigorous blowdown simulation tool. The remaining challenge was the more complex two-phase (or three-phase) flash problem and the non-equilibrium / partial equilibrium assumption. 
At all times the big inspiration has been the work done at Empirical College London, University College London and later on also Università Politecnica delle Marche on the codes BLOWDOWN, BLOWSIM and VBsim, respectively. The former was also acquired by AspenTech and made available in HYSYS. Further the motivation has also been to have a tool easily accessible as a supplement to commerical (and expensive) tools. Both to reduce load on license pools, but also to provide more efficient workflows. We wanted to make the tool open source and available to the public, but license limitations on the legacy flash calculation by Prof. Michelsen required an alternative flash calculation. Several tools are now available such as Python [*thermo*](https://github.com/CalebBell/thermo) [@thermo], [NeqSim](https://equinor.github.io/neqsimhome/) [@neqsim] and [*thermopack*](https://thermotools.github.io/thermopack/) . Handling vessel depressurisation, which is effectively an UV-flash problem (Internal Energy - Volume) requires extremely many flash calculations to be performed. Thus, a fast and stable flash calculation is required. In order to provide speed and stability the preliminary choice has been [*thermopack* from SINTEF](https://thermotools.github.io/thermopack/), although it may change in the future in order to provide a three-phase (VLLE) flash.   

## Limitations and implementation details 
A few choices has been made to keep things simple:

- [*thermopack*](https://thermotools.github.io/thermopack/) [@thermopack] is used as thermodynamic backend
- Thermodynamic and transport properties (enthalpy, entropy, internal energy, liquid density, thermal conductivity, viscosity) is provided via Python *thermo* / *chemicals*.
- No temperature stratification inside bulk phases
- No temperature gradient through vessel wall. 
- Heat transfer is modelled as constant or simplified using empirical correlations
- Only single and two-phase (VLE) is handled. Three-phase (VLLE) cannot currently be modelled in the open source version.
- Currently only the Peng-Robinson (PR) [@peng_new_1976] and Soave-Redlich-Kwong [@SOAVE19721197] cubic equations of state are made available.  
- *openthermo* and it's legacy code has been built by engineers NOT software developers.  

Ignoring temperature gradients in the vessel wall is an acceptable assumption for (not too thick) steel vessel walls. However, low thermal conductivity materials cannot accurately be modelled. In order to do so a 1-D heat transfer model shall be implemented. 

One limitation of *thermopack* is that including pseudo components is not easy if the internal property calculations shall be used. Further, *thermopack* does not provide transport properties. In order to allow an implementation for pseudo components and transport properties the *thermopack* library is only used for flash calculations and used as a plugin, the rest (and less computationally demanding tasks) is provided by Python *thermo*. 

While *openthermo* has been extensively validated and fidelity has been built, we have not had the work force or man-months/years available like the legacy university projects or the commerical software providers. This means that there likely will be situations where the code fails to provide results, cases that are not covered or simply the code will not work. This is the cost and risk of using open source software.   

## Getting the software
The source code can be obtained either from GitHub (via `git` or via the latest tar-ball release) or via **pip** . No packaged releases have currently been planned for **conda**.  

The main branch is located here:

[`https://github.com/ORS-Consulting/ORS-openthermo`](https://github.com/ORS-Consulting/ORS-openthermo)

Clone the repo by:

    git clone https://github.com/ORS-Consulting/ORS-openthermo.git

Running from source via `git` the dependencies must be installed manually from the repo root dir:

    pip install -r requirements.txt
    pip install -e .

Alternatively `pip` can install directly from github:

    pip install git+https://github.com/ORS-Consulting/ORS-openthermo.git

Installation of latest release via **pip** also installs dependencies automatically (still pending):

    pip install openthermo

## Requirements
- [Python](http://www.Python.org) 
- [thermopack](https://thermotools.github.io/thermopack/)
- [thermno](https://thermo.readthedocs.io/index.html)
- [chemicals](https://chemicals.readthedocs.io/index.html)
- [ht](https://ht.readthedocs.io/en/release/)
- [fluids](https://fluids.readthedocs.io/)
- [Scipy](https://www.scipy.org/)
- [Numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [SciencePlots](https://github.com/garrettj403/SciencePlots)
- [cerberus](https://docs.Python-cerberus.org/en/latest/)
- [PyYaml](https://pypi.org/project/PyYAML/)
- [pandas](https://pandas.pydata.org/)
- [tqdm](https://tqdm.github.io/)
- [openpyxl](https://openpyxl.readthedocs.io/en/stable/)

The code is running on Windows 10/11 x64, with stock Python installation from Python.org and packages installed using pip.
It should also run on Linux (it does on an Ubuntu image on GitHub) or in any conda environment as well, but this hasn't been checked.

## Testing
Although testing is mainly intended for automated testing (CI) during development using github actions, testing of the installed package can be done for source install by:

    python -m pytest

run from the root folder. 

## Units of measure
The SI unit system is adapted for this project.
The following common units are used in the present project and this also applies to the units used in the input files:

Property | Unit | Comment
----    | ----  | ----
Temperature | K | $^\circ$ C may be used in plots
Pressure | Pa    | bar is used in plots
Mass | kg |
Volume | m$^3$ |
Time | s |
Energy | J |
Duty/power | W
Length | m
Area | m$^2$
Heat flux | W/m$^2$
Heat transfer coefficient | W/(m$^2$ K)
Thermal conductivity | W/(m K)
Density | kg/m$^3$
Heat capacity | J/(kg K)

: Unit system {#tbl:units}

As will be noted when presenting the equations implemented in the code, some of the equations utilise different units than the ones listed in [@tbl:units].
However, it is important to note that unit conversions are built in to the methods implemented, so the user shall not worry about unit conversion.  

## Credit
In the making of this document a great deal of material has been sourced (and modified) from a former colleague's M.Sc. thesis [@iskov], from co-published papers [@Bjerre2017;@safety4010011] and from on-line material published under permissive licenses (with proper citation). In particular the HydDown manual [@Andreasen2024] has been heavily scavenged for usable material,  

Further, the making of this project would not have possible without the magnificent [*thermopack*](https://thermotools.github.io/thermopack/) library [@thermopack].
I am also thankful for enlightening discussions with former colleague Jacob Gram Iskov Eriksen (Ramboll Energy, Denmark)  and former and now present colleague Carsten Stegelmann (ORS Consulting).

The present document is typeset using Markdown + [pandoc](https://pandoc.org/) with the [Eisvogel](https://github.com/Wandmalfarbe/pandoc-latex-template) template.

## License

MIT License

Copyright (c) 2021-2025 Anders Andreasen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# Usage
## Basic usage
A minimal example for running *openthermo* is provided below. The example is for an isentropic blowdown i.e. with a rigorous energy and mass balance, taking into accout work done by the discharged fluid, but without any heat transfer. 

```
from openthermo.vessel.blowdown import Blowdown
from openthermo.flash.michelsen import get_flash_dry

input = {}
P = 12e5
T = 298.15

input["mode"] = "isentropic"
input["eos_model"] = "PR"
input["liquid_density"] = "eos"
input["max_time"] = 900
input["length"] = 10
input["diameter"] = 3
input["vessel_type"] = "Flat-end"
input["orientation"] = "horizontal"
input["liquid_level"] = 1.5
input["operating_temperature"] = T
input["operating_pressure"] = P
input["ambient_temperature"] = 273
input["back_pressure"] = 1.01e5
input["bdv_orifice_size"] = 0.03  # m
input["bdv_orifice_cd"] = 0.84

names = ["methane", "propane", "n-butane", "i-butane", "n-decane"]
molefracs = [0.8, 0.05, 0.01, 0.01, 0.10]

input["molefracs"] = molefracs

input["flash"] = get_flash_dry(
    names,
    molefracs,
    P=P,
    T=T,
    rho=input["liquid_density"],
    model=input["eos_model"],
)

segment = Blowdown(input)
segment.depressurize()
segment.plot()    
```

## Demos

## Calculation methods
The following methods are implemented:

- Isothermal: constant temperature of the fluid during depressurization 
- Adiabatic: the energy balance does not account for work done by the fluid during discharge. No heat transfer is included.  
- Isentropic: the energy balance accounts for work done by the fluid during discharge. This method can be run without or with heat transfer modelled. A variation of this method includes external backgorund heat load from a fire (pool or jet) using the Stefan-Boltzmann equation. 
- Fire: This method is a variant of the isentropic method where an API 521 pool-fire heat load is applied to the wetted surface inside the vessel. No other heat transfer is included.  

For calculations the minimal input required are:

- Initial conditions (pressure, temperature, liquid level)
- vessel dimensions (ID, length)
- orifice parameters (Cd, diameter, back-pressure)
- Calculation setup (end time)
- Fluid components and composition

For rigorous heat transfer additional input are required:

- Vessel thickness
- Vessel material heat capacity
- Vessel material density
- Ambient temperature or a Stefan-Boltzmann external fire heat load type. 

More elaborate description of the required input for the different calculation types are provided in [@Sec:input].

## Input file examples
When using *openthermo* a dictionary holding all relevant input in order for the program to do vessel calculations shall be provided when the class is initialized. 

An example of a minimal file for an isentropic vessel depressurization (no heat transfer) is shown below. 



## Input fields and hierarchy {#sec:input}
In the following the full hierarchy of input for the different calculation types is summarised.

At the top level the following fields are accepted, with the last being optional and the second last dependant on calculation type:


# Theory
In this chapter the basic theory and governing equations for the model implementation in *openthermo* is presented.
The following main topics are covered:

- thermodynamics
- mass transfer
- heat transfer
- Partial equilibrium 

## Thermodynamics

### Equation of state
Currently only the Peng-Robinson [@peng_new_1976] and Soave-Redlich-Kwong [SOAVE19721197] cubic equations of stare are available despite many being available in both *thermopack* and *thermo*. In the future more cubic equations of state may be made available. 

For more background information and theory regarding the equations of state, please refer to e.g. [*thermo* documentation](https://thermo.readthedocs.io/thermo.eos_mix.html), the [DWSIM user guide](https://github.com/DanWBR/dwsim/blob/windows/PlatformFiles/Common/docs/User_Guide.pdf), a [*thermopack* memo](https://github.com/thermotools/thermopack/blob/main/docs/memo/thermopack/thermopack2013.pdf).

### First law for flow process {#sec:firstlaw}
The control volume sketched in [@Fig:firstlaw], separated from the surrounding by a control surface, is used as a basis for the analysis of an open thermodynamic system with flowing streams (fs) in and out, according to [@sva]

A general mass balance or continuity equation can be written:

![Control volume with one entrance and one exit. The image has been sourced from [@firstlaw].](docs/img/First_law_open_system.png){#fig:firstlaw}

$$ \frac{m_{cv}}{dt} + \Delta \left( \dot{m} \right)_{fs}= 0 $$ {#eq:continuity}

The first term is the accumulation term i.e. the rate of change of the mass inside the control volume, $m_cv$, and the $\Delta$ in the second term represents the difference between the outflow and the inflow

$$ \Delta \left( \dot{m} \right) _{fs} = \dot{m}_2 - \dot{m}_1 $$

An energy balance for the control volume, with the first law of thermodynamics applied, needs to account for all the energy modes with can cross the control surface. Each stream has a total energy

$$ U + \frac{1}{2}u^2 + zg $$

where the first term is the specific internal energy, the second term is the kinetic energy and the last term is the potential energy.
The rate at which each stream transports energy in or out of the control volume is given by

$$ \dot{m} (U + \frac{1}{2}u^2 + zg) $$

and in total

$$  \Delta \left[ \dot{m} (U + \frac{1}{2}u^2 + zg) \right]_{fs}$$

Furthermore, work (not to be confused with shaft work) is also associated with each stream in order to move the stream into or out from the control volume (one can think of a hypothetical piston pushing the fluid at constant pressure), and the work is equal to $PV$ on the basis of the specific fluid volume.
The work rate for each stream is

$$ \dot{m}(PV) $$

and in total

$$ \Delta\left[ \dot{m}(PV) \right]_{fs} $$

Further, heat may be transferred to (or from) the control volume at a rate $\dot{Q}$ and shaft work may be applied, $\dot{W}_{shaft}$.
Combining all this with the accumulation term given by the change in total internal energy the following general energy balance can be written:

$$ \frac{d(mU)_{cv}}{dt} + \Delta \left[ \dot{m} (U + \frac{1}{2}u^2 + zg) \right]_{fs} + \Delta \left[ \dot{m}(PV) \right]_{fs} = \dot{Q} +\dot{W}_{shaft}   $$

Applying the relation $H = U + PV$, setting $\dot{W}_{shaft} = 0$ since no shaft work is applied to or extracted from the vessel, and assuming that kinetic and potential energy changes can be omitted the energy balance simplifies to:

$$ \frac{d(mU)_{cv}}{dt} + \Delta \left[ \dot{m} H \right]_{fs} = \dot{Q}  $$

The equation can be further simplified if only a single port acting as either inlet or outlet is present:

$$ \frac{d(mU)_{cv}}{dt} + \dot{m} H  = \dot{Q}  $$ {#eq:energybalance}

where the sign of $\dot{m}$ determines if the control volume is either emptied or filled.
The continuity equation [@Eq:continuity] and the energy balance [@Eq:energybalance] combined with the equation of state are the key equations that shall be solved/integrated in order to calculate the change in temperature and pressure as a function of time.

## Flow devices {#sec:flow}
### Restriction Orifice
When a fluid flows through a restriction or opening such as an orifice, the velocity will be affected by conditions upstream and downstream.
If the upstream pressure is high enough, relative to the downstream pressure, the velocity will reach the speed of sound (Ma = 1) and the flow rate obtained will be the critical flow rate.
This condition is referred to as choked flow.
The maximum downstream pressure for the flow to still be sonic (Ma = 1), is when $P_d = P_c$.
The ratio of the critical and upstream pressure is defined by equation [@Eq:P_critical].

$$ \frac{P_{c}}{P_u}=\left (\frac{2}{k+1} \right)^\frac{k}{k-1} $$ {#eq:P_critical}

- $P_c$ is the critical pressure. [kPa]
- $P_d$ is the downstream pressure. [kPa]
- $P_u$ is the upstream pressure. [kPa]
- k is the isentropic coefficient, approximated by the ideal gas heat capacity ratio $C_p/C_v$. [-]

In order to calculate the mass flow rate through an orifice equations are used based on literature from the Committee for the Prevention of Disasters [@yellowbook].

To account for the difference in choked and non-choked flow a set limit pressure is introduced as in equation [@Eq:plimit].
If the downstream pressure, $P_{down}$, is below the pressure limit, $ P_{limit}$, then the flow is choked, and the pressure used, $P_{used}$, in equation [@Eq:massfloworifice] should be the pressure limit, $P_{limit}$.
Otherwise if the downstream pressure, $P_{down}$, is greater than or equal to the pressure limit, $P_{limit}$, the flow is no longer choked and the pressure used should be the downstream pressure, $P_{down}$ [@yellowbook].

$$ P_{limit}=P_{up} \cdot \left ( \frac{2}{k+1} \right ) ^{\frac{k}{k-1}} $$ {#eq:plimit}

$$ \dot{m}_{flow}= C_d  \cdot A \cdot\sqrt{\left ( \frac{2 k}{k-1}\right )  \cdot  P_{up} \cdot \rho \cdot \left (  \frac{P_{used}}{P_{up}} \right )^{\frac{2}{k}} \left (1-\left ( \frac{P_{used}}{P_{up}} \right ) ^{\frac{k-1}{k}} \right )} $$ {#eq:massfloworifice}

- $\rho$ is the density of the gas upstream. $[kg/m^3]$
- $P_{limit}$ is the pressure limit of the upstream absolute pressure. $[bara]$
- $P_{up}$ is the absolute pressure upstream of the orifice. $[bara]$
- $k$ is the ratio of the heat capacities at constant pressure, $C_p$, and at constant volume, $C_v$.
- $P_{down}$ is the absolute pressure downstream of the orifice. $[bara]$
- $P_{used}$ is the pressure used in the mass flow equation based on choked or non-choked conditions. $[bara]$
- $\dot{m}_{flow}$ is the mass flow through the orifice. $[kg/s]$
- $C_d$ is the discharge coefficient of the orifice opening. $[-]$
- $A$ is the cross sectional area of the orifice. $[m^2]$

## Heat transfer {#sec:heat}

### Natural convection
Experiments have indicated that the internal heat transfer mechanism for a vessel subject to depressurization can be well approximated by that of natural convection as found from measured Nusselt numbers being well correlated with Rayleigh number, with no apparent improvement in model performance by inclusion of the Reynolds number in the model [@woodfield].

To determine the heat transfer for the gas-wall interface, the following is applied cf. equation [@Eq:newton]:

$$\frac{dQ}{dt} = \dot{Q} = h A ( T_{s} - T_{gas} ) $$  {#eq:newton}

- $d Q$ is the change in thermal energy due to convective heat transfer. [J]
- $d t$ is the change in time during the heat transfer. [s]
- $h$ is the convective heat transfer. [W/m$$^2$$K ]
- $A$ is the area normal to the direction of the heat transfer. [m$^2$]
- $T_{s}$ is the inner surface temperature of the geometry. [K]
- $T_{gas}$ is the temperature of the bulk gas inside the vessel. [K]

The convective heat transfer will need to be estimated for the the gas-wall interface, by the use of empirical relations for the Nusselt number.
The Nusselt number describes the ratio of convective heat transfer to conductive heat transfer, normal to a surface area, as given in equation [@eq:Nu].

$$ Nu=\frac{hL}{k} $$ {#eq:Nu}

- $Nu$ is the Nusselt number. [-]
- $h$ is the convective heat transfer. [W/m$^2$K]
- $L$ is a characteristic length of the geometry. [m]
- $k$ is the thermal conductivity of the gas. [W/m$\cdot$K]

The characteristic length $L$ used is the height of the gas volume i.e. the length of a vertical vessel or the diameter of a horizontal vessel.

The empirical correlations used to calculate the Nusselt number of the gas-wall interface is a function of the Rayleigh number, which can be defined by the Grashof number and Prandtl number, as in equation [@Eq:rayleigh_gas]:

$$ Ra=Gr \cdot Pr $$ {#eq:rayleigh_gas}

- $Ra$ is the Rayleigh number. [-]
- $Gr$ is the Grashof  number. [-]
- $Pr$ is the Prandtl number. [-]

The Grashof number is a dimensionless number which approximates the ratio of the buoyancy forces to viscous forces, as given in equation [@Eq:grashof_gas}:

$$ Gr=\frac{\beta g\rho^2 L^3 \Delta T }{\mu^2} $$ {#eq:grashof_gas}

The Prandtl number is a dimensionless number defined as the ratio of the momentum diffusivity to thermal diffusivity, as given in equation [@Eq:prandtl_gas]:

$$ Pr=\frac{c_p \mu}{k} $$ {#eq:prandtl_gas}

- $\beta$ is the coefficient of volume expansion. [1/K]
- $g$ is the standard acceleration of gravity. [m/s$^2$]
- $\rho$ is the gas density. [kg/m$^3$]
- $L$ is the characteristic length. [m]
- $\Delta T$ is the temperature difference of the surface and gas. [K]
- $\mu$ is the dynamic viscosity. [kg/m$\cdot$s]
- $c_p$ is the heat capacity of gas. [J/kg$\cdot$K]
- $k$ is the thermal conductivity of gas. [J/m$\cdot$K]

It is important to note that the properties in the above equations shall be evaluated at the fluid film temperature which can be approximated by the average of the the fluid bulk temperature and the vessel wall temperature [@geankoplis].

The mechanism of heat transfer on the outside of the vessel at ambient conditions is similar to the above.
In HydDown the external heat transfer is not modelled currently.
A heat transfer coefficient shall be provided.

### Conduction
For accurate prediction of the outer and especially the inner wall temperature for correct estimation of internal convective heat transfer and the average material temperature, the general equation of 1-D unsteady heat transfer shall be solved:

$$ \frac{\delta T}{\delta t} = \frac{k}{\rho c_p} \frac{\delta^2 T}{\delta x^2} $$

- T is temperature
- x is the spatial (1-D) coordinate
- k is the thermal conductivity
- $\rho$ is the material density
- $C_p$ is the specific heat capacity

To be solved, the initial values and boundary values must be specified.
In its default state (if thermal cobductivity is not applied for the vessel), HydDown does not include the unsteady heat transfer model, i.e. the assumption is that the temperature from outer to inner surface is uniform and equal to the average temperature.
This is obviously a crude approximation, but might be justified depending in the Biot number:

$$ Bi = \frac{hL}{k} $$

The Biot number is a simple measure of the ratio of the heat transfer resistances at the surface of a body to the inside of a body.
The ratio gives an indication to which extent the temperature will vary in space (gradient) when the body is subject to a displacement in temperature at the surface boundary layer.
Striednig *et al.* [@STRIEDNIG] concluded that for a type I (steel) cylinder the Biot number was approx. 0.03 and hence the error in assuming a uniform temperature in the vessel wall was low.

With a typical thermal conductivity of 45 $W/m K$ for steel and a heat transfer coefficient up to 600 $W/m^2 K$ [@woodfield] the Biot number for a vessel with a wall thickness of 2 cm is 0.27.
This is significantly higher than that approximated by [@STRIEDNIG].
However, the Biot number is significantly lower than 1, and the assumption of a uniform temperature is reasonable.
However, for increased wall thickness, and/or for different materials with lower thermal conductivity, the error may grow to an unacceptable level.

Especially for vessels with low conductivity materials (or very thick walls) accurate estimation of the vessel wall temperatures requires the 1-D transient heat transfer problem to be solved. 

### Fire heat loads
The heat transfer from the flame to the shell is modelled using the recommended approach from Scandpower [@scandpower] and API [@API521].
The heat transfer from the flame to the vessel shell is divided into radiation, convection, and reradiation as seen in equation [@Eq:flame]:

$$ q_f=\underbrace{{\alpha}_s \cdot {\varepsilon}_f \cdot \sigma \cdot T_f^4}_\text{Radiation}+\underbrace{h_f \cdot (T_f-T_s(t))}_\text{Convection}-\underbrace{{\varepsilon}_s \cdot \sigma \cdot T_s(t)^4 }_\text{Reradiation} $$ {#eq:flame}

- $q_f$ is the flame heat flux. [W/m$^2$]
- ${\alpha}_s$ is the vessel surface absorptivity. [-]
- ${\varepsilon}_f$ is the flame emissivity. [-]
- $\sigma$ is the Stefan-Boltzmann constant, $\sigma$ = $5.67 \cdot 10 ^ {-8}$  [W/m$^2 \cdot$ K$^4$]
- $T_f$ is the flame temperature. [K]
- $h_f$ is the convection heat transfer coefficient between the flame and the surface. [W/m$^2 \cdot$ K]
- $T_s(t)$ is the time dependent surface temperature. [K]
- ${\varepsilon}_s$ is the surface emissivity. [-]

This model assumes that the pressure vessel is fully engulfed by the flame.
This means that the view factor for the radiation is unity and is therefore not taken into consideration.
The convective heat transfer coefficients for a jet fire and a pool fire, and recommended values for the emissivity and absorptivity, are given by Scandpower as [@scandpower]:

- $h_{jet~fire}$ = 100 [W/m$^2 \cdot$K]
- $h_{pool~fire}$ = 30 [W/m$^2 \cdot$K]
- ${\alpha}_s$ = 0.85
- ${\varepsilon}_s$ = 0.85
- ${\varepsilon}_f$ = 1.0 (optical thick flames, thickness > 1 m)

The flame temperature is found by solving equation [@Eq:flame2] for the incident heat flux in relation to the ambient conditions.
The flame temperature is kept constant throughout the simulation:

$$ q_{total}=\sigma \cdot T_f^4 + h_f \cdot (T_f-T_{amb})$$ {#eq:flame2}

- $q_{total}$ is the incident flame heat flux as given in table [@Tbl:heatfluxes1]. [W/m$^2$]
- $T_{amb}$ is the ambient temperature $\approx$ 293 K (20$^\circ$ C)

The heat flux used to calculate the flame temperature is given in table [@tbl:heatfluxes1].

|                        | Small jet fire  [kW/m$^2$]  |  Large jet fire  [kW/m$^2$] |  Pool fire  [kW/m$^2$]
| ----                   |  ----           |  ----           | ----
| Peak heat load         |  250            |  350            | 150         
| Background heat load   |   0             |   100           |  100         

: Incident heat fluxes for various fire scenarios given by Scandpower [@scandpower] {#tbl:heatfluxes1}


## Model implementation
A simple (naive) explicit Euler scheme is implemented to integrate the mass balance over time, with the mass rate being calculated from an orifice/valve equation as described in [@Sec:flow]:

$$ m_{cv}(i+1) =  m_{cv}(i) + \dot{m}(i) \Delta t  $$ {#eq:euler_mass}

$$\dot{m}(i) = f(P,T,) $$

For each step, the mass relief/ left in the vessel is known.
Since the volume is fixed the mass density is directly given.

For the calculation methods (isentropic, isenthalpic, isenergetic, etc.), CoolProp allows specifying density and either H, S, or U directly - this is very handy and normally only TP, PH, or TS property pairs are implemented and you would need to code a second loop to make it into am UV, VH, or SV calculation.

In the following, the high-level steps in the solution procedure are outlined for the different calculation types.
To illustrate calls to the equation of state in CoolProp, the notation $EOS(D,T)$ corresponds to calculation of the state at known density and temperature, for example.

### Isothermal process
For an isothermal process, the solution procedure for each calculation step is the following with the mass in the control volume being calculated from [@Eq:euler_mass]:

### Isenthalpic process

### Isentropic process

### Isenergetic process

### Energy balance
The general first law applied to a flow process as outlined in [@Sec:firstlaw] subject to an explicit Euler scheme is:

$$ D(i+1) = \frac{m_{cv}(i+1)}{V} $$
$$ U_{cv}(i+1) = \frac{m_{cv}(i)U_{cv}(i) - \left( \dot{m}(i) H (i) +  \dot{Q}(i) \right) \Delta t}{m_{cv}(i+1)}  $$ {#eq:firstlaw_euler}

The above assumes that mass flow is positive when leaving the control volume and heat rate is positive when transferred to the control volume.
$H(i)$ is the specific enthalpy of the fluid in the control volume for a discharging process and it is equal to the energy of the entering stream for a filling process.
The heat rate is calculated as outlined in [#Sec:heat].

For the vessel wall a simple heat balance is also made:

$$ T_{wall}(i+1) = T_{wall}(i)  + \frac{\dot{Q}_{outer} - \dot{Q}_{inner} } {c_p m_{vessel}} \Delta t $$

where $\dot{Q}_{outer}$ is the convective heat transfer to or from the ambient surroundings from the outer surface of the vessel, with positive values indicating that heat is transferred from the surroundings to the vessel wall.
This is either a fixed heat transfer coefficient with a specified ambient temperature (outer surface) or a calculated fire heat load.
$\dot{Q}_{inner}$ is the internal heat transfer, either a fixed number or calculated as natural convection.

The remaining steps are update of temperature and pressure.


