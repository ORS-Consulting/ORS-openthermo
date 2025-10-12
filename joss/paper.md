---
title: 'openthermo: A Python package for equilibrium and partial equilibrium vessel blowdown'
tags:
  - Python
  - Depressurisation
  - Blowdown
  - Process safety
  - Loss prevention
authors:
  - name: Anders Andreasen
    orcid: 0000-0003-0475-323X
    affiliation: 1
affiliations:
 - name: ORS Consulting, Borgergade 66 st. th., DK-6700 Esbjerg, Denmark
   index: 1
date: 12 October 2025
bibliography: paper.bib

---

# Summary

ORS *openthermo* is a Python package for calculation of pressure vessel behaviour during depressurisation (blowdown). The software allows calculation of vessel pressure, fluid inventory temperature as well as vessel wall temperature as a function of time during depressurisation. In addition to a full equilibrium assumption in which liquid and gas (if two-phases are present) have the same temperature the program also includes a partial equilibrium approach which allows the gas and liquid phase to have different temperatures. The latter approach has been necessary in order to match experimental data. 

A typical system modelled using the partial equilibrium approach is shown in \autoref{fig:sketch} and it visualizes the key parameters and transport phenomena during vessel depressurisation. The thermodynamic state inside the vessel changes over time as seen from immediately observable variables temperature (T) and pressure (P). This is caused by a change in fluid inventory (density) due to flow of gas either in or out of the vessel. Furthermore, heat is transferred from or to the surroundings via convective heat transfer on the in- and outside of the vessel - with heat being conducted through the vessel wall. 

![Partial equilibrium visualised. \label{fig:sketch}](../joss/vessel_sketch.png){ width=100% }

In its essence, the code solves the mass and energy balances with gas thermodynamics calculated using a combination of [thermopack](https://github.com/thermotools/thermopack) [@thermopack] and [thermo](https://github.com/CalebBell/thermo) [@thermo] where the former are used for flash calculations and the latter for thermodynamic and transport properties (enthalpy, entropy, internal energy, density, viscosity, thermal conductivity). Either the Peng-Robinson or the Soave-Redlich-Kwong cubic equation of state can be applied. The energy balance is the first law of thermodynamics for an open system exchanging both heat and mass with the surroundings [@sva]. Heat transfer between gas inventory and vessel wall is accounted for using either natural convection (wall in contact with gas) or nucleate boiling heat transfer (wall in contact with liquid) [@geankoplis]. The mass balance is closed using an applicable flow equation [@yellowbook]. The partial equilibrium approach implemented is heavily inspired by that implemented in VBsim [@speranza_blowdown_2005;@ricci_unsteady_2015;@DALESSANDRO2015719].
The code also allows an external heat load to be applied using the API 521 equation for fire heat load from a pool fire [@API521] and leaks (gas, liquid or two-phase). 

A few choices have been made to keep things simple:

- Only gas, liquid or two-phase gas/liquid is modelled. Three phase (Vapour-Liquid-Liquid/VLLE) is not possible.
- Only real chemical components allowed. No pseudo-components / hypothetical components can be added.
- No temperature stratification in vessel inventory
- No temperature gradient through vessel wall (applicable for high heat conductivity / thin-walled vessels)

Still the code can manage a number of different assumption: isothermal depressurisation (very fast), isentropic (1. law) depressurisation with and without heat transfer, rigorous partial equilibrium depressurisation. Typical calculation output is shown in \autoref{fig:pres} and \autoref{fig:wall} with experimental data included for comparison [@Szczepanski;@WONG].

![Calculated pressure compared with experimental pressure. Blowdown of condensing/two-phase hydrocarbon mixture conducted at Spadeadam. \label{fig:pres}](../joss/condensable_gas_pressure_rig.png)

![Calculated vessel wall temperatures in contact with gas and liquid compared with experimental data.  \label{fig:wall}](../joss/condensable_gas_inner_wall_rig.png)

# Statement of need
The rapid depressurisation (blowdown) of pressure vessels containing hazardous, either toxic or flammable,  substances in a chemical process plant, is an essential part of the plant process safety measures. 
Blowdown has both preventive and mitigating effects in case of a potential hazardous event. In case of a leak, the depressurisation significantly shortens the duration of the leak as well as rate of the leak is reduced, thus decreasing the potential fire and explosion risk.
The total amount of released inventory is reduced thereby reducing potential overpressure and impulse from explosion as well is reduced length and radiation from a jet fire. The blowdown also effectively reduces the inventory in nearby sections, thereby also reducing the risk of escalation in case of loss of containment due to exposure to fire or explosion.

Despite the importance of being able to predict vessel response during blowdown, no free tool exists today which accomplishes the same tasks as *openthermo*. One reason is likely the complexity of the problem modelled and also the requirement of high quality thermodynamic models, especially a very robust and computationally efficient flash calculation is required. Developing such a tool requires a tremendous work effort counted in several years of full time devoted work. Existing or legacy tools are only available as a part of proprietary/commercial tools and comes at a significant license fee or via academic tools not available to the general public such as BLOWSIM [@blowsim;@mahgerefteh_numerical_1999] or VBsim [@DALESSANDRO2015719]. Even one of legacy tools (BLOWDOWN) [@haque_rapid_1990] was acquired and is now part of one of the commerical tools. For review of more codes, please refer to @SHAFIQ2020104. A few other tools are available as open source. This includes HydDown [@Andreasen2021] and HYRAM+ [@groth_hyram_2017] both relying on CoolProp [@CoolProp]. Both have limitations in handling multicomponent two-phase mixtures and the latter does not rigorously model heat transfer. None of the two employs a partial equilibrium approach.


# Acknowledgements
 The author is thankful for fruitfull discussions with former colleague Jacob Gram Iskov Eriksen (Ramboll Energy, Denmark) and  colleague Carsten Stegelmann (ORS Consulting) in relation to vessel depressurisation, nozzle flow and heat transfer considerations. Carsten Stegelmann have made a huge contribution of developing VBA code for modelling the blowdown process using the equillibrium assumption and this has formed the foundation for the present python code including the extension of the partial equilibrium approach.  

# References
