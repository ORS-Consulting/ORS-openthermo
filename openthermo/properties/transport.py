"""
Heat and mass transfer correlations for vessel blowdown calculations.

This module provides functions for calculating:
- Dimensionless numbers (Grashof, Prandtl, Nusselt, Rayleigh)
- Heat transfer coefficients for natural convection (gas and liquid)
- Nucleate boiling heat transfer using Rohsenow correlation
- Liquid density estimation using COSTALD correlation

The correlations are based on established literature including:
- Geankoplis, Transport Processes and Unit Operations (Natural convection)
- Rohsenow correlation for nucleate pool boiling (ht package)
- COSTALD correlation for saturated liquid densities

All functions use SI units unless otherwise specified.
The thermo package is used as the thermodynamic backend for fluid property calculations.
"""

import math
import copy
from scipy.constants import g
from ht import Rohsenow
from thermo.volume import COSTALD_mixture
import warnings

warnings.simplefilter("once", UserWarning)


def COSTALD_rho(phase):
    """
    Calculate liquid mass density using COSTALD correlation.

    The COSTALD (Corresponding States Liquid Density) correlation provides more
    accurate saturated liquid densities compared to cubic equations of state,
    especially for hydrocarbon mixtures.

    Parameters
    ----------
    phase : thermo.EquilibriumState
        Liquid phase object containing composition, temperature, and critical properties

    Returns
    -------
    float
        Liquid mass density [kg/m³]

    Notes
    -----
    Uses COSTALD_Vm to calculate molar volume, then converts to mass density
    using the phase molecular weight.

    References
    ----------
    Hankinson, R. W., & Thomson, G. H. (1979). A new correlation for saturated
    densities of liquids and their mixtures. AIChE Journal, 25(4), 653-663.
    """
    # Vm = COSTALD_mixture(phase.zs, phase.T, phase.Tcs, phase.Vcs, phase.omegas)
    Vm = COSTALD_Vm(phase)
    rho = 1 / Vm * phase.MW() / 1000
    return rho


def COSTALD_Vm(phase):
    """
    Calculate liquid molar volume using COSTALD correlation with corrected parameters.

    Uses component-specific characteristic volumes (Vchar) and SRK omega values
    for improved accuracy compared to using critical volumes directly. Falls back
    to standard parameters for components not in the special list.

    Parameters
    ----------
    phase : thermo.EquilibriumState
        Liquid phase object containing composition, temperature, and critical properties

    Returns
    -------
    float
        Liquid molar volume [m³/mol]

    Notes
    -----
    The function maintains a database of corrected parameters for common components:
    - Light gases: N₂, CO₂, CH₄
    - Hydrocarbons: C₂-C₁₀ (including isomers)
    - Other: H₂O, H₂S

    For components not in the database, uses critical volume and acentric factor
    from the phase object directly.

    References
    ----------
    Hankinson, R. W., & Thomson, G. H. (1979). A new correlation for saturated
    densities of liquids and their mixtures. AIChE Journal, 25(4), 653-663.
    """
    name = [
        "nitrogen",
        "hydrogen sulfide",
        "carbon dioxide",
        "water",
        "methane",
        "ethane",
        "propane",
        "i-butane",
        "n-butane",
        "i-pentane",
        "n-pentane",
        "hexane",
        "heptane",
        "octane",
        "nonane",
        "decane",
    ]
    Vchar = [
        9.01e-02,
        9.94e-02,
        9.38e-02,
        4.36e-02,
        9.94e-02,
        0.145750001,
        0.200080007,
        0.256830007,
        0.254390001,
        0.309590012,
        0.311320007,
        0.368200004,
        0.430440009,
        0.490420014,
        0.552900016,
        0.619220018,
    ]
    srk_omega = [
        3.58e-02,
        9.30e-02,
        0.23725,
        -0.654420018,
        7.40e-03,
        9.83e-02,
        0.153200001,
        0.182500005,
        0.200800002,
        0.239950001,
        0.252200007,
        0.300700009,
        0.350690007,
        0.399800003,
        0.44780001,
        0.491600007,
    ]
    names = phase.names
    V_char = copy.deepcopy(phase.Vcs)
    omega = copy.deepcopy(phase.omegas)
    for i in range(len(name)):
        if name[i] in names:
            idx = names.index(name[i])
            V_char[idx] = Vchar[i] / 1000
            omega[idx] = srk_omega[i]

    # Vm = COSTALD_mixture(phase.zs, phase.T, phase.Tcs, phase.Vcs, phase.omegas)
    Vm = COSTALD_mixture(phase.zs, phase.T, phase.Tcs, V_char, omega)
    return Vm


def h_inside(L, Tvessel, Tfluid, fluid):
    """
    Calculate internal natural convective heat transfer coefficient for gas phase.

    Uses empirical correlations for natural convection from vertical surfaces or
    horizontal cylinders. Properties are evaluated at the film temperature (average
    of bulk fluid and wall temperatures) for improved accuracy.

    Parameters
    ----------
    L : float
        Characteristic length [m]
        For vertical vessels: vessel height
        For horizontal vessels: vessel diameter
    Tvessel : float
        Temperature of the vessel wall [K]
    Tfluid : float
        Temperature of the bulk fluid [K]
    fluid : thermo.EquilibriumState
        Gas phase object equilibrated at film temperature [(Tfluid + Tvessel)/2]

    Returns
    -------
    float
        Internal heat transfer coefficient [W/(m²·K)]

    Notes
    -----
    The calculation follows these steps:
    1. Calculate Prandtl number: Pr = Cp·μ/k
    2. Calculate Grashof number: Gr = g·β·ΔT·L³/ν²
    3. Calculate Rayleigh number: Ra = Pr·Gr
    4. Calculate Nusselt number from Ra (see Nu function)
    5. Calculate h = Nu·k/L

    Includes fallback values for thermal conductivity (0.023 W/(m·K)) and
    viscosity (1×10⁻⁵ Pa·s) if property evaluation fails.

    References
    ----------
    Geankoplis, C. J. (1993). Transport Processes and Unit Operations,
    International Edition, Prentice-Hall. Equations 4.7-4 and Table 4.7-1.
    """
    cond = fluid.k()
    if math.isnan(cond):
        cond = 0.023
    visc = fluid.mu()
    if math.isnan(visc):
        visc = 1e-5
    nu = visc / fluid.rho_mass()

    cp = fluid.Cp_mass()
    Pr = cp * visc / cond
    beta = fluid.isobaric_expansion()
    Gr = 9.81 * beta * abs(Tvessel - Tfluid) * L**3 / nu**2
    Ra = Pr * Gr
    NNu = Nu(Ra, Pr)
    h_inner = NNu * cond / L
    return h_inner


def h_inside_liquid(L, Tvessel, Tfluid, fluid):
    """
    Calculate internal natural convective heat transfer coefficient for liquid phase.

    Similar to gas phase natural convection but uses liquid-specific properties.
    Thermal conductivity is calculated using mass-weighted mixing rule for
    multi-component liquids.

    Parameters
    ----------
    L : float
        Characteristic length [m]
    Tvessel : float
        Temperature of the vessel wall [K]
    Tfluid : float
        Temperature of the bulk liquid [K]
    fluid : thermo.EquilibriumState
        Liquid phase object

    Returns
    -------
    float
        Internal heat transfer coefficient [W/(m²·K)]

    Notes
    -----
    For mixtures, thermal conductivity is calculated using a mass-weighted
    average of pure component liquid thermal conductivities:
    k_mix = Σ(w_i · k_i)

    where w_i is the mass fraction and k_i is the thermal conductivity of
    component i.

    References
    ----------
    Geankoplis, C. J. (1993). Transport Processes and Unit Operations,
    International Edition, Prentice-Hall.
    """

    cond = sum([k * w for k, w in zip(fluid.ws(), fluid.kls())])
    visc = fluid.mu()
    cp = fluid.Cp_mass()
    Pr = cp * visc / cond
    beta = fluid.isobaric_expansion()
    Gr = 9.81 * beta * abs(Tvessel - Tfluid) * L**3 / fluid.nu() ** 2
    Ra = Pr * Gr
    NNu = Nu(Ra, Pr)
    h_inner = NNu * cond / L
    return h_inner


def h_inside_wetted(L, Tvessel, Tfluid, fluid):
    """
    Calculate internal heat transfer coefficient for boiling liquid using Rohsenow correlation.

    Combines nucleate pool boiling and natural convection, returning the maximum
    of the two mechanisms with practical bounds (1000-3000 W/(m²·K)).

    Parameters
    ----------
    L : float
        Characteristic length [m]
    Tvessel : float
        Temperature of the vessel wall [K]
    Tfluid : float
        Temperature of the bulk liquid [K]
    fluid : thermo.FlashPureVLS
        Two-phase flash result containing liquid and vapor phases

    Returns
    -------
    float
        Internal heat transfer coefficient [W/(m²·K)]
        Bounded between 1000 and 3000 W/(m²·K)

    Notes
    -----
    The Rohsenow correlation is used for nucleate pool boiling:

    h = μ_L · ΔH_vap · [g(ρ_L - ρ_V)/σ]^0.5 · [C_p,L · ΔT_e / (C_sf · ΔH_vap · Pr_L^n)]^3

    Parameters:
    - C_sf = 0.013 (surface-fluid combination factor)
    - n = 1.7 (exponent, typical for non-water fluids)
    - ΔT_e = excess wall temperature (capped at 30 K for stability)

    Falls back to natural convection if boiling calculation fails.
    Uses default surface tension (0.0175 N/m) if evaluation fails.

    References
    ----------
    Rohsenow, W. M. (1951). A method of correlating heat transfer data for
    surface boiling of liquids. Technical Report.

    Pioro, I. L., et al. (1999). Experimental evaluation of constants for the
    Rohsenow pool boiling correlation. International Journal of Heat and Mass
    Transfer, 42(11), 2003-2013.
    """
    liq = fluid.liquid_bulk
    # if liq.beta == 0:
    #    return 0

    kl = sum([k * w for k, w in zip(fluid.liquid0.ws(), fluid.liquid0.kls())])

    try:
        sigma = liq.sigma()
    except:
        warnings.warn(
            "Surface tension evaluation failed, using default 0.0175 N/m. May not be accurate.",
            UserWarning,
        )
        sigma = 0.0175

    h_boil = Rohsenow(
        rhol=liq.rho_mass(),
        rhog=fluid.gas.rho_mass(),
        mul=liq.mu(),
        kl=kl,
        Cpl=liq.Cp_mass(),
        Hvap=(fluid.gas.H_mass() - liq.H_mass()),
        sigma=sigma,
        Te=min(max((Tvessel - Tfluid), 0), 30),
        Csf=0.013,
        # n=1.3,
        # Csf=0.018,
        n=1.7,
    )

    if math.isnan(h_boil):
        h_boil = 0

    h_conv = h_inside_liquid(L, Tvessel, Tfluid, liq)
    # return 3000
    # return h_conv
    return min(max(h_boil, h_conv, 1000), 3000)


def Pr(cp, mu, k):
    """
    Calculate Prandtl number.

    The Prandtl number is a dimensionless number representing the ratio of
    momentum diffusivity (kinematic viscosity) to thermal diffusivity.

    Parameters
    ----------
    cp : float
        Specific heat capacity at constant pressure [J/(kg·K)]
    mu : float
        Dynamic viscosity [Pa·s]
    k : float
        Thermal conductivity [W/(m·K)]

    Returns
    -------
    float
        Prandtl number [-]

    Notes
    -----
    Pr = (μ · Cp) / k = ν / α

    where:
    - ν is kinematic viscosity
    - α is thermal diffusivity

    Typical values:
    - Gases: 0.7-1.0
    - Water: 2-13 (decreases with temperature)
    - Oils: 50-100,000

    References
    ----------
    Geankoplis, C. J. (1993). Transport Processes and Unit Operations,
    International Edition, Prentice-Hall. Equation 4.5-6.
    """
    return cp * mu / k


def Nu(Ra, Pr):
    """
    Calculate Nusselt number for natural convection from vertical surfaces.

    The Nusselt number represents the ratio of convective to conductive heat
    transfer normal to a boundary. Different correlations are used depending
    on the Rayleigh number regime.

    Parameters
    ----------
    Ra : float
        Rayleigh number [-]
    Pr : float
        Prandtl number [-] (not currently used in calculation)

    Returns
    -------
    float
        Nusselt number [-]

    Notes
    -----
    The correlation depends on the Rayleigh number regime:

    - Ra ≥ 10⁹ (turbulent): Nu = 0.13 · Ra^(1/3)
    - 10⁴ < Ra < 10⁹ (transition): Nu = 0.59 · Ra^(1/4)
    - Ra ≤ 10⁴ (laminar): Nu = 1.36 · Ra^(1/5)

    Once Nu is known, the heat transfer coefficient is calculated as:
    h = Nu · k / L

    where k is thermal conductivity and L is the characteristic length.

    References
    ----------
    Geankoplis, C. J. (1993). Transport Processes and Unit Operations,
    International Edition, Prentice-Hall. Equation 4.7-4 and Table 4.7-1.
    """
    if Ra >= 1e9:
        NNu = 0.13 * Ra**0.333
    elif Ra < 1e9 and Ra > 1e4:
        NNu = 0.59 * Ra**0.25
    else:
        NNu = 1.36 * Ra**0.20
    return NNu
