import math
import copy
from scipy.constants import g
from ht import Rohsenow
from thermo.volume import COSTALD_mixture


def COSTALD_rho(phase):
    # Vm = COSTALD_mixture(phase.zs, phase.T, phase.Tcs, phase.Vcs, phase.omegas)
    Vm = COSTALD_Vm(phase)
    rho = 1 / Vm * phase.MW() / 1000
    return rho


def COSTALD_Vm(phase):
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
    return Vm  # COSTALD_mixture(phase.zs, phase.T, phase.Tcs, phase.Vcs, phase.omegas)


def h_inside(L, Tvessel, Tfluid, fluid):
    """
    Calculation of internal natural convective heat transfer coefficient from Nusselt number

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    fluid : obj
            Gas object equilibrated at film temperature

    Returns
    ----------
    h_inner : float
        Heat transfer coefficient (W/m2 K)
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
    Calculation of internal natural convective heat transfer coefficient from Nusselt number

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    fluid : obj
            Gas object equilibrated at film temperature

    Returns
    ----------
    h_inner : float
        Heat transfer coefficient (W/m2 K)
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
    Calculation of internal heat transfer coefficient for boiling liquid

    Parameters
    ----------
    L : float
        Vessel length
    Tfluid : float
        Temperature of the bulk fluid inventory
    Tvessel : float
        Temperature of the vessel wall (bulk)
    fluid : obj
            Gas object equilibrated at film temperature

    Returns
    ----------
    h_inner : float
        Heat transfer coefficient (W/m2 K)
    """
    liq = fluid.liquid_bulk
    # if liq.beta == 0:
    #    return 0

    kl = sum([k * w for k, w in zip(fluid.liquid0.ws(), fluid.liquid0.kls())])

    h_boil = Rohsenow(
        rhol=liq.rho_mass(),
        rhog=fluid.gas.rho_mass(),
        mul=liq.mu(),
        kl=kl,
        Cpl=liq.Cp_mass(),
        Hvap=(fluid.gas.H_mass() - liq.H_mass()),
        sigma=liq.sigma(),
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
    Calculation of Prandtl number, eq. 4.5-6 in
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993

    Parameters
    ----------

    Returns
    ----------
    Pr : float
        Prantdl number
    """
    return cp * mu / k


def Nu(Ra, Pr):
    """
    Calculation of Nusselt number for natural convection. See eq. 4.7-4  and Table 4.7-1 in
    C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993

    Parameters
    ----------
    Ra : float
        Raleigh number
    Pr : float
        Prandtl number

    Returns
    ----------
    Nu : float
        Nusselt number
    """
    if Ra >= 1e9:
        NNu = 0.13 * Ra**0.333
    elif Ra < 1e9 and Ra > 1e4:
        NNu = 0.59 * Ra**0.25
    else:
        NNu = 1.36 * Ra**0.20
    return NNu
