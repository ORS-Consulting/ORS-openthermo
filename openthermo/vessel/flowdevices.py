import math
from scipy.constants import g


def hem_release_rate():
    pass


def liquid_release_bernouilli(P1, P2, rho, CD, area, H):
    """
    Liquid mass flow (kg/s) trough a hole or orifice.
    flow conditions. The formula is based on Yellow Book equation 2.194.

    Methods for the calculation of physical effects, CPR 14E, van den Bosch and Weterings (Eds.), 1996

    Parameters
    ----------
    P1 : float
        Upstream pressure
    P2 : float
        Downstream pressure
    rho : float
        Fluid density
    CD : float
        Coefficient of discharge
    area : float
        Orifice area
    H : float
        Liquid static height

    Returns
    ----------
        : float
        liquid  release rate / mass flow of discharge
    """
    P1 = P1 + H * rho * g
    if H > 0 and P1 > P2:
        return CD * area * math.sqrt(2 * (P1 - P2) * rho)
    else:
        return 0


def two_phase_release_fauske(P1, Pc, rho, CD, area):
    """
    Two-phase mass flow (kg/s) trough a hole or orifice,.
    flow conditions. The formula is based on Yellow Book equation 2.91.

    Methods for the calculation of physical effects, CPR 14E, van den Bosch and Weterings (Eds.), 1996
    World Bank Technical Paper Number 55, Techniques for Assessing Industrial Hazards, 2nd print, 1990

    Parameters
    ----------
    P1 : float
        Upstream pressure
    Pc : float
        Critical pressure or ambient which ever is the highest. Pc = 0.55 P1
    rho : float
        Fluid density
    CD : float
        Coefficient of discharge
    are : float
        Orifice area

    Returns
    ----------
        : float
        Two-phase   release rate / mass flow of discharge
    """

    return max(CD * area * math.sqrt(2 * (P1 - Pc) * rho), 0)


def gas_release_rate(P1, P2, rho, k, CD, area):
    """
    Gas massflow (kg/s) trough a hole at critical (sonic) or subcritical
    flow conditions. The formula is based on Yellow Book equation 2.22.

    Methods for the calculation of physical effects, CPR 14E,
    van den Bosch and Weterings (Eds.), 1996

    Parameters
    ----------
    P1 : float
        Upstream pressure
    P2 : float
        Downstream pressure
    rho : float
        Fluid density
    k : float
        Ideal gas k (Cp/Cv)
    CD : float
        Coefficient of discharge
    are : float
        Orifice area

    Returns
    ----------
        : float
        Gas release rate / mass flow of discharge
    """
    if P1 > P2:
        if P1 / P2 > ((k + 1) / 2) ** ((k) / (k - 1)):
            flow_coef = 1
        else:
            flow_coef = (
                2
                / (k - 1)
                * (((k + 1) / 2) ** ((k + 1) / (k - 1)))
                * ((P2 / P1) ** (2 / k))
                * (1 - (P2 / P1) ** ((k - 1) / k))
            )

        retval = (
            math.sqrt(flow_coef)
            * CD
            * area
            * math.sqrt(rho * P1 * k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
        )
    else:
        retval = 0

    return retval
