import math
from scipy.constants import g

psv_state = "closed"


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


def relief_valve(P1, Pback, Pset, blowdown, k, CD, T1, Z, MW, area):
    """
    Pop action relief valve model including hysteresis.
    The pressure shall rise above P_set to open and
    decrease below P_reseat (P_set*(1-blowdown)) to close

    Parameters
    ----------
    P1 : float
        Upstream pressure
    Pback : float
        Downstream / backpressure
    Pset : float
        Set pressure of the PSV / relief valve
    blowdown : float
        The percentage of the set pressure at which the valve reseats
    k : float
        Ideal gas k (Cp/Cv)
    CD : float
        Coefficient of discharge
    T1 : float
        Upstream temperature
    Z : float
        Compressibility
    MW : float
        Molecular weight of the gas relieved
    area : float
        PSV orifice area

    Returns
    ----------
        : float
        Relief rate / mass flow
    """

    global psv_state
    if P1 > Pset:
        eff_area = area
        psv_state = "open"
    elif P1 < Pset * (1 - blowdown):
        eff_area = 0
        psv_state = "closed"
    else:
        if psv_state == "open":
            eff_area = area
        elif psv_state == "closed":
            eff_area = 0
        else:
            raise ValueError("Unknown PSV open/close state.")

    if eff_area > 0:
        return api_psv_release_rate(P1, Pback, k, CD, T1, Z, MW, area)
    else:
        return 0.0


def api_psv_release_rate(P1, Pback, k, CD, T1, Z, MW, area):
    """
    PSV vapour relief rate calculated according to API 520 Part I 2014
    Eq. 5, 9, 15, 18

    Parameters
    ----------
    P1 : float
        Upstream pressure
    Pback : float
        Downstream / backpressure
    k : float
        Ideal gas k (Cp/Cv)
    CD : float
        Coefficient of discharge
    T1 : float
        Upstream temperature
    Z : float
        Compressibility
    MW : float
        Molecular weight of the gas relieved
    area : float
        PSV orifice area

    Returns
    ----------
        : float
        Relief rate / mass flow
    """

    # Convert units for API 520 equations
    P1 = P1 / 1000  # Pa to kPa
    Pback = Pback / 1000  # Pa to kPa
    area = area * 1e6  # m² to mm²
    MW = MW * 1000  # kg/mol to g/mol

    # API 520 critical flow coefficient
    C = 0.03948 * (k * (2 / (k + 1)) ** ((k + 1) / (k - 1))) ** 0.5

    # Check if flow is critical (choked) or subcritical
    if P1 / Pback > ((k + 1) / 2) ** ((k) / (k - 1)):
        # Critical flow (choked at throat)
        w = CD * area * C * P1 / math.sqrt(T1 * Z / MW)
    else:
        # Subcritical flow (not choked)
        r = Pback / P1
        # Subcritical flow correction factor f2
        f2 = ((k / (k - 1)) * r ** (2 / k) * (1 - r ** ((k - 1) / k)) / (1 - r)) ** 0.5
        w = CD * area * f2 / (T1 * Z / (MW * P1 * (P1 - Pback))) ** 0.5 / 17.9
    return w / 3600


def cv_vs_time(Cv_max, t, time_constant=0, characteristic="linear"):
    """
    Control valve flow coefficient vs time / actuator postion
    assuming a linear rate of actuator for the three archetypes of
    characteristics: linear, equal percentage and fast/quick opening.

    Parameters
    ----------
    Cv_max : float
        Valve flow coefficient at full open position
    t : float
        Time
    time_constant : float (optional)
        The time required for the actuator to fully open.
        Default to instant open
    characteristic : string (optional)
        Valve characteristic
        Default to linear.
    """

    if time_constant == 0:
        return Cv_max
    else:
        if characteristic == "linear":
            return Cv_max * min(t / time_constant, 1)
        elif characteristic == "eq":
            # https://www.spiraxsarco.com/learn-about-steam/control-hardware-electric-pneumatic-actuation/control-valve-characteristics
            tau = 50
            travel = min(t / time_constant, 1)
            return Cv_max * math.exp(travel * math.log(tau)) / tau
        elif characteristic == "fast":
            # square root function used
            return Cv_max * min(t / time_constant, 1) ** (0.5)
        else:
            return Cv_max


def control_valve(P1, P2, T, Z, MW, gamma, Cv, xT=0.75, FP=1):
    """
    Flow calculated from ANSI/ISA control valve equations for single phase gas flow.
    Equation 19 pp. 132 in
    Control Valves / Guy Borden, editor; Paul Friedmann, style editor

    Parameters
    ----------
    P1 : float
        Upstream pressure
    P2 : float
        Downstream / backpressure
    T : float
        Upstream temperature
    Z : float
        Upstream compressibility
    MW : float
        Molecular weight of the gas relieved
    gamma : float
        Upstream Ideal gas k (Cp/Cv)
    Cv : float
        Valve coefficient
    xT : float
        Value of xT for valve fitting assembly, default value
    FP : float
        Piping geometry factor

    Returns
    ----------
        : float
        Mass flow
    """

    P1 = P1 / 1e5
    P2 = P2 / 1e5
    MW = MW * 1000
    N8 = 94.8
    Fk = gamma / 1.4
    x = (P1 - P2) / P1
    if x < 0:
        x = 0
    Y = 1.0 - min(x, Fk * xT) / (3.0 * Fk * xT)
    mass_flow = N8 * FP * Cv * P1 * Y * (MW * min(x, xT * Fk) / T / Z) ** 0.5
    return mass_flow / 3600  # kg/s
