"""
Material property database for pressure vessel calculations.

This module provides temperature-dependent material properties for common
pressure vessel materials including:
- Stainless steels (SS316, Duplex, Super Duplex)
- Carbon steels (low temperature grade)

Properties available:
- Heat capacity (Cp) as function of temperature [J/(kg·K)]
- Ultimate tensile strength (UTS) as function of temperature [Pa]
- Allowable tensile stress (ATS) as function of temperature [Pa]
- von Mises equivalent stress calculations

Data sources:
- Scandpower Risk Management AS guidelines for fire scenarios
- EN standards for structural steel properties

All temperature data are stored in Kelvin [K].
Property values are interpolated using numpy.interp() for intermediate temperatures.

Available materials:
- CS_235LT: Carbon Steel 235 Low Temperature grade (ASTM A-333/A-671)
- CS_360LT: Carbon Steel 360 Low Temperature grade (ASTM A-333/A-671)
- SS316: Stainless Steel 316 (ASTM A-320, ASME A-358)
- Duplex: Duplex Stainless Steel 2205 (SA-790/ASTM A-790)
- 6Mo: Super Duplex SMo254 (ASTM B-677)

References
----------
Hekkelstrand, B., & Skulstad, P. (2004). Guidelines for the Protection of
Pressurised Systems Exposed to Fire. Scandpower Risk Management AS, Kjeller, Norway.

Andreasen, A., Borroni, F., Zan Nieto, M., Stegelmann, C., & Nielsen, R. P. (2018).
On the Adequacy of API 521 Relief-Valve Sizing Method for Gas-Filled Pressure Vessels
Exposed to Fire. Safety, 4(1), 11. https://doi.org/10.3390/safety4010011
"""

import math
import numpy as np

# ============================================================================
# HEAT CAPACITY DATA
# ============================================================================
# Temperature array for heat capacity data (converted from °C to K)
# Source: Scandpower Risk Management AS guideline
T_Cp = (
    np.array([20, 100, 200, 300, 400, 500, 600, 700, 750, 800, 900, 1000, 1100])
    + 273.15
)

# Stainless Steel 316 heat capacity [J/(kg·K)] at temperatures in T_Cp
SS316_Cp = np.array((472, 487, 503, 512, 520, 530, 541, 551, 555, 559, 565, 571, 577))

# Duplex stainless steel heat capacity [J/(kg·K)] at temperatures in T_Cp
# Missing points from Scandpower guideline were added manually based on interpolation
Duplex_Cp = np.array([480, 500, 530, 560, 600, 635, 670, 710, 730, 750, 790, 840, 840])

# Super Duplex (SMo254) stainless steel heat capacity [J/(kg·K)] at temperatures in T_Cp
SMo_Cp = np.array([500, 520, 540, 555, 570, 580, 590, 600, 610, 610, 610, 610, 610])

# Carbon Steel Low Temperature grade heat capacity [J/(kg·K)] at temperatures in T_Cp
# Note: Discontinuity at 750°C due to phase transformation (austenite formation)
CS_LT_Cp = np.array([450, 480, 510, 550, 600, 660, 750, 900, 1450, 820, 540, 540, 540])

# ============================================================================
# ULTIMATE TENSILE STRENGTH (UTS) DATA
# ============================================================================
# Temperature array for UTS and ATS data (converted from °C to K)
# Source: Scandpower Risk Management AS guideline
T = (
    np.array(
        [
            20,
            50,
            100,
            150,
            200,
            250,
            300,
            350,
            400,
            450,
            500,
            550,
            600,
            650,
            700,
            750,
            800,
            900,
            1000,
            1100,
        ]
    )
    + 273.17
)

# Duplex stainless steel ultimate tensile strength [Pa] at temperatures in T
# Converted from MPa to Pa by multiplying by 1e6
Duplex_UTS = (
    np.array(
        [
            730,
            701,
            668,
            650,
            640,
            630,
            621,
            606,
            591,
            540,
            482,
            423,
            358,
            299,
            234,
            164,
            124,
            88,
            69,
            58,
        ]
    )
    * 1e6
)

# Stainless steel 316 ultimate tensile strength [Pa] at temperatures in T
SS_UTS = (
    np.array(
        [
            575,
            549,
            523,
            503,
            489,
            477,
            472,
            466,
            463,
            460,
            449,
            431,
            397,
            357,
            299,
            242,
            184,
            98,
            70,
            60,
        ]
    )
    * 1e6
)

# Super Duplex (SMo254) stainless steel ultimate tensile strength [Pa] at temperatures in T
SMo_UTS = (
    np.array(
        [
            730,
            710,
            680,
            660,
            645,
            633,
            625,
            618,
            610,
            595,
            580,
            550,
            510,
            445,
            380,
            305,
            230,
            130,
            90,
            65,
        ]
    )
    * 1e6
)

# Carbon Steel 235 Low Temperature grade ultimate tensile strength [Pa] at temperatures in T
# Note: Last two values are manual hand-interpolations where data was unavailable
CS_235LT_UTS = (
    np.array(
        [
            420,
            414,
            407,
            403,
            397,
            391,
            382,
            378,
            370,
            353,
            308,
            252,
            189,
            139,
            92,
            59,
            46,
            38,
            30,  # Manual hand-interpolation
            22,  # Manual hand-interpolation
        ]
    )
    * 1e6
)

# Carbon Steel 360 Low Temperature grade ultimate tensile strength [Pa] at temperatures in T
CS_360LT_UTS = (
    np.array(
        [
            545,
            537,
            529,
            523,
            515,
            507,
            496,
            491,
            480,
            458,
            400,
            327,
            245,
            180,
            120,
            76,
            60,
            49,
            38,  # Manual hand-interpolation
            27,  # Manual hand-interpolation
        ]
    )
    * 1e6
)


def von_mises(p, d, wt, sigma_a=30e6):
    """
    Calculate von Mises equivalent stress for thin-walled pressure vessel.

    Computes the von Mises stress for a cylindrical pressure vessel subjected
    to internal pressure. This is used to assess vessel integrity under fire
    exposure conditions.

    Parameters
    ----------
    p : float
        Internal pressure [Pa]
    d : float
        Inner diameter [m]
    wt : float
        Wall thickness [m]
    sigma_a : float, optional
        Axial stress component [Pa]. Default is 30 MPa (typical for fire scenarios).

    Returns
    -------
    float
        von Mises equivalent stress [Pa]

    Notes
    -----
    The von Mises stress is calculated as:

    σ_e = √(3 · (p·D²/(D²-d²))² + σ_a)

    where D = d + 2·wt is the outer diameter.

    This formulation accounts for:
    - Hoop stress from internal pressure
    - Axial stress from end caps and thermal expansion
    - Thin-wall assumption (wt << d)

    The von Mises stress is compared against the allowable tensile strength
    (ATS) to determine if the vessel remains within safe operating limits.

    Examples
    --------
    >>> # 0.5 m ID vessel, 10 mm wall, 100 bar pressure
    >>> p = 100e5  # Pa
    >>> d = 0.5  # m
    >>> wt = 0.01  # m
    >>> stress = von_mises(p, d, wt)
    >>> print(f"von Mises stress: {stress/1e6:.1f} MPa")

    References
    ----------
    Hekkelstrand, B., & Skulstad, P. (2004). Guidelines for the Protection of
    Pressurised Systems Exposed to Fire. Scandpower Risk Management AS, Kjeller, Norway.

    Andreasen, A., Borroni, F., Zan Nieto, M., Stegelmann, C., & Nielsen, R. P. (2018).
    On the Adequacy of API 521 Relief-Valve Sizing Method for Gas-Filled Pressure Vessels
    Exposed to Fire. Safety, 4(1), 11. https://doi.org/10.3390/safety4010011
    """

    D = d + 2 * wt

    sigma_e = math.sqrt(3 * ((p * D**2) / (D**2 - d**2)) ** 2 + sigma_a)
    return sigma_e


def ATS(temperature, material, k_s=0.85, k_y=1):
    """
    Calculate Allowable Tensile Strength (ATS) with safety factors.

    The ATS is derived from the Ultimate Tensile Strength (UTS) by applying
    safety factors to account for material variability, data uncertainty, and
    design margins.

    Parameters
    ----------
    temperature : float
        Material temperature [K]
    material : str
        Material type. Valid options:
        - "CS_235LT": Carbon Steel 235 Low Temperature (ASTM A-333/A-671)
        - "CS_360LT": Carbon Steel 360 Low Temperature (ASTM A-333/A-671)
        - "SS316": Stainless Steel 316 (ASTM A-320, ASME A-358)
        - "Duplex": Duplex Stainless Steel 2205 (SA-790/ASTM A-790)
        - "6Mo": Super Duplex SMo254 (ASTM B-677)
    k_s : float, optional
        General safety factor [-]. Default is 0.85.
        - 0.85: Standard safety factor for typical materials
        - 1.0: Use for "guaranteed" minimum material properties
    k_y : float, optional
        Uncertainty factor for incomplete material data [-]. Default is 1.0.
        - 1.0: Complete material data available
        - <1.0: Apply additional margin for uncertain or incomplete data

    Returns
    -------
    float
        Allowable Tensile Strength [Pa]

    Notes
    -----
    The ATS is calculated as:

    ATS = UTS(T) · k_s · k_y

    The safety factor k_s accounts for:
    - Scatter in material test data
    - Fabrication quality variations
    - Service degradation
    - Design uncertainties

    The uncertainty factor k_y provides additional margin when material
    data is incomplete or uncertain.

    Examples
    --------
    >>> # SS316 at 400°C with standard safety factor
    >>> T = 400 + 273.15  # K
    >>> ats = ATS(T, "SS316", k_s=0.85)
    >>> print(f"Allowable stress: {ats/1e6:.1f} MPa")

    >>> # With additional uncertainty margin
    >>> ats_uncertain = ATS(T, "SS316", k_s=0.85, k_y=0.9)

    References
    ----------
    Hekkelstrand, B., & Skulstad, P. (2004). Guidelines for the Protection of
    Pressurised Systems Exposed to Fire. Scandpower Risk Management AS, Kjeller, Norway.
    """

    return UTS(temperature, material=material) * k_s * k_y


def UTS(temperature, material):
    """
    Retrieve Ultimate Tensile Strength (UTS) as function of temperature.

    Performs linear interpolation of tabulated UTS data to obtain the ultimate
    tensile strength at any temperature within the data range. UTS decreases
    with increasing temperature as materials lose strength.

    Parameters
    ----------
    temperature : float
        Material temperature [K]
    material : str
        Material type. Valid options:
        - "CS_235LT": Carbon Steel 235 Low Temperature (ASTM A-333/A-671)
        - "CS_360LT": Carbon Steel 360 Low Temperature (ASTM A-333/A-671)
        - "SS316": Stainless Steel 316 (ASTM A-320, ASME A-358)
        - "Duplex": Duplex Stainless Steel 2205 (SA-790/ASTM A-790)
        - "6Mo": Super Duplex SMo254 (ASTM B-677)

    Returns
    -------
    float
        Ultimate Tensile Strength [Pa]

    Notes
    -----
    Temperature range: 293 K (20°C) to 1373 K (1100°C)

    Typical room temperature values:
    - CS_235LT: 420 MPa
    - CS_360LT: 545 MPa
    - SS316: 575 MPa
    - Duplex: 730 MPa
    - 6Mo: 730 MPa

    At elevated temperatures (e.g., 500°C):
    - Carbon steels lose ~30% strength
    - Stainless steels lose ~20% strength

    UTS is used for:
    - Vessel integrity assessment during fire exposure
    - Relief valve sizing (API 521)
    - Failure pressure calculations
    - Material selection for high-temperature service

    Examples
    --------
    >>> # SS316 at room temperature
    >>> T_room = 293.15  # K
    >>> uts_room = UTS(T_room, "SS316")
    >>> print(f"UTS at 20°C: {uts_room/1e6:.0f} MPa")
    UTS at 20°C: 575 MPa

    >>> # SS316 at 500°C (fire exposure)
    >>> T_fire = 773.15  # K
    >>> uts_fire = UTS(T_fire, "SS316")
    >>> print(f"UTS at 500°C: {uts_fire/1e6:.0f} MPa")
    UTS at 500°C: 449 MPa

    References
    ----------
    Hekkelstrand, B., & Skulstad, P. (2004). Guidelines for the Protection of
    Pressurised Systems Exposed to Fire. Scandpower Risk Management AS, Kjeller, Norway.
    """

    if material == "CS_235LT":
        return np.interp(temperature, T, CS_235LT_UTS)
    elif material == "CS_360LT":
        return np.interp(temperature, T, CS_360LT_UTS)
    elif material == "SS316":
        return np.interp(temperature, T, SS_UTS)
    elif material == "Duplex":
        return np.interp(temperature, T, Duplex_UTS)
    elif material == "6Mo":
        return np.interp(temperature, T, SMo_UTS)


def steel_Cp(temperature, material):
    """
    Retrieve specific heat capacity (Cp) as function of temperature.

    Performs linear interpolation of tabulated heat capacity data to obtain
    the specific heat capacity at any temperature within the data range.
    Heat capacity generally increases with temperature.

    Parameters
    ----------
    temperature : float
        Material temperature [K]
    material : str
        Material type. Valid options:
        - "CS_235LT": Carbon Steel 235 Low Temperature (ASTM A-333/A-671)
        - "CS_360LT": Carbon Steel 360 Low Temperature (ASTM A-333/A-671)
        - "SS316": Stainless Steel 316 (ASTM A-320, ASME A-358)
        - "Duplex": Duplex Stainless Steel 2205 (SA-790/ASTM A-790)
        - "6Mo": Super Duplex SMo254 (ASTM B-677)

    Returns
    -------
    float
        Specific heat capacity [J/(kg·K)]

    Notes
    -----
    Temperature range: 293 K (20°C) to 1373 K (1100°C)

    Typical room temperature values:
    - Carbon steels: 450-480 J/(kg·K)
    - SS316: 472 J/(kg·K)
    - Duplex: 480 J/(kg·K)
    - 6Mo: 500 J/(kg·K)

    Important features:
    - Carbon steels show discontinuity at ~750°C due to ferrite-austenite
      phase transformation (peak at 1450 J/(kg·K))
    - Stainless steels show smooth increase with temperature
    - 6Mo plateaus at ~610 J/(kg·K) above 750°C

    Heat capacity is used for:
    - Thermal transient calculations during blowdown
    - Heat-up rates during fire exposure
    - Energy balance calculations
    - Thermal stress analysis

    Examples
    --------
    >>> # SS316 at room temperature
    >>> T_room = 293.15  # K
    >>> cp_room = steel_Cp(T_room, "SS316")
    >>> print(f"Cp at 20°C: {cp_room:.0f} J/(kg·K)")
    Cp at 20°C: 472 J/(kg·K)

    >>> # Carbon steel showing phase transformation
    >>> T_phase = 750 + 273.15  # K
    >>> cp_phase = steel_Cp(T_phase, "CS_235LT")
    >>> print(f"Cp at 750°C: {cp_phase:.0f} J/(kg·K)")
    Cp at 750°C: 1450 J/(kg·K)

    References
    ----------
    Hekkelstrand, B., & Skulstad, P. (2004). Guidelines for the Protection of
    Pressurised Systems Exposed to Fire. Scandpower Risk Management AS, Kjeller, Norway.
    """
    if material == "CS_235LT":
        return np.interp(temperature, T_Cp, CS_LT_Cp)
    elif material == "CS_360LT":
        return np.interp(temperature, T_Cp, CS_LT_Cp)
    elif material == "SS316":
        return np.interp(temperature, T_Cp, SS316_Cp)
    elif material == "Duplex":
        return np.interp(temperature, T_Cp, Duplex_Cp)
    elif material == "6Mo":
        return np.interp(temperature, T_Cp, SMo_Cp)


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    plt.figure(1)
    plt.plot(T_Cp, Duplex_Cp, "--", label="22Cr Duplex")
    plt.plot(T_Cp, SS316_Cp, "-.", label="SS316")
    plt.plot(T_Cp, SMo_Cp, "-", label="6Mo")
    plt.plot(T_Cp, CS_LT_Cp, "k--", label="CS")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Steel heat capacity (J/kg K)")
    plt.savefig("heat_capacity.png", dpi=600)
    plt.legend(loc="best")

    plt.figure(2)
    plt.plot(T, Duplex_UTS, "--", label="22Cr Duplex")
    plt.plot(T, SS_UTS, "-.", label="SS316")
    plt.plot(T, SMo_UTS, "-", label="6Mo")
    plt.plot(T, CS_235LT_UTS, "k--", label="235LT")
    plt.plot(T, CS_360LT_UTS, "k-.", label="360LT")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Ultimate Tensile Strength (Pa)")
    plt.legend(loc="best")
    plt.savefig("UTS.png", dpi=600)
    plt.show()
