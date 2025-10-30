import math
import numpy as np


# Material UTS data from Scandpower guideline
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
    von Mises stress calculated according to:
    Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised
    Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.

    Parameters
    ----------
    p : float
        Pressure (Pa)
    d : float
        Inner diameter (m)
    D : float
        Outer diameter (m)
    wt: float
        Wall thickness (m)
    sigma_a: float
        Default

    Returns
    ----------
    sigma_e : von Mises stress (Pa)

    """

    D = d + 2 * wt

    sigma_e = math.sqrt(3 * ((p * D**2) / (D**2 - d**2)) ** 2 + sigma_a)
    return sigma_e


def ATS(temperature, material, k_s=0.85, k_y=1):
    """
    Calculation of Allowable Tensile Strength according to:
    Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised
    Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.

    Parameters
    ----------
    temperature : float
        Temperature (K)
    material : string
        Material type: 235LT, 360LT (ASTM A-333/A-671), 2205 (SA-790/ASTM A-790),
        316 (ASTM A-320, ASME A-358), 6Mo (ASTM B-677)
    k_s : float
        General safety factor. For typical materials 0.85 is used. If "guaranteed" minimum
        values a factor 1.0 can be used.
    k_y : float
        Additional factor used for materials with missing or uncertain material data.
        Normally 1.0.

    Returns
    ----------
    ATS : Allowable Tensile Strength (Pa)

    """

    return UTS(temperature, material=material) * k_s * k_y


def UTS(temperature, material):
    """
    Tabulation look-up / interpolation to retrieve the Ultimate Tensile Strength
    as a function of temperature for various typical materials according to:

    Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised
    Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.


    Parameters
    ----------
    temperature : float
        Temperature (K)
    material : string
        Material type: 235LT, 360LT (ASTM A-333/A-671), Duplex (2205, SA-790/ASTM A-790),
        316 (ASTM A-320, ASME A-358), 6Mo (ASTM B-677)

    Return
    ----------
    UTS : Ultimate Tensile Strength (Pa)

    """

    if material == "CS_235LT":
        UTS = np.interp(temperature, T, CS_235LT_UTS)
        return UTS
    elif material == "CS_360LT":
        UTS = np.interp(temperature, T, CS_360LT_UTS)
        return UTS
    elif material == "SS316":
        UTS = np.interp(temperature, T, SS_UTS)
        return UTS
    elif material == "Duplex":
        UTS = np.interp(temperature, T, Duplex_UTS)
        return UTS
    elif material == "6Mo":
        UTS = np.interp(temperature, T, SMo_UTS)
        return UTS


if __name__ == "__main__":
    from matplotlib import pyplot as plt

    plt.plot(T, Duplex_UTS, "--", label="22Cr Duplex")
    plt.plot(T, SS_UTS, "-.", label="SS316")
    plt.plot(T, SMo_UTS, "-", label="6Mo")
    plt.plot(T, CS_235LT_UTS, "k--", label="235LT")
    plt.plot(T, CS_360LT_UTS, "k-.", label="360LT")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Ultimate Tensile Strength (Pa)")
    plt.legend(loc="best")
    plt.show()
