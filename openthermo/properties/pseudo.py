from math import log, exp
from scipy.constants import psi


def Tc_Kesler_Lee_SG_Tb(SG, Tb):
    r"""Estimates critical temperature of a hydrocarbon compound or petroleum
    fraction using only its specific gravity and boiling point, from
    [1]_ as presented in [2]_.

    .. math::
        T_c = 341.7 + 811.1SG + [0.4244 + 0.1174SG]T_b
        + \frac{[0.4669 - 3.26238SG]10^5}{T_b}

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]
    Tb : float
        Boiling point the fluid [K]

    Returns
    -------
    Tc : float
        Estimated critical temperature [K]

    Notes
    -----
    Model shows predictions for Tc, Pc, MW, and omega.
    Original units in degrees Rankine.

    Examples
    --------
    Example 2.2 from [2]_, but with K instead of R.

    >>> Tc_Kesler_Lee_SG_Tb(0.7365, 365.555)
    545.0124354151242

    References
    ----------
    .. [1] Kesler, M. G., and B. I. Lee. "Improve Prediction of Enthalpy of
       Fractions." Hydrocarbon Processing (March 1976): 153-158.
    .. [2] Ahmed, Tarek H. Equations of State and PVT Analysis: Applications
       for Improved Reservoir Modeling. Gulf Pub., 2007.
    """
    Tb = 9 / 5.0 * Tb  # K to R
    Tc = (
        341.7
        + 811.1 * SG
        + (0.4244 + 0.1174 * SG) * Tb
        + ((0.4669 - 3.26238 * SG) * 1e5) / Tb
    )
    Tc = 5 / 9.0 * Tc  # R to K
    return Tc


def Tc_Riazi_Daubert_SG_Tb(SG, Tb):
    r"""Estimates critical temperature of a hydrocarbon compound or petroleum
    fraction using only its specific gravity and boiling point, from
    [1]_ as presented in [2]_.

    .. math::
        T_c = 19.06232 T_b^{0.58848} SG^0.3596

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]
    Tb : float
        Boiling point the fluid [K]

    Returns
    -------
    Tc : float
        Estimated critical temperature [K]

    Examples
    --------
    Example 2.2 from [2]_

    >>> Tc_Riazi_Daubert_SG_Tb(0.7365, 365.555)
    550.3734796518954

    References
    ----------
    .. [1] Riazi, M. R. and Daubert, T. E.,
       Simplify Property Predictions, Hydrocarbon Processing, Vol. 59, No. 3, 1980, pp. 115-116.
    .. [2] Riazi, M.R, Characterization and Properties of Petroleum Fractions
       American Society for Testing and Materials, 2005.
    """
    Tc = 19.06232 * Tb**0.58848 * SG**0.3596
    return Tc


def Pc_Kesler_Lee_SG_Tb(SG, Tb):
    r"""Estimates critical pressure of a hydrocarbon compound or petroleum
    fraction using only its specific gravity and boiling point, from
    [1]_ as presented in [2]_.

    .. math::
        \ln(P_c) = 8.3634 - \frac{0.0566}{SG} - \left[0.24244 + \frac{2.2898}
        {SG} + \frac{0.11857}{SG^2}\right]10^{-3}T_b
        + \left[1.4685 + \frac{3.648}{SG} + \frac{0.47227}{SG^2}\right]
        10^{-7}T_b^2-\left[0.42019 + \frac{1.6977}{SG^2}\right]10^{-10}T_b^3

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]
    Tb : float
        Boiling point the fluid [K]

    Returns
    -------
    Pc : float
        Estimated critical pressure [Pa]

    Notes
    -----
    Model shows predictions for Tc, Pc, MW, and omega.
    Original units in degrees Rankine and psi.

    Examples
    --------
    Example 2.2 from [2]_, but with K instead of R and Pa instead of psi.

    >>> Pc_Kesler_Lee_SG_Tb(0.7365, 365.555)
    3238323.346840464

    References
    ----------
    .. [1] Kesler, M. G., and B. I. Lee. "Improve Prediction of Enthalpy of
       Fractions." Hydrocarbon Processing (March 1976): 153-158.
    .. [2] Ahmed, Tarek H. Equations of State and PVT Analysis: Applications
       for Improved Reservoir Modeling. Gulf Pub., 2007.
    """
    Tb = 9 / 5.0 * Tb  # K to R
    Pc = exp(
        8.3634
        - 0.0566 / SG
        - (0.24244 + 2.2898 / SG + 0.11857 / SG**2) * 1e-3 * Tb
        + (1.4685 + 3.648 / SG + 0.47227 / SG**2) * 1e-7 * Tb**2
        - (0.42019 + 1.6977 / SG**2) * 1e-10 * Tb**3
    )
    Pc = Pc * psi
    return Pc


def Pc_Riazi_Daubert_SG_Tb(SG, Tb):
    r"""Estimates critical pressure of a hydrocarbon compound or petroleum
    fraction using only its specific gravity and boiling point, from
    [1]_ as presented in [2]_.

    .. math::
        Pc = 5.53027\cdot 10^7 Tb^{-2.3125} SG^{2.3201}

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]
    Tb : float
        Boiling point the fluid [K]

    Returns
    -------
    Pc : float
        Estimated critical pressure [Pa]

    Notes
    -----
    Model shows predictions for Tc, Pc, MW, and omega.
    Original units in degrees Rankine and psi.

    Examples
    --------
    Example 2.2 from [2]_, but with K instead of R and Pa instead of psi.

    >>> Pc_Riazi_Daubert_SG_Tb(0.7365, 365.555)
    3219182.887436976

    References
    ----------
    .. [1] Riazi, M. R. and Daubert, T. E.,
       Simplify Property Predictions, Hydrocarbon Processing, Vol. 59, No. 3, 1980, pp. 115-116.
    .. [2] Riazi, M.R, Characterization and Properties of Petroleum Fractions
       American Society for Testing and Materials, 2005.
    """
    Pc = 5.53027e7 * Tb**-2.3125 * SG**2.3201
    return Pc * 1e5


def MW_Kesler_Lee_SG_Tb(SG, Tb):
    r"""Estimates molecular weight of a hydrocarbon compound or petroleum
    fraction using only its specific gravity and boiling point, from
    [1]_ as presented in [2]_.

    .. math::
        MW = -12272.6 + 9486.4SG + [4.6523 - 3.3287SG]T_b + [1-0.77084SG
        - 0.02058SG^2]\left[1.3437 - \frac{720.79}{T_b}\right]\frac{10^7}{T_b}
        + [1-0.80882SG + 0.02226SG^2][1.8828 - \frac{181.98}{T_b}]
        \frac{10^{12}}{T_b^3}

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]
    Tb : float
        Boiling point the fluid [K]

    Returns
    -------
    MW : float
        Estimated molecular weight [g/mol]

    Notes
    -----
    Model shows predictions for Tc, Pc, MW, and omega.
    Original units in degrees Rankine.

    Examples
    --------
    Example 2.2 from [2]_, but with K instead of R and Pa instead of psi.

    >>> MW_Kesler_Lee_SG_Tb(0.7365, 365.555)
    98.70887589833501

    References
    ----------
    .. [1] Kesler, M. G., and B. I. Lee. "Improve Prediction of Enthalpy of
       Fractions." Hydrocarbon Processing (March 1976): 153-158.
    .. [2] Ahmed, Tarek H. Equations of State and PVT Analysis: Applications
       for Improved Reservoir Modeling. Gulf Pub., 2007.
    """
    Tb = 9 / 5.0 * Tb  # K to R
    MW = (
        -12272.6
        + 9486.4 * SG
        + (4.6523 - 3.3287 * SG) * Tb
        + (1.0 - 0.77084 * SG - 0.02058 * SG**2) * (1.3437 - 720.79 / Tb) * 1e7 / Tb
        + (1.0 - 0.80882 * SG + 0.02226 * SG**2) * (1.8828 - 181.98 / Tb) * 1e12 / Tb**3
    )
    return MW


def omega_Kesler_Lee_SG_Tb_Tc_Pc(SG, Tb, Tc=None, Pc=None):
    r"""Estimates accentric factor of a hydrocarbon compound or petroleum
    fraction using only its specific gravity and boiling point, from
    [1]_ as presented in [2]_. If Tc and Pc are provided, the Kesler-Lee
    routines for estimating them are not used.

    For Tbr > 0.8:
    .. math::
        \omega = -7.904 + 0.1352K - 0.007465K^2 + 8.359T_{br}
        + ([1.408-0.01063K]/T_{br})

    Otherwise:
    .. math::
        \omega = \frac{-\ln\frac{P_c}{14.7} - 5.92714 + \frac{6.09648}{T_{br}}
        + 1.28862\ln T_{br} - 0.169347T_{br}^6}{15.2518 - \frac{15.6875}{T_{br}}
         - 13.4721\ln T_{br} + 0.43577T_{br}^6}

        K = \frac{T_b^{1/3}}{SG}

        T_{br} = \frac{T_b}{T_c}

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]
    Tb : float
        Boiling point the fluid [K]
    Tc : float, optional
        Estimated critical temperature [K]
    Pc : float, optional
        Estimated critical pressure [Pa]

    Returns
    -------
    omega : float
        Acentric factor [-]

    Notes
    -----
    Model shows predictions for Tc, Pc, MW, and omega.
    Original units in degrees Rankine and psi.

    Examples
    --------
    Example 2.2 from [2]_, but with K instead of R and Pa instead of psi.

    >>> omega_Kesler_Lee_SG_Tb_Tc_Pc(0.7365, 365.555, 545.012, 3238323.)
    0.306392118159797

    References
    ----------
    .. [1] Kesler, M. G., and B. I. Lee. "Improve Prediction of Enthalpy of
       Fractions." Hydrocarbon Processing (March 1976): 153-158.
    .. [2] Ahmed, Tarek H. Equations of State and PVT Analysis: Applications
       for Improved Reservoir Modeling. Gulf Pub., 2007.
    """
    if Tc is None:
        Tc = Tc_Kesler_Lee_SG_Tb(SG, Tb)
    if Pc is None:
        Pc = Pc_Kesler_Lee_SG_Tb(SG, Tb)
    Tb = 9 / 5.0 * Tb  # K to R
    Tc = 9 / 5.0 * Tc  # K to R
    K = Tb ** (1 / 3.0) / SG
    Tbr = Tb / Tc
    if Tbr > 0.8:
        omega = (
            -7.904
            + 0.1352 * K
            - 0.007465 * K**2
            + 8.359 * Tbr
            + ((1.408 - 0.01063 * K) / Tbr)
        )
    else:
        omega = (
            -log(Pc / 101325.0)
            - 5.92714
            + 6.09648 / Tbr
            + 1.28862 * log(Tbr)
            - 0.169347 * Tbr**6
        ) / (15.2518 - 15.6875 / Tbr - 13.4721 * log(Tbr) + 0.43577 * Tbr**6)
    return omega


def HC_atomic_ratio(SG, Tb):
    r"""Estimates the HC atomic ratio for petroleum fractions from specific gravity and
    boiling point, from [1]_ modified in [2]_. Applicable for C6-C50

    The CH weight ratio (Carbon-to-hydrogen ratio) is calculated from:
    .. math::
        CH = 8.7743\cdot10^{-10} \right[ \exp{7.176 \cdot 10^{-3}T_b + 30.06242 SG -7.35\cdot 10^{-3} Tb SG} \left] Tb^{-0.98445}SG^{-18.2753}

    The Hydrogen-to-Carbon ratio is calculated from:
    .. math::
        HC_atomic_ratio = 11.9147 / CH

    Parameters
    ----------
    SG : float
        Specific gravity of the fluid at 60 degrees Farenheight [-]
    Tb : float
        Boiling point the fluid [K]

    Returns
    -------
    HC_atomic ratio : float
        The ratio of moles H atoms to C atoms in the petroleum fraction [-]

    Examples
    --------
    Example from sec. 2.6.3 [2]_ for n-tetradecylbenzene (C20H24) with HC
    atimic ratio of 1.7

    >>> HC_atomic_ratio(0.8587,627)
    1.7022132444770803

    References
    ----------
    .. [1] Riazi, M. R. and Daubert, T. E., Prediction of Molecular Type
       Analysis of Petroleum Fractions and Coal Liquids, Industrial
       and Engineering Chemistry, Process Design and Development,
       Vol. 25, No. 4, 1986, pp. 1009-1015.
    .. [2] Riazi, M.R, Characterization and Properties of Petroleum Fractions
       American Society for Testing and Materials, 2005.
    """

    res = (
        8.7743e-10
        * (exp(7.176e-3 * Tb + 30.06242 * SG - 7.35e-3 * Tb * SG))
        * Tb**-0.98445
        * SG**-18.2753
    )
    res = 11.914683 / res
    return res


def Zc_pseudo(omega):
    """
    The function calculates critical compresibility [-] for pseudo components based on equation 21 7 in [1]_

    Parameters
    ----------
    omega : float
        Accentric factor

    Returns
    -------
    Zc : float
        Critical compressibility

    References
    ----------
    .. [1] B. I. Lee and M. G. Kessler, A Generalized thermodynamic correlation based on three-parameter
        corresponding states, AIChE Journal 1975,  https://doi.org/10.1002/aic.690210313
    """
    zc_pseudo = 0.2905 - 0.085 * omega
    return zc_pseudo


def Vc_pseudo(Zc, Tc, Pc):
    """
    The function calculates critical volume [m3/mol] for pseudo components based on definition of compresibility.

    Input: critical compresibility ZC_i [-], critical temperature TC_i [K], critical pressure PC_i [Pa].
    """
    Vc = Zc * Tc * 8.314 / Pc
    return Vc
