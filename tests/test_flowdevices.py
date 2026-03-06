from openthermo.vessel import flowdevices as tp
from CoolProp.CoolProp import PropsSI
import CoolProp.CoolProp as CP
import pytest


def test_orifice():
    P1 = 10.0e5
    P2 = 5.5e5
    D = PropsSI("D", "P", P1, "T", 298.15, "HEOS::N2")
    cpcv = PropsSI("CP0MOLAR", "T", 298.15, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", 298.15, "P", P1, "HEOS::N2"
    )
    assert tp.gas_release_rate(
        P1, P2, D, cpcv, 0.85, 0.01**2 / 4 * 3.1415
    ) == pytest.approx(9.2 / 60, rel=0.001)


def test_orifice1():
    P1 = 10.0e5
    P2 = 6.5e5
    D = PropsSI("D", "P", P1, "T", 298.15, "HEOS::N2")
    cpcv = PropsSI("CP0MOLAR", "T", 298.15, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", 298.15, "P", P1, "HEOS::N2"
    )
    assert tp.gas_release_rate(
        P1, P2, D, cpcv, 0.85, 0.01**2 / 4 * 3.1415
    ) == pytest.approx(9.2 / 60, rel=0.2)


def test_controlvalve():
    P1 = 10.0e5
    P2 = 5.5e5
    T1 = 20.0 + 273.15
    MW = PropsSI("M", "P", P1, "T", T1, "HEOS::N2")
    Z1 = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    assert tp.control_valve(P1, P2, T1, Z1, MW, gamma, 500) == pytest.approx(
        21.92, rel=0.05
    )


def test_cv_vs_time():
    assert tp.cv_vs_time(1, 0.5, time_constant=1, characteristic="linear") == 0.5
    assert tp.cv_vs_time(1, 0.5, time_constant=1, characteristic="eq") == pytest.approx(
        0.14, abs=0.002
    )
    assert tp.cv_vs_time(
        1, 0.5, time_constant=1, characteristic="fast"
    ) == pytest.approx(0.707, abs=0.002)
    assert tp.cv_vs_time(1, 0.5, time_constant=0) == 1.0


def test_psv3():
    Pback = 1e5
    Pset = 18.2e5
    blowdown = 0.1
    P1 = 0.99 * Pset * (1 - blowdown)
    T1 = 100.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M", "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(P1, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area) == 0
    psv_state = "closed"
    assert (
        tp.relief_valve(P1 * 1.01, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area)
        == 0
    )


def test_psv2():
    P1 = 21.0e5
    Pback = 1e5
    Pset = 20.99e5
    blowdown = 0.1
    T1 = 100.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M", "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(
        P1, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == pytest.approx(1046 / 3600, rel=0.03)
    psv_state = "open"
    assert tp.relief_valve(
        Pset * 0.99, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == pytest.approx(1046 / 3600, rel=0.03)


def test_psv():
    P1 = 100e5
    Pback = 1e5
    Pset = 99.2e5
    blowdown = 0.1
    T1 = 25.0 + 273.15
    Z = PropsSI("Z", "P", P1, "T", T1, "HEOS::N2")
    MW = PropsSI("M", "HEOS::N2")
    gamma = PropsSI("CP0MOLAR", "T", T1, "P", P1, "HEOS::N2") / PropsSI(
        "CVMOLAR", "T", T1, "P", P1, "HEOS::N2"
    )
    CD = 0.975
    area = 71e-6
    assert tp.relief_valve(
        P1, Pback, Pset, blowdown, gamma, CD, T1, Z, MW, area
    ) == pytest.approx(1.57, rel=0.02)


def test_api_psv_relief():
    assert tp.api_psv_release_rate(
        121.9e5, 71e5, 1.39, 0.975, 298.15, 1.01, 2 / 1e3, 71e-6
    ) == pytest.approx(1846 / 3600, rel=0.01)
    assert tp.api_psv_release_rate(
        121.9e5, 1e5, 1.39, 0.975, 298.15, 1.01, 2 / 1e3, 71e-6
    ) == pytest.approx(1860 / 3600, rel=0.01)
