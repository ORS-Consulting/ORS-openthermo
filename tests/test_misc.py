import pytest
import ht
import math
from openthermo.properties.transport import h_inside, h_inside_liquid
import openthermo.properties.transport as tp
from openthermo.flash.michelsen import get_flash_dry
from openthermo.vessel import fire
from openthermo.vessel.blowdown import (
    liquid_release_bernouilli,
    two_phase_release_fauske,
    gas_release_rate,
)
from openthermo.properties.pseudo import (
    Tc_Riazi_Daubert_SG_Tb,
    Tc_Kesler_Lee_SG_Tb,
    Pc_Kesler_Lee_SG_Tb,
    Pc_Riazi_Daubert_SG_Tb,
    MW_Kesler_Lee_SG_Tb,
    omega_Kesler_Lee_SG_Tb_Tc_Pc,
    HC_atomic_ratio,
)


def test_HC_atmic_ratio():
    assert 1.7022132444770803 == HC_atomic_ratio(0.8587, 627)


def test_omega_Lee_Kesler_SG_Tb():
    assert 0.306392118159797 == omega_Kesler_Lee_SG_Tb_Tc_Pc(
        0.7365, 365.555, 545.012, 3238323.0
    )


def test_MW_Kesler_Lee_SG_Tb():
    assert 98.70887589833501 == MW_Kesler_Lee_SG_Tb(0.7365, 365.555)


def test_Pc_Riazi_Daubert_SG_Tb():
    assert 3219182.887436976 == Pc_Riazi_Daubert_SG_Tb(0.7365, 365.555)


def test_Pc_Kesler_Lee_SG_Tb():
    assert 3238323.346840464 == Pc_Kesler_Lee_SG_Tb(0.7365, 365.555)


def test_Tc_Lee_Kessler_SG_Tb():
    assert 545.0124354151242 == Tc_Kesler_Lee_SG_Tb(0.7365, 365.555)


def test_Tc_Riazi_Daubert_SG_Tb():
    assert 550.3734796518954 == Tc_Riazi_Daubert_SG_Tb(0.7365, 365.555)


def test_NNu():
    assert tp.Nu(1.27e8, 0) == pytest.approx(0.59 * (1.27e8) ** 0.25, abs=1.0)


def test_NPr():
    assert tp.Pr(120, 7.1e-4, 13) == pytest.approx(0.00655, rel=0.01)


def test_liquid_release():
    """
    Liquid mass flow (kg/s) trough a hole or orifice.
    flow conditions. The formula is based on Yellow Book equation 2.194.
    Example from sec. 2.6.4.1

    Methods for the calculation of physical effects, CPR 14E, van den Bosch and Weterings (Eds.), 1996
    """
    d = 0.1
    A = math.pi * (d / 2) ** 2
    P1 = 1.013e5
    P2 = 1.013e5
    rho = 812.5
    Cd = 0.62
    H = 11.2

    assert 60.915 == pytest.approx(
        liquid_release_bernouilli(P1, P2, rho, Cd, A, H), abs=2.5
    )


def test_two_phase_relief():
    """
    Made up example to test the two phase relief function
    Fauske, H.K., 1982. Two-phase critical flow. In:
    Proceedings of the 10th International Conference on Nuclear Engineering,
    Arlington, VA, USA, April 19-23, 1982. American Society of Mechanical Engineers, New York, pp. 576-586.
    Example sec.

    """
    P1 = 50e5
    rho = 4.08
    k = 1.41
    Cd = 0.62
    d = 0.1
    A = math.pi * (d / 2) ** 2
    assert 20.86497695114019 == two_phase_release_fauske(P1, 0.55 * P1, rho, Cd, A)


def test_gas_release():
    """
    Gas massflow (kg/s) trough a hole at critical (sonic) or subcritical
    flow conditions. The formula is based on Yellow Book equation 2.22.

    Methods for the calculation of physical effects, CPR 14E,
    van den Bosch and Weterings (Eds.), 1996

    Example sec. 2.6.2.1 with hydrogen properties from CoolProp
    """
    P1 = 50e5
    P2 = 1.013e5
    rho = 4.08
    k = 1.41
    Cd = 0.62
    d = 0.1
    A = math.pi * (d / 2) ** 2

    assert 15.31 == pytest.approx(gas_release_rate(P1, P2, rho, k, Cd, A), abs=0.3)


def test_h_inside():
    """
    Example 4.7.1 from C. J. Geankoplis Transport Processes and Unit Operations, International Edition,
    Prentice-Hall, 1993
    """
    L = 0.3
    names = ["nitrogen", "oxygen"]
    molefracs = [0.79, 0.21]

    Tfluid = 311
    Tvessel = 505
    Tavg = (Tfluid + Tvessel) / 2
    P = 1.013e5
    flash = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=Tavg,
        rho="eos",
        model="PR",
    )

    gas = flash.flash(P=P, T=Tavg, zs=molefracs).gas
    print(flash.flash(P=P, T=Tavg, zs=molefracs).betas[2] == 0)
    print(flash.flash(P=P, T=Tavg, zs=molefracs).liquid0.zs)
    h_inner = h_inside(L, Tvessel, Tfluid, gas)

    assert h_inner == pytest.approx(7.03, rel=0.02)


def test_boiling_h():
    """
    Pseudo test comparing free convective heat transfer with nucleate boiling via
    the Rohsenow correaltion and that nucleate boiling exceeds convective heat
    transfer when the fluid-wall temperature difference approaches 5 C and above.
    """
    L = 0.3
    names = ["methane", "propane", "n-butane", "i-butane", "n-decane"]
    molefracs = [0.1, 0.05, 0.025, 0.025, 0.60]

    Tfluid = 300
    Tvessel = 304.85
    Tavg = (Tfluid + Tvessel) / 2
    P = 12.013e5
    flash = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=Tavg,
        rho="eos",
        model="PR",
    )
    res = flash.flash(P=P, T=Tavg, zs=molefracs)

    # something wrong with the thermal conductivity returned by thermo
    # Below is a simple mass weighted mixing rule for the individual
    # thermal conductivities.
    kl = sum([k * w for k, w in zip(res.liquid0.ws(), res.liquid0.kls())])

    h_boil = ht.Rohsenow(
        rhol=res.liquid0.rho_mass(),
        rhog=res.gas.rho_mass(),
        mul=res.liquid0.mu(),
        kl=kl,
        Cpl=res.liquid0.Cp_mass(),
        Hvap=(res.gas.H_mass() - res.liquid0.H_mass()),
        sigma=res.liquid0.sigma(),
        Te=(Tvessel - Tfluid),
        Csf=0.011,
        n=1.26,
    )
    h_conv = h_inside_liquid(L, Tvessel, Tfluid, res.liquid0)
    assert h_boil == pytest.approx(h_conv, rel=0.04)


def test_boiling_nucleic_Rohsenow():
    # Checked with 10.30 Problem set 8.
    # Copied from HT package tests
    h_calc = [
        ht.Rohsenow(
            Te=i,
            Cpl=4180,
            kl=0.688,
            mul=2.75e-4,
            sigma=0.0588,
            Hvap=2.25e6,
            rhol=958,
            rhog=0.597,
            Csf=0.013,
            n=1,
        )
        for i in [4.3, 9.1, 13]
    ]
    h_values = [2860.6242230238613, 12811.697777642301, 26146.321995188344]
    assert h_calc == h_values
    q_test = (
        ht.Rohsenow(
            Te=4.9,
            Cpl=4217.0,
            kl=0.680,
            mul=2.79e-4,
            sigma=0.0589,
            Hvap=2.257e6,
            rhol=957.854,
            rhog=0.595593,
            Csf=0.011,
            n=1.26,
        )
        * 4.9
    )
    assert 18245.91080863059 == pytest.approx(q_test, rel=0.01)

    h_Te = 1316.2269561541964
    h_q = ht.Rohsenow(
        q=5 * h_Te,
        Cpl=4180,
        kl=0.688,
        mul=2.75e-4,
        sigma=0.0588,
        Hvap=2.25e6,
        rhol=958,
        rhog=0.597,
    )
    assert h_Te == pytest.approx(h_q, rel=0.01)

    with pytest.raises(Exception):
        ht.Rohsenow(
            Cpl=4180,
            kl=0.688,
            mul=2.75e-4,
            sigma=0.0588,
            Hvap=2.25e6,
            rhol=958,
            rhog=0.597,
        )


def test_stefan_boltzmann():
    alpha = 1
    e_flame = 1
    e_surface = 0
    h = 100
    Tflame = 635 + 273.15
    Tradiative = 635 + 273.15
    assert (
        fire.stefan_boltzmann(
            alpha, e_flame, e_surface, h, Tflame, Tradiative, 20 + 273.15
        )
    ) == pytest.approx(1e5, abs=100)


def test_pool_fire_api521():
    assert fire.pool_fire_api521(273 + 50) == pytest.approx(45.5e3, abs=100)


def test_jet_fire_api521():
    assert fire.jet_fire_api521(273 + 50) == pytest.approx(83.5e3, abs=500)


def test_jet_fire_scandpower():
    assert fire.jet_fire_scandpower(273 + 20) == pytest.approx(94.5e3, abs=1000)


def test_pool_fire_scandpower():
    assert fire.pool_fire_scandpower(273 + 20) == pytest.approx(88.5e3, abs=500)


def test_peak_large_fire_scandpower():
    assert fire.jet_fire_peak_large_scandpower(273 + 25) == pytest.approx(
        313.6e3, abs=500
    )


def test_peak_small_fire_scandpower():
    assert fire.jet_fire_peak_small_scandpower(273 + 25) == pytest.approx(
        226.7e3, abs=500
    )


def test_peak_pool_scandpower():
    assert fire.pool_fire_peak_scandpower(273 + 25) == pytest.approx(131.2e3, abs=500)


def test_sb():
    assert fire.sb_fire(273 + 50, "api_jet") == pytest.approx(83.5e3, abs=500)
    assert fire.sb_fire(273 + 50, "api_pool") == pytest.approx(45.5e3, abs=100)
    assert fire.sb_fire(273 + 20, "scandpower_pool") == pytest.approx(88.5e3, abs=500)
    assert fire.sb_fire(273 + 20, "scandpower_jet") == pytest.approx(94.5e3, abs=1000)
    assert fire.sb_fire(273 + 25, "scandpower_jet_peak_large") == pytest.approx(
        313.6e3, abs=1000
    )
    assert fire.sb_fire(273 + 25, "scandpower_jet_peak_small") == pytest.approx(
        226.7e3, abs=1000
    )
    assert fire.sb_fire(273 + 25, "scandpower_pool_peak") == pytest.approx(
        131.2e3, abs=1000
    )

    try:
        Q = fire.sb_fire(273 + 20, "scand_jet") == pytest.approx(94.5e3, abs=1000)
    except ValueError:
        pass


if __name__ == "__main__":
    pass
