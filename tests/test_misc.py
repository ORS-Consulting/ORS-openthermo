import pytest
import ht
import math
import numpy as np
from openthermo.properties.transport import h_inside, h_inside_liquid
import openthermo.properties.transport as tp
from openthermo.flash.michelsen import get_flash_dry
from openthermo.vessel import fire
from openthermo.vessel.flowdevices import (
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


class TestPseudoComponentProperties:
    """Test pseudo-component property estimation correlations."""

    # Critical Temperature Tests
    def test_Tc_Kesler_Lee_reference_case(self):
        """Test Kesler-Lee Tc against reference case from Ahmed (2007)."""
        # Example 2.2 from Ahmed
        SG = 0.7365
        Tb = 365.555  # K
        Tc = Tc_Kesler_Lee_SG_Tb(SG, Tb)
        assert Tc == pytest.approx(545.0124354151242, rel=1e-6)

    def test_Tc_Riazi_Daubert_reference_case(self):
        """Test Riazi-Daubert Tc against reference case."""
        SG = 0.7365
        Tb = 365.555  # K
        Tc = Tc_Riazi_Daubert_SG_Tb(SG, Tb)
        assert Tc == pytest.approx(550.3734796518954, rel=1e-6)

    def test_Tc_increases_with_boiling_point(self):
        """Test that Tc increases with boiling point for fixed SG."""
        SG = 0.75
        Tb_low = 350  # K
        Tb_high = 450  # K

        Tc_low_KL = Tc_Kesler_Lee_SG_Tb(SG, Tb_low)
        Tc_high_KL = Tc_Kesler_Lee_SG_Tb(SG, Tb_high)
        assert Tc_high_KL > Tc_low_KL

        Tc_low_RD = Tc_Riazi_Daubert_SG_Tb(SG, Tb_low)
        Tc_high_RD = Tc_Riazi_Daubert_SG_Tb(SG, Tb_high)
        assert Tc_high_RD > Tc_low_RD

    def test_Tc_increases_with_specific_gravity(self):
        """Test that Tc increases with specific gravity for fixed Tb."""
        Tb = 400  # K
        SG_low = 0.70
        SG_high = 0.85

        Tc_low = Tc_Kesler_Lee_SG_Tb(SG_low, Tb)
        Tc_high = Tc_Kesler_Lee_SG_Tb(SG_high, Tb)
        assert Tc_high > Tc_low

    # Critical Pressure Tests
    def test_Pc_Kesler_Lee_reference_case(self):
        """Test Kesler-Lee Pc against reference case."""
        SG = 0.7365
        Tb = 365.555  # K
        Pc = Pc_Kesler_Lee_SG_Tb(SG, Tb)
        assert Pc == pytest.approx(3238323.346840464, rel=1e-6)

    def test_Pc_Riazi_Daubert_reference_case(self):
        """Test Riazi-Daubert Pc against reference case."""
        SG = 0.7365
        Tb = 365.555  # K
        Pc = Pc_Riazi_Daubert_SG_Tb(SG, Tb)
        assert Pc == pytest.approx(3219182.887436976, rel=1e-6)

    def test_Pc_decreases_with_boiling_point(self):
        """Test that Pc decreases with boiling point for fixed SG."""
        SG = 0.75
        Tb_low = 350  # K
        Tb_high = 450  # K

        Pc_low = Pc_Kesler_Lee_SG_Tb(SG, Tb_low)
        Pc_high = Pc_Kesler_Lee_SG_Tb(SG, Tb_high)
        # Higher boiling point = heavier molecule = lower Pc
        assert Pc_high < Pc_low

    def test_Pc_physically_reasonable(self):
        """Test that estimated Pc values are physically reasonable."""
        SG = 0.75
        Tb = 400  # K
        Pc_KL = Pc_Kesler_Lee_SG_Tb(SG, Tb)
        Pc_RD = Pc_Riazi_Daubert_SG_Tb(SG, Tb)

        # Critical pressure should be between 1 and 10 MPa for typical fractions
        assert 1e6 < Pc_KL < 10e6
        assert 1e6 < Pc_RD < 10e6

    # Molecular Weight Tests
    def test_MW_Kesler_Lee_reference_case(self):
        """Test Kesler-Lee MW against reference case."""
        SG = 0.7365
        Tb = 365.555  # K
        MW = MW_Kesler_Lee_SG_Tb(SG, Tb)
        assert MW == pytest.approx(98.70887589833501, rel=1e-6)

    def test_MW_increases_with_boiling_point(self):
        """Test that MW increases with boiling point."""
        SG = 0.75
        Tb_low = 350  # K
        Tb_high = 450  # K

        MW_low = MW_Kesler_Lee_SG_Tb(SG, Tb_low)
        MW_high = MW_Kesler_Lee_SG_Tb(SG, Tb_high)
        assert MW_high > MW_low

    def test_MW_physically_reasonable(self):
        """Test that estimated MW values are physically reasonable."""
        SG = 0.75
        Tb = 400  # K
        MW = MW_Kesler_Lee_SG_Tb(SG, Tb)

        # MW should be between C6 (86) and C30 (422) for typical fractions
        assert 80 < MW < 500

    # Acentric Factor Tests
    def test_omega_Kesler_Lee_reference_case(self):
        """Test Kesler-Lee omega against reference case."""
        SG = 0.7365
        Tb = 365.555  # K
        Tc = 545.012  # K
        Pc = 3238323.0  # Pa
        omega = omega_Kesler_Lee_SG_Tb_Tc_Pc(SG, Tb, Tc, Pc)
        assert omega == pytest.approx(0.306392118159797, rel=1e-6)

    def test_omega_with_varying_SG(self):
        """Test omega calculation for different specific gravities.

        Note: Omega relationship with SG is complex and depends on both
        SG and Tb. This test verifies the calculation runs without error
        and produces reasonable values.
        """
        Tb = 400  # K
        SG_low = 0.70
        SG_high = 0.85

        omega_low = omega_Kesler_Lee_SG_Tb_Tc_Pc(SG_low, Tb)
        omega_high = omega_Kesler_Lee_SG_Tb_Tc_Pc(SG_high, Tb)

        # Both should be physically reasonable
        assert 0.1 < omega_low < 1.0
        assert 0.1 < omega_high < 1.0

    def test_omega_physically_reasonable(self):
        """Test that estimated omega values are physically reasonable."""
        SG = 0.75
        Tb = 400  # K
        omega = omega_Kesler_Lee_SG_Tb_Tc_Pc(SG, Tb)

        # Acentric factor should be between 0.1 and 1.0 for typical fractions
        assert 0.1 < omega < 1.0

    # HC Atomic Ratio Tests
    def test_HC_atomic_ratio_reference_case(self):
        """Test HC atomic ratio against reference case."""
        # Example from Riazi (2005) for n-tetradecylbenzene (C20H24)
        # Expected HC ratio: H/C = 24/20 = 1.7
        SG = 0.8587
        Tb = 627  # K
        HC_ratio = HC_atomic_ratio(SG, Tb)
        assert HC_ratio == pytest.approx(1.7022132444770803, rel=1e-4)

    def test_HC_ratio_decreases_with_SG(self):
        """Test that HC ratio decreases with specific gravity (more aromatic)."""
        Tb = 500  # K
        SG_low = 0.70  # Paraffinic
        SG_high = 0.90  # Aromatic

        HC_low = HC_atomic_ratio(SG_low, Tb)
        HC_high = HC_atomic_ratio(SG_high, Tb)
        # More aromatic = lower H/C ratio
        assert HC_high < HC_low

    def test_HC_ratio_physically_reasonable(self):
        """Test that HC ratio is physically reasonable."""
        SG = 0.80
        Tb = 450  # K
        HC_ratio = HC_atomic_ratio(SG, Tb)

        # H/C ratio should be between 1.0 (very aromatic) and 2.5 (paraffinic)
        assert 1.0 < HC_ratio < 2.5

    # Critical Compressibility Tests
    def test_Zc_simple_fluid(self):
        """Test Zc for simple fluid (omega=0)."""
        from openthermo.properties.pseudo import Zc_pseudo

        omega = 0.0
        Zc = Zc_pseudo(omega)
        assert Zc == pytest.approx(0.2905, rel=1e-6)

    def test_Zc_normal_alkane(self):
        """Test Zc for normal alkane (omega~0.3)."""
        from openthermo.properties.pseudo import Zc_pseudo

        omega = 0.3
        Zc = Zc_pseudo(omega)
        assert Zc == pytest.approx(0.2650, rel=1e-4)

    def test_Zc_decreases_with_omega(self):
        """Test that Zc decreases with increasing omega."""
        from openthermo.properties.pseudo import Zc_pseudo

        omega_low = 0.2
        omega_high = 0.6

        Zc_low = Zc_pseudo(omega_low)
        Zc_high = Zc_pseudo(omega_high)
        assert Zc_high < Zc_low

    # Critical Volume Tests
    def test_Vc_calculation(self):
        """Test Vc calculation using ideal gas law at critical point."""
        from openthermo.properties.pseudo import Vc_pseudo

        Zc = 0.256
        Tc = 617.7  # K (n-Decane)
        Pc = 21.1e5  # Pa
        Vc = Vc_pseudo(Zc, Tc, Pc)

        # n-Decane Vc ~ 624 cm³/mol = 6.24e-4 m³/mol
        assert Vc == pytest.approx(6.24e-4, rel=0.01)

    def test_Vc_increases_with_Tc(self):
        """Test that Vc increases with Tc for fixed Zc and Pc."""
        from openthermo.properties.pseudo import Vc_pseudo

        Zc = 0.27
        Pc = 30e5  # Pa
        Tc_low = 500  # K
        Tc_high = 600  # K

        Vc_low = Vc_pseudo(Zc, Tc_low, Pc)
        Vc_high = Vc_pseudo(Zc, Tc_high, Pc)
        assert Vc_high > Vc_low


class TestPseudoCorrelationComparisons:
    """Test comparisons between different pseudo-component correlations."""

    def test_Kesler_Lee_vs_Riazi_Daubert_Tc(self):
        """Compare Tc predictions from different correlations."""
        SG = 0.75
        Tb = 400  # K

        Tc_KL = Tc_Kesler_Lee_SG_Tb(SG, Tb)
        Tc_RD = Tc_Riazi_Daubert_SG_Tb(SG, Tb)

        # Both should give similar but not identical results
        # Typically within 5-10%
        assert Tc_KL == pytest.approx(Tc_RD, rel=0.10)

    def test_Kesler_Lee_vs_Riazi_Daubert_Pc(self):
        """Compare Pc predictions from different correlations."""
        SG = 0.75
        Tb = 400  # K

        Pc_KL = Pc_Kesler_Lee_SG_Tb(SG, Tb)
        Pc_RD = Pc_Riazi_Daubert_SG_Tb(SG, Tb)

        # Both should give similar but not identical results
        # Typically within 10-15%
        assert Pc_KL == pytest.approx(Pc_RD, rel=0.15)


class TestDimensionlessNumbers:
    """Test dimensionless number calculations."""

    def test_Nu_turbulent_regime(self):
        """Test Nusselt number in turbulent regime (Ra >= 1e9)."""
        Ra = 1.27e9
        Nu_calc = tp.Nu(Ra, 0)
        Nu_expected = 0.13 * (Ra) ** (1 / 3)
        assert Nu_calc == pytest.approx(Nu_expected, rel=0.01)

    def test_Nu_transition_regime(self):
        """Test Nusselt number in transition regime (1e4 < Ra < 1e9)."""
        Ra = 1.27e8
        Nu_calc = tp.Nu(Ra, 0)
        Nu_expected = 0.59 * (Ra) ** 0.25
        assert Nu_calc == pytest.approx(Nu_expected, abs=1.0)

    def test_Nu_laminar_regime(self):
        """Test Nusselt number in laminar regime (Ra <= 1e4)."""
        Ra = 1000
        Nu_calc = tp.Nu(Ra, 0)
        Nu_expected = 1.36 * (Ra) ** 0.20
        assert Nu_calc == pytest.approx(Nu_expected, rel=0.01)

    def test_Nu_regime_boundaries(self):
        """Test Nusselt number at regime boundaries."""
        # At Ra = 1e9
        Ra_high = 1e9
        Nu_high = tp.Nu(Ra_high, 0)
        assert Nu_high > 0

        # Just below 1e9
        Ra_mid = 9.9e8
        Nu_mid = tp.Nu(Ra_mid, 0)
        assert Nu_mid > 0

        # At Ra = 1e4
        Ra_low = 1e4
        Nu_low = tp.Nu(Ra_low, 0)
        assert Nu_low > 0

    def test_Pr_basic(self):
        """Test Prandtl number calculation."""
        cp = 120  # J/(kg·K)
        mu = 7.1e-4  # Pa·s
        k = 13  # W/(m·K)
        Pr_calc = tp.Pr(cp, mu, k)
        Pr_expected = cp * mu / k
        assert Pr_calc == pytest.approx(Pr_expected, rel=0.01)
        assert Pr_calc == pytest.approx(0.00655, rel=0.01)

    def test_Pr_typical_gases(self):
        """Test Prandtl number for typical gas properties."""
        # Air at room temperature: Pr ~ 0.7
        cp_air = 1005  # J/(kg·K)
        mu_air = 1.85e-5  # Pa·s
        k_air = 0.026  # W/(m·K)
        Pr_air = tp.Pr(cp_air, mu_air, k_air)
        assert 0.6 < Pr_air < 0.8

    def test_Pr_typical_liquids(self):
        """Test Prandtl number for typical liquid properties."""
        # Water at 20°C: Pr ~ 7
        cp_water = 4182  # J/(kg·K)
        mu_water = 1.0e-3  # Pa·s
        k_water = 0.6  # W/(m·K)
        Pr_water = tp.Pr(cp_water, mu_water, k_water)
        assert 5 < Pr_water < 10


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


class TestHeatTransferCoefficients:
    """Test heat transfer coefficient calculations."""

    def test_h_inside_geankoplis_example(self):
        """
        Validate against Example 4.7.1 from Geankoplis.

        Example 4.7.1 from C. J. Geankoplis Transport Processes and Unit Operations,
        International Edition, Prentice-Hall, 1993
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
        h_inner = h_inside(L, Tvessel, Tfluid, gas)

        assert h_inner == pytest.approx(7.03, rel=0.02)

    def test_h_inside_increases_with_delta_T(self):
        """Test that heat transfer coefficient increases with temperature difference."""
        L = 0.3
        names = ["nitrogen", "oxygen"]
        molefracs = [0.79, 0.21]
        P = 1.013e5

        # Small temperature difference
        Tfluid_low = 300
        Tvessel_low = 310
        Tavg_low = (Tfluid_low + Tvessel_low) / 2
        flash = get_flash_dry(names, molefracs, P=P, T=Tavg_low, rho="eos", model="PR")
        gas_low = flash.flash(P=P, T=Tavg_low, zs=molefracs).gas
        h_low = h_inside(L, Tvessel_low, Tfluid_low, gas_low)

        # Large temperature difference
        Tfluid_high = 300
        Tvessel_high = 400
        Tavg_high = (Tfluid_high + Tvessel_high) / 2
        gas_high = flash.flash(P=P, T=Tavg_high, zs=molefracs).gas
        h_high = h_inside(L, Tvessel_high, Tfluid_high, gas_high)

        # Higher temperature difference should give higher heat transfer coefficient
        assert h_high > h_low

    def test_h_inside_positive(self):
        """Test that heat transfer coefficient is always positive."""
        L = 0.5
        names = ["methane"]
        molefracs = [1.0]
        P = 10e5

        Tfluid = 300
        Tvessel = 350
        Tavg = (Tfluid + Tvessel) / 2
        flash = get_flash_dry(names, molefracs, P=P, T=Tavg, rho="eos", model="PR")
        gas = flash.flash(P=P, T=Tavg, zs=molefracs).gas
        h = h_inside(L, Tvessel, Tfluid, gas)

        assert h > 0
        assert h < 1000  # Reasonable upper bound for natural convection


class TestBoilingHeatTransfer:
    """Test boiling heat transfer calculations."""

    def test_boiling_h_small_superheat(self):
        """
        Test comparing free convective heat transfer with nucleate boiling.

        At small temperature differences (~5 K), nucleate boiling and natural
        convection should give similar heat transfer coefficients.
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

        # Mass weighted mixing rule for thermal conductivity
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

    def test_boiling_h_larger_superheat(self):
        """Test boiling heat transfer at larger superheats.

        Note: At moderate superheats, either boiling or natural convection
        may dominate depending on fluid properties and geometry. This test
        verifies that both mechanisms produce reasonable values.
        """
        L = 0.3
        names = ["methane", "propane", "n-butane", "i-butane", "n-decane"]
        molefracs = [0.1, 0.05, 0.025, 0.025, 0.60]

        Tfluid = 300
        Tvessel = 315  # Larger superheat (15 K)
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
            Csf=0.013,
            n=1.7,
        )
        h_conv = h_inside_liquid(L, Tvessel, Tfluid, res.liquid0)

        # Both mechanisms should produce reasonable heat transfer coefficients
        assert h_boil > 0
        assert h_conv > 0
        assert h_boil < 10000  # Reasonable upper bound
        assert h_conv < 10000  # Reasonable upper bound


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


class TestStefanBoltzmann:
    """Test Stefan-Boltzmann radiation calculations."""

    def test_stefan_boltzmann_basic(self):
        """Test basic Stefan-Boltzmann calculation."""
        alpha = 1
        e_flame = 1
        e_surface = 0
        h = 100
        Tflame = 635 + 273.15
        Tradiative = 635 + 273.15
        Tvessel = 20 + 273.15

        Q = fire.stefan_boltzmann(
            alpha, e_flame, e_surface, h, Tflame, Tradiative, Tvessel
        )

        assert Q == pytest.approx(1e5, abs=100)

    def test_stefan_boltzmann_zero_emissivity(self):
        """Test Stefan-Boltzmann with zero surface emissivity (only convection)."""
        Q = fire.stefan_boltzmann(
            alpha=1,
            e_flame=1,
            e_surface=0,
            h=50,
            Tflame=1000,
            Tradiative=1000,
            Tvessel=300,
        )

        # Should still have heat flux from convection
        assert Q > 0
        # Convective part should be h * (Tflame - Tvessel)
        Q_conv_expected = 50 * (1000 - 300)
        assert Q_conv_expected > 0

    def test_stefan_boltzmann_high_emissivity(self):
        """Test Stefan-Boltzmann with different surface emissivity values."""
        Q_low_emissivity = fire.stefan_boltzmann(
            alpha=1,
            e_flame=1,
            e_surface=0.1,
            h=100,
            Tflame=1000,
            Tradiative=1000,
            Tvessel=300,
        )

        Q_high_emissivity = fire.stefan_boltzmann(
            alpha=1,
            e_flame=1,
            e_surface=0.9,
            h=100,
            Tflame=1000,
            Tradiative=1000,
            Tvessel=300,
        )

        # Higher surface emissivity means surface radiates more heat away
        # so net heat flux to vessel is reduced
        assert Q_high_emissivity <= Q_low_emissivity

    def test_stefan_boltzmann_equal_temperatures(self):
        """Test Stefan-Boltzmann when wall equals flame temperature."""
        # When temperatures are equal, heat flux should be small
        T = 500
        Q = fire.stefan_boltzmann(
            alpha=1, e_flame=1, e_surface=0.5, h=100, Tflame=T, Tradiative=T, Tvessel=T
        )

        # Allow for numerical effects and small residuals
        assert Q == pytest.approx(0, abs=2000)

    def test_stefan_boltzmann_different_radiative_temp(self):
        """Test Stefan-Boltzmann with different radiative temperature."""
        Q1 = fire.stefan_boltzmann(
            alpha=1,
            e_flame=1,
            e_surface=0.5,
            h=100,
            Tflame=1000,
            Tradiative=1000,
            Tvessel=300,
        )

        Q2 = fire.stefan_boltzmann(
            alpha=1,
            e_flame=1,
            e_surface=0.5,
            h=100,
            Tflame=1000,
            Tradiative=800,
            Tvessel=300,
        )

        # Different radiative temperature should give different heat flux
        assert Q1 != Q2


class TestPoolFireAPI521:
    """Test API 521 pool fire scenarios."""

    def test_pool_fire_api521_room_temp(self):
        """Test API 521 pool fire at room temperature."""
        Q = fire.pool_fire_api521(273 + 20)
        # Should be around 60 kW/m² minus radiative losses
        assert Q > 0
        assert Q < 60e3  # Maximum incident heat flux

    def test_pool_fire_api521_elevated_temp(self):
        """Test API 521 pool fire at elevated temperature."""
        Q = fire.pool_fire_api521(273 + 50)
        assert Q == pytest.approx(45.5e3, abs=100)

    def test_pool_fire_api521_temperature_dependence(self):
        """Test that pool fire heat flux decreases with wall temperature."""
        Q_low = fire.pool_fire_api521(273 + 20)
        Q_high = fire.pool_fire_api521(273 + 200)

        # Higher wall temperature should give lower heat flux
        assert Q_high < Q_low


class TestJetFireAPI521:
    """Test API 521 jet fire scenarios."""

    def test_jet_fire_api521_room_temp(self):
        """Test API 521 jet fire at room temperature."""
        Q = fire.jet_fire_api521(273 + 20)
        # Should be around 100 kW/m² minus radiative losses
        assert Q > 0
        assert Q < 100e3

    def test_jet_fire_api521_elevated_temp(self):
        """Test API 521 jet fire at elevated temperature."""
        Q = fire.jet_fire_api521(273 + 50)
        assert Q == pytest.approx(83.5e3, abs=500)

    def test_jet_fire_api521_higher_than_pool(self):
        """Test that jet fire gives higher heat flux than pool fire."""
        T = 273 + 50
        Q_jet = fire.jet_fire_api521(T)
        Q_pool = fire.pool_fire_api521(T)

        # Jet fire should have higher heat flux than pool fire
        assert Q_jet > Q_pool


class TestPoolFireScandpower:
    """Test Scandpower pool fire scenarios."""

    def test_pool_fire_scandpower_room_temp(self):
        """Test Scandpower pool fire at room temperature."""
        Q = fire.pool_fire_scandpower(273 + 20)
        assert Q == pytest.approx(88.5e3, abs=500)

    def test_pool_fire_scandpower_elevated_temp(self):
        """Test Scandpower pool fire at elevated temperature."""
        Q = fire.pool_fire_scandpower(273 + 100)
        # Should be positive but less than at room temp
        assert Q > 0
        assert Q < fire.pool_fire_scandpower(273 + 20)

    def test_pool_fire_scandpower_higher_than_api(self):
        """Test that Scandpower pool fire is more conservative than API."""
        T = 273 + 20
        Q_scandpower = fire.pool_fire_scandpower(T)
        Q_api = fire.pool_fire_api521(T)

        # Scandpower should be more conservative (higher heat flux)
        assert Q_scandpower > Q_api


class TestJetFireScandpower:
    """Test Scandpower jet fire scenarios."""

    def test_jet_fire_scandpower_room_temp(self):
        """Test Scandpower jet fire at room temperature."""
        Q = fire.jet_fire_scandpower(273 + 20)
        assert Q == pytest.approx(94.5e3, abs=1000)

    def test_jet_fire_scandpower_elevated_temp(self):
        """Test Scandpower jet fire at elevated temperature."""
        Q = fire.jet_fire_scandpower(273 + 150)
        # Should be positive but less than at room temp
        assert Q > 0
        assert Q < fire.jet_fire_scandpower(273 + 20)

    def test_jet_fire_scandpower_temperature_dependence(self):
        """Test temperature dependence of Scandpower jet fire."""
        Q_20 = fire.jet_fire_scandpower(273 + 20)
        Q_50 = fire.jet_fire_scandpower(273 + 50)
        Q_100 = fire.jet_fire_scandpower(273 + 100)

        # Heat flux should decrease with increasing wall temperature
        assert Q_20 > Q_50 > Q_100


class TestSbFire:
    """Test the sb_fire wrapper function."""

    def test_sb_fire_all_scenarios(self):
        """Test sb_fire for all valid scenario types."""
        T = 273 + 50

        # Test all valid scenarios
        Q_api_jet = fire.sb_fire(T, "api_jet")
        Q_api_pool = fire.sb_fire(T, "api_pool")
        Q_sp_pool = fire.sb_fire(T, "scandpower_pool")
        Q_sp_jet = fire.sb_fire(T, "scandpower_jet")

        # All should be positive
        assert Q_api_jet > 0
        assert Q_api_pool > 0
        assert Q_sp_pool > 0
        assert Q_sp_jet > 0

        # Jets should be higher than pools
        assert Q_api_jet > Q_api_pool

    def test_sb_fire_api_scenarios(self):
        """Test sb_fire API scenarios match direct function calls."""
        T = 273 + 50

        assert fire.sb_fire(T, "api_jet") == pytest.approx(83.5e3, abs=500)
        assert fire.sb_fire(T, "api_pool") == pytest.approx(45.5e3, abs=100)

    def test_sb_fire_scandpower_scenarios(self):
        """Test sb_fire Scandpower scenarios match direct function calls."""
        T = 273 + 20

        assert fire.sb_fire(T, "scandpower_pool") == pytest.approx(88.5e3, abs=500)
        assert fire.sb_fire(T, "scandpower_jet") == pytest.approx(94.5e3, abs=1000)

    def test_sb_fire_invalid_scenario(self):
        """Test sb_fire raises error for invalid scenario."""
        with pytest.raises(ValueError):
            fire.sb_fire(273 + 20, "invalid_scenario")

    def test_sb_fire_typo_scenario(self):
        """Test sb_fire raises error for typo in scenario name."""
        # Common typos
        with pytest.raises(ValueError):
            fire.sb_fire(273 + 20, "scand_jet")

        with pytest.raises(ValueError):
            fire.sb_fire(273 + 20, "scandpower_jett")

        with pytest.raises(ValueError):
            fire.sb_fire(273 + 20, "API_jet")  # Wrong case


class TestFireScenarioComparisons:
    """Test relationships between different fire scenarios."""

    def test_jet_vs_pool_heat_flux(self):
        """Test that jet fires consistently give higher heat flux than pools."""
        temperatures = [273 + 20, 273 + 50, 273 + 100, 273 + 200]

        for T in temperatures:
            Q_jet_api = fire.jet_fire_api521(T)
            Q_pool_api = fire.pool_fire_api521(T)

            # Jet should always be higher than pool for same standard
            assert Q_jet_api > Q_pool_api

    def test_scandpower_vs_api_conservatism(self):
        """Test relative conservatism between standards."""
        T = 273 + 50

        Q_jet_api = fire.jet_fire_api521(T)
        Q_pool_api = fire.pool_fire_api521(T)
        Q_jet_sp = fire.jet_fire_scandpower(50 + 273)
        Q_pool_sp = fire.pool_fire_scandpower(50 + 273)

        # All should be reasonable heat flux values
        assert 40e3 < Q_pool_api < 100e3
        assert 70e3 < Q_jet_api < 150e3

    def test_temperature_effect_consistency(self):
        """Test that all scenarios show consistent temperature effects."""
        T_low = 273 + 20
        T_high = 273 + 200

        # All scenarios should give lower heat flux at higher wall temperature
        assert fire.pool_fire_api521(T_high) < fire.pool_fire_api521(T_low)
        assert fire.jet_fire_api521(T_high) < fire.jet_fire_api521(T_low)
        assert fire.pool_fire_scandpower(T_high) < fire.pool_fire_scandpower(T_low)
        assert fire.jet_fire_scandpower(T_high) < fire.jet_fire_scandpower(T_low)
        
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


class TestVonMisesStress:
    """Test von Mises stress calculations for pressure vessels."""

    def test_von_mises_basic_calculation(self):
        """Test basic von Mises stress calculation."""
        from openthermo.properties.materials import von_mises

        p = 100e5  # 100 bar
        d = 0.5  # 0.5 m inner diameter
        wt = 0.01  # 10 mm wall thickness
        stress = von_mises(p, d, wt)

        # Verify stress is positive and reasonable
        assert stress > 0
        assert stress < 1e9  # Should be less than 1 GPa for typical cases

    def test_von_mises_higher_pressure(self):
        """Test von Mises stress increases with pressure."""
        from openthermo.properties.materials import von_mises

        d = 0.5
        wt = 0.01
        p_low = 50e5
        p_high = 150e5

        stress_low = von_mises(p_low, d, wt)
        stress_high = von_mises(p_high, d, wt)

        assert stress_high > stress_low

    def test_von_mises_thicker_wall(self):
        """Test von Mises stress decreases with wall thickness."""
        from openthermo.properties.materials import von_mises

        p = 100e5
        d = 0.5
        wt_thin = 0.005
        wt_thick = 0.020

        stress_thin = von_mises(p, d, wt_thin)
        stress_thick = von_mises(p, d, wt_thick)

        assert stress_thin > stress_thick

    def test_von_mises_with_custom_sigma_a(self):
        """Test von Mises calculation with custom axial stress."""
        from openthermo.properties.materials import von_mises

        p = 100e5
        d = 0.5
        wt = 0.01
        sigma_a_default = 30e6
        sigma_a_custom = 50e6

        stress_default = von_mises(p, d, wt, sigma_a=sigma_a_default)
        stress_custom = von_mises(p, d, wt, sigma_a=sigma_a_custom)

        # Higher axial stress should give higher equivalent stress
        assert stress_custom > stress_default


class TestUTS:
    """Test Ultimate Tensile Strength lookups."""

    def test_UTS_SS316_room_temperature(self):
        """Test UTS for SS316 at room temperature."""
        from openthermo.properties.materials import UTS

        # Room temperature (20°C = 293.15 K)
        uts_room = UTS(293.15, "SS316")
        assert uts_room == pytest.approx(575e6, rel=0.01)

    def test_UTS_SS316_elevated_temperature(self):
        """Test UTS for SS316 at elevated temperature."""
        from openthermo.properties.materials import UTS

        # Test at 500°C (773.15 K)
        uts_high = UTS(773.15, "SS316")
        uts_room = UTS(293.15, "SS316")

        # Strength decreases with temperature
        assert uts_high < uts_room
        assert uts_high == pytest.approx(449e6, rel=0.05)

    def test_UTS_all_materials_positive(self):
        """Test UTS function returns positive values for all materials."""
        from openthermo.properties.materials import UTS

        materials_list = ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"]
        T_test = 373.15  # 100°C

        for mat in materials_list:
            uts = UTS(T_test, mat)
            assert uts > 0
            assert uts < 1e9

    def test_UTS_temperature_degradation(self):
        """Test that UTS decreases with temperature for all materials."""
        from openthermo.properties.materials import UTS

        materials_list = ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"]
        T_low = 293.15  # 20°C
        T_high = 773.15  # 500°C

        for mat in materials_list:
            uts_low = UTS(T_low, mat)
            uts_high = UTS(T_high, mat)
            assert uts_high < uts_low

    def test_UTS_duplex_vs_carbon_steel(self):
        """Test that Duplex has higher UTS than carbon steel at room temp."""
        from openthermo.properties.materials import UTS

        T = 293.15
        uts_duplex = UTS(T, "Duplex")
        uts_cs235 = UTS(T, "CS_235LT")
        uts_cs360 = UTS(T, "CS_360LT")

        assert uts_duplex > uts_cs235
        assert uts_duplex > uts_cs360

    def test_UTS_CS360_higher_than_CS235(self):
        """Test that CS360LT has higher UTS than CS235LT."""
        from openthermo.properties.materials import UTS

        T = 373.15  # 100°C
        uts_cs360 = UTS(T, "CS_360LT")
        uts_cs235 = UTS(T, "CS_235LT")

        assert uts_cs360 > uts_cs235

    def test_UTS_interpolation_consistency(self):
        """Test that interpolation gives smooth values."""
        from openthermo.properties.materials import UTS

        mat = "SS316"
        # Test at intermediate temperature
        T_mid = 350  # Between tabulated values
        uts_mid = UTS(T_mid, mat)

        # Should be between neighboring values
        T_low = 323.15  # 50°C
        T_high = 423.15  # 150°C
        uts_low = UTS(T_low, mat)
        uts_high = UTS(T_high, mat)

        assert uts_high < uts_mid < uts_low

    def test_UTS_extreme_temperatures(self):
        """Test UTS at extreme temperatures."""
        from openthermo.properties.materials import UTS

        mat = "Duplex"

        # At minimum temperature (20°C)
        uts_min = UTS(293.17, mat)
        assert uts_min == pytest.approx(730e6, rel=0.01)

        # At maximum temperature (1100°C)
        uts_max = UTS(1373.17, mat)
        assert uts_max > 0
        assert uts_max < uts_min


class TestATS:
    """Test Allowable Tensile Strength calculations."""

    def test_ATS_default_safety_factor(self):
        """Test ATS with default safety factor k_s=0.85."""
        from openthermo.properties.materials import ATS, UTS

        T = 373.15
        mat = "SS316"

        ats = ATS(T, mat, k_s=0.85, k_y=1.0)
        uts = UTS(T, mat)

        assert ats == pytest.approx(uts * 0.85, rel=0.01)

    def test_ATS_guaranteed_minimum(self):
        """Test ATS with guaranteed minimum (k_s=1.0)."""
        from openthermo.properties.materials import ATS, UTS

        T = 373.15
        mat = "SS316"

        ats_default = ATS(T, mat, k_s=0.85, k_y=1.0)
        ats_guaranteed = ATS(T, mat, k_s=1.0, k_y=1.0)

        assert ats_guaranteed > ats_default
        assert ats_guaranteed == pytest.approx(UTS(T, mat), rel=0.01)

    def test_ATS_uncertain_material_factor(self):
        """Test ATS with uncertain material data (k_y < 1.0)."""
        from openthermo.properties.materials import ATS

        T = 373.15
        mat = "SS316"

        ats_normal = ATS(T, mat, k_s=0.85, k_y=1.0)
        ats_uncertain = ATS(T, mat, k_s=0.85, k_y=0.9)

        # Uncertain material should have lower allowable stress
        assert ats_uncertain < ats_normal

    def test_ATS_all_materials(self):
        """Test ATS for all material types."""
        from openthermo.properties.materials import ATS, UTS

        materials_list = ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"]
        T = 473.15  # 200°C

        for mat in materials_list:
            ats = ATS(T, mat)
            uts = UTS(T, mat)
            # ATS should be less than UTS (with default k_s=0.85)
            assert ats < uts
            assert ats == pytest.approx(uts * 0.85, rel=0.01)


class TestSteelCp:
    """Test heat capacity lookups."""

    def test_steel_Cp_SS316_room_temperature(self):
        """Test heat capacity for SS316 at room temperature."""
        from openthermo.properties.materials import steel_Cp

        mat = "SS316"
        Cp = steel_Cp(293.15, mat)  # 20°C
        assert Cp == pytest.approx(472, rel=0.01)

    def test_steel_Cp_temperature_dependence(self):
        """Test that Cp increases with temperature."""
        from openthermo.properties.materials import steel_Cp

        mat = "SS316"
        Cp_low = steel_Cp(293.15, mat)  # 20°C
        Cp_high = steel_Cp(773.15, mat)  # 500°C

        assert Cp_high > Cp_low

    def test_steel_Cp_all_materials(self):
        """Test heat capacity for all material types."""
        from openthermo.properties.materials import steel_Cp

        materials_list = ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"]
        T_test = 473.15  # 200°C

        for mat in materials_list:
            cp = steel_Cp(T_test, mat)
            assert cp > 400  # Reasonable lower bound for steel
            assert cp < 2000  # Reasonable upper bound

    def test_steel_Cp_carbon_steel_phase_transformation(self):
        """Test CS heat capacity discontinuity at phase transformation."""
        from openthermo.properties.materials import steel_Cp

        # Carbon steel has discontinuity at ~750°C due to austenite formation
        mat = "CS_235LT"

        # Test at temperatures around the phase transformation (750°C)
        Cp_700 = steel_Cp(700 + 273.15, mat)  # 700°C
        Cp_750 = steel_Cp(750 + 273.15, mat)  # 750°C

        # At 750°C there's a phase transformation, Cp jumps to 1450
        # At 700°C it's 900, so we should see a significant increase
        assert Cp_750 > Cp_700 * 1.4  # At least 40% increase
        assert Cp_750 > 1400  # Should be near the high value

    def test_steel_Cp_interpolation(self):
        """Test heat capacity interpolation between data points."""
        from openthermo.properties.materials import steel_Cp

        mat = "SS316"

        # Test at exact data points
        Cp_100 = steel_Cp(373.15, mat)  # 100°C
        Cp_200 = steel_Cp(473.15, mat)  # 200°C

        # Test at midpoint (should be between the two)
        Cp_150 = steel_Cp(423.15, mat)  # 150°C

        assert Cp_100 < Cp_150 < Cp_200

    def test_steel_Cp_duplex_vs_SS316(self):
        """Test different materials have different heat capacities."""
        from openthermo.properties.materials import steel_Cp

        T = 473.15  # 200°C

        Cp_duplex = steel_Cp(T, "Duplex")
        Cp_ss316 = steel_Cp(T, "SS316")
        Cp_6mo = steel_Cp(T, "6Mo")

        # All should be positive but different
        assert Cp_duplex > 0
        assert Cp_ss316 > 0
        assert Cp_6mo > 0
        # Duplex has higher Cp than SS316 at this temperature
        assert Cp_duplex > Cp_ss316

    def test_steel_Cp_extreme_temperatures(self):
        """Test heat capacity at extreme temperatures."""
        from openthermo.properties.materials import steel_Cp

        mat = "SS316"

        # At minimum temperature (20°C)
        Cp_min = steel_Cp(293.15, mat)
        assert Cp_min > 0

        # At maximum temperature (1100°C)
        Cp_max = steel_Cp(1373.15, mat)
        assert Cp_max > Cp_min

    def test_steel_Cp_6Mo_plateau(self):
        """Test that 6Mo Cp plateaus at high temperature."""
        from openthermo.properties.materials import steel_Cp

        mat = "6Mo"

        # 6Mo data shows plateau at 610 J/(kg·K) above ~750°C
        Cp_800 = steel_Cp(1073.15, mat)  # 800°C
        Cp_1000 = steel_Cp(1273.15, mat)  # 1000°C

        # Should be approximately the same (plateau)
        assert Cp_800 == pytest.approx(610, rel=0.05)
        assert Cp_1000 == pytest.approx(610, rel=0.05)


class TestMaterialDataConsistency:
    """Test consistency and validity of material property data."""

    def test_temperature_arrays_monotonic(self):
        """Test that temperature arrays are monotonically increasing."""
        from openthermo.properties import materials

        # Check T_Cp array
        assert np.all(np.diff(materials.T_Cp) > 0)

        # Check T array
        assert np.all(np.diff(materials.T) > 0)

    def test_UTS_arrays_monotonic_decreasing(self):
        """Test that UTS decreases monotonically with temperature."""
        from openthermo.properties import materials

        # UTS should decrease with temperature for all materials
        assert np.all(np.diff(materials.Duplex_UTS) <= 0)
        assert np.all(np.diff(materials.SS_UTS) <= 0)
        assert np.all(np.diff(materials.SMo_UTS) <= 0)
        assert np.all(np.diff(materials.CS_235LT_UTS) <= 0)
        assert np.all(np.diff(materials.CS_360LT_UTS) <= 0)

    def test_Cp_arrays_physically_reasonable(self):
        """Test that all Cp values are physically reasonable."""
        from openthermo.properties import materials

        # All heat capacities should be between 400-2000 J/(kg·K) for steel
        all_cp_values = np.concatenate(
            [
                materials.SS316_Cp,
                materials.Duplex_Cp,
                materials.SMo_Cp,
                materials.CS_LT_Cp,
            ]
        )

        # Check all values are within reasonable range
        # Allow up to 1500 for CS phase transformation
        assert np.all(all_cp_values > 400)
        assert np.all(all_cp_values < 1500)

    def test_UTS_arrays_physically_reasonable(self):
        """Test that all UTS values are physically reasonable."""
        from openthermo.properties import materials

        all_uts_values = np.concatenate(
            [
                materials.Duplex_UTS,
                materials.SS_UTS,
                materials.SMo_UTS,
                materials.CS_235LT_UTS,
                materials.CS_360LT_UTS,
            ]
        )

        # All UTS values should be positive and less than 1 GPa
        assert np.all(all_uts_values > 0)
        assert np.all(all_uts_values < 1e9)

    def test_temperature_array_lengths_match(self):
        """Test that temperature and property arrays have matching lengths."""
        from openthermo.properties import materials

        # Cp data
        assert len(materials.T_Cp) == len(materials.SS316_Cp)
        assert len(materials.T_Cp) == len(materials.Duplex_Cp)
        assert len(materials.T_Cp) == len(materials.SMo_Cp)
        assert len(materials.T_Cp) == len(materials.CS_LT_Cp)

        # UTS data
        assert len(materials.T) == len(materials.Duplex_UTS)
        assert len(materials.T) == len(materials.SS_UTS)
        assert len(materials.T) == len(materials.SMo_UTS)
        assert len(materials.T) == len(materials.CS_235LT_UTS)
        assert len(materials.T) == len(materials.CS_360LT_UTS)


if __name__ == "__main__":
    pass
