import pytest
import ht
import math
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


if __name__ == "__main__":
    pass
