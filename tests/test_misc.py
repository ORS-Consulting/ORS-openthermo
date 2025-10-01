import pytest
import ht
from openthermo.properties.transport import h_inside, h_inside_liquid
from openthermo.flash.michelsen import get_flash_dry


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
    assert h_boil == pytest.approx(h_conv, rel=0.02)


if __name__ == "__main__":
    # names = ["water", "methane", "propane", "n-butane", "i-butane", "n-decane"]
    # molefracs = [0.1, 0.1, 0.05, 0.025, 0.025, 0.60]
    # Tfluid = 310
    # Tvessel = 300
    # Tavg = (Tfluid + Tvessel) / 2
    # P = 12.013e5
    # flash = get_flash(
    #     names,
    #     molefracs,
    #     P=P,
    #     T=Tavg,
    #     rho="eos",
    #     model="PR",
    # )
    # res = flash.flash(P=P, T=Tavg, zs=molefracs)

    # kls = 0
    # for i in range(len(res.liquid0.ws())):
    #     kls += res.liquid0.ws()[i] * res.liquid0.kls()[i]

    # print(kls, res.liquid0.k())

    test_h_inside()
