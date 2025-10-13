from openthermo.flash.michelsen import get_flash_dry
import pytest


def test_flash_dry():
    P = 5.0e6
    T = 300.0
    names = ["methane", "ethane", "propane", "n-butane"]
    molefracs = [0.64, 0.06, 0.28, 0.02]
    flash = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=T,
        rho="eos",
        model="PR",
    )
    res = flash.flash(P=P, T=T, zs=molefracs)

    # running sanity check on various flashes
    res2 = flash.flash(P=P * 0.1, H=res.H(), zs=molefracs)
    assert res.H() == pytest.approx(res2.H(), rel=1e-5)
    res2 = flash.flash(P=P * 1.1, S=res.S(), zs=molefracs)
    assert res.S() == pytest.approx(res2.S(), rel=1e-5)
    res2 = flash.flash(
        U=res.U(),
        V=res.V(),
        zs=molefracs,
        Pguess=res.P * 1.1,
        Tguess=res.T * 1.1,
    )
    assert res.U() == pytest.approx(res2.U(), rel=1e-5)
    assert res.V() == pytest.approx(res2.V(), rel=1e-5)
    res2 = flash.flash(
        P=res.P * 1.1,
        V=res.V(),
        zs=molefracs,
    )
    assert res.V() == pytest.approx(res2.V(), rel=1e-5)
    res2 = flash.flash(
        T=res.T * 1.1,
        V=res.V(),
        zs=molefracs,
    )
    assert res.V() == pytest.approx(res2.V(), rel=1e-5)
    res2 = flash.flash(
        P=res.P * 1.1,
        U=res.U(),
        zs=molefracs,
        Tguess=res.T * 1.1,
    )
    assert res.U() == pytest.approx(res2.U(), rel=1e-5)
    res2 = flash.flash(
        T=res.T * 1.1,
        U=res.U(),
        zs=molefracs,
        Pguess=res.P * 1.1,
    )
    assert res.U() == pytest.approx(res2.U(), rel=1e-5)


def test_flash_SRK():
    P = 5.0e6
    T = 300.0
    names = ["methane", "ethane", "propane", "n-butane"]
    molefracs = [0.64, 0.06, 0.28, 0.02]
    flash = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=T,
        rho="eos",
        model="SRK",
    )
    res = flash.flash(P=P, T=T, zs=molefracs)

    # running sanity check on various flashes
    res2 = flash.flash(P=P * 0.1, H=res.H(), zs=molefracs)
    assert res.H() == pytest.approx(res2.H(), rel=1e-5)


def test_flash_COCO():
    "Testing som basic results from COCO"
    P = 5.0e6
    T = 300.0
    names = ["methane", "ethane", "propane", "n-butane"]
    molefracs = [0.64, 0.06, 0.28, 0.02]
    flash = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=T,
        rho="eos",
        model="PR",
    )
    res = flash.flash(P=P, T=T, zs=molefracs)

    assert res.gas.beta == pytest.approx(0.932555, rel=1e-3)
    assert res.gas.rho_mass() == pytest.approx(67.92, rel=0.006)
    assert res.gas.V() == pytest.approx(3.6e-4, rel=0.01)

    assert res.liquid0.rho_mass() == pytest.approx(439.92, rel=0.002)
    assert res.liquid0.V() == pytest.approx(8.5792e-5, rel=0.001)

    assert res.gas.zs[0] == pytest.approx(0.66926, rel=1e-3)
    assert res.gas.zs[1] == pytest.approx(0.059714, rel=1e-3)
    assert res.gas.zs[2] == pytest.approx(0.25553, rel=2e-3)
    assert res.gas.zs[3] == pytest.approx(0.015499, rel=9e-3)
    assert res.liquid0.zs[0] == pytest.approx(0.23548, rel=1e-3)
    assert res.liquid0.zs[1] == pytest.approx(0.063952, rel=1e-3)
    assert res.liquid0.zs[2] == pytest.approx(0.61833, rel=1e-3)
    assert res.liquid0.zs[3] == pytest.approx(0.08224, rel=4e-3)
    assert res.H() == pytest.approx(-2867.08, rel=2e-3)
    assert res.gas.H() == pytest.approx(-2163.36, rel=4e-3)
    assert res.liquid0.H() == pytest.approx(-12597.3, rel=2e-3)
    assert res.gas.S() == pytest.approx(-30.3583, rel=5e-4)
    assert res.liquid0.S() == pytest.approx(-57.4518, rel=1e-3)
    assert res.gas.U() == pytest.approx(-3981.73, rel=2e-3)
    assert res.liquid0.U() == pytest.approx(-13026.3, rel=2e-3)
    assert res.gas.Cp() == pytest.approx(66.524, rel=0.02)
    assert res.liquid0.Cp() == pytest.approx(123.90, rel=0.02)
    assert res.gas.mu() == pytest.approx(1.0092e-5, rel=2e-2)


def test_flash_PH_COCO():
    "Testing som basic results from COCO"
    P = 5.0e6
    T = 300.0
    names = ["methane", "ethane", "propane", "n-butane"]
    molefracs = [0.64, 0.06, 0.28, 0.02]
    flash = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=T,
        rho="eos",
        model="PR",
    )
    res = flash.flash(P=P, T=T, zs=molefracs)

    res2 = flash.flash(P=1e5, S=res.S(), zs=molefracs)
    assert res2.T == pytest.approx(194.24, rel=5e-4)

    res2 = flash.flash(P=1e5, H=res.H(), zs=molefracs)
    assert res2.T == pytest.approx(237.53, rel=4e-3)
