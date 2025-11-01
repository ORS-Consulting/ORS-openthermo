import pytest
from openthermo.properties import materials


def test_Cp():
    assert materials.steel_Cp(500.1 + 273.15, "Duplex") == pytest.approx(635, rel=0.01)
    assert materials.steel_Cp(500.1 + 273.15, "SS316") == pytest.approx(530, rel=0.01)
    assert materials.steel_Cp(500.1 + 273.15, "6Mo") == pytest.approx(580, rel=0.01)
    assert materials.steel_Cp(500.1 + 273.15, "CS_235LT") == pytest.approx(
        660, rel=0.01
    )
    assert materials.steel_Cp(500.1 + 273.15, "CS_360LT") == pytest.approx(
        660, rel=0.01
    )


def test_Cp_2():
    assert materials.steel_Cp(1200 + 273.15, "Duplex") == pytest.approx(840, rel=0.01)
    assert materials.steel_Cp(1200 + 273.15, "SS316") == pytest.approx(577, rel=0.01)
    assert materials.steel_Cp(1200 + 273.15, "6Mo") == pytest.approx(610, rel=0.01)
    assert materials.steel_Cp(1200 + 273.15, "CS_235LT") == pytest.approx(540, rel=0.01)
    assert materials.steel_Cp(1200 + 273.15, "CS_360LT") == pytest.approx(540, rel=0.01)


def test_UTS():
    assert materials.UTS(1100 + 273.15, "Duplex") == pytest.approx(58e6, rel=0.01)
    assert materials.UTS(1100 + 273.15, "SS316") == pytest.approx(60e6, rel=0.01)
    assert materials.UTS(1100 + 273.15, "6Mo") == pytest.approx(65e6, rel=0.01)
    assert materials.UTS(1100 + 273.15, "CS_235LT") == pytest.approx(22e6, rel=0.01)
    assert materials.UTS(1100 + 273.15, "CS_360LT") == pytest.approx(27e6, rel=0.01)


def test_UTS_1():
    assert materials.UTS(-100 + 273.15, "Duplex") == pytest.approx(730e6, rel=0.01)
    assert materials.UTS(-100 + 273.15, "SS316") == pytest.approx(575e6, rel=0.01)
    assert materials.UTS(-100 + 273.15, "6Mo") == pytest.approx(730e6, rel=0.01)
    assert materials.UTS(-100 + 273.15, "CS_235LT") == pytest.approx(420e6, rel=0.01)
    assert materials.UTS(-100 + 273.15, "CS_360LT") == pytest.approx(545e6, rel=0.01)


def test_UTS_2():
    assert materials.UTS(525 + 273.15, "Duplex") == pytest.approx(450e6, rel=0.01)
    assert materials.UTS(525 + 273.15, "SS316") == pytest.approx(440e6, rel=0.01)
    assert materials.UTS(525 + 273.15, "6Mo") == pytest.approx(575e6, rel=0.02)
    assert materials.UTS(525 + 273.15, "CS_235LT") == pytest.approx(277e6, rel=0.02)
    assert materials.UTS(525 + 273.15, "CS_360LT") == pytest.approx(365e6, rel=0.02)


def test_ATS():
    assert materials.ATS(1100 + 273.15, "Duplex") == pytest.approx(
        58e6 * 0.85, rel=0.01
    )
    assert materials.ATS(1100 + 273.15, "SS316") == pytest.approx(60e6 * 0.85, rel=0.01)
    assert materials.ATS(1100 + 273.15, "6Mo") == pytest.approx(65e6 * 0.85, rel=0.01)
    assert materials.ATS(1100 + 273.15, "CS_235LT") == pytest.approx(
        22e6 * 0.85, rel=0.01
    )
    assert materials.ATS(1100 + 273.15, "CS_360LT") == pytest.approx(
        27e6 * 0.85, rel=0.01
    )


def test_von_mises():
    """
    Test using manual readings of figure C.2 in:
    Hekkelstrand, B.; Skulstad, P. Guidelines for the Protection of Pressurised
    Systems Exposed to Fire; Scandpower Risk Management AS: Kjeller, Norway, 2004.
    """

    assert materials.von_mises(100e5, 1.0, 0.125) == pytest.approx(49e6, rel=0.02)
    assert materials.von_mises(100e5, 1.0 * 2.5, 0.1) == pytest.approx(121e6, rel=0.02)
