from openthermo.flash.michelsen import get_flash_dry

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
        "eos",
        "PR",
    )
    print(flash)
    pass