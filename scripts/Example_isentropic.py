from openthermo.vessel.blowdown import Blowdown
from openthermo.flash.michelsen import get_flash_dry

input = {}
P = 12e5
T = 298.15

input["mode"] = "isentropic"
input["eos_model"] = "PR"
input["liquid_density"] = "eos"
input["max_time"] = 900
input["length"] = 10
input["diameter"] = 3
input["vessel_type"] = "Flat-end"
input["orientation"] = "horizontal"
input["liquid_level"] = 1.5
input["operating_temperature"] = T
input["operating_pressure"] = P
input["ambient_temperature"] = 273
input["back_pressure"] = 1.01e5
input["bdv_orifice_size"] = 0.03  # m
input["bdv_orifice_cd"] = 0.84

names = ["methane", "propane", "n-butane", "i-butane", "n-decane"]
molefracs = [0.8, 0.05, 0.01, 0.01, 0.10]

input["molefracs"] = molefracs

input["flash"] = get_flash_dry(
    names,
    molefracs,
    P=P,
    T=T,
    rho=input["liquid_density"],
    model=input["eos_model"],
)

segment = Blowdown(input)
segment.depressurize()
