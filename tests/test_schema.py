from openthermo.validation import validate_mandatory_ruleset


def test_validate_mandatory_ruleset():
    input = {}
    P = 116e5
    T = 293.0

    input["mode"] = "isentropic"
    input["heat_transfer"] = "rigorous"
    input["wall_thickness"] = 0.059  # m
    input["eos_model"] = "PR"
    input["liquid_density"] = "eos"
    input["max_time"] = 1500
    input["delay"] = 0
    input["time_step"] = 10
    input["length"] = 2.25
    input["diameter"] = 1.13
    input["vessel_type"] = "ASME F&D"
    input["orientation"] = "vertical"
    input["liquid_level"] = 0
    input["water_level"] = 0.0
    input["operating_temperature"] = T
    input["operating_pressure"] = P
    input["ambient_temperature"] = 293.0
    input["back_pressure"] = 1.01e5
    input["bdv_orifice_size"] = 0.01  # m
    input["bdv_orifice_cd"] = 0.8

    input["leak_active"] = 0
    input["leak_size"] = 0.01  # m

    input["leak_cd"] = 0.65
    input["leak_type"] = "liquid"

    names = ["methane", "ethane", "propane", "n-butane"]
    molefracs = [0.64, 0.06, 0.28, 0.02]

    input["molefracs"] = molefracs

    assert validate_mandatory_ruleset(input)


if __name__ == "__main__":
    test_validate_mandatory_ruleset()
    print("test_validate_mandatory_ruleset passed")
