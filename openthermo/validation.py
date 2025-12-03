from cerberus import Validator


def validate_mandatory_ruleset(input):
    """
    Validate input file using cerberus

    Parameters
    ----------
    input : dict
        Structure holding input

    Return
    ----------
    retval : bool
        True for success, False for failure
    """
    max_components = 20
    schema_general = {
        "eos_model": {"required": True, "type": "string", "allowed": ["SRK", "PR"]},
        "liquid_density": {
            "required": True,
            "type": "string",
            "allowed": ["eos", "costald"],
        },
        "molefracs": {
            "required": True,
            "type": "list",
            "minlength": 1,
            "maxlength": max_components,
            "schema": {"type": "number", "min": 1e-5, "max": 1},
        },
        "component_names": {
            "required": True,
            "type": "list",
            "minlength": 1,
            "maxlength": max_components,
            "schema": {"type": "string"},
        },
        "pseudo_names": {
            "required": False,
            "type": "list",
            "minlength": 1,
            "maxlength": max_components,
            "schema": {"type": "string"},
        },
        "pseudo_molefracs": {
            "required": False,
            "type": "list",
            "minlength": 1,
            "maxlength": max_components,
            "schema": {"type": "number", "min": 1e-5, "max": 1},
            "dependencies": ["pseudo_names", "pseudo_SGs", "pseudo_Tbs"],
        },
        "pseudo_SGs": {
            "required": False,
            "type": "list",
            "minlength": 1,
            "maxlength": max_components,
            "schema": {"type": "number", "min": 1e-5, "max": 1},
            "dependencies": ["pseudo_names", "pseudo_molefracs", "pseudo_Tbs"],
        },
        "pseudo_Tbs": {
            "required": False,
            "type": "list",
            "minlength": 1,
            "maxlength": max_components,
            "schema": {"type": "number", "min": 1e-5, "max": 1},
            "dependencies": ["pseudo_names", "pseudo_SGs", "pseudo_molefracs"],
        },
        "mode": {
            "required": True,
            "type": "string",
            "allowed": ["isothermal", "adiabatic", "isentropic", "fire"],
        },
        "fire_type": {
            "required": False,
            "type": "string",
            "allowed": ["API521", "API521_CONFINED"],
            "dependencies": {"mode": "fire"},
        },
        "drain_fire_fighting": {
            "required": False,
            "type": "string",
            "allowed": ["Inadequate", "Adequate"],
        },
        "environment_factor": {
            "required": False,
            "type": "number",
            "min": 0.01,
            "max": 1.0,
        },
        "leak_type": {
            "required": False,
            "type": "string",
            "allowed": ["liquid", "gas", "two-phase"],
        },
        "heat_transfer": {
            "required": False,
            "type": "string",
            "allowed": ["rigorous", "rigorous_sb_fire"],
        },
        "external_heat_transfer_coefficient": {
            "required": False,
            "type": "number",
            "min": 0.0,
            "max": 1000,
        },
        "sb_fire_type": {
            "required": False,
            "type": "string",
            "allowed": ["scandpower_pool", "scandpower_jet", "api_pool", "api_jet"],
            "dependencies": {"heat_transfer": "rigorous_sb_fire"},
        },
        "sb_peak_fire_type": {
            "required": False,
            "type": "string",
            "allowed": [
                "scandpower_jet_peak_large",
                "scandpower_jet_peak_small",
                "scandpower_pool_peak",
            ],
            "dependencies": {"heat_transfer": ["rigorous", "rigorous_sb_fire"]},
        },
        "vessel_material": {
            "required": False,
            "type": "string",
            "allowed": ["CS_235LT", "CS_360LT", "SS316", "Duplex", "6Mo"],
            "dependencies": ["sb_peak_fire_type"],
        },
        "orientation": {
            "required": True,
            "type": "string",
            "allowed": ["horizontal", "vertical"],
        },
        "vessel_type": {
            "required": True,
            "type": "string",
            "allowed": [
                "Flat-end",
                "ASME F&D",
                "DIN",
                "2:1 Semi-elliptical",
                "Hemispherical",
            ],
        },
        "leak_active": {
            "required": False,
            "type": "number",
            "allowed": [0, 1],
            "dependencies": ["leak_size", "leak_cd", "leak_type"],
        },
        "leak_size": {
            "required": False,
            "type": "number",
            "dependencies": ["leak_active", "leak_cd", "leak_type"],
        },
        "leak_cd": {
            "required": False,
            "type": "number",
            "min": 0.01,
            "max": 1.0,
            "dependencies": ["leak_active", "leak_size", "leak_type"],
        },
        "liquid_level": {
            "required": True,
            "type": "number",
            "min": 0,
        },
        "water_level": {
            "required": False,
            "type": "number",
            "min": 0,
        },
        "length": {
            "required": True,
            "type": "number",
            "min": 0.1,
        },
        "diameter": {
            "required": True,
            "type": "number",
            "min": 0.1,
        },
        "bdv_orifice_size": {
            "required": True,
            "type": "number",
            "min": 1e-4,
        },
        "bdv_orifice_cd": {
            "required": True,
            "type": "number",
            "min": 0.1,
            "max": 1.0,
        },
        "wall_thickness": {
            "required": False,
            "type": "number",
            "min": 1e-3,
        },
        "operating_temperature": {
            "required": True,
            "type": "number",
            "min": 10,
        },
        "max_time": {"required": True, "type": "number", "min": 1},
        "time_step": {
            "required": False,
            "type": "number",
            "min": 0.01,
        },
        "delay": {
            "required": False,
            "type": "number",
            "min": 0,
        },
        "ambient_temperature": {
            "required": False,
            "type": "number",
            "min": 250,
        },
        "operating_pressure": {
            "required": True,
            "type": "number",
            "min": 1e5,
        },
        "back_pressure": {
            "required": True,
            "type": "number",
            "min": 1e5,
        },
    }

    v = Validator(schema_general)
    retval = v.validate(input)
    if v.errors:
        raise ValueError("Input validation error", v.errors)

    return retval
