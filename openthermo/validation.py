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
            "maxlength": 20,
            "schema": {"type": "number", "min": 1e-5, "max": 1},
        },
        "mode": {
            "required": True,
            "type": "string",
            "allowed": ["isothermal", "adiabatic", "isentropic", "fire"],
        },
        "drain_fire_fighting": {
            "required": False,
            "type": "string",
            "allowed": ["Inadequate", "Adequate"],
            "dependencies": {"mode": "fire"},
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
        "sb_fire_type": {
            "required": False,
            "type": "string",
            "allowed": ["scandpower_pool", "scandpower_jet", "api_pool", "api_jet"],
            "dependencies": {"heat_transfer": "rigorous_sb_fire"},
        },
        "orientation": {
            "required": True,
            "type": "string",
            "allowed": ["horizontal", "vertical"],
        },
        "vessel_type": {
            "required": True,
            "type": "string",
            "allowed": ["Flat-end", "ASME F&D", "DIN"],
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
            "min": 100,
        },
        "max_time": {
            "required": True,
            "type": "number",
        },
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
