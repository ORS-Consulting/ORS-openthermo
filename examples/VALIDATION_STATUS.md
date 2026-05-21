# Validation Data Status

This document shows which example files have validation data included.

## Examples with Validation Data

### From Tab-Delimited TXT Files
1. **condensable_gas.yml** - ✓ Complete
   - Pressure (time, pres)
   - Gas temperature (high, low)
   - Wall temperature (high, low)
   - Source: Haque et al. (1992) experiments

2. **condensable_gas_rig.yml** - ✓ Complete
   - Pressure (time, pres)
   - Gas temperature (high, low)
   - Liquid temperature (high, low)
   - Wall temperature (high, low)
   - Source: Haque et al. (1992) experiments

3. **non_condensable.yml** - ✓ Complete
   - Pressure (time, pres)
   - Gas temperature (high, low)
   - Wall temperature (mean)
   - Source: Haque et al. (1992), Wong (1998)

4. **nitrogen_co2.yml** - ✓ Complete
   - Gas temperature (high, low)
   - Liquid temperature (high, low)
   - HYSYS gas temperature
   - HYSYS liquid temperature
   - Source: Experimental + HYSYS V9

### From HYSYS/UniSim CSV Files
5. **isentropic.yml** - ✓ Complete
   - Pressure (gauge, bar)
   - Fluid temperature (mean, converted C to K)
   - Mass flow (kg/s)
   - Source: HYSYS Depressurisation utility

6. **sbfire_multiphase.yml** - ✓ Complete
   - Pressure (time, pres)
   - Fluid temperature (mean)
   - Wall temperature (mean)
   - Mass flow
   - Source: UniSim simulation

7. **api_dry_inadequate.yml** - ✓ Complete
   - Pressure (gauge, bar)
   - Fluid temperature (mean, converted C to K)
   - Mass flow (kg/s)
   - Source: HYSYS Depressurisation - API 521

### From Experimental Publications (Embedded in YAML)
8. **nitrogen_orifice.yml** - ✓ Complete
   - Pressure (time, pres)
   - Gas temperature (high, low)
   - Wall temperature (high, low)
   - Source: N2_I1.yml validation data

9. **nitrogen_control_valve.yml** - ✓ Complete
   - Pressure (time, pres)
   - Gas temperature (high, low)
   - Wall temperature (high, low)
   - Source: N2_I1.yml validation data

10. **byrnes_run7.yml** - ✓ Complete
    - Pressure (experimental, hyddown)
    - Gas temperature (experimental, hyddown)
    - Wall temperature (experimental, hyddown)
    - Source: Byrnes et al. (1964)

11. **byrnes_run8.yml** - ✓ Complete
    - Pressure (experimental, hyddown)
    - Gas temperature (experimental, hyddown)
    - Wall temperature (experimental, hyddown)
    - Source: Byrnes et al. (1964)

12. **byrnes_run9.yml** - ✓ Complete
    - Pressure (experimental, hyddown)
    - Gas temperature (experimental, hyddown)
    - Wall temperature (experimental, hyddown)
    - Source: Byrnes et al. (1964)

13. **woodfield_discharge.yml** - ✓ Complete
    - Pressure (experimental, hyddown)
    - Gas temperature (experimental, hyddown)
    - Wall temperature (experimental, hyddown)
    - Source: Woodfield et al. (2007)

## Examples without Validation Data

The following examples don't have validation data as they are comparison/demonstration cases:

14. isothermal.yml
15. adiabatic.yml
16. adiabatic_cold.yml
17. co2.yml
18. ineris_exp16.yml
19. sbfire_n2.yml
20. sbfire_n2_rupture.yml

## Validation Data Formats

### Standard Format
```yaml
validation:
  pressure:
    time: [0, 10, 20, ...]
    pres: [150, 100, 50, ...]  # bar
  temperature:
    gas_mean:  # or gas_high/gas_low for range
      time: [0, 10, 20, ...]
      temp: [289, 250, 240, ...]  # K
    liquid_mean:  # or liquid_high/liquid_low for multiphase
      time: [0, 10, 20, ...]
      temp: [289, 250, 240, ...]  # K
    wall_mean:  # or wall_high/wall_low
      time: [0, 10, 20, ...]
      temp: [289, 287, 285, ...]  # K
```

**Important**: For multiphase cases with separate wall measurements:
- `wall_high` = unwetted wall (gas side) temperature
- `wall_low` = wetted wall (liquid side) temperature

### Byrnes/HydDown Format (Nested)
```yaml
validation:
  pressure:
    experimental:
      time: [...]
      pres: [...]
    hyddown:
      time: [...]
      pres: [...]
  temperature:
    gas_experimental:
      time: [...]
      temp: [...]
    gas_hyddown:
      time: [...]
      temp: [...]
```

The `openthermo_main.py` script automatically handles both formats.

## Usage

To run a simulation with validation comparison:

```bash
# Run simulation and show validation comparison
python3 scripts/openthermo_main.py examples/byrnes_run7.yml

# Run with plotting to visualize validation data
python3 scripts/openthermo_main.py examples/nitrogen_orifice.yml --plot
```

The script will automatically:
1. Display final value comparisons in the console
2. Plot experimental/simulation data together (if --plot is used)
3. Handle both flat and nested validation structures
