# OpenThermo Examples

This directory contains example YAML input files for running OpenThermo blowdown simulations. Each example corresponds to a test case from `tests/test_blowdown.py` and demonstrates different simulation capabilities.

## Running Examples

Use the `openthermo_main.py` script to run any example:

```bash
# Basic usage
python3 scripts/openthermo_main.py examples/nitrogen_orifice.yml

# With plotting
python3 scripts/openthermo_main.py examples/nitrogen_orifice.yml --plot

# Without progress bar
python3 scripts/openthermo_main.py examples/nitrogen_orifice.yml --no-pbar
```

## Example Files

All 20 test cases from `tests/test_blowdown.py` have been converted to YAML examples.

### Single-Component Gas Examples

#### 1. nitrogen_orifice.yml
**Test:** `test_blowdown_nitrogen`
**Description:** Nitrogen blowdown through a 6.35mm orifice with rigorous heat transfer

**Key Features:**
- Single component (N2) gas phase blowdown
- Isentropic expansion with wall heat transfer
- Small vessel (1.524m x 0.273m)
- High pressure (150 bar → 1 bar)
- Includes experimental validation data from N2_I1.yml

**Expected Results:**
- Final pressure: ~1.08 bar
- Final gas temperature: ~239 K
- Final wall temperature: ~285 K

### 2. nitrogen_control_valve.yml
**Test:** `test_blowdown_nitrogen_control_valve`
**Description:** Same as nitrogen_orifice but using a control valve (Cv=1.25)

**Key Features:**
- Demonstrates control valve discharge vs orifice
- Cv calibrated to match orifice behavior
- Includes experimental validation data

**Expected Results:**
- Final pressure: ~1.14 bar
- Final gas temperature: ~240 K
- Final wall temperature: ~285 K

### 3. isothermal.yml
**Test:** `test_isothermal`
**Description:** Isothermal blowdown of methane-rich mixture

**Key Features:**
- Isothermal mode (constant temperature)
- Multi-component mixture (CH4/C3H8/nC4/iC4/nC10)
- Two-phase flow (liquid level 1.5m)
- Large vessel (10m x 3m)

**Expected Results:**
- Temperature remains constant at 298.15 K
- Final pressure: ~2 bar

### 4. adiabatic.yml
**Test:** `test_adiabatic`
**Description:** Adiabatic blowdown of methane-rich mixture

**Key Features:**
- Adiabatic mode (no heat transfer)
- Multi-component mixture
- Two-phase flow
- Demonstrates Joule-Thomson cooling effect

**Expected Results:**
- Significant temperature decrease due to expansion cooling
- Final pressure: ~1 bar

### 5. isentropic.yml
**Test:** `test_isentropic`
**Description:** Isentropic blowdown of methane-rich mixture

**Key Features:**
- Isentropic expansion
- Multi-component mixture
- Two-phase flow
- Validated against HYSYS Depressurisation utility

**Expected Results:**
- Final pressure: ~1.99 bar
- Final temperature: ~295.5 K (22.3°C)

### 6. adiabatic_cold.yml
**Test:** `test_adiabatic_cold`
**Description:** Cold blowdown - vessel cools to ambient before depressurization

**Key Features:**
- Cold blowdown mode (vessel cools from 350K to 273K before blowdown)
- Multi-component mixture
- Two-phase flow
- Demonstrates pre-cooling effect

**Expected Results:**
- Initial cooling phase 350K → 273K
- Then adiabatic expansion cooling

### Hydrogen Blowdown Examples

#### 7. byrnes_run7.yml
**Test:** `test_byrnes_run7`
**Description:** Fast hydrogen discharge from high pressure vessel

**Key Features:**
- Pure hydrogen, 138 bar → 15 bar
- Fast discharge (2.7mm orifice, 30s duration)
- Vertical vessel (1.394m x 0.217m)
- Rigorous heat transfer with no external HTC
- Validation: Byrnes et al. (1964)

**Expected Results:**
- Final P: ~14 bar
- Final T: ~227-233 K
- Significant Joule-Thomson cooling

#### 8. byrnes_run8.yml
**Test:** `test_byrnes_run8`
**Description:** Slow hydrogen discharge with elevated back pressure

**Key Features:**
- Pure hydrogen, 138 bar → 13.5 bar back pressure
- Slow discharge (0.7mm orifice, 480s duration)
- External HTC: 10 W/(m²·K) natural convection
- Validation: Byrnes et al. (1964)

#### 9. byrnes_run9.yml
**Test:** `test_byrnes_run9`
**Description:** Medium-rate hydrogen discharge

**Key Features:**
- Pure hydrogen, 138 bar
- Medium discharge (4mm orifice, 15s duration)
- No external heat transfer
- Validation: Byrnes et al. (1964)

**Expected Results:**
- Final P: ~11-15 bar
- Final T: ~193-195 K

#### 10. woodfield_discharge.yml
**Test:** `test_woodfield_discharge`
**Description:** Small vessel hydrogen discharge

**Key Features:**
- 99.99% H2, 100 bar
- Very small vessel (0.212m x 0.075m)
- Very slow discharge (0.5mm orifice, 50s)
- External HTC: 5 W/(m²·K)
- Validation: Woodfield et al. (2007)

**Note:** May have convergence issues on some platforms

### CO2 and Mixed Gas Examples

#### 11. co2.yml
**Test:** `test_blowdown_co2`
**Description:** Large CO2 storage vessel blowdown

**Key Features:**
- 95% CO2, 5% N2, 16 bar, 248 K (cold storage)
- Large vessel (35m x 7.8m hemispherical)
- Solver: depressurize_euler (multiphase tracking)
- Duration: 3600s (1 hour)
- Adiabatic walls (external HTC = 0)

#### 12. nitrogen_co2.yml
**Test:** `test_blowdown_nitrogen_co2`
**Description:** N2-CO2 mixture blowdown

**Key Features:**
- 70% N2, 30% CO2, 150 bar, 290 K
- Vertical DIN vessel (1.455m x 0.273m)
- Solver: depressurize_euler (multiphase)
- Duration: 59s
- Validation: Experimental + HYSYS V9

#### 13. ineris_exp16.yml
**Test:** `test_blowdown_ineris_exp16`
**Description:** CO2 pipeline blowdown

**Key Features:**
- 96% CO2, 1.9% N2, 2.1% CH4
- Long pipeline (37m x 0.050m, L/D = 740)
- Solver: depressurize_euler
- Small time step (0.3s) for stability
- Adiabatic walls

### Hydrocarbon Mixture Examples

#### 14. condensable_gas.yml
**Test:** `test_blowdown_condensable_gas`
**Description:** Condensable natural gas blowdown (Haque et al. 1992)

**Key Features:**
- 64% C1, 6% C2, 28% C3, 2% n-C4
- 117.48 bar (116 atm), 293 K
- Horizontal ASME F&D (2.25m x 1.13m, wall=59mm)
- Duration: 1500s (25 minutes)
- Validation: Haque et al. (1992)

#### 15. condensable_gas_rig.yml
**Test:** `test_blowdown_condensable_gas_rig`
**Description:** Vertical condensable gas rig with separate phase tracking

**Key Features:**
- Same composition as condensable_gas.yml
- **Vertical** orientation (tracks liquid/gas separately)
- Solver: depressurize_euler
- Tracks gas_high and gas_low temperatures
- Validation: Experimental bounds

#### 16. non_condensable.yml
**Test:** `test_blowdown_non_condensable`
**Description:** Non-condensable gas blowdown (Haque et al. 1992)

**Key Features:**
- 91% C1, 9% C2 (stays gaseous)
- 121.56 bar (120 atm), 303 K
- Duration: 2000s (~33 minutes)
- Validation: Haque et al. (1992), Wong (1998)

### Fire Scenario Examples

#### 17. sbfire_multiphase.yml
**Test:** `test_blowdown_sbfire_multiphase`
**Description:** Scandpower-Bjordal jet fire with multiphase flow

**Key Features:**
- Heat transfer: rigorous_sb_fire
- Fire type: scandpower_jet
- Hydrocarbon mixture (C1/C3/nC4/iC4/nC10)
- Validation: HYSYS Depressurisation utility

#### 18. sbfire_n2.yml
**Test:** `test_blowdown_sbfire_n2`
**Description:** Nitrogen vessel with API jet fire and relief valve

**Key Features:**
- Fire type: api_jet (conservative API fire)
- Flow device: relief_valve (PSV)
- PSV: 110 bar set, 10% blowdown
- Solver: depressurize_euler
- Duration: 700s

**Expected Results:**
- PSV cycles open/closed
- Temperature rises due to fire heat input

#### 19. sbfire_n2_rupture.yml
**Test:** `test_blowdown_sbfire_n2_rupture`
**Description:** Vessel rupture time evaluation under fire

**Key Features:**
- Fire types: api_jet + scandpower_jet_peak_large
- Material: CS_360LT carbon steel
- Thick wall (136.3mm)
- Undersized relief (8mm orifice)
- Duration: 3600s or until rupture

**Expected Results:**
- Rupture time: ~1585 seconds
- Reference: Andreasen et al. (2018)

#### 20. api_dry_inadequate.yml
**Test:** `test_blowdown_api_dry_inadequate_costald`
**Description:** API 521 fire relief with inadequate drainage

**Key Features:**
- Mode: fire (API 521 sizing)
- Drain/fire fighting: Inadequate
- Liquid density: COSTALD correlation
- Validation: HYSYS

## Complete Example Index

| # | File | Test Type | Fluid | Solver | Duration |
|---|------|-----------|-------|--------|----------|
| 1 | nitrogen_orifice.yml | Validation | N2 | depressurize | 100s |
| 2 | nitrogen_control_valve.yml | Control Valve | N2 | depressurize | 100s |
| 3 | isothermal.yml | Thermodynamic | HC mix | depressurize | 900s |
| 4 | adiabatic.yml | Thermodynamic | HC mix | depressurize | 900s |
| 5 | isentropic.yml | Thermodynamic | HC mix | depressurize | 900s |
| 6 | adiabatic_cold.yml | Cold Blowdown | HC mix | depressurize | 900s |
| 7 | byrnes_run7.yml | Hydrogen | H2 | depressurize | 30s |
| 8 | byrnes_run8.yml | Hydrogen | H2 | depressurize | 480s |
| 9 | byrnes_run9.yml | Hydrogen | H2 | depressurize | 15s |
| 10 | woodfield_discharge.yml | Hydrogen | H2 | depressurize | 50s |
| 11 | co2.yml | CO2 Storage | CO2/N2 | euler | 3600s |
| 12 | nitrogen_co2.yml | Mixed Gas | N2/CO2 | euler | 59s |
| 13 | ineris_exp16.yml | Pipeline | CO2 mix | euler | 100s |
| 14 | condensable_gas.yml | Hydrocarbon | Natural gas | depressurize | 1500s |
| 15 | condensable_gas_rig.yml | Hydrocarbon | Natural gas | euler | 1500s |
| 16 | non_condensable.yml | Hydrocarbon | CH4/C2H6 | depressurize | 2000s |
| 17 | sbfire_multiphase.yml | Fire | HC mix | depressurize | varies |
| 18 | sbfire_n2.yml | Fire + PSV | N2 | euler | 700s |
| 19 | sbfire_n2_rupture.yml | Fire Rupture | N2 | euler | 3600s |
| 20 | api_dry_inadequate.yml | API Fire | HC mix | depressurize | varies |

## YAML File Structure

Each YAML file follows this structure:

```yaml
# Metadata (optional but recommended)
metadata:
  description: "Brief description of the test case"
  reference: "Reference to validation data or source"
  test_name: "Original test function name"

# Solver configuration (required)
solver_method: depressurize  # or 'depressurize_euler'

# Simulation parameters (required)
eos_model: PR              # Equation of state (PR or SRK)
liquid_density: eos        # Liquid density model (eos or costald)
mode: isentropic           # Thermodynamic mode (isothermal, adiabatic, isentropic, fire)

# Vessel geometry (required)
length: 1.524              # m
diameter: 0.273            # m
vessel_type: Flat-end      # Flat-end, ASME F&D, DIN, 2:1 Semi-elliptical, Hemispherical
orientation: horizontal    # horizontal or vertical

# Initial conditions (required)
operating_pressure: 15000000   # Pa
operating_temperature: 289.0   # K
ambient_temperature: 288       # K
back_pressure: 101000          # Pa

# Fluid levels (required)
liquid_level: 0.0          # m (0 for gas-only)
water_level: 0.0           # m

# Flow device (required - orifice OR control_valve OR relief_valve)
# For orifice:
bdv_orifice_size: 0.00635  # m
bdv_orifice_cd: 0.8        # discharge coefficient

# For control valve:
flow_device: control_valve
bdv_Cv: 1.25               # valve flow coefficient
bdv_xT: 0.75               # pressure recovery factor

# For relief valve/PSV:
flow_device: relief_valve
psv_set_pressure: 10000000 # Pa
psv_blowdown: 0.1          # fraction (10%)

# Heat transfer (optional)
heat_transfer: rigorous    # rigorous, rigorous_sb_fire, or omit for simple
wall_thickness: 0.025      # m (required if heat_transfer specified)

# Leak configuration (optional)
leak_active: 0             # 0 or 1
leak_size: 0.01            # m
leak_cd: 0.65              # discharge coefficient
leak_type: liquid          # liquid, gas, or two-phase

# Fluid composition (required)
component_names: [nitrogen]
molefracs: [1.0]

# Simulation control (required)
max_time: 100              # s
delay: 0                   # s (optional)
time_step: 1.0             # s (optional)

# Validation data (optional)
validation:
  temperature:
    gas_mean:              # Use gas_mean for single measurement
      time: [0, 10, 20]
      temp: [289, 250, 240]
    gas_high:              # Or gas_high/gas_low for range
      time: [0, 10, 20]
      temp: [289, 260, 250]
    gas_low:
      time: [0, 10, 20]
      temp: [289, 240, 230]
    wall_mean:             # Use wall_mean for single measurement
      time: [0, 10, 20]
      temp: [289, 287, 285]
    wall_high:             # Or wall_high/wall_low for range
      time: [0, 10, 20]
      temp: [289, 287, 286]
    wall_low:
      time: [0, 10, 20]
      temp: [289, 286, 284]
  pressure:
    time: [0, 10, 20]
    pres: [150, 100, 50]   # bar
```

## Solver Methods

- **depressurize**: Default solver using adaptive time stepping (recommended)
- **depressurize_euler**: Euler method solver (for debugging or special cases)

## Validation Data Conventions

Temperature data uses consistent naming:
- **Gas temperature**: `gas_mean` (single measurement) OR `gas_high`/`gas_low` (range)
- **Wall temperature**: `wall_mean` (single measurement) OR `wall_high`/`wall_low` (range)

The script will:
1. Prefer `*_mean` values if available
2. Fall back to `*_high` values for comparison
3. Plot both high/low ranges when available

## Creating Custom Examples

To create a new example:

1. Copy an existing example file
2. Modify the parameters as needed
3. Update the metadata section
4. Optionally add validation data if experimental results are available
5. Test with: `python3 scripts/openthermo_main.py your_example.yml`

## Validation Data Sources

Validation data included in these examples comes from:
- **N2_I1.yml**: Nitrogen blowdown experimental data
- **HYSYS**: Aspen HYSYS Depressurisation utility results
- Internal test benchmarks

## Notes

- All pressures in Pa (1 bar = 100000 Pa)
- All temperatures in K
- All lengths in m
- Component names must match those in the openthermo database
- Validation data is optional and only used for comparison/plotting
