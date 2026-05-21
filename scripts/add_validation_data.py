#!/usr/bin/env python3
"""
Script to add validation data to example YAML files from validation sources
"""

import os
import yaml
import numpy as np
import pandas as pd

# Paths
validation_dir = "validation"
examples_dir = "examples"

def read_txt_data(filepath):
    """Read tab-delimited .txt file"""
    data = np.loadtxt(filepath, delimiter="\t")
    return {"time": data[:, 0].tolist(), "values": data[:, 1].tolist()}

def read_csv_data(filepath, skiprows=None):
    """Read CSV file with optional row skipping"""
    if skiprows is None:
        # Try to detect header rows automatically
        with open(filepath, 'r') as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                # Look for the data separator line (all dashes)
                if '-----' in line:
                    skiprows = i + 2  # Skip separator and column header rows
                    break

    data = pd.read_csv(filepath, skiprows=skiprows)
    return data

def update_yaml_with_validation(yaml_file, validation_data):
    """Update YAML file with validation data"""
    with open(yaml_file, 'r') as f:
        content = yaml.safe_load(f)

    # Add or update validation section
    content['validation'] = validation_data

    with open(yaml_file, 'w') as f:
        yaml.dump(content, f, default_flow_style=False, sort_keys=False, width=120)

    print(f"Updated {yaml_file}")

# Example mappings of test cases to validation files
validation_mappings = {
    'isentropic.yml': {
        'csv_file': 'Water_dry_isentropic_history.csv',
        'extract': lambda df: {
            'pressure': {
                'comment': 'Gauge pressure in bar',
                'pres_gauge': df.iloc[:, 2].tolist(),  # Column 3
                'comment_total': 'Add 1.013 bar for absolute pressure'
            },
            'temperature': {
                'fluid_mean': {
                    'temp': (df.iloc[:, 1] + 273.15).tolist(),  # Column 2, convert C to K
                    'comment': 'Converted from Celsius to Kelvin'
                }
            },
            'mass_flow': {
                'mdot': (df.iloc[:, 3] / -3600).tolist(),  # Column 4, convert kg/h to kg/s
                'comment': 'Converted from kg/h to kg/s'
            }
        }
    },

    'condensable_gas.yml': {
        'txt_files': {
            'pressure': 'condensable_gas_pressure.txt',
            'gas_high': 'condensable_gas_gas_temp_higher.txt',
            'gas_low': 'condensable_gas_gas_temp_lower.txt',
            'wall_high': 'condensable_gas_inner_wall_temp.txt',
            'wall_low': 'condensable_gas_liquid_inner_wall_temp.txt',
        }
    },

    'non_condensable.yml': {
        'txt_files': {
            'pressure': 'non_condensable_pressure.txt',
            'gas_high': 'non_condensable_gas_temp_higher.txt',
            'gas_low': 'non_condensable_gas_temp_lower.txt',
            'wall': 'non_condensable_wall_temp.txt',
        }
    },

    'nitrogen_co2.yml': {
        'txt_files': {
            'gas_high': 'co2_n2_gas_higher.txt',
            'gas_low': 'co2_n2_gas_lower.txt',
            'liquid_high': 'co2_n2_liquid_higher.txt',
            'liquid_low': 'co2_n2_liquid_lower.txt',
            'HYSYS_gas': 'co2_n2_HYSYS_gas.txt',
            'HYSYS_liquid': 'co2_n2_HYSYS_liquid.txt',
            'HYSYS_gas_wall': 'co2_n2_HYSYS_gas_wall.txt',
            'HYSYS_liquid_wall': 'co2_n2_HYSYS_liquid_wall.txt',
        }
    },
}

def main():
    """Main function to process all validation mappings"""

    # Process condensable_gas.yml
    print("Processing condensable_gas.yml...")
    try:
        validation_data = {'temperature': {}, 'pressure': {}}

        # Read pressure
        pres_data = read_txt_data(os.path.join(validation_dir, 'condensable_gas_pressure.txt'))
        validation_data['pressure'] = {
            'time': pres_data['time'],
            'pres': pres_data['values'],
            'comment': 'Pressure in bar (absolute)'
        }

        # Read gas temperatures
        gas_high = read_txt_data(os.path.join(validation_dir, 'condensable_gas_gas_temp_higher.txt'))
        validation_data['temperature']['gas_high'] = {
            'time': gas_high['time'],
            'temp': gas_high['values']
        }

        gas_low = read_txt_data(os.path.join(validation_dir, 'condensable_gas_gas_temp_lower.txt'))
        validation_data['temperature']['gas_low'] = {
            'time': gas_low['time'],
            'temp': gas_low['values']
        }

        # Read wall temperatures
        wall_high = read_txt_data(os.path.join(validation_dir, 'condensable_gas_inner_wall_higher.txt'))
        validation_data['temperature']['wall_high'] = {
            'time': wall_high['time'],
            'temp': wall_high['values']
        }

        wall_low = read_txt_data(os.path.join(validation_dir, 'condensable_gas_liquid_inner_wall_lower.txt'))
        validation_data['temperature']['wall_low'] = {
            'time': wall_low['time'],
            'temp': wall_low['values']
        }

        update_yaml_with_validation(os.path.join(examples_dir, 'condensable_gas.yml'), validation_data)
    except Exception as e:
        print(f"Error processing condensable_gas.yml: {e}")

    # Process condensable_gas_rig.yml (same as condensable_gas but with liquid temperatures)
    print("Processing condensable_gas_rig.yml...")
    try:
        validation_data = {'temperature': {}, 'pressure': {}}

        # Read pressure
        pres_data = read_txt_data(os.path.join(validation_dir, 'condensable_gas_pressure.txt'))
        validation_data['pressure'] = {
            'time': pres_data['time'],
            'pres': pres_data['values'],
            'comment': 'Pressure in bar (absolute)'
        }

        # Read gas temperatures
        gas_high = read_txt_data(os.path.join(validation_dir, 'condensable_gas_gas_temp_higher.txt'))
        validation_data['temperature']['gas_high'] = {
            'time': gas_high['time'],
            'temp': gas_high['values']
        }

        gas_low = read_txt_data(os.path.join(validation_dir, 'condensable_gas_gas_temp_lower.txt'))
        validation_data['temperature']['gas_low'] = {
            'time': gas_low['time'],
            'temp': gas_low['values']
        }

        # Read liquid temperatures
        liquid_high = read_txt_data(os.path.join(validation_dir, 'condensable_gas_liq_temp_higher.txt'))
        validation_data['temperature']['liquid_high'] = {
            'time': liquid_high['time'],
            'temp': liquid_high['values']
        }

        liquid_low = read_txt_data(os.path.join(validation_dir, 'condensable_gas_liq_temp_lower.txt'))
        validation_data['temperature']['liquid_low'] = {
            'time': liquid_low['time'],
            'temp': liquid_low['values']
        }

        # Read wall temperatures
        wall_high = read_txt_data(os.path.join(validation_dir, 'condensable_gas_inner_wall_higher.txt'))
        validation_data['temperature']['wall_high'] = {
            'time': wall_high['time'],
            'temp': wall_high['values']
        }

        wall_low = read_txt_data(os.path.join(validation_dir, 'condensable_gas_liquid_inner_wall_lower.txt'))
        validation_data['temperature']['wall_low'] = {
            'time': wall_low['time'],
            'temp': wall_low['values']
        }

        update_yaml_with_validation(os.path.join(examples_dir, 'condensable_gas_rig.yml'), validation_data)
    except Exception as e:
        print(f"Error processing condensable_gas_rig.yml: {e}")

    # Process non_condensable.yml
    print("Processing non_condensable.yml...")
    try:
        validation_data = {'temperature': {}, 'pressure': {}}

        pres_data = read_txt_data(os.path.join(validation_dir, 'non-condensable_pressure.txt'))
        validation_data['pressure'] = {
            'time': pres_data['time'],
            'pres': pres_data['values']
        }

        gas_high = read_txt_data(os.path.join(validation_dir, 'non-condensable_gas_temp_higher.txt'))
        validation_data['temperature']['gas_high'] = {
            'time': gas_high['time'],
            'temp': gas_high['values']
        }

        gas_low = read_txt_data(os.path.join(validation_dir, 'non-condensable_gas_temp_lower.txt'))
        validation_data['temperature']['gas_low'] = {
            'time': gas_low['time'],
            'temp': gas_low['values']
        }

        wall = read_txt_data(os.path.join(validation_dir, 'non-condensable_wall_temp.txt'))
        validation_data['temperature']['wall_mean'] = {
            'time': wall['time'],
            'temp': wall['values']
        }

        update_yaml_with_validation(os.path.join(examples_dir, 'non_condensable.yml'), validation_data)
    except Exception as e:
        print(f"Error processing non_condensable.yml: {e}")

    # Process nitrogen_co2.yml
    print("Processing nitrogen_co2.yml...")
    try:
        validation_data = {'temperature': {}}

        # Read gas temperatures (in Celsius, convert to Kelvin)
        gas_high = read_txt_data(os.path.join(validation_dir, 'co2_n2_gas_higher.txt'))
        validation_data['temperature']['gas_high'] = {
            'time': gas_high['time'],
            'temp': [t + 273.15 for t in gas_high['values']],
            'comment': 'Converted from Celsius to Kelvin'
        }

        gas_low = read_txt_data(os.path.join(validation_dir, 'co2_n2_gas_lower.txt'))
        validation_data['temperature']['gas_low'] = {
            'time': gas_low['time'],
            'temp': [t + 273.15 for t in gas_low['values']],
            'comment': 'Converted from Celsius to Kelvin'
        }

        # Read liquid temperatures (in Celsius, convert to Kelvin)
        liquid_high = read_txt_data(os.path.join(validation_dir, 'co2_n2_liquid_higher.txt'))
        validation_data['temperature']['liquid_high'] = {
            'time': liquid_high['time'],
            'temp': [t + 273.15 for t in liquid_high['values']],
            'comment': 'Converted from Celsius to Kelvin'
        }

        liquid_low = read_txt_data(os.path.join(validation_dir, 'co2_n2_liquid_lower.txt'))
        validation_data['temperature']['liquid_low'] = {
            'time': liquid_low['time'],
            'temp': [t + 273.15 for t in liquid_low['values']],
            'comment': 'Converted from Celsius to Kelvin'
        }

        # Add HYSYS data (already in Kelvin, no conversion needed)
        hysys_gas = read_txt_data(os.path.join(validation_dir, 'co2_n2_HYSYS_gas.txt'))
        validation_data['temperature']['HYSYS_gas'] = {
            'time': hysys_gas['time'],
            'temp': hysys_gas['values'],
            'comment': 'HYSYS simulation data (already in Kelvin)'
        }

        hysys_liquid = read_txt_data(os.path.join(validation_dir, 'co2_n2_HYSYS_liquid.txt'))
        validation_data['temperature']['HYSYS_liquid'] = {
            'time': hysys_liquid['time'],
            'temp': hysys_liquid['values'],
            'comment': 'HYSYS simulation data (already in Kelvin)'
        }

        update_yaml_with_validation(os.path.join(examples_dir, 'nitrogen_co2.yml'), validation_data)
    except Exception as e:
        print(f"Error processing nitrogen_co2.yml: {e}")

    # Process isentropic.yml (CSV file)
    print("Processing isentropic.yml...")
    try:
        csv_data = read_csv_data(os.path.join(validation_dir, 'Water_dry_isentropic_history.csv'))

        # Downsample HYSYS data - take every 100th point for cleaner plots
        downsample_factor = 100
        csv_data_downsampled = csv_data.iloc[::downsample_factor, :]

        validation_data = {
            'pressure': {
                'time': csv_data_downsampled.iloc[:, 0].tolist(),  # Column 1 (time)
                'pres': (csv_data_downsampled.iloc[:, 2] + 1.013).tolist(),  # Column 3, convert gauge to absolute
                'comment': 'Pressure in bar (absolute), converted from gauge by adding 1.013 bar, downsampled (every 100th point)'
            },
            'temperature': {
                'fluid_mean': {
                    'time': csv_data_downsampled.iloc[:, 0].tolist(),  # Column 1 (time)
                    'temp': (csv_data_downsampled.iloc[:, 1] + 273.15).tolist(),  # Column 2, convert C to K
                    'comment': 'Vapour temperature, converted from Celsius to Kelvin, downsampled (every 100th point)'
                }
            },
            'mass_flow': {
                'time': csv_data_downsampled.iloc[:, 0].tolist(),  # Column 1 (time)
                'mdot': (csv_data_downsampled.iloc[:, 3] / -3600).tolist(),  # Column 4, convert kg/h to kg/s
                'comment': 'Converted from kg/h to kg/s, negative indicates discharge, downsampled (every 100th point)'
            }
        }

        update_yaml_with_validation(os.path.join(examples_dir, 'isentropic.yml'), validation_data)
    except Exception as e:
        print(f"Error processing isentropic.yml: {e}")

    # Process sbfire_multiphase.yml (CSV files from unisim)
    print("Processing sbfire_multiphase.yml...")
    try:
        validation_data = {'temperature': {}, 'pressure': {}, 'mass_flow': {}}

        # Read pressure data (in kPa, convert to bar)
        pres_csv = read_csv_data(os.path.join(validation_dir, 'unisim_sb_fire_pressure.csv'))
        validation_data['pressure'] = {
            'time': pres_csv.iloc[:, 0].tolist(),
            'pres': (pres_csv.iloc[:, 1] / 100.0).tolist(),  # Convert kPa to bar
            'comment': 'Pressure from Unisim EO Blowdown, converted from kPa to bar'
        }

        # Read temperature data (in Celsius, convert to Kelvin)
        temp_csv = read_csv_data(os.path.join(validation_dir, 'unisim_sb_fire_fluid_temperatures.csv'))
        validation_data['temperature']['gas_mean'] = {
            'time': temp_csv.iloc[:, 0].tolist(),
            'temp': (temp_csv.iloc[:, 1] + 273.15).tolist(),  # Vapour zone temp, C to K
            'comment': 'Vapour zone temperature from Unisim EO Blowdown, converted from Celsius to Kelvin'
        }
        validation_data['temperature']['liquid_mean'] = {
            'time': temp_csv.iloc[:, 0].tolist(),
            'temp': (temp_csv.iloc[:, 2] + 273.15).tolist(),  # Liquid zone temp, C to K
            'comment': 'Liquid zone temperature from Unisim EO Blowdown, converted from Celsius to Kelvin'
        }

        # Read wall temperature data (in Celsius, convert to Kelvin)
        wall_csv = read_csv_data(os.path.join(validation_dir, 'unisim_wall_temperature.csv'))
        # Unwetted wall (vapour/gas zone)
        validation_data['temperature']['wall_unwetted_inner'] = {
            'time': wall_csv.iloc[:, 0].tolist(),
            'temp': (wall_csv.iloc[:, 1] + 273.15).tolist(),  # InsideWallTemp[Vapour], C to K
            'comment': 'Unwetted (vapour zone) inner wall temperature from Unisim EO Blowdown, converted from Celsius to Kelvin'
        }
        validation_data['temperature']['wall_unwetted_outer'] = {
            'time': wall_csv.iloc[:, 0].tolist(),
            'temp': (wall_csv.iloc[:, 2] + 273.15).tolist(),  # OutsideWallTemp[Vapour], C to K
            'comment': 'Unwetted (vapour zone) outer wall temperature from Unisim EO Blowdown, converted from Celsius to Kelvin'
        }
        # Wetted wall (liquid zone)
        validation_data['temperature']['wall_wetted_inner'] = {
            'time': wall_csv.iloc[:, 0].tolist(),
            'temp': (wall_csv.iloc[:, 3] + 273.15).tolist(),  # InsideWallTemp[Liquid], C to K
            'comment': 'Wetted (liquid zone) inner wall temperature from Unisim EO Blowdown, converted from Celsius to Kelvin'
        }
        validation_data['temperature']['wall_wetted_outer'] = {
            'time': wall_csv.iloc[:, 0].tolist(),
            'temp': (wall_csv.iloc[:, 4] + 273.15).tolist(),  # OutsideWallTemp[Liquid], C to K
            'comment': 'Wetted (liquid zone) outer wall temperature from Unisim EO Blowdown, converted from Celsius to Kelvin'
        }

        # Read mass flow data
        flow_csv = read_csv_data(os.path.join(validation_dir, 'unisim_sb_mass_flow.csv'))
        validation_data['mass_flow'] = {
            'time': flow_csv.iloc[:, 0].tolist(),
            'mdot': flow_csv.iloc[:, 1].tolist(),
            'comment': 'Mass flow from Unisim EO Blowdown'
        }

        update_yaml_with_validation(os.path.join(examples_dir, 'sbfire_multiphase.yml'), validation_data)
    except Exception as e:
        print(f"Error processing sbfire_multiphase.yml: {e}")

    # Process api_dry_inadequate.yml (CSV file)
    print("Processing api_dry_inadequate.yml...")
    try:
        csv_data = read_csv_data(os.path.join(validation_dir, 'Water_dry_API_inadequate_costald_history.csv'))

        # Downsample HYSYS data - take every 100th point for cleaner plots
        downsample_factor = 100
        csv_data_downsampled = csv_data.iloc[::downsample_factor, :]

        validation_data = {
            'pressure': {
                'time': csv_data_downsampled.iloc[:, 0].tolist(),  # Column 1 (time)
                'pres': (csv_data_downsampled.iloc[:, 2] + 1.013).tolist(),  # Column 3, convert gauge to absolute
                'comment': 'Pressure in bar (absolute), converted from gauge by adding 1.013 bar, downsampled (every 100th point)'
            },
            'temperature': {
                'fluid_mean': {
                    'time': csv_data_downsampled.iloc[:, 0].tolist(),  # Column 1 (time)
                    'temp': (csv_data_downsampled.iloc[:, 1] + 273.15).tolist(),  # Column 2, convert C to K
                    'comment': 'Converted from Celsius to Kelvin, downsampled (every 100th point)'
                }
            },
            'mass_flow': {
                'time': csv_data_downsampled.iloc[:, 0].tolist(),  # Column 1 (time)
                'mdot': (csv_data_downsampled.iloc[:, 3] / -3600).tolist(),  # Column 4, convert kg/h to kg/s
                'comment': 'Converted from kg/h to kg/s, negative indicates discharge, downsampled (every 100th point)'
            }
        }

        update_yaml_with_validation(os.path.join(examples_dir, 'api_dry_inadequate.yml'), validation_data)
    except Exception as e:
        print(f"Error processing api_dry_inadequate.yml: {e}")

    print("\nValidation data extraction complete!")

if __name__ == "__main__":
    main()
