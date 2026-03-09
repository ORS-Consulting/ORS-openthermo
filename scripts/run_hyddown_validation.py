"""
Run HydDown validation cases and save results for comparison with openthermo.
"""
import sys
sys.path.insert(0, '/home/anra/github/HydDown/src')

from hyddown import HydDown
import numpy as np
import yaml
import pickle

# Run Byrnes Run 7
print("Running HydDown - Byrnes Run 7...")
with open('/home/anra/github/HydDown/validation/Byrnes_run7.yml', 'r') as f:
    config7 = yaml.safe_load(f)

hd7 = HydDown(config7)
hd7.run()

print(f"Byrnes Run 7 - HydDown Results:")
print(f"Final pressure: {hd7.P[-1]/1e5:.1f} bar")
print(f"Final temperature: {hd7.T_fluid[-1]:.1f} K")
print(f"Final time: {hd7.time_array[-1]:.1f} s")

# Run Byrnes Run 8
print("\nRunning HydDown - Byrnes Run 8...")
with open('/home/anra/github/HydDown/validation/Byrnes_run8.yml', 'r') as f:
    config8 = yaml.safe_load(f)

hd8 = HydDown(config8)
hd8.run()

print(f"Byrnes Run 8 - HydDown Results:")
print(f"Final pressure: {hd8.P[-1]/1e5:.1f} bar")
print(f"Final temperature: {hd8.T_fluid[-1]:.1f} K")
print(f"Final time: {hd8.time_array[-1]:.1f} s")

# Run Byrnes Run 9
print("\nRunning HydDown - Byrnes Run 9...")
with open('/home/anra/github/HydDown/validation/Byrnes_run9.yml', 'r') as f:
    config9 = yaml.safe_load(f)

hd9 = HydDown(config9)
hd9.run()

print(f"Byrnes Run 9 - HydDown Results:")
print(f"Final pressure: {hd9.P[-1]/1e5:.1f} bar")
print(f"Final temperature: {hd9.T_fluid[-1]:.1f} K")
print(f"Final time: {hd9.time_array[-1]:.1f} s")

# Run Woodfield Discharge
print("\nRunning HydDown - Woodfield Discharge...")
with open('/home/anra/github/HydDown/validation/dischargeH2_woodfield.yml', 'r') as f:
    config_wood = yaml.safe_load(f)

hd_wood = HydDown(config_wood)
hd_wood.run()

print(f"Woodfield Discharge - HydDown Results:")
print(f"Final pressure: {hd_wood.P[-1]/1e5:.1f} bar")
print(f"Final temperature: {hd_wood.T_fluid[-1]:.1f} K")
print(f"Final time: {hd_wood.time_array[-1]:.1f} s")

# Save results for use in tests
results = {
    'run7': {
        'time': np.array(hd7.time_array),
        'pressure': np.array(hd7.P),
        'temperature': np.array(hd7.T_fluid),
        'wall_temp': np.array(hd7.T_outer_wall),
    },
    'run8': {
        'time': np.array(hd8.time_array),
        'pressure': np.array(hd8.P),
        'temperature': np.array(hd8.T_fluid),
        'wall_temp': np.array(hd8.T_outer_wall),
    },
    'run9': {
        'time': np.array(hd9.time_array),
        'pressure': np.array(hd9.P),
        'temperature': np.array(hd9.T_fluid),
        'wall_temp': np.array(hd9.T_outer_wall),
    },
    'woodfield': {
        'time': np.array(hd_wood.time_array),
        'pressure': np.array(hd_wood.P),
        'temperature': np.array(hd_wood.T_fluid),
        'wall_temp': np.array(hd_wood.T_outer_wall),
    }
}

# Save to pickle file
with open('validation/hyddown_results.pkl', 'wb') as f:
    pickle.dump(results, f)

print("\nHydDown results saved to validation/hyddown_results.pkl")
