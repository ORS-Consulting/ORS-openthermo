#!/usr/bin/env python3
# OpenThermo vessel blowdown simulation
# Main script to run blowdown simulations from YAML input files

import yaml
import sys
import os
import argparse
import numpy as np

try:
    from openthermo.vessel.blowdown import Blowdown
except ImportError:
    # If openthermo is not installed, try to import from parent directory
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from openthermo.vessel.blowdown import Blowdown


def run_blowdown(input_file, plot=False, disable_pbar=False):
    """
    Run a blowdown simulation from a YAML input file.

    Parameters
    ----------
    input_file : str
        Path to the YAML input file
    plot : bool, optional
        Whether to plot results (default: False)
    disable_pbar : bool, optional
        Whether to disable progress bar (default: False)

    Returns
    -------
    segment : Blowdown
        The blowdown simulation object with results
    """
    # Load input from YAML file
    with open(input_file, mode="r") as infile:
        config = yaml.load(infile, Loader=yaml.FullLoader)

    # Extract solver method (default to depressurize)
    solver_method = config.pop("solver_method", "depressurize")

    # Extract validation data if present (not used for simulation, just for reference)
    validation_data = config.pop("validation", None)

    # Extract metadata if present
    metadata = config.pop("metadata", None)

    # Create Blowdown instance
    print(f"Running blowdown simulation from: {input_file}")
    if metadata:
        print(f"Description: {metadata.get('description', 'N/A')}")

    segment = Blowdown(config)

    # Run simulation with specified solver method
    if solver_method == "depressurize":
        segment.depressurize()
    elif solver_method == "depressurize_euler":
        segment.depressurize_euler()
    else:
        raise ValueError(f"Unknown solver_method: {solver_method}. Use 'depressurize' or 'depressurize_euler'")

    # Print summary results
    print("\n" + "=" * 70)
    print("SIMULATION RESULTS")
    print("=" * 70)
    print(f"Solver method: {solver_method}")
    print(f"Initial pressure: {segment.pressure[0]/1e5:.2f} bar")
    print(f"Final pressure:   {segment.pressure[-1]/1e5:.2f} bar")
    print(f"Initial temperature: {segment.temperature[0]:.2f} K")
    print(f"Final temperature:   {segment.temperature[-1]:.2f} K")
    print(f"Simulation time: {segment.times[-1]:.2f} s")
    print("=" * 70)

    # Compare with validation data if present
    if validation_data:
        print("\nVALIDATION COMPARISON")
        print("-" * 70)

        if "pressure" in validation_data:
            pres_val = validation_data["pressure"]
            # Check for nested structure (experimental/hyddown) or flat structure
            pres_data = None
            pres_label = "Experimental"
            if "experimental" in pres_val:
                pres_data = pres_val["experimental"]
            elif "pres" in pres_val:
                pres_data = pres_val

            if pres_data and "pres" in pres_data and len(pres_data["pres"]) > 0:
                final_pres_exp = pres_data["pres"][-1]
                final_pres_sim = segment.pressure[-1] / 1e5
                print(f"Final pressure:")
                print(f"  Simulated:    {final_pres_sim:6.3f} bar")
                print(f"  {pres_label}: {final_pres_exp:6.3f} bar")
                print(f"  Difference:   {abs(final_pres_sim - final_pres_exp):6.3f} bar")

        if "temperature" in validation_data:
            temp_val = validation_data["temperature"]

            # Check if we have separate gas/liquid temperatures (multiphase euler solver)
            has_separate_phases = hasattr(segment, 'gas_temperature') and segment.gas_temperature is not None

            # === GAS/FLUID TEMPERATURE COMPARISON ===
            gas_data = None
            gas_label = "Experimental"
            if "fluid_mean" in temp_val and "temp" in temp_val["fluid_mean"]:
                # Single-phase fluid temperature (HYSYS isentropic cases)
                gas_data = temp_val["fluid_mean"]
                gas_label = "HYSYS"
            elif "gas_mean" in temp_val and "temp" in temp_val["gas_mean"]:
                gas_data = temp_val["gas_mean"]
            elif "gas_experimental" in temp_val and "temp" in temp_val["gas_experimental"]:
                gas_data = temp_val["gas_experimental"]
            elif "gas_high_experimental" in temp_val and "temp" in temp_val["gas_high_experimental"]:
                gas_data = temp_val["gas_high_experimental"]
                gas_label = "Experimental (high)"
            elif "gas_high" in temp_val and "temp" in temp_val["gas_high"]:
                gas_data = temp_val["gas_high"]
                gas_label = "Experimental (high)"

            if gas_data and len(gas_data["temp"]) > 0:
                final_temp_exp = gas_data["temp"][-1]
                # Use gas_temperature if available (euler), otherwise use general temperature
                if has_separate_phases:
                    final_temp_sim = segment.gas_temperature[-1]
                else:
                    final_temp_sim = segment.temperature[-1]
                print(f"\nFinal gas temperature:")
                print(f"  Simulated:    {final_temp_sim:6.2f} K")
                print(f"  {gas_label}: {final_temp_exp:6.2f} K")
                print(f"  Difference:   {abs(final_temp_sim - final_temp_exp):6.2f} K")

            # === LIQUID TEMPERATURE COMPARISON ===
            # Only if we have separate phases and liquid validation data
            if has_separate_phases and hasattr(segment, 'liquid_temperature') and segment.liquid_temperature is not None:
                liquid_data = None
                liquid_label = "Experimental"
                if "liquid_mean" in temp_val and "temp" in temp_val["liquid_mean"]:
                    liquid_data = temp_val["liquid_mean"]
                elif "liquid_high" in temp_val and "temp" in temp_val["liquid_high"]:
                    liquid_data = temp_val["liquid_high"]
                    liquid_label = "Experimental (high)"

                if liquid_data and len(liquid_data["temp"]) > 0:
                    final_liq_exp = liquid_data["temp"][-1]
                    final_liq_sim = segment.liquid_temperature[-1]
                    print(f"\nFinal liquid temperature:")
                    print(f"  Simulated:    {final_liq_sim:6.2f} K")
                    print(f"  {liquid_label}: {final_liq_exp:6.2f} K")
                    print(f"  Difference:   {abs(final_liq_sim - final_liq_exp):6.2f} K")

            # === WALL TEMPERATURE COMPARISON ===
            # For multiphase cases: wall_high = unwetted (gas), wall_low = wetted (liquid)
            # For single-phase gas: wall_high and wall_low are both unwetted at different positions
            # For Unisim cases: wall_unwetted_inner/outer and wall_wetted_inner/outer

            # Check unwetted wall (always present for gas phase)
            wall_data = None
            wall_label = "Experimental"
            if "wall_mean" in temp_val and "temp" in temp_val["wall_mean"]:
                wall_data = temp_val["wall_mean"]
            elif "wall_experimental" in temp_val and "temp" in temp_val["wall_experimental"]:
                wall_data = temp_val["wall_experimental"]
            elif "wall_unwetted_inner" in temp_val and "temp" in temp_val["wall_unwetted_inner"]:
                wall_data = temp_val["wall_unwetted_inner"]
                wall_label = "Unisim (unwetted inner)"
            elif "wall_high" in temp_val and "temp" in temp_val["wall_high"]:
                wall_data = temp_val["wall_high"]
                wall_label = "Experimental (high)"
            elif "wall_outer" in temp_val and "temp" in temp_val["wall_outer"]:
                wall_data = temp_val["wall_outer"]
                wall_label = "Unisim (outer)"

            if wall_data and len(wall_data["temp"]) > 0:
                final_wall_exp = wall_data["temp"][-1]
                final_wall_sim = segment.unwetted_wall_temp[-1]
                print(f"\nFinal unwetted wall temperature:")
                print(f"  Simulated:    {final_wall_sim:6.2f} K")
                print(f"  {wall_label}: {final_wall_exp:6.2f} K")
                print(f"  Difference:   {abs(final_wall_sim - final_wall_exp):6.2f} K")

            # Check wetted wall (for multiphase cases)
            if has_separate_phases and hasattr(segment, 'wetted_wall_temp') and segment.wetted_wall_temp is not None:
                wetted_data = None
                wetted_label = "Experimental"

                if "wall_wetted_inner" in temp_val and "temp" in temp_val["wall_wetted_inner"]:
                    wetted_data = temp_val["wall_wetted_inner"]
                    wetted_label = "Unisim (wetted inner)"
                elif "wall_low" in temp_val and "temp" in temp_val["wall_low"]:
                    wetted_data = temp_val["wall_low"]
                    wetted_label = "Experimental (wetted/liquid)"

                if wetted_data and len(wetted_data["temp"]) > 0:
                    final_wetted_exp = wetted_data["temp"][-1]
                    final_wetted_sim = segment.wetted_wall_temp[-1]
                    print(f"\nFinal wetted wall temperature:")
                    print(f"  Simulated:    {final_wetted_sim:6.2f} K")
                    print(f"  {wetted_label}: {final_wetted_exp:6.2f} K")
                    print(f"  Difference:   {abs(final_wetted_sim - final_wetted_exp):6.2f} K")

            # For single-phase gas with wall_low (unwetted at different position)
            elif "wall_low" in temp_val and "temp" in temp_val["wall_low"]:
                wall_low_data = temp_val["wall_low"]
                final_wall_low_exp = wall_low_data["temp"][-1]
                final_wall_low_sim = segment.unwetted_wall_temp[-1]
                print(f"\nFinal wall temperature (low/unwetted):")
                print(f"  Simulated:    {final_wall_low_sim:6.2f} K")
                print(f"  Experimental (low): {final_wall_low_exp:6.2f} K")
                print(f"  Difference:   {abs(final_wall_low_sim - final_wall_low_exp):6.2f} K")

        print("-" * 70)

    # Plot if requested
    if plot:
        try:
            from matplotlib import pyplot as plt

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

            # Pressure plot
            ax1.plot(segment.times, np.asarray(segment.pressure) / 1e5, "-", label="Simulated")
            if validation_data and "pressure" in validation_data:
                pres_val = validation_data["pressure"]
                # Check for nested structure (experimental/hyddown) or flat structure
                if "experimental" in pres_val:
                    if "time" in pres_val["experimental"] and "pres" in pres_val["experimental"]:
                        ax1.plot(pres_val["experimental"]["time"], pres_val["experimental"]["pres"],
                                "x", label="Experimental")
                    if "hyddown" in pres_val and "time" in pres_val["hyddown"] and "pres" in pres_val["hyddown"]:
                        ax1.plot(pres_val["hyddown"]["time"], pres_val["hyddown"]["pres"],
                                "o", fillstyle="none", label="HydDown")
                elif "time" in pres_val and "pres" in pres_val:
                    # Infer label based on temperature data structure
                    # HYSYS cases use "fluid_mean", experimental cases use "gas_high/gas_low"
                    pres_label = "HYSYS"
                    if "temperature" in validation_data:
                        temp_val = validation_data["temperature"]
                        if "gas_high" in temp_val or "gas_low" in temp_val:
                            pres_label = "Experimental"
                    ax1.plot(pres_val["time"], pres_val["pres"], "x", label=pres_label)
            ax1.set_xlabel("Time (s)")
            ax1.set_ylabel("Pressure (bar)")
            ax1.legend()
            ax1.grid(True)

            # Temperature plot
            # Check if we have separate gas/liquid temperatures (euler solver) or single fluid temp
            if hasattr(segment, 'gas_temperature') and segment.gas_temperature is not None:
                # depressurize_euler: separate gas and liquid temperatures
                ax2.plot(segment.times, segment.gas_temperature, "-", label="Simulated Gas")
                if hasattr(segment, 'liquid_temperature') and segment.liquid_temperature is not None:
                    ax2.plot(segment.times, segment.liquid_temperature, "-", label="Simulated Liquid")
            else:
                # depressurize: single fluid temperature
                ax2.plot(segment.times, segment.temperature, "-", label="Simulated Fluid")

            ax2.plot(segment.times, segment.unwetted_wall_temp, "-", label="Simulated Unwetted Wall")
            # Plot wetted wall if available
            if hasattr(segment, 'wetted_wall_temp') and segment.wetted_wall_temp is not None:
                ax2.plot(segment.times, segment.wetted_wall_temp, "--", label="Simulated Wetted Wall")

            if validation_data and "temperature" in validation_data:
                temp_val = validation_data["temperature"]

                # Plot fluid/gas temperature data
                if "fluid_mean" in temp_val:
                    # Single-phase fluid temperature (e.g., HYSYS isentropic cases)
                    if "time" in temp_val["fluid_mean"] and "temp" in temp_val["fluid_mean"]:
                        ax2.plot(temp_val["fluid_mean"]["time"], temp_val["fluid_mean"]["temp"],
                                "x", label="HYSYS Fluid")
                elif "gas_mean" in temp_val:
                    if "time" in temp_val["gas_mean"] and "temp" in temp_val["gas_mean"]:
                        ax2.plot(temp_val["gas_mean"]["time"], temp_val["gas_mean"]["temp"],
                                "x", label="Experimental Gas (mean)")
                elif "gas_experimental" in temp_val:
                    # Byrnes-style nested structure
                    if "time" in temp_val["gas_experimental"] and "temp" in temp_val["gas_experimental"]:
                        ax2.plot(temp_val["gas_experimental"]["time"], temp_val["gas_experimental"]["temp"],
                                "x", label="Experimental Gas")
                    if "gas_hyddown" in temp_val and "time" in temp_val["gas_hyddown"] and "temp" in temp_val["gas_hyddown"]:
                        ax2.plot(temp_val["gas_hyddown"]["time"], temp_val["gas_hyddown"]["temp"],
                                "o", fillstyle="none", label="HydDown Gas")
                else:
                    # Check for _experimental suffix first (Woodfield naming)
                    if "gas_high_experimental" in temp_val:
                        if "time" in temp_val["gas_high_experimental"] and "temp" in temp_val["gas_high_experimental"]:
                            ax2.plot(temp_val["gas_high_experimental"]["time"], temp_val["gas_high_experimental"]["temp"],
                                    "x", label="Experimental Gas (high)")
                    elif "gas_high" in temp_val:
                        if "time" in temp_val["gas_high"] and "temp" in temp_val["gas_high"]:
                            ax2.plot(temp_val["gas_high"]["time"], temp_val["gas_high"]["temp"],
                                    "x", label="Experimental Gas (high)")

                    if "gas_low_experimental" in temp_val:
                        if "time" in temp_val["gas_low_experimental"] and "temp" in temp_val["gas_low_experimental"]:
                            ax2.plot(temp_val["gas_low_experimental"]["time"], temp_val["gas_low_experimental"]["temp"],
                                    "o", fillstyle="none", label="Experimental Gas (low)")
                    elif "gas_low" in temp_val:
                        if "time" in temp_val["gas_low"] and "temp" in temp_val["gas_low"]:
                            ax2.plot(temp_val["gas_low"]["time"], temp_val["gas_low"]["temp"],
                                    "o", fillstyle="none", label="Experimental Gas (low)")

                    # Also plot HydDown gas temperature if available
                    if "gas_hyddown" in temp_val and "time" in temp_val["gas_hyddown"] and "temp" in temp_val["gas_hyddown"]:
                        ax2.plot(temp_val["gas_hyddown"]["time"], temp_val["gas_hyddown"]["temp"],
                                "s", fillstyle="none", label="HydDown Gas")

                # Plot liquid temperature data (for multiphase cases)
                if "liquid_high" in temp_val:
                    if "time" in temp_val["liquid_high"] and "temp" in temp_val["liquid_high"]:
                        ax2.plot(temp_val["liquid_high"]["time"], temp_val["liquid_high"]["temp"],
                                "^", label="Experimental Liquid (high)")
                if "liquid_low" in temp_val:
                    if "time" in temp_val["liquid_low"] and "temp" in temp_val["liquid_low"]:
                        ax2.plot(temp_val["liquid_low"]["time"], temp_val["liquid_low"]["temp"],
                                "v", fillstyle="none", label="Experimental Liquid (low)")

                # Plot wall temperature data
                # Determine if this is a multiphase case
                is_multiphase = hasattr(segment, 'gas_temperature') and segment.gas_temperature is not None

                if "wall_mean" in temp_val:
                    if "time" in temp_val["wall_mean"] and "temp" in temp_val["wall_mean"]:
                        ax2.plot(temp_val["wall_mean"]["time"], temp_val["wall_mean"]["temp"],
                                "+", label="Experimental Wall (mean)")
                elif "wall_experimental" in temp_val:
                    # Byrnes-style nested structure
                    if "time" in temp_val["wall_experimental"] and "temp" in temp_val["wall_experimental"]:
                        ax2.plot(temp_val["wall_experimental"]["time"], temp_val["wall_experimental"]["temp"],
                                "+", label="Experimental Wall")
                    if "wall_hyddown" in temp_val and "time" in temp_val["wall_hyddown"] and "temp" in temp_val["wall_hyddown"]:
                        ax2.plot(temp_val["wall_hyddown"]["time"], temp_val["wall_hyddown"]["temp"],
                                "s", fillstyle="none", label="HydDown Wall")
                else:
                    # Plot unwetted wall temperatures (gas zone)
                    if "wall_unwetted_inner" in temp_val:
                        if "time" in temp_val["wall_unwetted_inner"] and "temp" in temp_val["wall_unwetted_inner"]:
                            ax2.plot(temp_val["wall_unwetted_inner"]["time"], temp_val["wall_unwetted_inner"]["temp"],
                                    "+", label="Unisim Wall (unwetted inner)")
                    elif "wall_high" in temp_val:
                        if "time" in temp_val["wall_high"] and "temp" in temp_val["wall_high"]:
                            ax2.plot(temp_val["wall_high"]["time"], temp_val["wall_high"]["temp"],
                                    "+", label="Experimental Wall (high)")

                    if "wall_unwetted_outer" in temp_val:
                        if "time" in temp_val["wall_unwetted_outer"] and "temp" in temp_val["wall_unwetted_outer"]:
                            ax2.plot(temp_val["wall_unwetted_outer"]["time"], temp_val["wall_unwetted_outer"]["temp"],
                                    "x", label="Unisim Wall (unwetted outer)")

                    # Plot wetted wall temperatures (liquid zone) for multiphase
                    if "wall_wetted_inner" in temp_val:
                        if "time" in temp_val["wall_wetted_inner"] and "temp" in temp_val["wall_wetted_inner"]:
                            ax2.plot(temp_val["wall_wetted_inner"]["time"], temp_val["wall_wetted_inner"]["temp"],
                                    "s", fillstyle="none", label="Unisim Wall (wetted inner)")
                    elif "wall_low" in temp_val:
                        if "time" in temp_val["wall_low"] and "temp" in temp_val["wall_low"]:
                            # For multiphase: wall_low is wetted (liquid side)
                            # For single-phase: wall_low is unwetted at low position
                            wall_low_label = "Experimental Wall (wetted)" if is_multiphase else "Experimental Wall (low)"
                            ax2.plot(temp_val["wall_low"]["time"], temp_val["wall_low"]["temp"],
                                    "s", fillstyle="none", label=wall_low_label)

                    if "wall_wetted_outer" in temp_val:
                        if "time" in temp_val["wall_wetted_outer"] and "temp" in temp_val["wall_wetted_outer"]:
                            ax2.plot(temp_val["wall_wetted_outer"]["time"], temp_val["wall_wetted_outer"]["temp"],
                                    "o", fillstyle="none", label="Unisim Wall (wetted outer)")
            ax2.set_xlabel("Time (s)")
            ax2.set_ylabel("Temperature (K)")
            ax2.legend()
            ax2.grid(True)

            plt.tight_layout()
            plt.show()

        except ImportError:
            print("\nWarning: matplotlib not available, skipping plot")

    return segment


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run OpenThermo blowdown simulation from YAML input file"
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default="input.yml",
        help="Path to YAML input file (default: input.yml)"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Plot results after simulation"
    )
    parser.add_argument(
        "--no-pbar",
        action="store_true",
        help="Disable progress bar"
    )

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found")
        sys.exit(1)

    try:
        segment = run_blowdown(args.input_file, plot=args.plot, disable_pbar=args.no_pbar)
    except Exception as e:
        print(f"\nError running simulation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
