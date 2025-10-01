import os
import numpy as np
import pytest
from scipy.constants import atm
from openthermo.vessel.blowdown import Blowdown
from openthermo.flash.michelsen import get_flash_dry

validation_path = os.path.join(
    os.path.abspath(os.path.dirname(__file__)), "..", "validation"
)


def test_blowdown_api_dry_inadequate_costald(plot=False):
    """
    Test against HYSYS Depressurisation utility for API 521 pool fire
    exposed the wetted area of the vessel. No water present.
    Inadeequate drainage and fire fighting available. COSTAÆD liquid
    density is applied.
    """
    file_name = "Water_dry_API_inadequate_costald_history.csv"

    path = os.path.join(validation_path, file_name)

    data = np.loadtxt(path, skiprows=11, delimiter=",", usecols=(range(25)))

    input = {}
    P = 12e5
    T = 298.15

    input["mode"] = "fire"
    input["drain_fire_fighting"] = "Inadequate"
    input["eos_model"] = "PR"
    input["liquid_density"] = "costald"
    input["max_time"] = 900
    input["delay"] = 0
    input["length"] = 10
    input["diameter"] = 3
    input["vessel_type"] = "Flat-end"
    input["orientation"] = "horizontal"
    input["liquid_level"] = 1.5
    input["water_level"] = 0.0
    input["operating_temperature"] = T
    input["operating_pressure"] = P
    input["ambient_temperature"] = 273
    input["back_pressure"] = 1.01e5
    input["bdv_orifice_size"] = 0.03  # m
    input["bdv_orifice_cd"] = 0.84

    input["leak_active"] = 0
    input["leak_size"] = 0.01  # m
    input["leak_cd"] = 0.65
    input["leak_type"] = "liquid"

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
    r = segment.depressurize()
    assert segment.pressure[-1] == pytest.approx(data[:, 2][-1] * 1e5 + atm, rel=0.03)
    assert segment.temperature[-1] == pytest.approx(data[:, 1][-1] + 273.15, rel=0.01)
    assert segment.mdot[-1] == pytest.approx(data[:, 3][-1] / -3600, rel=0.03)

    import matplotlib.pyplot as plt

    name = "plots\\API521_inadequate_costald_water_dry"

    if plot:
        plt.figure(1)
        plt.plot(
            segment.times,
            np.asarray(segment.pressure) / 1e5,
            "bo",
            label="openthermo/python",
        )
        plt.plot(data[:, 0], data[:, 2] + 1.013, "r-", label="HYSYS Depressurisation")
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (bar)")
        plt.savefig(name + "_pressure.png", dpi=300)

        plt.figure(2)
        plt.plot(
            segment.times,
            np.asarray(segment.temperature) - 273.15,
            "bo",
            label="openthermo/python",
        )
        plt.plot(data[:, 0], data[:, 1], "r-", label="HYSYS Depressurisation")
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel(r"Temperature ($^\circ$C)")
        plt.savefig(name + "_temperature.png", dpi=300)

        plt.figure(3)
        plt.plot(
            segment.times,
            np.asarray(segment.mdot) * -3600,
            "bo",
            label="openthermo/python",
        )
        plt.plot(data[:, 0], data[:, 3], "r-", label="HYSYS Depressurisation")
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel("Mass flow (kg/h)")
        plt.savefig(name + "_mdot.png", dpi=300)
        plt.show()


def test_blowdown_condensable_gas(plot=False):
    """
    Test against experiment by Haque et al., data revealed by Szczepanski and
    sourced from Wong.


    Haque, M. A.,  Richardson, S. M.,  and Saville,  G., Chamberlain,  G.,  and Shirvill, L.,
    Blowdown  of Pressure Vessels. II.  Experimental  Validation  of  Computer Model  and
    Case Studies, Trans IChemE  Part B:  Process Safety Environmental  Protection.,  70
    (BI),  10 -  17 (1992b

    Szczepanski, R., 1994. Simulation programs for blowdown of
    pressure vessels. Proceedings of the IChemE SONG Meeting 12.

    Wong, M.A., 1998. Development of a Mathematical Model for
    Blowdown of Vessels Containing Multi-Component
    Hydrocarbon Mixtures. University of London.
    """
    file_name = "condensable_gas_pressure.txt"
    pres = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "condensable_gas_inner_wall_lower.txt"
    inner_wall_low = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    file_name = "condensable_gas_inner_wall_higher.txt"
    inner_wall_high = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    t = np.linspace(0, 1500, 20)
    iwh = np.interp(t, inner_wall_high[:, 0], inner_wall_high[:, 1])
    iwl = np.interp(t, inner_wall_low[:, 0], inner_wall_low[:, 1])

    file_name = "condensable_gas_liquid_inner_wall_lower.txt"
    inner_wall_low = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    file_name = "condensable_gas_liquid_inner_wall_higher.txt"
    inner_wall_high = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    t2 = np.linspace(min(inner_wall_low[:, 0]), max(inner_wall_low[:, 0]), 20)
    liwh = np.interp(t2, inner_wall_high[:, 0], inner_wall_high[:, 1], left=None)
    liwl = np.interp(t2, inner_wall_low[:, 0], inner_wall_low[:, 1], left=None)

    input = {}
    P = 116 * atm
    T = 293.0

    input["mode"] = "isentropic"
    input["heat_transfer"] = "rigorous"
    input["wall_thickness"] = 0.059  # m
    input["eos_model"] = "PR"
    input["liquid_density"] = "eos"
    input["max_time"] = 1500
    input["delay"] = 0

    input["length"] = 2.25
    # input["diameter"] = 1.130
    input["diameter"] = 1.13
    input["vessel_type"] = "ASME F&D"
    # input["length"] = 2.75
    # # input["diameter"] = 1.130
    # input["diameter"] = 1.09
    # input["vessel_type"] = "ASME F&D"
    input["orientation"] = "horizontal"
    input["liquid_level"] = 0.0
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

    names = ["water", "methane", "ethane", "propane", "n-butane"]
    molefracs = [0.0, 0.64, 0.06, 0.28, 0.02]

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
    r = segment.depressurize()

    if plot:
        from matplotlib import pyplot as plt
        import scienceplots

        plt.style.use(["science", "nature", "scatter"])

        plt.figure(1)
        # plt.plot(t2, liwl, "-")
        # plt.plot(t2, liwh, "-")
        plt.fill_between(t2, liwh, liwl, alpha=0.2, label="Exp. wetted")
        # plt.plot(t, iwl, "-")
        # plt.plot(t, iwh, "-")
        plt.fill_between(t, iwh, iwl, alpha=0.2, label="Exp. unwetted")
        plt.plot(segment.times, segment.wetted_wall_temp, label="Model wetted")
        plt.plot(segment.times, segment.unwetted_wall_temp, label="Model unwetted")

        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel(r"Inner wall temperature (K)")
        plt.savefig("plots\condensable_gas_inner_wal.png", dpi=300)

        plt.figure(2)
        plt.plot(segment.times, np.asarray(segment.pressure) / 1e5, label="Model")
        plt.plot(pres[:, 0], pres[:, 1] * 1.013, "-", label="Experiment")
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (bar)")
        plt.savefig("plots\condensable_gas_pressure.png", dpi=300)
        plt.show()


def test_blowdown_nitrogen(plot=False):
    import yaml

    file_name = "N2_I1.yml"
    input_file = os.path.join(validation_path, file_name)

    with open(input_file) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)

    P = 150e5
    T = 289.0

    input["mode"] = "isentropic"
    input["heat_transfer"] = "rigorous"
    input["wall_thickness"] = 0.025  # m
    input["eos_model"] = "PR"
    input["liquid_density"] = "eos"
    input["max_time"] = 100
    input["delay"] = 0
    input["length"] = 1.524
    # input["diameter"] = 1.130
    input["diameter"] = 0.273
    input["vessel_type"] = "Flat-end"
    input["orientation"] = "horizontal"
    input["liquid_level"] = 0.0
    input["water_level"] = 0.0
    input["operating_temperature"] = T
    input["operating_pressure"] = P
    input["ambient_temperature"] = 288
    input["back_pressure"] = 1.01e5
    input["bdv_orifice_size"] = 0.00635  # m
    input["bdv_orifice_cd"] = 0.8

    input["leak_active"] = 0
    input["leak_size"] = 0.01  # m

    input["leak_cd"] = 0.65
    input["leak_type"] = "liquid"

    names = ["water", "nitrogen"]
    molefracs = [0.0, 1.0]

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
    r = segment.depressurize()
    # segment.plot()

    if plot:
        from matplotlib import pyplot as plt
        import scienceplots

        plt.style.use(["science", "nature", "scatter"])

        t_gh = input["validation"]["temperature"]["gas_high"]["time"]
        gh = input["validation"]["temperature"]["gas_high"]["temp"]

        t_gl = input["validation"]["temperature"]["gas_low"]["time"]
        gl = input["validation"]["temperature"]["gas_low"]["temp"]

        t = np.linspace(0, 100, 20)
        gh = np.interp(t, t_gh, gh)
        gl = np.interp(t, t_gl, gl)

        t_wh = input["validation"]["temperature"]["wall_outer"]["time"]
        wh = input["validation"]["temperature"]["wall_outer"]["temp"]

        t_wl = input["validation"]["temperature"]["wall_inner"]["time"]
        wl = input["validation"]["temperature"]["wall_inner"]["temp"]
        t = np.linspace(0, 100, 20)
        wh = np.interp(t, t_wh, wh)
        wl = np.interp(t, t_wl, wl)

        plt.figure(1)
        plt.plot(t, gl, "-")
        plt.plot(t, gh, "-")

        plt.fill_between(t, gh, gl, alpha=0.2, label="Exp. gas")
        plt.plot(t, wl, "-")
        plt.plot(t, wh, "-")

        plt.fill_between(t, wh, wl, alpha=0.2, label="Exp. wall")
        plt.plot(segment.times, segment.unwetted_wall_temp, label="Model wall")
        plt.plot(segment.times, segment.temperature, label="Model gas")
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel(r"Temperature (K)")
        plt.savefig("plots\\nitrogen_wall.png", dpi=300)

        plt.figure(2)
        plt.plot(segment.times, np.asarray(segment.pressure) / 1e5, "-", label="Model")
        plt.plot(
            input["validation"]["pressure"]["time"],
            input["validation"]["pressure"]["pres"],
            label="Experiment",
        )
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (bar)")
        plt.savefig("plots\\nitorgen_pres.png", dpi=300)
        plt.show()


def test_blowdown_nitrogen_co2(plot=False):

    file_name = "co2_n2_gas_lower.txt"
    gas_low = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "co2_n2_gas_higher.txt"
    gas_high = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    t2 = np.linspace(min(gas_low[:, 0]), max(gas_low[:, 0]), 20)
    gh = np.interp(t2, gas_high[:, 0], gas_high[:, 1], left=None) + 273.15
    gl = np.interp(t2, gas_low[:, 0], gas_low[:, 1], left=None) + 273.15

    file_name = "co2_n2_liquid_lower.txt"
    liquid_low = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "co2_n2_liquid_higher.txt"
    liquid_high = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    t1 = np.linspace(min(liquid_low[:, 0]), max(liquid_low[:, 0]), 20)
    lh = np.interp(t1, liquid_high[:, 0], liquid_high[:, 1], left=None) + 273.15
    ll = np.interp(t1, liquid_low[:, 0], liquid_low[:, 1], left=None) + 273.15

    file_name = "co2_n2_HYSYS_gas.txt"
    HYSYS_gas = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "co2_n2_HYSYS_liquid.txt"
    HYSYS_liquid = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "co2_n2_HYSYS_gas_wall.txt"
    HYSYS_gas_wall = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    file_name = "co2_n2_HYSYS_liquid_wall.txt"
    HYSYS_liquid_wall = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    P = 150e5
    T = 290.0
    input = {}
    input["mode"] = "isentropic"
    input["heat_transfer"] = "rigorous"
    input["wall_thickness"] = 0.025  # m
    input["eos_model"] = "PR"
    input["liquid_density"] = "eos"
    input["max_time"] = 59
    input["delay"] = 0
    input["length"] = 1.4551999494616477
    # input["diameter"] = 1.130
    input["diameter"] = 0.273
    input["vessel_type"] = "DIN"
    input["orientation"] = "vertical"
    input["liquid_level"] = 0.0
    input["water_level"] = 0.0
    input["operating_temperature"] = T
    input["operating_pressure"] = P
    input["ambient_temperature"] = 293
    input["back_pressure"] = 1.01e5
    input["bdv_orifice_size"] = 0.00635  # m
    input["bdv_orifice_cd"] = 0.80

    input["leak_active"] = 0
    input["leak_size"] = 0.01  # m

    input["leak_cd"] = 0.65
    input["leak_type"] = "liquid"

    names = ["water", "nitrogen", "carbon dioxide"]
    molefracs = [0.0, 0.70, 0.30]

    names = ["nitrogen", "carbon dioxide"]
    molefracs = [0.70, 0.30]

    input["molefracs"] = molefracs

    input["flash"] = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=T,
        # rho=input["liquid_density"],
        # model=input["eos_model"],
    )

    segment = Blowdown(input)
    r = segment.depressurize_euler()
    # segment.plot()

    if plot:
        from matplotlib import pyplot as plt
        import scienceplots

        plt.style.use(["science", "nature", "scatter"])

        plt.fill_between(t1, lh, ll, alpha=0.2, label="Exp. liquid")
        plt.plot(segment.times, segment.liquid_temperature, label="Model liquid")
        plt.fill_between(t2, gh, gl, alpha=0.2, label="Exp. gas")
        plt.plot(segment.times, segment.gas_temperatrure, label="Model gas")

        plt.plot(segment.times, segment.wetted_wall_temp, label="Model wet wall")
        plt.plot(segment.times, segment.unwetted_wall_temp, label="Model dry wall")

        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel(r"Temperature (K)")
        plt.savefig("plots\\n2_co2_temp.png", dpi=300)

        plt.figure(2)
        plt.plot(segment.times, np.asarray(segment.pressure) / 1e5, "-", label="Model")
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (bar)")
        plt.savefig("plots\\n2_co2_pres.png", dpi=300)
        print(
            segment.wetted_wall_temp[-1] - 273.15,
            segment.unwetted_wall_temp[-1] - 273.15,
        )

        plt.figure(3)
        plt.plot(
            segment.times,
            segment.liquid_temperature,
            color="tab:blue",
            label="Model liquid",
        )
        plt.plot(
            HYSYS_liquid[:, 0],
            HYSYS_liquid[:, 1],
            "--",
            color="tab:blue",
            label="HYSYS V9 Blowdown liquid",
        )
        plt.plot(
            segment.times,
            segment.gas_temperatrure,
            label="Model gas",
            color="tab:green",
        )
        plt.plot(
            HYSYS_gas[:, 0],
            HYSYS_gas[:, 1],
            "--",
            color="tab:green",
            label="HYSYS V9 Blowdown gas",
        )
        plt.xlim((0, 60))
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel(r"Temperature (K)")
        plt.savefig("plots\\n2_co2_temp_HYSYS.png", dpi=300)

        plt.figure(4)
        plt.plot(
            segment.times,
            segment.wetted_wall_temp,
            color="tab:blue",
            label="Model wetted wall",
        )
        plt.plot(
            HYSYS_liquid_wall[:, 0],
            HYSYS_liquid_wall[:, 1],
            "--",
            color="tab:blue",
            label="HYSYS V9 Blowdown wetted wall",
        )
        plt.plot(
            segment.times,
            segment.unwetted_wall_temp,
            color="tab:green",
            label="Model dry wall",
        )
        plt.plot(
            HYSYS_gas_wall[:, 0],
            HYSYS_gas_wall[:, 1],
            "--",
            color="tab:green",
            label="HYSYS V9 Blowdown dry wall",
        )
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel(r"Temperature (K)")
        plt.xlim((0, 60))
        plt.savefig("plots\\n2_co2_wall_temp_HYSYS.png", dpi=300)
        plt.show()


def test_blowdown_non_condensable(plot=False):
    """
    Test against experiment by Haque et al., data sourced from Wong.


    Haque, M. A.,  Richardson, S. M.,  and Saville,  G., Chamberlain,  G.,  and Shirvill, L.,
    Blowdown  of Pressure Vessels. II.  Experimental  Validation  of  Computer Model  and
    Case Studies, Trans IChemE  Part B:  Process Safety Environmental  Protection.,  70
    (BI),  10 -  17 (1992b

    Wong, M.A., 1998. Development of a Mathematical Model for
    Blowdown of Vessels Containing Multi-Component
    Hydrocarbon Mixtures. University of London.
    """
    file_name = "non-condensable_pressure.txt"
    pres = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "non-condensable_gas_temp_higher.txt"
    gas_temp_high = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "non-condensable_gas_temp_lower.txt"
    gas_temp_low = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    t = np.linspace(0, 2000, 20)
    gth = np.interp(t, gas_temp_high[:, 0], gas_temp_high[:, 1])
    gtl = np.interp(t, gas_temp_low[:, 0], gas_temp_low[:, 1])

    file_name = "non-condensable_wall_temp.txt"
    wt = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    wt = np.interp(t, wt[:, 0], wt[:, 1], left=None)

    input = {}
    P = 120 * atm
    T = 303.0

    input["mode"] = "isentropic"
    input["heat_transfer"] = "rigorous"
    input["wall_thickness"] = 0.059  # m
    input["eos_model"] = "PR"
    input["liquid_density"] = "eos"
    input["max_time"] = 2000
    input["delay"] = 0
    input["length"] = 2.25
    # input["diameter"] = 1.130
    input["diameter"] = 1.13
    input["vessel_type"] = "ASME F&D"
    input["orientation"] = "horizontal"
    input["liquid_level"] = 0.0
    input["water_level"] = 0.0
    input["operating_temperature"] = T
    input["operating_pressure"] = P
    input["ambient_temperature"] = 293.0
    input["back_pressure"] = 1.01e5
    input["bdv_orifice_size"] = 0.00635  # m
    input["bdv_orifice_cd"] = 0.85

    input["leak_active"] = 0
    input["leak_size"] = 0.01  # m

    input["leak_cd"] = 0.65
    input["leak_type"] = "liquid"

    names = ["water", "methane", "ethane"]
    molefracs = [0.0, 0.91, 0.09]
    names = ["methane", "ethane"]
    molefracs = [0.91, 0.09]

    input["molefracs"] = molefracs

    input["flash"] = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=T,
        # rho=input["liquid_density"],
        # model=input["eos_model"],
    )

    segment = Blowdown(input)
    r = segment.depressurize()
    if plot:
        from matplotlib import pyplot as plt
        import scienceplots

        plt.style.use(["science", "nature", "scatter"])

        plt.figure(1)
        plt.plot(t, gtl, "-")
        plt.plot(t, gth, "-")
        plt.fill_between(t, gtl, gth, alpha=0.2, label="Exp. gas temp.")
        plt.plot(segment.times, segment.temperature, "-", label="Model gas temp.")
        plt.plot(t, wt, label="Exp. wall temp.")
        plt.plot(
            segment.times, segment.unwetted_wall_temp, "-", label="Model wall temp."
        )

        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel(r"Temperature (K)")
        plt.savefig("plots\\non_condensable_gas_temp.png", dpi=300)

        plt.figure(2)
        plt.plot(segment.times, np.asarray(segment.pressure) / 1e5, "-", label="Model")
        plt.plot(pres[:, 0], pres[:, 1] * 1.013, label="Experiment")
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (bar)")
        plt.savefig("plots\\non-condensable_gas_pressure.png", dpi=300)
        plt.show()


def test_blowdown_condensable_gas_rig(plot=False):
    """
    Test against experiment by Haque et al., data revealed by Szczepanski and
    sourced from Wong.


    Haque, M. A.,  Richardson, S. M.,  and Saville,  G., Chamberlain,  G.,  and Shirvill, L.,
    Blowdown  of Pressure Vessels. II.  Experimental  Validation  of  Computer Model  and
    Case Studies, Trans IChemE  Part B:  Process Safety Environmental  Protection.,  70
    (BI),  10 -  17 (1992b

    Szczepanski, R., 1994. Simulation programs for blowdown of
    pressure vessels. Proceedings of the IChemE SONG Meeting 12.

    Wong, M.A., 1998. Development of a Mathematical Model for
    Blowdown of Vessels Containing Multi-Component
    Hydrocarbon Mixtures. University of London.
    """
    file_name = "condensable_gas_pressure.txt"
    pres = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "condensable_gas_inner_wall_lower.txt"
    inner_wall_low = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    file_name = "condensable_gas_inner_wall_higher.txt"
    inner_wall_high = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    t = np.linspace(0, 1500, 20)
    iwh = np.interp(t, inner_wall_high[:, 0], inner_wall_high[:, 1])
    iwl = np.interp(t, inner_wall_low[:, 0], inner_wall_low[:, 1])

    file_name = "condensable_gas_liquid_inner_wall_lower.txt"
    inner_wall_low = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    file_name = "condensable_gas_liquid_inner_wall_higher.txt"
    inner_wall_high = np.loadtxt(
        os.path.join(validation_path, file_name), delimiter="\t"
    )

    t2 = np.linspace(min(inner_wall_low[:, 0]), max(inner_wall_low[:, 0]), 20)
    liwh = np.interp(t2, inner_wall_high[:, 0], inner_wall_high[:, 1], left=None)
    liwl = np.interp(t2, inner_wall_low[:, 0], inner_wall_low[:, 1], left=None)

    file_name = "condensable_gas_gas_temp_higher.txt"
    gas_high = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "condensable_gas_gas_temp_lower.txt"
    gas_low = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    t3 = np.linspace(min(gas_low[:, 0]), max(gas_low[:, 0]), 20)
    gh = np.interp(t3, gas_high[:, 0], gas_high[:, 1], left=None)
    gl = np.interp(t3, gas_low[:, 0], gas_low[:, 1], left=None)

    file_name = "condensable_gas_liq_temp_higher.txt"
    liq_high = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    file_name = "condensable_gas_liq_temp_lower.txt"
    liq_low = np.loadtxt(os.path.join(validation_path, file_name), delimiter="\t")

    t4 = np.linspace(min(liq_low[:, 0]), max(liq_low[:, 0]), 20)
    lh = np.interp(t3, liq_high[:, 0], liq_high[:, 1], left=None)
    ll = np.interp(t3, liq_low[:, 0], liq_low[:, 1], left=None)

    input = {}
    P = 116 * atm
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
    # input["diameter"] = 1.130
    input["diameter"] = 1.13
    input["vessel_type"] = "ASME F&D"
    # input["length"] = 2.75
    # # input["diameter"] = 1.130
    # input["diameter"] = 1.09
    # input["vessel_type"] = "ASME F&D"
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

    # names = ["water", "methane", "ethane", "propane", "n-butane"]
    names = ["methane", "ethane", "propane", "n-butane"]
    # molefracs = [0.855, 0.0545, 0.10, 0.0]
    molefracs = [0.64, 0.06, 0.28, 0.02]

    input["molefracs"] = molefracs

    input["flash"] = get_flash_dry(
        names,
        molefracs,
        P=P,
        T=T,
        # rho=input["liquid_density"],
        # model=input["eos_model"],
    )
    import time

    time1 = time.time()
    segment = Blowdown(input)
    r = segment.depressurize_euler()  # _euler()
    time2 = time.time()
    print(f"Elapsed time {time2-time1} sec.")
    if True:
        from matplotlib import pyplot as plt
        import scienceplots

        plt.style.use(["science", "nature", "scatter"])
        import pandas as pd

        BD = pd.read_excel("..//validation//HYSYS_BLOWDOWN//BLOWDOWN.xlsx")

        plt.figure(1)
        # plt.plot(t2, liwl, "-")
        # plt.plot(t2, liwh, "-")
        plt.fill_between(t2, liwh, liwl, alpha=0.2, label="Exp. wetted")
        # plt.plot(t, iwl, "-")
        # plt.plot(t, iwh, "-")
        plt.fill_between(t, iwh, iwl, alpha=0.2, label="Exp. unwetted")
        plt.plot(segment.times, segment.wetted_wall_temp, label="Model wetted")
        plt.plot(segment.times, segment.unwetted_wall_temp, label="Model unwetted")

        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel(r"Inner wall temperature (K)")
        plt.savefig("plots\condensable_gas_inner_wall_rig.png", dpi=300)

        plt.figure(2)
        plt.plot(segment.times, np.asarray(segment.pressure) / 1e5, label="Model")
        plt.plot(pres[:, 0], pres[:, 1] * 1.013, "-", label="Experiment")
        # plt.plot(BD["Time"], BD["P (barg)"] + 1.013, label="HYSYS")
        plt.legend(loc="best")
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (bar)")
        plt.savefig("plots\condensable_gas_pressure_rig.png", dpi=300)

        plt.figure(3)
        plt.plot(
            segment.times,
            np.asarray(segment.gas_temperatrure),
            label="Model gas temperature",
        )
        plt.plot(
            segment.times,
            np.asarray(segment.liquid_temperature),
            label="Model liquid temperature",
        )
        # plt.plot(BD["Time"], BD["Gas wall"], label="HYSYS Gas")
        plt.fill_between(t3, gh, gl, alpha=0.2, label="Exp. gas temperature")
        plt.fill_between(t4, lh, ll, alpha=0.2, label="Exp. liquid temperature")
        plt.legend(loc="best")
        plt.xlabel("Time (sec)")
        plt.ylabel("Temperature (K)")
        plt.savefig("plots\condensable_gas_bulk_rig.png", dpi=300)
        plt.show()


if __name__ == "__main__":
    test_blowdown_condensable_gas(plot=True)
    # test_blowdown_condensable_gas_rig(plot=True)
    # test_blowdown_non_condensable(plot=True)
    # test_blowdown_api_dry_inadequate_costald(plot=True)
    # test_blowdown_nitrogen(plot=True)
    # test_blowdown_nitrogen_co2(plot=True)
