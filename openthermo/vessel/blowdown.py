from fluids.geometry import TANK
from chemicals import normalize
import numpy as np
from scipy.integrate import ode
from scipy.optimize import minimize
import math
from scipy.constants import g
from thermo.volume import COSTALD_mixture
from openthermo.properties.transport import h_inside, h_inside_wetted
from openthermo import errors
from openthermo.flash.michelsen import get_flash_dry
from openthermo.vessel.fire import sb_fire


def hem_release_rate():
    pass


def liquid_release_bernouilli(P1, P2, rho, CD, area, H):
    """
    Liquid mass flow (kg/s) trough a hole or orifice.
    flow conditions. The formula is based on Yellow Book equation 2.194.

    Methods for the calculation of physical effects, CPR 14E, van den Bosch and Weterings (Eds.), 1996

    Parameters
    ----------
    P1 : float
        Upstream pressure
    P2 : float
        Downstream pressure
    rho : float
        Fluid density
    CD : float
        Coefficient of discharge
    area : float
        Orifice area
    H : float
        Liquid static height

    Returns
    ----------
        : float
        liquid  release rate / mass flow of discharge
    """
    P1 = P1 + H * rho * g
    if H > 0 and P1 > P2:
        return CD * area * math.sqrt(2 * (P1 - P2) * rho)
    else:
        return 0


def two_phase_release_fauske(P1, Pc, rho, CD, area):
    """
    Two-phase mass flow (kg/s) trough a hole or orifice,.
    flow conditions. The formula is based on Yellow Book equation 2.91.

    Methods for the calculation of physical effects, CPR 14E, van den Bosch and Weterings (Eds.), 1996
    World Bank Technical Paper Number 55, Techniques for Assessing Industrial Hazards, 2nd print, 1990

    Parameters
    ----------
    P1 : float
        Upstream pressure
    Pc : float
        Critical pressure or ambient which ever is the highest. Pc = 0.55 P1
    rho : float
        Fluid density
    CD : float
        Coefficient of discharge
    are : float
        Orifice area

    Returns
    ----------
        : float
        Two-phase   release rate / mass flow of discharge
    """

    return max(CD * area * math.sqrt(2 * (P1 - Pc) * rho), 0)


def gas_release_rate(P1, P2, rho, k, CD, area):
    """
    Gas massflow (kg/s) trough a hole at critical (sonic) or subcritical
    flow conditions. The formula is based on Yellow Book equation 2.22.

    Methods for the calculation of physical effects, CPR 14E,
    van den Bosch and Weterings (Eds.), 1996

    Parameters
    ----------
    P1 : float
        Upstream pressure
    P2 : float
        Downstream pressure
    rho : float
        Fluid density
    k : float
        Ideal gas k (Cp/Cv)
    CD : float
        Coefficient of discharge
    are : float
        Orifice area

    Returns
    ----------
        : float
        Gas release rate / mass flow of discharge
    """
    if P1 > P2:
        if P1 / P2 > ((k + 1) / 2) ** ((k) / (k - 1)):
            flow_coef = 1
        else:
            flow_coef = (
                2
                / (k - 1)
                * (((k + 1) / 2) ** ((k + 1) / (k - 1)))
                * ((P2 / P1) ** (2 / k))
                * (1 - (P2 / P1) ** ((k - 1) / k))
            )

        retval = (
            math.sqrt(flow_coef)
            * CD
            * area
            * math.sqrt(rho * P1 * k * (2 / (k + 1)) ** ((k + 1) / (k - 1)))
        )
    else:
        retval = 0

    return retval


class Blowdown:
    def __init__(self, input):
        self.mode = input["mode"]
        if "heat_transfer" in input:
            self.ambient_temperature = input["ambient_temperature"]
            self.heat_transfer = input["heat_transfer"]
            self.wall_thickness = input["wall_thickness"]
            if "wall_heat capacity" in input:
                self.wall_heat_capacity = input["wall_heat capacity"]
            else:
                self.wall_heat_capacity = 477
            if "wall_density" in input:
                self.wall_density = input["wall_heat capacity"]
            else:
                self.wall_density = 7800
            if input["heat_transfer"] == "rigorous_sb_fire":
                self.sb_fire_type = input["sb_fire_type"]
        else:
            self.heat_transfer = None
        self.liq_density = input["liquid_density"]
        self.delay = input["delay"]
        self.leak_active = input["leak_active"]
        self.max_time = input["max_time"]
        self.shutdown = 0
        self.length = input["length"]
        self.diameter = input["diameter"]
        self.vessel_type = input["vessel_type"]
        self.vessel_orientation = input["orientation"]
        self.liquid_level = input["liquid_level"]
        self.water_level = input["water_level"]
        self.operating_temperature = input["operating_temperature"]
        self.ambient_temperature = input["ambient_temperature"]
        self.operating_pressure = input["operating_pressure"]
        self.back_pressure = input["back_pressure"]
        self.orifice_size = input["bdv_orifice_size"]
        self.orifice_cd = input["bdv_orifice_cd"]
        self.leak_type = input["leak_type"]
        self.leak_size = input["leak_size"]
        self.leak_cd = input["leak_cd"]
        self.orifice_area = self.orifice_size**2 / 4 * math.pi
        self.leak_area = self.leak_size**2 / 4 * math.pi
        if "time_step" in input:
            self.dt = input["time_step"]
        else:
            self.dt = 1

        if self.mode == "fire":
            self._setup_fire(input)

        self._setup_tank()
        if "flash" in input:
            self.flash = input["flash"]
            # self.zi = self.flash.zs
            self.zi = (np.array(input["molefracs"]) / sum(input["molefracs"])).tolist()

        self._adjust_composition()

    def _liq_density(self, phase):
        if self.liq_density == "eos":
            return phase.rho_mass()
        elif self.liq_density == "costald":
            Vm = COSTALD_mixture(phase.zs, phase.T, phase.Tcs, phase.Vcs, phase.omegas)
            rho = 1 / Vm * phase.MW() / 1000
            return rho

    def _density(self, phases):
        if self.liq_density == "eos":
            return phases.rho_mass()
        elif self.liq_density == "costald":
            Vmn = 0
            for phase in phases:
                Vmn += (
                    COSTALD_mixture(
                        phase.zs, phase.T, phase.Tcs, phase.Vcs, phase.omegas
                    )
                    * phase.beta
                )
            return 1 / Vmn * (phases.MW() / 1000)

    def _Vm_liq(self, phase):
        if self.liq_density == "eos":
            return phase.V()
        elif self.liq_density == "costald":
            return COSTALD_mixture(
                phase.zs, phase.T, phase.Tcs, phase.Vcs, phase.omegas
            )

    def _setup_fire(self, input):
        if "fire_type" in input:
            self.fire_type = input["fire_type"]
        else:
            self.fire_type = "API521"

        if "drain_fire_fighting" in input:
            self.drainage_and_fire_fighting = input["drain_fire_fighting"]
        else:
            self.drainage_and_fire_fighting = "Inadequate"

        if "exposed_area" in input:
            self.exposed_area_type = input["exposed_area"]
        else:
            self.exposed_area_type = "Wetted"

        if "manual_heat" in input:
            self.manual_heat = input["manual_heat"]
            if "manual_area" in input:
                self.user_defined_area = input["manual_area"]
            else:
                raise ValueError(
                    "Missing input for manual/user defined area for fire exposure."
                )
        if "environmental_factor" in input:
            self.F = input["environmental_factor"]
        else:
            self.F = 1

        self.C = None

    def _setup_tank(self):
        if self.vessel_orientation == "horizontal":
            horizontal = True
        else:
            horizontal = False

        if self.vessel_type == "Flat-end":
            self.vessel = TANK(D=self.diameter, L=self.length, horizontal=horizontal)
        elif self.vessel_type == "ASME F&D":
            self.vessel = TANK(
                D=self.diameter,
                L=self.length,
                sideA="torispherical",
                sideB="torispherical",
                horizontal=horizontal,
            )
        elif self.vessel_type == "DIN":
            self.vessel = TANK(
                D=self.diameter,
                L=self.length,
                sideA="torispherical",
                sideB="torispherical",
                sideA_f=1,
                sideA_k=0.1,
                sideB_f=1,
                sideB_k=0.1,
                horizontal=horizontal,
            )

    def _get_heat_load_W(self):
        if self.fire_type == "API521":
            if self.drainage_and_fire_fighting == "Adequate":
                self.C = 43200
                self.n = 0.82
            elif self.drainage_and_fire_fighting == "Inadequate":
                self.C = 70900
                self.n = 0.82
        elif self.fire_type == "API521_CONFINED":
            self.C = 108900
            self.n = 1.0
        else:
            raise errors.InputError("Unsupported heat load")

        exposed_area = self.exposed_area()
        return self.C * self.F * exposed_area**self.n  # J/s

    def exposed_area(self):
        if self.exposed_area_type == "Manual Input":
            self.exposed_area_m2 = self.user_defined_area
        elif self.exposed_area_type == "Wetted":
            self.exposed_area_m2 = self.wetted_area()
        elif self.exposed_area_type == "Total":
            self.exposed_area_m2 = self.total_area()
        else:
            raise errors.InputError("Area calculation not supported yet!")

        return self.exposed_area_m2

    def wetted_area(self, level=None):
        if level is None:
            level = self.liquid_level
        return self.vessel.SA_from_h(level)

    def total_area(self):
        return self.vessel.A

    def _adjust_composition(self):
        # adjust the gas / liquid / water phase composition to give the user
        # supplied liquid levels
        res = self.flash.flash(
            T=self.operating_temperature, P=self.operating_pressure, zs=self.zi
        )

        for i in range(4):
            # print(res.betas)
            liq_rho = self._liq_density(res.liquid0)
            water_rho = self._liq_density(res.liquid1)
            V_gas = self.vessel.V_total - self.vessel.V_from_h(self.liquid_level)
            N_gas = V_gas * res.gas.rho_mass() / (res.gas.MW() / 1000)
            V_liquid = self.vessel.V_from_h(self.liquid_level) - self.vessel.V_from_h(
                self.water_level
            )
            N_liquid = V_liquid * liq_rho / (res.lightest_liquid.MW() / 1000)
            V_water = self.vessel.V_from_h(self.water_level)
            N_water = V_water * water_rho / (res.heaviest_liquid.MW() / 1000)
            V_tot = V_gas + V_liquid + V_water
            N_tot = N_gas + N_liquid + N_water
            zs = (
                (
                    N_gas * np.array(res.gas.zs)
                    + N_liquid * np.array(res.lightest_liquid.zs)
                    + N_water * np.array(res.heaviest_liquid.zs)
                )
                / N_tot
            ).tolist()
            res = self.flash.flash(
                T=self.operating_temperature, P=self.operating_pressure, zs=zs
            )

        self.m0 = V_gas * res.gas.rho_mass() + V_liquid * liq_rho + V_water * water_rho
        self.z_adjust = zs

    def _blowdown_gov_eqns(
        self, t, y
    ):  # , times):#, sol, mdot, pressure, temperature, gas_mass, molefracs, enthalpy):
        self.stamp.append(hash((t, (yi for yi in y))))
        ##############################################################################
        # The thermodynamic part of the code, all at constant volume
        ##############################################################################
        if self.mode == "isothermal":
            m, N = y[0], y[1]
            z = normalize(y[2:] / N)
            V_molar = self.vessel.V_total / N

            if len(self.pressure) > 0:
                Pguess = self.pressure[-1]
            else:
                Pguess = self.operating_pressure
            res = self.flash.flash(
                T=self.operating_temperature, V=V_molar, zs=z, Pguess=Pguess
            )

        elif (
            self.heat_transfer == "rigorous" or self.heat_transfer == "rigorous_sb_fire"
        ):
            m, N, U = y[0], y[1], y[2]
            Tuw, Tw = y[3], y[4]
            z = normalize(y[5:] / N)
            V_molar = self.vessel.V_total / N
            U_molar = U / N
            if len(self.pressure) > 0:
                Pguess = self.pressure[-1]
                Tguess = self.temperature[-1]
            else:
                Pguess = self.operating_pressure
                Tguess = self.operating_temperature

            res = self.flash.flash(
                U=U_molar, V=V_molar, zs=z, Pguess=Pguess, Tguess=Tguess
            )
        else:
            m, N, U = y[0], y[1], y[2]
            z = normalize(y[3:] / N)
            V_molar = self.vessel.V_total / N
            U_molar = U / N
            if len(self.pressure) > 0:
                Pguess = self.pressure[-1]
                Tguess = self.temperature[-1]
            else:
                Pguess = self.operating_pressure
                Tguess = self.operating_temperature

            res = self.flash.flash(
                U=U_molar, V=V_molar, zs=z, Pguess=Pguess, Tguess=Tguess
            )

        ##############################################################################
        # Mass mole/balance part (blowdown, leaks and inflows()
        ##############################################################################
        self.water_level = self.vessel.h_from_V(
            res.liquid1.beta * self._Vm_liq(res.liquid1) * N
        )
        self.liquid_level = self.vessel.h_from_V(
            res.liquid0.beta * self._Vm_liq(res.liquid0) * N
            + res.liquid1.beta * self._Vm_liq(res.liquid1) * N
        )
        static_height = self.liquid_level - self.water_level
        k = res.gas.Cp_ideal_gas() / res.gas.Cv_ideal_gas()

        dm_dt_bdv = self.shutdown * -gas_release_rate(
            res.P,
            self.back_pressure,
            res.gas.rho_mass(),
            k,
            self.orifice_cd,
            self.orifice_area,
        )

        if self.leak_type == "gas":
            dm_dt_leak = self.leak_active * -gas_release_rate(
                res.P,
                self.back_pressure,
                res.gas.rho_mass(),
                k,
                self.leak_cd,
                self.leak_area,
            )
            leak_MW = res.gas.MW()
            leak_zs = res.gas.zs
            leak_U = res.gas.U()
            leak_H = res.gas.H()
        elif self.leak_type == "liquid":
            dm_dt_leak = self.leak_active * -liquid_release_bernouilli(
                res.P,
                self.back_pressure,
                self._liq_density(res.liquid0),
                self.leak_cd,
                self.leak_area,
                static_height,
            )
            leak_MW = res.liquid0.MW()
            leak_zs = res.liquid0.zs
            leak_U = res.liquid0.U()
            leak_H = res.liquid0.H()
        elif self.leak_type == "two_phase":
            if res.P * 0.55 > self.back_pressure:
                Pamb = res.P * 0.55
            else:
                Pamb = self.back_pressure

            dm_dt_leak = self.leak_active * -two_phase_release_fauske(
                res.P, Pamb, self._density(res), self.leak_cd, self.leak_area
            )
            leak_MW = res.MW()
            leak_zs = res.zs
            leak_U = res.U()
            leak_H = res.H()
        else:
            dm_dt_leak = 0
            leak_MW = res.MW()
            leak_zs = res.zs
            leak_U = res.U()
            leak_H = res.H()

        dN_dt_bdv = dm_dt_bdv / (res.gas.MW() / 1000)
        dN_dt_leak = dm_dt_leak / (leak_MW / 1000)

        dm_dt = dm_dt_bdv + dm_dt_leak
        dN_dt = dN_dt_bdv + dN_dt_leak

        ##############################################################################
        # Wall temperature balances
        ##############################################################################
        if self.heat_transfer == "rigorous" or self.heat_transfer == "rigorous_sb_fire":
            h_amb = 8
            h_inner_uw, h_inner_w = h_inside(
                self.length, Tuw, res.T, res.gas
            ), h_inside_wetted(self.length, Tw, res.T, res)
            Auw = self.vessel.A - self.wetted_area()

            Aw = self.wetted_area()
            if self.heat_transfer == "rigorous":
                Quw = Auw * (
                    h_amb * (self.ambient_temperature - Tuw)
                    - h_inner_uw * (Tuw - res.T)
                )
                Qw = Aw * (
                    h_amb * (self.ambient_temperature - Tw) - h_inner_w * (Tw - res.T)
                )
            elif self.heat_transfer == "rigorous_sb_fire":
                Quw = Auw * (
                    sb_fire(Tuw, self.sb_fire_type) - h_inner_uw * (Tuw - res.T)
                )
                Qw = Aw * (sb_fire(Tw, self.sb_fire_type) - h_inner_w * (Tw - res.T))

            Cp = self.wall_heat_capacity
            m_uw = Auw * self.wall_thickness * self.wall_density
            m_w = Aw * self.wall_thickness * self.wall_density
            if m_uw > 0:
                dTuw_dt = Quw / (m_uw * Cp)
            else:
                dTuw_dt = 0

            if m_w > 0:
                dTw_dt = Qw / (m_w * Cp)
            else:
                dTw_dt = 0  # dTuw_dt
        else:
            dTuw_dt = 0
            dTw_dt = 0
            Quw = Qw = 0
        ##############################################################################
        # Energy balance
        ##############################################################################

        if self.mode == "adiabatic":
            dU_dt = res.gas.U() * dN_dt_bdv + leak_U * dN_dt_leak
        elif self.mode == "isentropic":
            dU_dt = res.gas.H() * dN_dt_bdv + leak_H * dN_dt_leak - (Quw + Qw)
        elif self.mode == "fire":
            dU_dt = (
                res.gas.H() * dN_dt_bdv + leak_H * dN_dt_leak + self._get_heat_load_W()
            )
        else:
            dU_dt = 0

        ##############################################################################
        # Component mole balance
        ##############################################################################

        dNs_dt = [
            dN_dt_bdv * zg + dN_dt_leak * zl for zg, zl in zip(res.gas.zs, leak_zs)
        ]

        ##############################################################################
        # Storing out at each integration step
        ##############################################################################

        if True:
            self.times.append(t)
            if self.mode == "isothermal":
                self.sol.append([y[0], y[1]])  # , y[2]])
            elif self.heat_transfer == "rigorous":
                self.sol.append([y[0], y[1], y[2], y[3], y[4]])
            else:
                self.sol.append([y[0], y[1], y[2]])
            self.mdot.append(dm_dt)
            self.mdot_bdv.append(dm_dt_bdv)
            self.mdot_leak.append(dm_dt_leak)
            self.gas_mass.append(N * res.gas.beta * (res.gas.MW() / 1000))
            self.liquid_mass.append(N * res.liquid0.beta * (res.liquid0.MW() / 1000))
            self.water_mass.append(N * res.liquid1.beta * (res.liquid1.MW() / 1000))
            self.pressure.append(res.P)
            self.temperature.append(res.T)
            self.molefracs.append(z)
            self.enthalpy.append(res.U() * N)
            self.liquid_dyn_level.append(self.liquid_level)
            if self.heat_transfer == "rigorous":
                self.wetted_wall_temp.append(y[4])
                self.unwetted_wall_temp.append(y[3])
            else:
                self.wetted_wall_temp.append(0)
                self.unwetted_wall_temp.append(0)

        if self.mode == "isothermal":
            return np.array([dm_dt, dN_dt] + dNs_dt)
        elif self.heat_transfer == "rigorous":
            return np.array([dm_dt, dN_dt, dU_dt, dTuw_dt, dTw_dt] + dNs_dt)
        else:
            return np.array([dm_dt, dN_dt, dU_dt] + dNs_dt)

    def _solout(self, t, y):
        self.stamp_success.append(hash((t, (yi for yi in y))))
        return None

    def depressurize(self):
        # find vessel mass
        # and start mole
        res = self.flash.flash(
            P=self.operating_pressure, T=self.operating_temperature, zs=self.z_adjust
        )

        N0 = self.m0 / (res.MW() / 1000)  # self.vessel.V_total / res.V()
        H0 = res.H() * N0
        U0 = res.U() * N0
        S0 = res.S() * N0
        Ns = N0 * np.asarray(self.z_adjust)
        T = res.T
        args = (
            self.times,
            self.sol,
            self.mdot,
            self.pressure,
            self.temperature,
            self.gas_mass,
            self.molefracs,
            self.enthalpy,
            self.mdot_bdv,
            self.mdot_leak,
            self.liquid_mass,
            self.water_mass,
            self.wetted_wall_temp,
            self.unwetted_wall_temp,
            self.liquid_dyn_level,
        ) = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )
        self.stamp, self.stamp_success = [], []

        fun = self._blowdown_gov_eqns
        if self.mode == "isothermal":
            y0 = [self.m0, N0] + Ns.tolist()
        elif self.heat_transfer == "rigorous":
            y0 = [self.m0, N0, U0] + [T, T] + Ns.tolist()
        else:
            y0 = [self.m0, N0, U0] + Ns.tolist()

        r = ode(fun).set_integrator("dopri5", rtol=1e-4, first_step=5)  # ('dopri5')
        r.set_initial_value(y0)
        r.set_solout(self._solout)
        dt = min(50, self.max_time / 20)
        step = 0

        if self.delay > 0:
            if (dt) >= self.delay:
                first_dt = 0.99 * dt
            else:
                first_dt = dt
            while r.successful() and (r.t + first_dt) < self.delay:
                step += 1
                r.integrate(r.t + dt)
            r.integrate(self.delay)
            self.shutdown = 1
        else:
            self.shutdown = 1

        while r.successful() and (r.t + dt) <= self.max_time:
            step += 1
            r.integrate(r.t + dt)
        r.integrate(self.max_time)

        # Identifying unsuccessfull cases
        fail_list = []
        for i in range(len(self.stamp)):
            if (
                self.stamp[i] not in self.stamp_success
                and self.times[i] != 0
                and self.times[i] != self.delay
            ):
                fail_list.append(i)

        # Reversing the matched successfull cases and deleting unsuccessfull items from the back and forth
        for i in fail_list[::-1]:
            [arg.pop(i) for arg in args]

        return r

    def _blowdown_gov_eqns_euler(
        self, t, y, dt
    ):  # , times):#, sol, mdot, pressure, temperature, gas_mass, molefracs, enthalpy):
        m, Ngas, Nliq, Ugas, Uliq = y[0], y[1], y[2], y[3], y[4]
        Tuw, Tw = y[5], y[6]
        N = Ngas + Nliq
        z = normalize(y[7 : 7 + len(self.z_adjust)] / (Ngas + Nliq))
        z_liq = normalize(y[7 + len(self.z_adjust) :] / (Ngas + Nliq))
        V_molar = self.vessel.V_total / (Ngas + Nliq)
        if Ngas == 0:
            V_molar_gas = 0
            U_molar_gas = 0
        else:
            V_molar_gas = (
                self.vessel.V_total - self.vessel.V_from_h(self.liquid_level)
            ) / Ngas
            U_molar_gas = Ugas / Ngas
        if Nliq == 0:
            V_molar_liq = 0
            U_molar_liq = 0
        else:
            V_molar_liq = self.vessel.V_from_h(self.liquid_level) / Nliq
            U_molar_liq = Uliq / Nliq

        U_molar = (Ugas + Uliq) / (Ngas + Nliq)
        if len(self.pressure) > 0:
            Pguess = self.pressure[-1]
            Tguess = self.temperature[-1]
        else:
            Pguess = self.operating_pressure
            Tguess = self.operating_temperature

        if Ngas == 0 or Nliq == 0:
            gas = liq = self.flash.flash(
                U=U_molar, V=V_molar, zs=z, Pguess=Pguess, Tguess=Tguess
            )
            print(gas.P, liq.P, gas.T, liq.T)
        else:

            def uv_mult(x):
                v_split = x[0]
                flash = self.flash
                Vm_gas = v_split * self.vessel.V_total / Ngas
                Vm_liq = (1 - v_split) * self.vessel.V_total / Nliq
                gres = flash.flash(
                    V=Vm_gas,
                    U=Ugas / Ngas,
                    zs=z,
                    Pguess=self.pressure[-1],
                    Tguess=self.gas_temperatrure[-1],
                )

                lres = flash.flash(
                    V=Vm_liq,
                    U=Uliq / Nliq,
                    zs=z_liq,
                    Pguess=self.pressure[-1],
                    Tguess=self.liquid_temperature[-1],
                )

                self.lresT = lres.T
                self.gresT = gres.T
                self.resP = (lres.P + gres.P) / 2
                return ((gres.P - lres.P) / 1e6) ** 2

            def tot_vol(x):
                Tg = x[0]
                Tl = x[1]
                P = x[2]
                flash = self.flash
                Vg = flash.flash(T=Tg, P=P, zs=z).V() * Ngas
                Vl = flash.flash(T=Tl, P=P, zs=z_liq).V() * Nliq
                V = Vg + Vl
                return (self.vessel.V_total - V) ** 4

            def tot_vol2(x):
                Tg = x[0]
                Tl = x[1]
                P = x[2]
                flash = self.flash
                g = flash.flash(T=Tg, P=P, zs=z)
                l = flash.flash(T=Tl, P=P, zs=z_liq)
                Vg = g.V() * Ngas
                Vl = l.V() * Nliq
                V = Vg + Vl
                retval = (
                    ((self.vessel.V_total - V) / self.vessel.V_total) ** 4
                    + ((Ugas - g.U() * Ngas) / Ugas) ** 4
                    + ((Uliq - l.U() * Nliq) / Uliq) ** 4
                )
                return retval

            def energy_cons_gas(x):
                Tg = x[0]
                Tl = x[1]
                P = x[2]
                flash = self.flash

                Ug = flash.flash(T=Tg, P=P, zs=z).U() * Ngas
                return Ug - Ugas

            def energy_cons_liq(x):
                Tg = x[0]
                Tl = x[1]
                P = x[2]
                flash = self.flash

                Ul = flash.flash(T=Tl, P=P, zs=z_liq).U() * Nliq
                return Ul - Uliq

            x0 = [
                self.gas_temperatrure[-1],
                self.liquid_temperature[-1],
                self.pressure[-1],
            ]

            bounds = (
                (self.gas_temperatrure[-1] - 20, self.gas_temperatrure[-1] + 20),
                (self.liquid_temperature[-1] - 20, self.liquid_temperature[-1] + 20),
                (max(self.pressure[-1] - 10e5, 1.013e5), self.pressure[-1] + 10e5),
            )
            cons = (
                {"type": "eq", "fun": energy_cons_gas},
                {"type": "eq", "fun": energy_cons_liq},
            )

            try:
                # print("Running double UV flash")
                # )
                ret = minimize(
                    fun=tot_vol2,
                    x0=x0,
                    # constraints=cons,
                    bounds=bounds,
                    method="Nelder-Mead",
                    options={"maxiter": 5000},
                )
                if not ret.success:
                    ret = minimize(
                        fun=tot_vol2,
                        x0=x0,
                        # constraints=cons,
                        bounds=bounds,
                        method="SLSQP",
                        # options={"maxiter": 2000},
                    )

                    gas = self.flash.flash(T=ret.x[0], P=ret.x[2], zs=z)
                    liq = self.flash.flash(T=ret.x[1], P=ret.x[2], zs=z_liq)

                else:
                    gas = self.flash.flash(T=ret.x[0], P=ret.x[2], zs=z)
                    liq = self.flash.flash(T=ret.x[1], P=ret.x[2], zs=z_liq)
            except:
                raise

            print(
                t,
                gas.P,
                liq.P,
                gas.T,
                liq.T,
                Ugas,
                Uliq,
                # Ngas,
                # Nliq,
                # z,
                # gas.betas,
                # liq.betas,
            )
            # if not ret.success:
            #    print(ret)
            #    assert 0 == 1
        # ##############################################################################
        # Mass mole/balance part (blowdown, leaks and inflows()
        ##############################################################################

        self.liquid_level = self.vessel.h_from_V(liq.liquid0.V() * Nliq)

        k = gas.gas.Cp_ideal_gas() / gas.gas.Cv_ideal_gas()

        dm_dt_bdv = self.shutdown * -gas_release_rate(
            gas.P,
            self.back_pressure,
            gas.gas.rho_mass(),
            k,
            self.orifice_cd,
            self.orifice_area,
        )

        dN_dt_bdv = dm_dt_bdv / (gas.MW() / 1000)

        dm_dt = dm_dt_bdv
        if gas.gas.beta > 0 and gas.gas.beta != 1:
            dn_gas_to_liq_dt = (1 - gas.gas.beta) * Ngas / (dt * 1.1)
        else:
            dn_gas_to_liq_dt = 0

        dn_liq_to_gas_dt = liq.gas.beta * Nliq / (dt * 1.1)
        dNgas_dt = dN_dt_bdv - dn_gas_to_liq_dt + dn_liq_to_gas_dt
        dNliq_dt = dn_gas_to_liq_dt - dn_liq_to_gas_dt

        # print(t, dn_gas_to_liq_dt, dn_liq_to_gas_dt)

        ##############################################################################
        # Wall temperature balances
        ##############################################################################

        if self.heat_transfer == "rigorous":
            h_amb = 8

            if self.vessel_orientation == "vertical":
                L = self.length - self.liquid_level
            else:
                L = self.diameter - self.liquid_level

            h_inner_uw, h_inner_w = h_inside(
                self.length, Tuw, gas.T, gas
            ), h_inside_wetted(self.length, Tw, gas.T, liq)

            self.h_inner_uw = h_inner_uw

            Auw = self.vessel.A - self.wetted_area()

            Aw = self.wetted_area()
            Quw = Auw * (
                h_amb * (self.ambient_temperature - Tuw)
                - self.h_inner_uw * (Tuw - gas.T)
            )

            Qw = Aw * (
                h_amb * (self.ambient_temperature - Tw) - h_inner_w * (Tw - liq.T)
            )

            if self.vessel_orientation == "vertical":
                L = self.diameter
            else:
                L = self.length
            Qlg = h_inside(L, liq.T, gas.T, gas) * self.vessel.A_cross_sectional(
                h=self.liquid_level
            )

            if liq.T > gas.T:
                Qlg = abs(Qlg)
            else:
                Qlg = -1 * abs(Qlg)

            Cp = self.wall_heat_capacity
            m_uw = Auw * self.wall_thickness * self.wall_density
            m_w = Aw * self.wall_thickness * self.wall_density
            if m_uw > 0:
                dTuw_dt = Quw / (m_uw * Cp)
            else:
                dTuw_dt = 0

            if m_w > 0:
                dTw_dt = Qw / (m_w * Cp)
            else:
                dTw_dt = 0  # dTuw_dt
        ##############################################################################
        # Energy balance
        ##############################################################################
        if Ngas == 0 or Nliq == 0:
            dUgas_dt = (
                gas.H() * dN_dt_bdv
                - Quw
                + Qlg
                - (gas.liquid0.U() + gas.liquid0.U()) / 2 * dn_gas_to_liq_dt
                + (liq.gas.U() + liq.gas.U()) / 2 * dn_liq_to_gas_dt
            )

            dUliq_dt = (
                -Qw
                - Qlg
                + (gas.liquid0.U() + gas.liquid0.U()) / 2 * dn_gas_to_liq_dt
                - (liq.gas.U() + liq.gas.U()) / 2 * dn_liq_to_gas_dt
                # res.gas.H() * dN_dt_bdv - (Quw + Qw)
            )
        else:
            dUgas_dt = (
                gas.H() * dN_dt_bdv
                - Quw
                + Qlg
                - (gas.liquid0.U() + gas.liquid0.H()) / 2 * dn_gas_to_liq_dt
                + (liq.gas.U() + liq.gas.H()) / 2 * dn_liq_to_gas_dt
            )

            dUliq_dt = (
                -Qw
                - Qlg
                + (gas.liquid0.U() + gas.liquid0.H()) / 2 * dn_gas_to_liq_dt
                - (liq.gas.U() + liq.gas.H()) / 2 * dn_liq_to_gas_dt
                # res.gas.H() * dN_dt_bdv - (Quw + Qw)
            )
        ##############################################################################
        # Component mole balance
        ##############################################################################

        dNs_gas_dt = (
            np.asarray([dN_dt_bdv * zg for zg in gas.zs])
            - np.asarray([dn_gas_to_liq_dt * zl for zl in gas.liquid0.zs])
            + np.asarray([dn_liq_to_gas_dt * zg for zg in liq.gas.zs])
        ).tolist()
        dNs_liq_dt = (
            np.asarray([dn_gas_to_liq_dt * zl for zl in gas.liquid0.zs])
            - np.asarray([dn_liq_to_gas_dt * zg for zg in liq.gas.zs])
        ).tolist()
        # (np.asarray([0 * zg for zg in gas.zs])).tolist()
        ##############################################################################
        # Storing out at each integration step
        ##############################################################################

        self.times.append(t)
        self.sol.append([y[0], y[1], y[2], y[3], y[4], y[5], y[6]])
        self.mdot.append(dm_dt)
        self.mdot_bdv.append(dm_dt_bdv)
        # self.gas_mass.append(N * res.gas.beta * (res.gas.MW() / 1000))
        # self.liquid_mass.append(N * res.liquid0.beta * (res.liquid0.MW() / 1000))
        # self.water_mass.append(N * res.liquid1.beta * (res.liquid1.MW() / 1000))
        self.pressure.append(gas.P)
        self.temperature.append(gas.T)
        self.molefracs.append(z)
        # self.enthalpy.append(res.U() * N)
        self.wetted_wall_temp.append(y[6])
        self.unwetted_wall_temp.append(y[5])
        self.gas_temperatrure.append(gas.T)
        self.liquid_temperature.append(liq.T)
        self.liquid_level = self.vessel.h_from_V(liq.V() * (Nliq + dNliq_dt * dt))
        return np.array(
            [dm_dt, dNgas_dt, dNliq_dt, dUgas_dt, dUliq_dt, dTuw_dt, dTw_dt]
            + dNs_gas_dt
            + dNs_liq_dt
        )

    def depressurize_euler(self):
        # find vessel mass
        # and start mole
        res = self.flash.flash(
            P=self.operating_pressure, T=self.operating_temperature, zs=self.z_adjust
        )

        N0 = self.m0 / (res.MW() / 1000)  # self.vessel.V_total / res.V()
        H0 = res.H() * N0
        U0 = res.U() * N0
        S0 = res.S() * N0
        Ns = N0 * np.asarray(self.z_adjust)
        T = res.T
        args = (
            self.times,
            self.sol,
            self.mdot,
            self.pressure,
            self.temperature,
            self.gas_mass,
            self.molefracs,
            self.enthalpy,
            self.mdot_bdv,
            self.mdot_leak,
            self.liquid_mass,
            self.water_mass,
            self.wetted_wall_temp,
            self.unwetted_wall_temp,
            self.gas_temperatrure,
            self.liquid_temperature,
        ) = (
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
            [],
        )
        self.stamp, self.stamp_success = [], []

        y0 = [self.m0, N0, 0, U0, 0] + [T, T] + Ns.tolist() + (Ns * 0).tolist()

        # r = ode(fun).set_integrator("dopri5", rtol=1e-4, first_step=5)  # ('dopri5')
        # r.set_initial_value(y0)
        # r.set_solout(self._solout)
        dt = self.dt  # min(50, self.max_time / 20)
        step = 0

        t = 0
        self.shutdown = 1

        while t < self.max_time:
            # if t < 10:
            #     dt = 0.1
            # elif t < 100:
            #     dt = 1
            # else:
            #     dt = 10
            dy = self._blowdown_gov_eqns_euler(t, y0, dt)
            y0 = y0 + dy * dt
            # print(f"Time {t}")
            t += dt
        r = None

        return r

    def plot(self):
        # print(self.pressure)
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(self.times, np.array(self.pressure) / 1e5, "bo", label="Pressure")
        plt.xlabel("Time (sec)")
        plt.ylabel("Pressure (bar)")

        plt.legend(shadow=True)
        plt.figure()
        plt.plot(self.times, -np.array(self.mdot) * 3600, label="Mass inventory")
        plt.xlabel("Time (sec)")
        plt.ylabel("Blowdown mass flow(kg/h)")

        # plt.figure()
        # plt.plot(self.times, np.array(self.enthalpy) / 1e6, label="Enthalpy")
        # plt.xlabel("Time (sec)")
        # plt.ylabel("Enthalpy (MJ)")
        # plt.show()

        plt.figure()
        plt.plot(
            self.times, np.asarray(self.temperature) - 273.15, label="Fluid temperature"
        )
        plt.plot(
            self.times,
            np.asarray(self.unwetted_wall_temp) - 273.15,
            label="Unwetted wall",
        )
        plt.plot(
            self.times, np.asarray(self.wetted_wall_temp) - 273.15, label="Wetted wall"
        )
        plt.legend(loc="best")
        plt.xlabel("Time (sec)")
        plt.ylabel("Temperature (K)")

        plt.show()


if __name__ == "__main__":
    import time

    t1 = time.time()
    input = {}
    P = 12e5
    T = 298.15

    input["mode"] = "isentropic"
    input["heat_transfer"] = "rigorous"
    input["wall_thickness"] = 0.010  # m
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
    input["ambient_temperature"] = 273.15
    input["back_pressure"] = 1.01e5
    input["bdv_orifice_size"] = 0.03  # m
    input["bdv_orifice_cd"] = 0.9

    input["leak_active"] = 0
    input["leak_size"] = 0.01  # m

    input["leak_cd"] = 0.65
    input["leak_type"] = "liquid"

    names = ["methane", "propane", "n-butane", "i-butane", "n-decane"]
    molefracs = [0.8, 0.05, 0.01, 0.01, 0.20]

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

    t2 = time.time()
    print(f"Calculation time {t2 - t1} s")
