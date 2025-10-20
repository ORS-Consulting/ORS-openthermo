from thermopack.cubic import SoaveRedlichKwong, PengRobinson
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
from scipy.optimize import root
from scipy.optimize import root_scalar
from ctypes import *
import os
from chemicals import normalize
from thermo.equilibrium import EquilibriumState
from thermo import (
    ChemicalConstantsPackage,
    PropertyCorrelationsPackage,
    CEOSGas,
    CEOSLiquid,
    PRMIX,
    SRKMIX,
    FlashVLN,
    FlashVL,
)
from openthermo.properties.transport import COSTALD_rho, COSTALD_Vm

# name_map = np.loadtxt("name_mapping.csv", dtype=str, delimiter=",", usecols=(0, 1, 2))
name_map = np.loadtxt(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "name_mapping.csv"),
    dtype=str,
    delimiter=",",
    usecols=(0, 1, 2),
)


class FlashVL:
    def __init__(self, zs, constants, properties, kijs, model="PR", rho=None):
        # names = []
        # for CAS in constants.CASs:
        #    if CAS not in name_map[:, 1]:
        #        raise ValueError("CAS %s not in name mapping" % CAS)
        #    else:
        #        names.append(name_map[np.where(name_map[:, 1] == CAS)[0][0], 2])
        #        print(names[-1])

        self.vectorized = False
        self.constants = constants
        self.properties = properties
        # self.zs = normalize(zs)
        self.model = 1

        if model == 2 or model == "SRK":
            self.model = 2
            # self.eos = SoaveRedlichKwong(",".join(names))
            self.eos = SoaveRedlichKwong(",".join(len(constants.names) * ["PSEUDO"]))
        else:
            # self.eos = PengRobinson(",".join(names))
            self.eos = PengRobinson(",".join(len(constants.names) * ["PSEUDO"]))
        cindices = range(1, self.eos.nc + 1)
        self.eos.init_pseudo(
            comps=",".join(constants.names),
            Tclist=constants.Tcs,
            Pclist=constants.Pcs,
            acflist=constants.omegas,
            Mwlist=constants.MWs,
        )
        _ = [
            [self.eos.set_kij(i, j, kijs[i - 1][j - 1]) for i in cindices]
            for j in cindices
        ]

        self.kijs = kijs
        self.zs = list(np.asarray(zs) / sum(zs))
        water_mole_fracs = np.zeros(len(constants.CASs))
        self.rho_model = rho

        if "7732-18-5" in constants.CASs:
            water_mole_fracs[constants.CASs.index("7732-18-5")] = 1.0
        self.water_mole_fracs = list(water_mole_fracs)
        self.eos_kwargs = dict(
            Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas, kijs=kijs
        )

        if self.model == 1:
            self.gas = CEOSGas(
                PRMIX,
                self.eos_kwargs,
                HeatCapacityGases=self.properties.HeatCapacityGases,
            )
            self.liq = CEOSLiquid(
                PRMIX,
                self.eos_kwargs,
                HeatCapacityGases=self.properties.HeatCapacityGases,
            )
            self.liq2 = CEOSLiquid(
                PRMIX,
                self.eos_kwargs,
                HeatCapacityGases=self.properties.HeatCapacityGases,
            )
        else:
            self.gas = CEOSGas(
                SRKMIX,
                self.eos_kwargs,
                HeatCapacityGases=self.properties.HeatCapacityGases,
            )
            self.liq = CEOSLiquid(
                SRKMIX,
                self.eos_kwargs,
                HeatCapacityGases=self.properties.HeatCapacityGases,
            )
            self.liq2 = CEOSLiquid(
                SRKMIX,
                self.eos_kwargs,
                HeatCapacityGases=self.properties.HeatCapacityGases,
            )
        self.fall_back = FlashVLN(
            self.constants, self.properties, gas=self.gas, liquids=[self.liq, self.liq2]
        )
        self.dest = EquilibriumState

    def flash(
        self,
        zs=None,
        T=None,
        P=None,
        V=None,
        H=None,
        S=None,
        G=None,
        U=None,
        A=None,
        D=None,
        Tguess=None,
        Pguess=None,
    ):
        r"""Method to perform a flash calculation and return the result as an
        :obj:`EquilibriumState <thermo.equilibrium.EquilibriumState>` object.

        Parameters
        ----------
        zs : list[float], optional
            Mole fractions of each component, required unless there is only
            one component, [-]
        T : float, optional
            Temperature, [K]
        P : float, optional
            Pressure, [Pa]
        V : float, optional
            Molar volume of the overall bulk, [m^3/mol]
        H : float, optional
            Molar enthalpy of the overall bulk, [J/mol]
        S : float, optional
            Molar entropy of the overall bulk, [J/(mol*K)]
        G : float, optional
            Molar Gibbs free energy of the overall bulk, [J/mol]
        U : float, optional
            Molar internal energy of the overall bulk, [J/mol]
        A : float, optional
            Molar Helmholtz energy of the overall bulk, [J/mol]
        """
        if zs is None:
            if self.zs is not None:
                zs = self.zs
            else:
                raise ValueError("Composition missing for flash")
        else:
            self.zs = normalize(zs)

        T_spec = T is not None
        P_spec = P is not None
        V_spec = V is not None
        H_spec = H is not None
        S_spec = S is not None
        U_spec = U is not None
        A_spec = A is not None
        G_spec = G is not None
        D_spec = D is not None

        HSGUA_spec_count = H_spec + S_spec + G_spec + U_spec + A_spec

        flash_specs = {"zs": zs}
        if T_spec:
            flash_specs["T"] = T
            if T <= 0.0:
                raise ValueError("Specified temperature (%s K) is unphysical" % (T,))
        if P_spec:
            flash_specs["P"] = P
            if P <= 0.0:
                raise ValueError("Specified pressure (%s Pa) is unphysical" % (P,))
        if V_spec:
            flash_specs["V"] = V
            if V <= 0.0:
                raise ValueError(
                    "Specified molar volume (%s m^3/mol) is unphysical" % (V,)
                )
        if H_spec:
            flash_specs["H"] = H
        if S_spec:
            flash_specs["S"] = S
        if U_spec:
            flash_specs["U"] = U
        if G_spec:
            flash_specs["G"] = G
        if A_spec:
            flash_specs["A"] = A

        self.flash_specs = flash_specs

        if T_spec and P_spec:
            res = self.PT_flash(T=T, P=P)
            return res

        if P_spec and H_spec:
            dH = lambda T2: (self.PT_flash(T=T2, P=P).H() - H)
            dH_min = lambda T2: ((self.PT_flash(T=T2, P=P).H() - H) ** 2)
            if Tguess == None:
                Tguess = 298
            res = root_scalar(dH, x0=Tguess, method="secant", xtol=1e-3)

            x = res.root

            if abs((dH(x)) / H) > 1e-6:
                res = root(dH, Tguess, method="broyden1", options={"ftol": 1e-2})
                x = res.x
                if abs((dH(x)) / H) > 1e-6:
                    print("Trying Nelser-Mead as last resort")
                    res = minimize(
                        dH_min,
                        Tguess,
                        method="Nelder-Mead",
                        # bounds=(Tguess - 10, Tguess + 10),
                    )
                    x = res.x
                    if abs((dH(x)) / H) > 1e-6:
                        raise ValueError("PH-flash failed to converge")
                x = res.x
            res = self.PT_flash(T=float(x), P=P)
            return res
        elif P_spec and U_spec:
            dU = lambda T2: (self.PT_flash(T=T2, P=P).U() - U)
            if Tguess == None:
                Tguess = 298
            res = root(
                dU,
                Tguess,
                method="hybr",
                jac=None,
                tol=None,
                callback=None,
                options=None,
            )
            if abs((dU(res.x[0])) / U) > 1e-6:
                raise ValueError("PU-flash failed to converge")
            res = self.PT_flash(T=float(res.x[0]), P=P)
            return res
        elif P_spec and S_spec:
            dS = lambda T2: (self.PT_flash(T=T2, P=P).S() - S)
            dS_min = lambda T2: ((self.PT_flash(T=T2, P=P).S() - S) ** 2)
            if Tguess == None:
                Tguess = 298

            # res = minimize_scalar(
            #     dS_min,
            #     Tguess,
            #     bounds=(0, 3000),
            # )  # ,bounds=((1,2000),(1,2000e5)))
            # res = root(dS, Tguess, method="krylov")
            res = root_scalar(dS, x0=Tguess, method="secant", bracket=[10, 3000])
            x = res.root
            if abs((dS(x)) / S) > 1e-5:
                res = root(dS, Tguess, method="broyden1", options={"ftol": 1e-4})
                x = res.x
                if abs((dS(x)) / S) > 1e-5:
                    res = minimize(
                        dS_min,
                        Tguess,
                        method="Nelder-Mead",
                    )
                    x = res.x
                    if abs((dS(x)) / S) > 1e-5:
                        raise ValueError("PS-flash failed to converge")
            res = self.PT_flash(T=float(x), P=P)
            return res
        elif P_spec and V_spec:
            if Tguess == None:
                Tguess = 298

            if self.rho_model == "costald":
                dV = lambda T2: (self.V_molar(self.PT_flash(T=T2, P=P)) - V)
                dVmin = lambda T2: (self.V_molar(self.PT_flash(T=T2, P=P)) - V) ** 2
            else:
                dV = lambda T2: (self.PT_flash(T=T2, P=P).V() - V)
                dVmin = lambda T2: (self.PT_flash(T=T2, P=P).V() - V) ** 2

            res = root(
                dV,
                Tguess,
                method="hybr",
                jac=None,
                tol=None,
                callback=None,
                options=None,
            )
            # print(res)
            if abs((dV(res.x)) / V) > 1e-4:
                raise ValueError("PV-flash failed to converge")
            res = self.PT_flash(T=float(res.x[0]), P=P)
            return res
        elif T_spec and V_spec:
            if Pguess == None:
                Pguess = 10e5

            if self.rho_model == "costald":
                dV = lambda P2: (self.V_molar(self.PT_flash(T=T, P=P2)) - V)
                dVmin = lambda P2: (self.V_molar(self.PT_flash(T=T, P=P2)) - V) ** 2
                res = root_scalar(
                    dV, x0=Pguess, bracket=[0, 1000e5], method="secant"
                )  # , method="broyden1", options={"ftol": 1e-4})
                x = res.root
            else:
                dV = lambda P2: (self.PT_flash(T=T, P=P2).V() - V)
                dVmin = lambda P2: (self.PT_flash(T=T, P=P2).V() - V) ** 2
                res = minimize_scalar(
                    dVmin,
                    Pguess,
                    bounds=(100, 1000e5),
                )
                x = res.x
            if abs(dV(x) / V) > 1e-4:
                raise ValueError("VT-flash failed to converge")
            if x < 0:
                raise ("Zero pressure error")
            res = self.PT_flash(T=float(T), P=float(x))
            return res
        elif T_spec and U_spec:
            dU = lambda P2: (self.PT_flash(T=T, P=P2).U() - U)
            if Pguess == None:
                Pguess = 10e5
            res = root(dU, Pguess, method="broyden1", options={"ftol": 1e-4})
            # print(res)
            if abs(dU(res.x) / U) > 1e-4:
                raise ValueError("UT-flash failed to converge")
            res = self.PT_flash(T=float(T), P=float(res.x))
            return res
        elif U_spec and V_spec:
            if Pguess == None:
                Pguess = self.flash(V=V, T=298).P
            if Tguess == None:
                Tguess = self.flash(P=Pguess, U=U).T
            x0 = [Tguess, Pguess]
            res = minimize(
                self.UVres,
                x0,
                args=(U, V),
                method="Nelder-Mead",
                bounds=((1, 2000), (1, 1000e5)),
            )  # ,bounds=((1,2000),(1,2000e5)))

            T2 = res.x[0]
            P2 = res.x[1]
            if self.rho_model == "costald":
                err = abs((self.PT_flash(T=res.x[0], P=res.x[1]).U() - U) / U)
                +abs((self.V_molar(self.PT_flash(T=res.x[0], P=res.x[1])) - V) / V)
            else:
                err = abs((self.PT_flash(T=res.x[0], P=res.x[1]).U() - U) / U)
                +abs((self.PT_flash(T=res.x[0], P=res.x[1]).V() - V) / V)

            if err > 1e-6:
                res = minimize(
                    self.UVres,
                    x0,
                    args=(U, V),
                    method="Powell",
                    bounds=((200, 300), (1e5, 100e5)),
                )  # ,bounds=((1,2000),(1,2000e5)))
                T2 = res.x[0]
                P2 = res.x[1]
                if self.rho_model == "costald":
                    err = abs((self.PT_flash(T=res.x[0], P=res.x[1]).U() - U) / U)
                    +abs((self.V_molar(self.PT_flash(T=res.x[0], P=res.x[1])) - V) / V)
                else:
                    err = abs((self.PT_flash(T=res.x[0], P=res.x[1]).U() - U) / U)
                    +abs((self.PT_flash(T=res.x[0], P=res.x[1]).V() - V) / V)
                if err > 1e-6:
                    raise ValueError("UV-flash failed to converge")

            res = self.PT_flash(T=float(T2), P=float(P2))
            return res
        else:
            raise Exception("Flash inputs unsupported")

    def UVres(self, x, U, V):
        if self.rho_model == "costald":
            res = ((U - self.PT_flash(T=x[0], P=x[1]).U()) / U) ** 2 + (
                (V - self.V_molar(self.PT_flash(T=x[0], P=x[1]))) / V
            ) ** 2
        else:
            res = ((U - self.PT_flash(T=x[0], P=x[1]).U()) / U) ** 2 + (
                (V - self.PT_flash(T=x[0], P=x[1]).V()) / V
            ) ** 2
        return res

    def V_molar(self, res):
        nVm_gas = res.gas.V() * res.gas.beta
        nVm_liquid = COSTALD_Vm(res.liquid0) * res.liquid0.beta
        nVm_water = COSTALD_Vm(res.liquid1) * res.liquid1.beta
        return nVm_gas + nVm_liquid + nVm_water

    def PT_flash(self, T=None, P=None):

        # res = self.flashvln.flash(T=T, P=P, zs=self.zs)
        flash = self.eos.two_phase_tpflash(temp=T, press=P, z=self.zs)

        if flash.betaL > 0 and flash.betaL < 1:
            liq_beta = flash.betaL
            x = flash.x
            gas_beta = flash.betaV
            y = flash.y
        elif flash.betaL == 0:
            gas_beta = 1
            liq_beta = 0
            x = flash.y
            y = flash.y
        elif flash.betaL == 1:
            gas_beta = 1
            liq_beta = 0
            x = flash.x
            y = flash.x
        else:
            gas_beta = 1
            liq_beta = 0
            x = flash.x
            y = flash.x

        # if len(res.betas) :
        #    wat_beta = res.liquid1.beta
        # else:
        wat_beta = 0

        betas = [gas_beta, liq_beta, wat_beta]
        gas, liq, liq2 = None, None, None
        liquids = []
        w = self.water_mole_fracs

        #################################
        # Check this, could give problems
        # Consider not settign to global composition
        ##################################

        w = w

        # gas = CEOSGas(PRMIX, self.eos_kwargs, HeatCapacityGases=self.properties.HeatCapacityGases, T=T, P=P, zs=list(y))
        # liq = CEOSLiquid(PRMIX, self.eos_kwargs, HeatCapacityGases=self.properties.HeatCapacityGases, T=T, P=P, zs=list(x))
        # liq2 = CEOSLiquid(PRMIX, self.eos_kwargs, HeatCapacityGases=self.properties.HeatCapacityGases, T=T, P=P, zs=w)

        gas = self.gas.to(
            T=T, P=P, zs=list(y)
        )  # CEOSGas(PRMIX, eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases, T=T, P=P, zs=z)

        liq = self.liq.to(
            T=T, P=P, zs=list(x)
        )  # CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases, T=T, P=P, zs=z)
        liq2 = self.liq2.to(
            T=T, P=P, zs=list(x)
        )  # CEOSLiquid(PRMIX, eos_kwargs, HeatCapacityGases=properties.HeatCapacityGases, T=T, P=P, zs=z)

        return self.dest(
            T,
            P,
            self.zs,
            gas=gas,
            liquids=[liq, liq2],
            solids=[],
            betas=betas,
            flash_specs=self.flash_specs,
            constants=self.constants,
            correlations=self.properties,
            flasher=self,
        )


def get_flash_dry(
    pure_comp_names,
    pure_comp_molefracs,
    pseudo_names=None,
    pseudo_molefracs=None,
    pseudo_SGs=None,
    pseudo_Tbs=None,
    P=None,
    T=None,
    rho=None,
    model="PR",
):
    from chemicals import similarity_variable
    from openthermo.data import interaction_parameters as ip
    import openthermo.properties.pseudo as lk

    pure_constants, properties = ChemicalConstantsPackage.from_IDs(pure_comp_names)
    kijs = ip.get_interaction_parameters(pure_constants.CASs)

    if pseudo_names is None and pseudo_molefracs is None:
        zs = normalize(
            pure_comp_molefracs
        )  # pure_comp_molefracs/sum(pure_comp_molefracs)
        # print(kijs)
        flash = FlashVL(zs, pure_constants, properties, kijs)
    else:
        pseudo_Tcs = [
            lk.Tc_Kesler_Lee_SG_Tb(SG, Tb) for SG, Tb in zip(pseudo_SGs, pseudo_Tbs)
        ]
        pseudo_Pcs = [
            lk.Pc_Kesler_Lee_SG_Tb(SG, Tb) for SG, Tb in zip(pseudo_SGs, pseudo_Tbs)
        ]
        pseudo_MWs = [
            lk.MW_Kesler_Lee_SG_Tb(SG, Tb) for SG, Tb in zip(pseudo_SGs, pseudo_Tbs)
        ]
        pseudo_omegas = [
            lk.omega_Kesler_Lee_SG_Tb_Tc_Pc(SG, Tb)
            for SG, Tb in zip(pseudo_SGs, pseudo_Tbs)
        ]
        pseudo_Zcs = [lk.Zc_pseudo(omega) for omega in pseudo_omegas]
        pseudo_Vcs = [
            lk.Vc_pseudo(Zc, Tc, Pc)
            for Zc, Tc, Pc in zip(pseudo_Zcs, pseudo_Tcs, pseudo_Pcs)
        ]
        atomic_ratio = [
            lk.HC_atomic_ratio(SG, Tb) for SG, Tb in zip(pseudo_SGs, pseudo_Tbs)
        ]
        monomer_MWs = [12.01 + HC_ratio * 1.008 for HC_ratio in atomic_ratio]
        pseudo_carbon_numbers = [
            MW / monomer_MW for MW, monomer_MW in zip(pseudo_MWs, monomer_MWs)
        ]
        hydrogen_counts = [c * HC for c, HC in zip(pseudo_carbon_numbers, atomic_ratio)]
        pseudo_atoms = [
            {"C": C, "H": H} for C, H in zip(pseudo_carbon_numbers, hydrogen_counts)
        ]
        similarity_variables = [
            similarity_variable(atoms=atoms) for atoms in pseudo_atoms
        ]

        pseudos = ChemicalConstantsPackage(
            names=pseudo_names,
            MWs=pseudo_MWs,
            Tbs=pseudo_Tbs,
            atomss=pseudo_atoms,
            Tcs=pseudo_Tcs,
            Pcs=pseudo_Pcs,
            omegas=pseudo_omegas,
            similarity_variables=similarity_variables,
            Zcs=pseudo_Zcs,
            Vcs=pseudo_Vcs,
        )
        # Add the pure components and the pseudocomponents to create a new package of constant values
        # which will be used by the phase and flash objects
        constants = pure_constants + pseudos
        # Obtain the temperature and pressure dependent objects
        properties = PropertyCorrelationsPackage(constants=constants)

        zs = normalize(pure_comp_molefracs + pseudo_molefracs)
        pseudo_kijs = np.zeros((len(zs), len(zs)))
        n_pure = len(pure_comp_names)
        pseudo_kijs[:n_pure, :n_pure] = np.asarray(kijs)
        kijs = pseudo_kijs.tolist()

        flash = FlashVL(zs, constants, properties, kijs, model=model, rho=rho)
    return flash


def flash_results(res):
    """
    Make list of lists for excel output to range/cell. Input is a flash result object
    """
    out = []
    out.append(
        ["Phase", "", "", "Overall", "Gas", "HC liquid", "Aqueous"]
    )  # Make this more robust
    out.append(
        ["Phase fraction (mole)"] + ["", "", ""] + [phase.beta for phase in res.phases]
    )
    out.append(
        ["Phase fraction (mass)"]
        + ["", "", ""]
        + [phase.beta_mass for phase in res.phases]
    )
    out.append(
        ["Compressibility (--)", "", ""]
        + [res.Z()]
        + [phase.Z() for phase in res.phases]
    )
    out.append(
        ["MW (kg/kmole)", "", ""] + [res.MW()] + [phase.MW() for phase in res.phases]
    )
    if res.flasher.rho_model == "costald":
        nVm_gas = res.gas.V() * res.gas.beta
        nVm_liquid = COSTALD_Vm(res.liquid0) * res.liquid0.beta
        nVm_water = COSTALD_Vm(res.liquid1) * res.liquid1.beta
        Vm = nVm_gas + nVm_liquid + nVm_water
        rho = 1 / Vm * res.MW() / 1000
        rho_liq = COSTALD_rho(res.liquid0)
        rho_water = COSTALD_rho(res.liquid1)

        out.append(
            ["Density (kg/m3)", "", ""]
            + [rho]
            + [res.gas.rho_mass(), rho_liq, rho_water]
        )

    else:
        out.append(
            ["Density (kg/m3)", "", ""]
            + [res.rho_mass()]
            + [phase.rho_mass() for phase in res.phases]
        )
    out.append(
        ["Enthalpy (kJ/kg)", "", ""]
        + [res.H_mass() / 1000]
        + [phase.H_mass() / 1000 for phase in res.phases]
    )
    out.append(
        ["Entropy (kJ/kg K)", "", ""]
        + [res.S_mass() / 1000]
        + [phase.S_mass() / 1000 for phase in res.phases]
    )
    out.append(
        ["Gibbs free energy (kJ/kg)", "", ""]
        + [res.G_mass() / 1000]
        + [phase.G_mass() / 1000 for phase in res.phases]
    )
    out.append(
        ["Cp (kJ/kg K)", "", ""]
        + [res.Cp_mass() / 1000]
        + [phase.Cp_mass() / 1000 for phase in res.phases]
    )
    out.append(
        ["Cv (kJ/kg K)", "", ""]
        + [res.Cv_mass() / 1000]
        + [phase.Cp_mass() / 1000 for phase in res.phases]
    )
    out.append(
        ["k (real gas)", "", ""]
        + [res.Cp() / res.Cv()]
        + [phase.Cp() / phase.Cv() for phase in res.phases]
    )
    out.append(
        ["k (ideal gas)", "", ""]
        + [res.Cp_ideal_gas() / res.Cv_ideal_gas()]
        + [phase.Cp_ideal_gas() / phase.Cv_ideal_gas() for phase in res.phases]
    )
    out.append(
        ["Viscosity (cP)", "", ""]
        + [res.mu() * 1000]
        + [phase.mu() * 1000 for phase in res.phases]
    )

    out.append(
        ["Thermal conductivity (W/m K)", "", ""]
        + [res.k()]
        + [phase.k() for phase in res.phases]
    )
    out.append(
        ["Component mole fractions: ", "", "", "Overall", "Gas", "HC liquid", "Aqueous"]
    )

    comps = np.array(
        [
            res.names,
            ["" for i in range(len(res.names))],
            ["" for i in range(len(res.names))],
            res.zs,
            res.gas.zs,
            res.lightest_liquid.zs,
            res.heaviest_liquid.zs,
        ]
    ).T.tolist()
    out = out + comps
    return out


def phase_string(res):
    phase_str = ""
    for phase in res.phases:
        if phase.beta > 0.0001:
            if phase.phase == "g":
                phase_str += str("V")
            else:
                phase_str += str("L")
    return phase_str


if __name__ == "__main__":
    from scipy.constants import atm
    from chemicals import *
    from thermo import *
    from openthermo.data import interaction_parameters as ip
    import openthermo.properties.pseudo as lk

    pure_constants, properties = ChemicalConstantsPackage.from_IDs(
        ["water", "methane", "propane", "n-butane", "i-butane", "n-decane"]
    )
    pseudo_carbon_numbers = [10, 16]
    MWs = [140.0, 325.0]
    hydrogen_counts = [
        (MW - C * periodic_table.C.MW) / periodic_table.H.MW
        for C, MW in zip(pseudo_carbon_numbers, MWs)
    ]
    pseudo_atoms = [
        {"C": C, "H": H} for C, H in zip(pseudo_carbon_numbers, hydrogen_counts)
    ]
    # Calculate the similarity variable of each species
    similarity_variables = [similarity_variable(atoms=atoms) for atoms in pseudo_atoms]
    Tcs = [606.28, 825.67]
    Pcs = [25.42 * atm, 14.39 * atm]
    omegas = [0.4019, 0.7987]
    MWs = [140.0, 325.0]
    Zcs = [lk.Zc_pseudo(omega) for omega in [0.4019, 0.7987]]

    Vcs = [lk.Vc_pseudo(Zc, Tc, Pc) for Zc, Tc, Pc in zip(Zcs, Tcs, Pcs)]

    pseudos = ChemicalConstantsPackage(
        names=["PS1", "PS2"],
        similarity_variables=similarity_variables,
        atomss=pseudo_atoms,
        Tcs=[606.28, 825.67],
        Pcs=[25.42 * atm, 14.39 * atm],
        omegas=[0.4019, 0.7987],
        MWs=[140.0, 325.0],
        Zcs=[lk.Zc_pseudo(omega) for omega in [0.4019, 0.7987]],
        Vcs=Vcs,
    )
    constants = pure_constants + pseudos
    properties = PropertyCorrelationsPackage(constants=constants)

    model = 1
    T = 94 + 273
    P = 25e5
    zs = [0.2, 0.2, 0.1, 0.1, 0.1, 0.3, 0.1, 0.1]
    zs = (np.asarray(zs) / np.asarray(zs).sum()).tolist()
    z = np.asarray(zs) / np.asarray(zs).sum()
    tc = constants.Tcs
    pc = constants.Pcs
    omega = constants.omegas
    kijs = [
        [0.0, 0.5, 0.5, 0.5, 0.5, 0.5],
        [0.5, 0.0, 0.006829, 0.013110, 0.012300, 0.043610],
        [0.5, 0.006829, 0.0, 0.001041, 0.000819, 0.016630],
        [0.5, 0.1311, 0.001041, 0.0, 1.33e-5, 0.009448],
        [0.5, 0.0123, 0.000819, 1.33e-5, 0, 0.01016],
        [0.5, 0.04361, 0.01663, 0.009448, 0.01016, 0.0],
    ]

    kijs = ip.get_interaction_parameters(pure_constants.CASs)
    kij = np.zeros((len(zs), len(zs)))
    kij[: len(kijs), : len(kijs)] = kijs
    kijs = kij.tolist()
    x, y, phi_x, phi_y, compfact, beta, ier = michelsen_flash(
        model, T, P, zs, tc, pc, omega, kij=kijs
    )
    eos_kwargs = dict(
        Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas, kijs=kijs
    )

    gas = CEOSGas(
        PRMIX,
        eos_kwargs,
        HeatCapacityGases=properties.HeatCapacityGases,
        T=T,
        P=P,
        zs=zs,
    )
    liq = CEOSLiquid(
        PRMIX,
        eos_kwargs,
        HeatCapacityGases=properties.HeatCapacityGases,
        T=T,
        P=P,
        zs=zs,
    )
    liq2 = CEOSLiquid(
        PRMIX,
        eos_kwargs,
        HeatCapacityGases=properties.HeatCapacityGases,
        T=T,
        P=P,
        zs=zs,
    )
    phase_list = [gas, liq, liq]

    flashN = FlashVLN(constants, properties, gas=gas, liquids=[liq, liq2])
    # flashN.PT_SS_TOL = 1e-18
    res1 = flashN.flash(T=T, P=P, zs=zs)
    print("There are %s phases present" % (res1.phase_count))

    flash = FlashVLW(z, constants, properties, kijs)
    res = flash.flash(T=T, P=P)
    res2 = flash.flash(P=1e5, H=res.H())
    V = res.V()
    U = res.U()
    D = res.rho_mass()
    # res = flash.flash(U=U * 1.1, V=0.8 * V, Pguess=25e5, Tguess=367)
    # res = flash.flash(P=P, V=.8*V, Tguess=320)
    # zs[0] = 0
    # zs[-1] = 0.5
    # flash = FlashVLW(zs, constants, properties, kijs)
    # res3 = flash.flash(T=T, P=P)
