from thermopack.cubic import PengRobinson
from thermopack.tcPR import tcPR
import scienceplots
import os


# Importing Numpy (math, arrays, etc...)
import numpy as np

# Importing Matplotlib (plotting)
import matplotlib.pyplot as plt

plt.style.use(["science", "nature"])

tc_pr = PengRobinson("N2")

# Plot phase envelope
z = np.array([1.0])
Tc, vc, Pc = tc_pr.critical(z)
file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "n2.txt")
n2 = np.loadtxt(file)
plt.plot(n2[:, 0], n2[:, 1], color="blue")
plt.plot([Tc], [Pc * 1.0e-5], "o", color="blue")
plt.plot([16 + 273.15], [150], "^", color="blue", label="Experiment A")


tc_pr = PengRobinson("C1,C2")

# Plot phase envelope
z = np.array([0.91, 0.09])
T, P, v = tc_pr.get_envelope_twophase(1.0e4, z, maximum_pressure=1.5e7, calc_v=True)
Tc, vc, Pc = tc_pr.critical(z)
plt.plot(T, P * 1.0e-5, "--", color="green")
plt.plot([Tc], [Pc * 1.0e-5], "o", color="green")
plt.plot([30 + 273.15], [121.6], "v", color="green", label="Experiment B")

tc_pr = PengRobinson("C1,C2,C3,iC4")

# Plot phase envelope
z = np.array([0.64, 0.06, 0.28, 0.02])
T, P, v = tc_pr.get_envelope_twophase(1.0e4, z, maximum_pressure=1.5e7, calc_v=True)
Tc, vc, Pc = tc_pr.critical(z)
plt.plot(T, P * 1.0e-5, "--", color="violet")
plt.plot([Tc], [Pc * 1.0e-5], "o", color="violet")
plt.plot([20 + 293.15], [117.5], ">", color="violet", label="Experiment C")

tc_pr = tcPR("CO2,N2")

# Plot phase envelope
z = np.array([0.70, 0.30])
T, P, v = tc_pr.get_envelope_twophase(1.0e4, z, maximum_pressure=2e7, calc_v=True)
Tc, vc, Pc = tc_pr.critical(z)
plt.plot(T, P * 1.0e-5, "-.", color="red")
plt.plot([Tc], [Pc * 1.0e-5], "o", color="red")
plt.plot([30 + 293.15], [150], "<", color="red", label="Experiment D")

plt.ylabel("Pressure (bar)")
plt.xlabel("Temperature (K)")
# plt.title("tcPR phase diagram")
plt.legend(loc="best")
plt.savefig("phase_diagram.png", dpi=600)
plt.show()

plt.clf()
