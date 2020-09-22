# SI
Ns = 20000  # Number of cells in series
A = 0.49 * Ns  # Area of the solar cell in cm²
I_0 = (1 / (10 ** (13))) * Ns
Isc_ref = (17.2 / 1000) * Ns

rs = 1.4 / Ns
rsh = 7800 / Ns
eg = (1.24 * 1.602176634) / (10 ** 19)  # eV
n = 1
temp_ref = 25  # °C
alpha = -0.0037  # 1/K  temperature coefficient of Jsc


p_mp = 0.00446 * Ns
EQE_filename = "CHEN_2020_EQE_curve_si_corrected.csv"
