#SI
Ns = 30  # Number of cells in series
A = 0.049 * Ns  # Area of the solar cell in cm²
I_0 = 1/10**(13)*A
Isc_ref=(17.2/1000) * A

rs =1.4/A
rsh =6230/A
eg =1.24 *1.602176634 / 10 ** 19 # eV
n = 1
temp_ref = 25  # °C
alpha = -0.0037  # 1/K  temperature coefficient of Jsc


p_mp_si = 0.000446
EQE_filename= "CHEN_2020_EQE_curve_si_corrected.csv"