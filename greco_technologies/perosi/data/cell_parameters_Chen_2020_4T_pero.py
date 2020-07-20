#PERO
Ns = 30  # Number of cells in series
A = 0.049 * Ns  # Area of the solar cell in cm²
I_0 = 7.6/10**(15)*A
Isc_ref=(22.3/1000) * A # A for 700 nm thickness

rs = 3.2/A
rsh =6230/A
eg =1.636 *1.602176634 / 10 ** 19  # eV
n = 1.5
temp_ref = 25  # °C
alpha = -0.0017 # 1/K  temperature coefficient of Jsc

p_mp_pero = 0.0009
p_mp = 0.013 + 0.009 #Watt/Zelle perovskite + silicone
EQE_filename= "CHEN_2020_EQE_curve_pero_corrected.csv"