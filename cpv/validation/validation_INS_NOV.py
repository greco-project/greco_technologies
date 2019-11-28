import pandas as pd
import pvlib.atmosphere
import matplotlib.pyplot as plt
import numpy as np
import math
from sklearn.metrics import mean_squared_error
import numpy.polynomial.polynomial as poly

import sys
sys.path.append('/home/local/RL-INSTITUT/inia.steinbach/Dokumente/greco_technologies_to_pvlib/CPV/')
from pvlib_CPVsystem import CPVSystem
import pvlib_CPVsystem as cpv



df= pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/INS/INS_Nov_twoaxistracking/INS_november_prepared_dataset.csv', sep=',', index_col=0)
# Converting the index as date

location= pvlib.location.Location(latitude=40.453,longitude=-3.727,
                                         tz=1, altitude=658)


module_params = {'gamma_ref' : 5.524, 'mu_gamma' : 0.003, 'I_L_ref' : 0.96,
                 'I_o_ref' : 0.00000000017, 'R_sh_ref' : 5226,
                 'R_sh_0': 21000, 'R_sh_exp' : 5.50, 'R_s' : 0.01,
                 'alpha_sc' : 0.00, 'EgRef' : 3.91, 'irrad_ref' : 1000,
                 'temp_ref' : 25, 'cells_in_series' : 12, 'eta_m' : 0.32,
                 'alpha_absorption' : 0.9}

csys = CPVSystem(module=None, module_parameters=module_params,
                     modules_per_string=1, strings_per_inverter=1,
                     inverter=None, inverter_parameters=None,
                     racking_model='freestanding',
                     losses_parameters=None, name=None)


celltemp = csys.pvsyst_celltemp(df['GNI'],                                      #todo: get wind speed
                                df['temp'],
                                df['wind_speed'])

(photocurrent, saturation_current, resistance_series,
 resistance_shunt, nNsVth) = (csys.calcparams_pvsyst(df['DNI'],
                                                     celltemp))

csys.diode_params = (photocurrent, saturation_current, resistance_series,
                     resistance_shunt, nNsVth)

csys.dc = csys.singlediode(photocurrent, saturation_current,
                           resistance_series,
                           resistance_shunt, nNsVth)

real_power = df['Pmp']
estimation = csys.dc['p_mp']



airmass = location.get_airmass(df.index)
relative_airmass= airmass['airmass_relative'].fillna(0)

# calculate single utilization factors

IscDNI_top = 0.96/1000

thld_am = 4.574231933073185
m_low_am = 3.906372068620377e-06
m_high_am = -3.0335768119184845e-05
thld_temp = 50
m_low_temp = 4.6781224141650075e-06
m_high_temp = 0




uf_am = []
for i, v in relative_airmass.items():
    uf_am.append(cpv.get_single_util_factor(v, thld_am,
                                            m_low_am/IscDNI_top,
                                            m_high_am/IscDNI_top))


uf_temp = []
for i, v in df['temp'].items():
    uf_temp.append(cpv.get_single_util_factor(v, thld_temp,
                                            m_low_temp/IscDNI_top,
                                              m_high_temp/IscDNI_top))


weight_am_final = 1.0
rmsd = 10000
rmsd_list=[]

for weight_am in np.arange(0, 1, 0.05):
    weight_temp = 1.0 - weight_am

    modeled_power = estimation * (np.multiply(weight_am, uf_am) +
                                  np.multiply(weight_temp, uf_temp))
    rmsd_temp = math.sqrt(mean_squared_error(real_power, modeled_power))
    rmsd_list.append(rmsd_temp)
    if rmsd_temp < rmsd:
        weight_am_final = weight_am
        weight_temp_final = weight_temp
        rmsd = rmsd_temp

print(weight_am_final, weight_temp_final)


modeled_power = estimation * (np.multiply(weight_am_final, uf_am) +
                                  np.multiply(weight_temp_final, uf_temp))

residualUF = modeled_power - real_power
residualwithoutUF= estimation - real_power


plt.plot(df['temp'], uf_temp, 'b.', label='UF(temp)')
plt.xlabel("Temperature in C")
plt.ylabel("Utilization Factor")
plt.legend()
plt.show()

plt.plot(relative_airmass, uf_am, 'r.', label='UF(AM)')
plt.xlabel("Airmass")
plt.ylabel("Utilization Factor")
plt.legend()
plt.show()




plt.plot(real_power, 'r--', label='real_power')
plt.plot(modeled_power, 'b', label='modeled_power')

plt.plot(estimation, 'g', label='estimation')
plt.xlabel("Time in Days")
plt.ylabel("Power in W")
plt.legend()
plt.show()

p1=poly.polyfit(real_power,modeled_power,1)

plt.plot(real_power,modeled_power, 'bo', markersize=1, label='modeled_power with UF over measured power')
plt.plot(real_power, estimation,'ro', markersize=1, label='modeled_power_without UF over measured power')
plt.plot(real_power, poly.polyval(real_power, p1), 'y-', label='model_power_fit')
plt.plot(real_power, real_power, 'g', label='real power')
plt.xlabel("Power in W")
plt.ylabel("Power in W")
plt.legend()
plt.show()


plt.plot(airmass['airmass_relative'], residualwithoutUF, 'go', markersize=1, label='Airmass residual without UF')
plt.plot(airmass['airmass_relative'].fillna(0), residualUF, 'ro', markersize=1, label='Airmass residual with UF')
plt.xlabel("Airmass")
plt.ylabel("Residual Pmpp in %")
plt.legend()
plt.show()


plt.plot(df['temp'], residualUF, 'ro', markersize=1, label='Temperature residual with UF')
plt.plot(df['temp'], residualwithoutUF, 'go', markersize=1, label='Temperature residual without UF')
plt.xlabel("Air Temperature in T")
plt.ylabel("Residual Pmpp in %")
plt.legend()
plt.show()