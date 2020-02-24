import pandas as pd
import numpy as np
import math
import pvlib.atmosphere
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import datetime
import numpy.polynomial.polynomial as poly

from cpvtopvlib import cpvsystem as cpv



df= pd.read_csv('../inputs/m300_data_filtered.txt', header=None)
columnnames= pd.read_csv('../inputs/m300_measurements_headers.csv', sep=',', dtype={'Daytime':str})
n=columnnames.columns.tolist()
df.columns = n

datetimestring = np.genfromtxt(
        '../inputs/m300_datetime.txt',
        dtype='str', delimiter = '\n')

datetimeobject = []
for i in range(len(datetimestring)):
    datetimeobject.append(datetime.datetime.strptime(datetimestring[i],
                                                     '%d-%b-%Y %H:%M:%S'))
df['Datetime'] = pd.to_datetime(datetimeobject)
df = df.set_index(pd.DatetimeIndex(df['Datetime']))

location = pvlib.location.Location(latitude=45.641603,longitude=5.875387,
                                         tz=1, altitude=234)

module_params = {'gamma_ref': 4.456, 'mu_gamma': 0.0012, 'I_L_ref': 3.346,
                 'I_o_ref': 0.000000000004, 'R_sh_ref': 4400,
                 'R_sh_0': 17500, 'R_sh_exp': 5.50, 'R_s': 0.736,
                 'alpha_sc': 0.00, 'irrad_ref': 1000, 'temp_ref': 25,
                 'cells_in_series': 42}

csys = cpv.StaticCPVSystem(module=None, module_parameters=module_params,
                     modules_per_string=1, strings_per_inverter=1,
                     inverter=None, inverter_parameters=None,
                     racking_model='freestanding',
                     losses_parameters=None, name=None)

celltemp = csys.pvsyst_celltemp(df['GNI'],
                                df['AirTemperature'],
                                df['WindSpeedAvg'])

(photocurrent, saturation_current, resistance_series,
 resistance_shunt, nNsVth) = (csys.calcparams_pvsyst(df['DNI'],
                                                     celltemp))

csys.diode_params = (photocurrent, saturation_current, resistance_series,
                     resistance_shunt, nNsVth)

csys.dc = csys.singlediode(photocurrent, saturation_current,
                           resistance_series,
                           resistance_shunt, nNsVth)

real_power = df['Pmpp']
estimation = csys.dc['p_mp']



airmass = location.get_airmass(df['Datetime'])
relative_airmass= airmass['airmass_relative'].fillna(0)

# calculate single utilization factors

thld_am = 2.022411098853249
m_low_am = 0.0423037910485609
m_high_am = -0.0210539236615148

uf_am = []
for i, v in relative_airmass.items():
    uf_am.append(cpv.get_simple_util_factor(v, thld_am,
                                            m_low_am, m_high_am))


thld_temp = 200
m_low_temp = 0.000923828521724516
m_high_temp = 0.0

uf_temp = []
for i, v in df['AirTemperature'].items():
    uf_temp.append(cpv.get_simple_util_factor(v, thld_temp,
                                            m_low_temp, m_high_temp))


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

modeled_power = estimation * (np.multiply(weight_am_final, uf_am) +
                                  np.multiply(weight_temp_final, uf_temp))

residualUF = ((real_power - modeled_power) / real_power)*100
residualwithoutUF= ((real_power - estimation) /real_power) *100

plt.plot(real_power, 'r', label='real_power')
#plt.plot(estimation, 'g', label='estimation')
plt.plot(modeled_power, 'b', label='modeled_power')
plt.xlabel("Time in Days")
plt.ylabel("Power in W")
plt.legend(loc='upper right')
plt.show()

p1=poly.polyfit(real_power,modeled_power,1)

plt.plot(real_power,modeled_power, 'bo', markersize=1, label='with UF')
plt.plot(real_power, estimation,'r+', markersize=1, label='without UF')
plt.plot(real_power, poly.polyval(real_power, p1), 'y-', label='modeled power fit')
plt.plot(real_power, real_power, 'g', label='optimal match')
plt.xlabel("Modeled power in W")
plt.ylabel("Real power in W")
plt.legend(loc='upper right')
plt.show()


plt.plot(airmass['airmass_relative'], residualwithoutUF, 'go', markersize=1, label='Airmass residual without UF')
plt.plot(airmass['airmass_relative'].fillna(0), residualUF, 'r+', markersize=1, label='Airmass residual with UF')
plt.xlabel("Airmass")
plt.ylabel("Residual Pmpp")
plt.legend(loc='upper right')
plt.show()



plt.plot(df['AirTemperature'], residualwithoutUF, 'go', markersize=1, label='Temperature residual without UF')
plt.plot(df['AirTemperature'], residualUF, 'r+', markersize=1, label='Temperature residual with UF')
plt.xlabel("Air Temperature in T")
plt.ylabel("Residual Pmpp")
plt.legend(loc='upper right')
plt.show()