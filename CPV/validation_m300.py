import pandas as pd
import numpy as np
import math

import pvlib.atmosphere

import pvlib_CPVsystem as cpv
from pvlib_CPVsystem import StaticCPVSystem
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import datetime

df= pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/Dokumente/CPV/data/m300_data_filtered.txt', header=None)
columnnames= pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/datasheets_commercial_CPV/m300_measurements_headers.csv', sep=',', dtype={'Daytime':str})
n=columnnames.columns.tolist()
df.columns = n

datetimestring = np.genfromtxt(
        '/home/local/RL-INSTITUT/inia.steinbach/Dokumente/CPV/data/m300_datetime.txt',
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

csys = StaticCPVSystem(module=None, module_parameters=module_params,
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

real_power = df['Isc']
estimation = csys.dc['p_mp']



airmass = location.get_airmass(df['Datetime'])
relative_airmass= airmass['airmass_relative'].fillna(0)

uf_am=cpv.ufam(relative_airmass)
uf_tamb=cpv.ufam(df['AirTemperature'])

weight_am_final = 1.0
rmsd = 10000

for weight_am in np.arange(0, 5, 0.1):
    for weight_tamb in np.arange(0,5,0.1):

        modeled_power = estimation * (np.multiply(weight_am, uf_am) +
                                  np.multiply(weight_tamb, uf_tamb))
        rmsd_temp = math.sqrt(mean_squared_error(real_power, modeled_power))

        if rmsd_temp < rmsd:
            weight_am_final = weight_am
            weight_tamb_final=weight_tamb
            rmsd = rmsd_temp

modeled_power = estimation * (np.multiply(weight_am_final, uf_am) +
                                  np.multiply(weight_tamb_final, uf_tamb))


plt.plot(real_power, 'r', label='real_power')
plt.plot(modeled_power, 'b', label='modeled_power')
plt.legend()
plt.show()