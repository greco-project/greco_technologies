
import pandas as pd
import os
from pvlib.pvsystem import PVSystem
from pvlib.location import Location
import matplotlib.pyplot as plt

import pvlib.atmosphere
import pvlib.pvsystem as pvsystem

import INS_Si
import visualizing_data


df= pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/INS/MAY/InsolightMay2019_filtered.csv', sep=',', index_col=0)
df.index= pd.to_datetime(df.index)

df['dhi']=df['ghi']-df['dni']

sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
module_ref= sandia_modules['Canadian_Solar_CS5P_220M___2009_']
module_ref['Cells_in_Series']=1
module_ref['Cells_in_Parallel']=1


system_ref = PVSystem(surface_tilt=30, surface_azimuth=180,
                  module_parameters=module_ref,
                  inverter_parameters=None)

# prepare dictionary with poa_global, poa_direct_poa_diffuse, absolute_airmass, aoi
prepared_poas= INS_Si.INS_Si_prepare_weather(df, lat=45.641603, lon=-3.727, surface_tilt=30, surface_azimuth=180)


#todo: adapt function for effective irradiance to increase diffuse fraction
# todo: adjust sapm_spectral_loss(airmass_absolute, module) and sapm_aoi_loss(aoi, module) that are included in the effective_irradiance


effective_irradiance = pvsystem.sapm_effective_irradiance(poa_direct=0,
                                                    poa_diffuse=prepared_poas['poa_diffuse'],
                                                    airmass_absolute = prepared_poas['airmass'],
                                                    aoi= prepared_poas['aoi'],
                                                    module=module_ref)

df['effective_irrandiance']=effective_irradiance
temp_cell= pvsystem.sapm_celltemp(poa_global=prepared_poas['poa_global'],
                                  wind_speed=df['wind_speed'],
                                  temp_air=df['temp_air'],
                                  model='open_rack_cell_glassback')             # todo: adjust model
temp_cell=temp_cell.fillna(0)
output= pvsystem.sapm(effective_irradiance=effective_irradiance,
                 temp_cell=temp_cell['temp_cell'],
                 module=module_ref)

modeled_power=output.p_mp
modeled_isc=output.i_sc

real_power=df['Pmp_Si']
real_isc=df['Isc_Si']

plt.plot(modeled_power, 'b')
plt.plot(real_power, 'r')
plt.show()

plt.plot(modeled_isc, 'b')
plt.plot(real_isc, 'r')
plt.show()