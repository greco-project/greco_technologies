import pandas as pd
import pvlib.atmosphere
from pvlib_CPVsystem import StaticCPVSystem
import matplotlib.pyplot as plt
import pvlib_CPVsystem as cpv
import numpy as np
import math
from sklearn.metrics import mean_squared_error
import numpy.polynomial.polynomial as poly


df= pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/INS/MAY/InsolightMay2019_filtered.csv', sep=',', index_col=0)
# Converting the index as date

location = pvlib.location.Location(latitude=40.453,longitude=-3.727,
                                         tz=1, altitude=658)

module_params = {'gamma_ref' : 5.524, 'mu_gamma' : 0.003, 'I_L_ref' : 0.96,
                 'I_o_ref' : 0.00000000017, 'R_sh_ref' : 5226,
                 'R_sh_0': 21000, 'R_sh_exp' : 5.50, 'R_s' : 0.01,
                 'alpha_sc' : 0.00, 'EgRef' : 3.91, 'irrad_ref' : 1000,
                 'temp_ref' : 25, 'cells_in_series' : 12, 'eta_m' : 0.32,
                 'alpha_absorption' : 0.9}

csys = StaticCPVSystem(surface_tilt=30, surface_azimuth=180,
                       module=None, module_parameters=module_params,
                     modules_per_string=1, strings_per_inverter=1,
                     inverter=None, inverter_parameters=None,
                     racking_model='insulated',
                     losses_parameters=None, name=None)


spa=pvlib.solarposition.spa_python(time=df.index, latitude=40.453,
                                   longitude=-3.72)


airmass = location.get_airmass(df.index)
relative_airmass= airmass['airmass_relative'].fillna(0)

#calculate AOI
aoi_list = pd.Series(name='aoi')
ot_list = pd.Series(name='ot')
gt_list = pd.Series(name='gt')
for index, row in spa.iterrows(): #todo: correct surface_tilt and surface_azimuth
    aoi = pvlib.irradiance.aoi(surface_tilt=30,
                               surface_azimuth=180,
                               solar_zenith=row['zenith'],
                               solar_azimuth=row['azimuth'])
    # calculate optical losses
    aoi_list[index] = aoi
    ot_list[index]=cpv.uf_aoi(aoi=aoi)
    gt_list[index]=cpv.glass_transmission_losses(aoi=aoi)

alignement_transmission = 0.95 #emperical parameter for Insolight module
df['aoi']=aoi_list
#weather_loc['glass_transmission']=gt_list
df['DII_new'] = df['dii'] * alignement_transmission *ot_list * gt_list



celltemp = csys.pvsyst_celltemp(df['gii'],
                                df['temp_air'],
                                df['wind_speed'])

(photocurrent, saturation_current, resistance_series,
 resistance_shunt, nNsVth) = (csys.calcparams_pvsyst(df['dii'],
                                                     celltemp))

csys.diode_params = (photocurrent, saturation_current, resistance_series,
                     resistance_shunt, nNsVth)

csys.dc = csys.singlediode(photocurrent, saturation_current,
                           resistance_series,
                           resistance_shunt, nNsVth)

real_power = df['Pmp']
estimation = csys.dc['p_mp']


# calculate single utilization factors

thld_am =  4.125860936553121
m_low_am =  0.0634016984325695
m_high_am =  -0.21442236571732923

thld_temp =  50
m_low_temp =  0.019546056846064516
m_high_temp =  0.0

thld_aoi =  62.74476733755362
m_low_aoi =  -0.0004954994407100744
m_high_aoi =  -0.02069421162410764



weight_am=0.85
weight_temp=0.15


IscDNI_top=0.96/1000

uf_am = []
for i, v in relative_airmass.items():
    uf_am.append(cpv.get_single_util_factor(v, thld_am,
                                            m_low_am,
                                            m_high_am))



uf_temp = []
for i, v in df['temp_air'].items():
    uf_temp.append(cpv.get_single_util_factor(v, thld_temp,
                                            m_low_temp,
                                              m_high_temp))


uf_aoi = []
for i,v in df['aoi'].items():
    uf_aoi.append(cpv.get_single_util_factor(aoi, thld_aoi, m_low_aoi,
                                    m_high_aoi))


uf_aoi_ast = cpv.get_single_util_factor(0, thld_aoi, m_low_aoi,
                                    m_high_aoi)

uf_aoi_norm = np.divide(uf_aoi, uf_aoi_ast)


uf_am_temp = np.multiply(weight_am, uf_am) + np.multiply(weight_temp, uf_temp)
UF_global = np.multiply(uf_am_temp, uf_aoi_norm)


modeled_power = estimation * uf_aoi

residualUF = modeled_power - real_power
residualwithoutUF= estimation - real_power

plt.plot(modeled_power, 'b', label='modeled power with UF_aoi')
plt.plot(real_power, 'r', label='real_power')
plt.plot(estimation, 'g', label='pvlib-calculated power')
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


plt.plot(df['temp_air'], residualUF, 'ro', markersize=1, label='Temperature residual with UF')
plt.plot(df['temp_air'], residualwithoutUF, 'go', markersize=1, label='Temperature residual without UF')
plt.xlabel("Air Temperature in T")
plt.ylabel("Residual Pmpp in %")
plt.legend()
plt.show()