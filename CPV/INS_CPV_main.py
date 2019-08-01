import pandas as pd
import os
from pvlib.pvsystem import PVSystem
from pvlib.location import Location


import pvlib.atmosphere
import pvlib.pvsystem as pvsystem
import pvlib.irradiance as irrad

import pvlib_CPVsystem as cpv
from pvlib_CPVsystem import StaticCPVSystem
import visualizing_data
import matplotlib.pyplot as plt



#GENERAL PARAMETERS
surface_tilt= 20
surface_azimuth= 180


#LOAD WEATHER DATA
filename = os.path.abspath(
            "/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/163_Open_FRED/03-Projektinhalte/AP2 Wetterdaten/open_FRED_TestWetterdaten_csv/fred_data_test_2016.csv")
weather_df = pd.read_csv(filename, skiprows=range(1, 50), nrows=(5000),
                                 index_col=0,
                                 date_parser=lambda idx: pd.to_datetime(idx,
                                                                        utc=True))
weather_df.index = pd.to_datetime(weather_df.index).tz_convert(
            'Europe/Berlin')
# calculate ghi
weather_df['ghi'] = weather_df.dirhi + weather_df.dhi
#weather = weather_df[['wind_speed', 'temp_air', 'P', 'dhi', 'dirhi', 'ghi']]
coordinates = weather_df.loc[:, ["lat", "lon"]].values


for lat, lon in coordinates:

    location = Location(latitude=lat, longitude=lon)

    # FILTER WEATHER
    weather_loc = weather_df.loc[(weather_df['lat'] == lat)  # kann das weg?
                                 & (weather_df['lon'] == lon)]
    location = Location(latitude=lat, longitude=lon)
    times = weather_loc.index  # kann das weg?


    #set up CPV-module

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

    celltemp = csys.pvsyst_celltemp(weather_loc['ghi'],
                                    weather_loc['temp_air'],
                                    weather_loc['wind_speed'])

#    DNIseries=weather_loc['dni']
#    DNI=DNIseries.to_numpy()


    # calculate airmass
    airmass = location.get_airmass(times)
    relative_airmass= airmass['airmass_relative'].fillna(0)

    # calculate aoi, optical transmission losses and glass transmission losses
    spa_python = pvlib.solarposition.spa_python(time=times, latitude=lat,
                                                longitude=lon)
    aoi_list = pd.Series(name='aoi')
    ot_list = pd.Series(name='ot')
    gt_list = pd.Series(name='gt')
    for index, row in spa_python.iterrows():
        aoi = pvlib.irradiance.aoi(surface_tilt=surface_tilt,
                                   surface_azimuth=surface_azimuth,
                                   solar_zenith=row['zenith'],
                                   solar_azimuth=row['azimuth'])
        # calculate optical losses
        aoi_list[index] = aoi
        ot_list[index]=cpv.optical_transmission_losses(aoi=aoi)
        gt_list[index]=cpv.glass_transmission_losses(aoi=aoi)

    alignement_transmission = 0.95 #emperical parameter for Insolight module
    weather_loc['aoi']=aoi_list
    weather_loc['glass_transmission']=gt_list
    weather_loc['dni_after_losses'] = weather_loc[
                        'dni'] * ot_list * gt_list * alignement_transmission

    # calculate DNI after the utilization factor
    weather_loc['dni_uf']= cpv.calculate_utilization_factor(
                               am=relative_airmass,
                               t_ambient=weather_loc['temp_air'],
                               dni=weather_loc['dni_after_losses'],
                                calculate_ufdni=False) #todo: check coefficients


    (photocurrent, saturation_current, resistance_series,
     resistance_shunt, nNsVth) = (csys.calcparams_pvsyst(weather_loc['dni_uf'],
                                                         celltemp))

    csys.diode_params = (photocurrent, saturation_current, resistance_series,
                         resistance_shunt, nNsVth)

    csys.dc = csys.singlediode(photocurrent, saturation_current,
                               resistance_series,
                               resistance_shunt, nNsVth)


    plt.plot(csys.dc)
    plt.show()
    break



