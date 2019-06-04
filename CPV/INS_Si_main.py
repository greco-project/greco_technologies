

import pandas as pd
import os
from pvlib.pvsystem import PVSystem
from pvlib.location import Location


import pvlib.atmosphere
import pvlib.pvsystem as pvsystem

import INS_Si
#import flatplate_SAPM_effective_irradiance as eff
import visualizing_data

# this scipt loads weather data and prepares it regarding the special setup of
# the Insolight-SI module for AOI <> 60°



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

    #set up module
    cpv_module = pd.read_csv(
        '/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/Parameters/CPV_parameters2.csv',
        sep=',', encoding='utf-7', converters={"Name": str, "CPV": float})
    sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
    sandia_module = sandia_modules['Canadian_Solar_CS5P_220M___2009_']
    cec_inverters = pvlib.pvsystem.retrieve_sam('cecinverter')
    inverter_parameters = cec_inverters[
        'ABB__MICRO_0_25_I_OUTD_US_208_208V__CEC_2014_']
    df = cpv_module
    df = df.set_index('Name')
    df = df['CPV']
    system = PVSystem(surface_tilt=20, surface_azimuth=200,
                      module_parameters=sandia_module,
                      inverter_parameters=inverter_parameters)
    location = Location(latitude=lat, longitude=lon)

    weather_loc = weather_df.loc[(weather_df['lat'] == lat) # kann das weg?
                                  & (weather_df['lon'] == lon)]
    location=Location(latitude=lat, longitude=lon)
    times=weather_loc.index #kann das weg?

    # prepare dictionary with poa_global, poa_direct_poa_diffuse, absolute_airmass, aoi
    prepared_poas= INS_Si.INS_Si_prepare_weather(weather_loc, lat, lon, surface_tilt, surface_azimuth)


    #todo: adapt function for effective irradiance to increase diffuse fraction
    # todo: adjust sapm_spectral_loss(airmass_absolute, module) and sapm_aoi_loss(aoi, module) that are included in the effective_irradiance
    effective_irradiance = pvsystem.sapm_effective_irradiance(poa_direct=prepared_poas['poa_direct'],
                                                        poa_diffuse=prepared_poas['poa_diffuse'],
                                                        airmass_absolute = prepared_poas['airmass'],
                                                        aoi= prepared_poas['aoi'],
                                                        module=sandia_module,
                                                        reference_irradiance=1000)
    temp_cell= pvsystem.sapm_celltemp(poa_global=prepared_poas['poa_global'],
                                      wind_speed=weather_loc['wind_speed'],
                                      temp_air=weather_loc['temp_air'],
                                      model='open_rack_cell_glassback') # todo: adjust model

    output= pvsystem.sapm(effective_irradiance=effective_irradiance,
                     temp_cell=temp_cell['temp_cell'],
                     module=sandia_module)


    visualizing_data.plot_sapm(output, effective_irradiance)
    break
#print(output)