import pandas as pd
import os
from pvlib.pvsystem import PVSystem
from pvlib.location import Location


import pvlib.atmosphere
import pvlib.pvsystem as pvsystem
import pvlib.irradiance as irrad

import pvlib_CPVsystem as cpv
import visualizing_data



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

    #set up CPV-module
    cpv_module = pd.read_csv(
        '/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/Parameters/CPV_parameters2.csv',
        sep=',', encoding='utf-7', converters={"Name": str, "CPV": float})
    sandia_modules = pvsystem.retrieve_sam('SandiaMod')
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


    #FILTER WEATHER
    weather_loc = weather_df.loc[(weather_df['lat'] == lat) # kann das weg?
                                  & (weather_df['lon'] == lon)]
    location=Location(latitude=lat, longitude=lon)
    times=weather_loc.index #kann das weg?

    # CALCULATE AIRMASS
    spa_python = pvlib.solarposition.spa_python(time=times, latitude=lat,
                                                longitude=lon)
    zenith = spa_python[['zenith']].fillna(0)
    airmass_rel = pvlib.atmosphere.get_relative_airmass(zenith)
    airmass_relative = pd.DataFrame(data=airmass_rel, index=times,
                               columns=['airmass']).fillna(0)
    absolute_airmass = pvlib.atmosphere.get_absolute_airmass(airmass_relative,
                                                             pressure=101325.)

    # CALCULATE AOI AND OPTICAL TRANSMISSION LOSSES RELATED TO AOI AND OPTICAL
    # TRANSMISSION LOSSES
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

    alignement_transmission = 0.95
    weather_loc['dni_losses']= weather_loc['dni']*ot_list*gt_list * alignement_transmission
    weather_loc['aoi']=aoi_list
    weather_loc['glass_transmission']=gt_list

    # CALCULATE THE DNI AFTER THE UTILIZATION FACTOR
#    weather_loc['dni_uf']= cpv.UF_corrected_DNI(am=absolute_airmass, t_ambient=weather_loc['temp_air'], dni=weather_loc['dni_losses']) #todo: check coefficients


    #CALCULATE POA COMPONENTS
    prepared_poas= irrad.poa_components(aoi=weather_loc['aoi'], dni=weather_loc['dni_losses'], poa_sky_diffuse=0, poa_ground_diffuse=0)

    effective_irradiance = pvsystem.sapm_effective_irradiance(poa_direct=prepared_poas['poa_direct'],
                                                        poa_diffuse=prepared_poas['poa_diffuse'],
                                                        airmass_absolute = absolute_airmass['airmass'],
                                                        aoi= weather_loc['aoi'],
                                                        module=sandia_module,
                                                        reference_irradiance=1000)

    temp_cell= pvsystem.sapm_celltemp(poa_global=prepared_poas['poa_global'],
                                      wind_speed=weather_loc['wind_speed'],
                                      temp_air=weather_loc['temp_air'],
                                      model='open_rack_cell_glassback') # todo: adjust model
    # CALCULATE OUTPUT
    output= pvsystem.sapm(effective_irradiance=effective_irradiance,
                     temp_cell=temp_cell['temp_cell'],
                     module=sandia_module)


    visualizing_data.plot_sapm(output, effective_irradiance)
    break



