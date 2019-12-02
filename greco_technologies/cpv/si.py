
import pandas as pd

import pvlib.solarposition
import pvlib.atmosphere
import pvlib.irradiance
from pvlib.pvsystem import PVSystem
import pvlib.pvsystem as pvsystem
import os

import cpvsystem as cpv
import matplotlib.pyplot as plt


def prepare_weather(weather_loc, lat, lon, surface_tilt, surface_azimuth):
    '''

    This function prepares the given weather data (in the form of PVLib) to the
     special case of INS_Si flatplate module. It includes galass transmission
     losses due to the fresnel lense and DNI-dependency of AOI:

    For AOI < 60° -> DNI=0, DHI is calculating using the perez-approach

    For API > 60° -> DHI is calculating using the perez-approach

    Parameters
    -----------
    weather_loc: pd.DataFrame
        with columns defined for PVLib
    lat: float
        latitude
    lon: float
        longitude
    surface_tilt: int
    surface_azimuth: integer

    Returns
    ---------
    pd.DataFrame
        the dataframe includes:
            poa_global
            poa_direct
            poa_diffuse
            aoi
            absulute_airmass
    '''



    times = weather_loc.index
    spa_python = pvlib.solarposition.spa_python(time=times, latitude=lat,
                                                longitude=lon)
    #calculate extraterrestrial radiation for diffuse model
    extra = pvlib.irradiance.get_extra_radiation(weather_loc.index)
    extra = extra.to_frame('extra')
    zenith = spa_python[['zenith']].fillna(0)
    azimuth = spa_python[['azimuth']].fillna(0)
    airmass_rel = pvlib.atmosphere.get_relative_airmass(zenith)
    airmass_relative = pd.DataFrame(data=airmass_rel, index=times,
                               columns=['airmass']).fillna(0)
    #pressure = weather_loc[['P']]
    absulute_airmass=pvlib.atmosphere.get_absolute_airmass(airmass_relative,
                                                           pressure=101325.) # todo: correct pressure

    parameters_diffuse = pd.concat(
        [zenith, extra, azimuth, airmass_relative, weather_loc[['dni']],
         weather_loc[['dhi']], weather_loc[['ghi']]], axis=1)

    # calculate DHI
    pz = pd.Series()
    for index, row in parameters_diffuse.iterrows():

        perez_diffuse = pvlib.irradiance.perez(surface_tilt, surface_azimuth,
                                    solar_zenith=row['zenith'],
                                    solar_azimuth=row['azimuth'],
                                    dni=row['dni'],
                                    dhi=row['dhi'],
                                    dni_extra=row['extra'],
                                    airmass=row['airmass'])
        pz[index] = perez_diffuse

    perez = pd.DataFrame(pz, columns=["perez"])
    weather_loc['perez_diffuse']=perez

    # calculate DNI
    # prepare weather data for AOI >< 60° and change DNI values accordingly
    aoi_list = pd.Series(name='aoi')
    gt_list = pd.Series(name='gt')
    for index, row in spa_python.iterrows():
        aoi = pvlib.irradiance.aoi(surface_tilt=surface_tilt,
                                   surface_azimuth=surface_azimuth,
                                   solar_zenith=row['zenith'],
                                   solar_azimuth=row['azimuth'])
        # calculate optical losses
        gt=cpv.glass_transmission_losses(aoi)

        aoi_list[index] = aoi
        gt_list[index] = gt

    alignement_transmission = 0.95
    weather_loc['dni']= weather_loc['dni']*gt*alignement_transmission
    weather_loc['perez_diffuse'] = weather_loc['perez_diffuse'] * \
                                   gt*alignement_transmission

    smallaoi = aoi_list[aoi_list < 60]
    bigaoi = aoi_list[aoi_list > 60]

    weather_smallaoi = weather_loc[weather_loc.index.isin(smallaoi.index)]
    weather_smallaoi['dni'] = 0
    weather_bigaoi = weather_loc[weather_loc.index.isin(bigaoi.index)] #todo: adjust DNI for AOI>60°: Durchnitt über Fläche
    weather_data=pd.concat([weather_smallaoi, weather_bigaoi])
    weather_data.sort_index(inplace=True)

    dni = weather_data['dni']
    ground_diffuse = pvlib.irradiance.get_ground_diffuse(surface_tilt=surface_tilt,
                                        ghi=weather_data['ghi'],
                                        albedo=.25, surface_type=None)

    poa_components=pvlib.irradiance.poa_components(aoi= aoi_list, dni=dni,        #todo: check if that is needed? DNI set to 0!!
                                poa_sky_diffuse=weather_loc['perez_diffuse'],
                                poa_ground_diffuse=ground_diffuse)
    poa_components['poa_diffuse']=weather_loc['perez_diffuse']
    prepared_dict= pd.concat(
        [poa_components['poa_global'], poa_components['poa_direct'],
         poa_components['poa_diffuse'], aoi_list, absulute_airmass], axis=1)

    plt.plot(aoi_list, dni, 'b.', label='DNI')
    plt.plot(aoi_list, weather_loc['perez_diffuse'], 'g.', label='DHI')
    plt.xlabel('AOI')
    plt.ylabel('Irradiance')
    plt.legend()
    plt.show()
    return prepared_dict

def create_si_timeseries(lat, lon, weather, surface_azimuth, surface_tilt):


    sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
    module_ref = sandia_modules['Canadian_Solar_CS5P_220M___2009_']
    module_ref['Cells_in_Series'] = 4
    module_ref['Cells_in_Parallel'] = 1
    module_ref['Isco'] = 1
    module_ref['Voco'] = 4
    module_ref['Impo'] = 6.52
    module_ref['Vmpo'] = 5

    system_ref = PVSystem(surface_tilt=surface_tilt,
                          surface_azimuth=surface_azimuth,
                          module_parameters=module_ref,
                          inverter_parameters=None)

    # prepare dictionary with poa_global, poa_direct_poa_diffuse, absolute_airmass, aoi
    prepared_poas = prepare_weather(weather, lat=lat,
                                        lon=lon, surface_tilt=surface_tilt,
                                        surface_azimuth=surface_azimuth)

    # todo: adapt function for effective irradiance to increase diffuse fraction
    # todo: adjust sapm_spectral_loss(airmass_absolute, module) and sapm_aoi_loss(aoi, module) that are included in the effective_irradiance

    effective_irradiance = pvsystem.sapm_effective_irradiance(
        poa_direct=prepared_poas['poa_direct'],
        poa_diffuse=prepared_poas['poa_diffuse'],
        airmass_absolute=prepared_poas['airmass'],
        aoi=prepared_poas['aoi'],
        module=module_ref)

    weather['effective_irrandiance'] = effective_irradiance
    temp_cell = pvsystem.sapm_celltemp(poa_global=prepared_poas['poa_global'],
                                       wind_speed=weather['wind_speed'],
                                       temp_air=weather['temp_air'],
                                       model='open_rack_cell_glassback')  # todo: adjust model
    temp_cell = temp_cell.fillna(0)
    output = pvsystem.sapm(effective_irradiance=effective_irradiance,
                           temp_cell=temp_cell['temp_cell'],
                           module=module_ref)

    return output.p_mp

if __name__ == '__main__':

    filename = os.path.abspath(
        "/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/163_Open_FRED/03-Projektinhalte/AP2 Wetterdaten/open_FRED_TestWetterdaten_csv/fred_data_test_2016.csv")
    weather_df = pd.read_csv(filename, skiprows=range(1, 50), nrows=(5000),
                             index_col=0,
                             date_parser=lambda idx: pd.to_datetime(idx,
                                                                    utc=True))
    weather_df.index = pd.to_datetime(weather_df.index).tz_convert(
        'Europe/Berlin')


    coordinates = weather_df.loc[:, ["lat", "lon"]].values

    for lat, lon in coordinates:
        weather_loc = weather_df.loc[
            (weather_df['lat'] == lat)  # kann das weg?
            & (weather_df['lon'] == lon)]

        create_si_timeseries(lat=52.11113, lon=12.48062, weather=weather_loc,
                             surface_azimuth=180, surface_tilt=30)
        break