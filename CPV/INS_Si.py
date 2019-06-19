
import pandas as pd

import pvlib.solarposition
import pvlib.atmosphere
import pvlib.irradiance

import pvlib_CPVsystem as cpv



def INS_Si_prepare_weather(weather_loc, lat, lon, surface_tilt, surface_azimuth):
    '''

    This function prepares the given weather data (in the form of PVLib) to the
     special case of INS_Si flatplate module. It includes galass transmission
     losses due to the fresnel lense and DNI-dependency of AOI:

    For AOI < 60° -> DNI=0, DHI is calculating using the perez-approach

    For API > 60° -> DHI is calculating using the perez-approach





    :param weather_loc: pd.DataFrame with columns defined for the PVLib
    :param lat: numeric, latitude
    :param lon: numeric, longitude
    :param surface_tilt: integer
    :param surface_azimuth: integer

    :return: pd.DataFrame that includes:
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
    extra = pvlib.irradiance.get_extra_radiation(times, solar_constant=1366.1,
                                      method='spencer', epoch_year=2014)
    extra = extra.to_frame('extra')
    zenith = spa_python[['zenith']].fillna(0)
    azimuth = spa_python[['azimuth']].fillna(0)
    airmass_rel = pvlib.atmosphere.get_relative_airmass(zenith)
    airmass_relative = pd.DataFrame(data=airmass_rel, index=times,
                               columns=['airmass']).fillna(0)
    #pressure = weather_loc[['P']]
    absulute_airmass=pvlib.atmosphere.get_absolute_airmass(airmass_relative, pressure=101325.) # todo: correct pressure

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

    weather_loc['dni-gt']= weather_loc['dni']*gt


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

    poa_components=pvlib.irradiance.poa_components(aoi= aoi_list, dni=dni, #todo: check if that is needed?
                                       poa_sky_diffuse=perez['perez'],
                                        poa_ground_diffuse=ground_diffuse)

    prepared_dict= pd.concat(
        [poa_components['poa_global'], poa_components['poa_direct'], poa_components['poa_diffuse'], aoi_list, absulute_airmass], axis=1)

    return prepared_dict
