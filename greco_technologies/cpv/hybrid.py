import greco_technologies.cpv.cpv as cpv
import cpvtopvlib.cpvsystem as cpvsystem
import pandas as pd
import os
from pvlib.temperature import pvsyst_cell, TEMPERATURE_MODEL_PARAMETERS
import matplotlib.pyplot as plt
import pvlib
import pvlib.pvsystem as pvsystem


def create_hybrid_time_series(lat, lon, weather, surface_tilt, surface_azimuth):

    """
    This function adds up the output power of the cpv and the backplane si
    module.
    :param lat: num
        latitude
    :param lon: num
    :param weather:
    :param surface_tilt:
    :param surface_azimuth:
    :return:
    """

    power_cpv = cpv.create_cpv_time_series(
        lat=lat,
        lon=lon,
        weather=weather,
        surface_azimuth=surface_azimuth,
        surface_tilt=surface_tilt,
        type="ins",
    )
    power_si = create_si_time_series(
        lat=lat,
        lon=lon,
        weather=weather,
        surface_azimuth=surface_azimuth,
        surface_tilt=surface_tilt,
    )
    hybrid_power = power_cpv + power_si.fillna(0)
    return hybrid_power


def hybrid_weather_data(weather_loc, lat, lon, surface_tilt, surface_azimuth):
    """

    This function prepares the given weather data (in the form of PVLib) to the
     special case of INS_Si flatplate module. It includes galass transmission
     losses due to the fresnel lense and DNI-dependency of AOI:

    For AOI < 60° -> DNI=0, DHI is calculated using the perez-approach

    For API > 60° -> DHI is calculated using the perez-approach

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
            absolute_airmass
    """

    times = weather_loc.index
    spa_python = pvlib.solarposition.spa_python(time=times, latitude=lat, longitude=lon)
    # calculate extraterrestrial radiation for diffuse model
    extra = pvlib.irradiance.get_extra_radiation(weather_loc.index).to_frame("extra")
    zenith = spa_python[["zenith"]].fillna(0)
    azimuth = spa_python[["azimuth"]].fillna(0)
    # calculate relative and absolute airmass
    airmass_relative = pvlib.atmosphere.get_relative_airmass(zenith)
    airmass_relative = pd.DataFrame(
        data=airmass_relative, index=times, columns=["airmass"]
    ).fillna(0)
    #    pressure = weather_loc[['P']]
    absulute_airmass = pvlib.atmosphere.get_absolute_airmass(
        airmass_relative, pressure=101325.0
    )  # todo: correct pressure
    parameters_diffuse = pd.concat(
        [
            zenith,
            extra,
            azimuth,
            airmass_relative,
            weather_loc[["dni"]],
            weather_loc[["dhi"]],
            weather_loc[["ghi"]],
        ],
        axis=1,
    )

    # calculate DHI with perez_diffuse approach
    pz = pd.Series()
    for index, row in parameters_diffuse.iterrows():

        perez_diffuse = pvlib.irradiance.perez(
            surface_tilt,
            surface_azimuth,
            solar_zenith=row["zenith"],
            solar_azimuth=row["azimuth"],
            dni=row["dni"],
            dhi=row["dhi"],
            dni_extra=row["extra"],
            airmass=row["airmass"],
        )
        pz[index] = perez_diffuse

    perez = pd.DataFrame(pz, columns=["perez"])
    weather_loc["perez_diffuse"] = perez
    # calculate DNI
    # prepare weather data for AOI >< 60° and change DNI values accordingly
    aoi_list = pd.Series(name="aoi")
    gt_list = pd.Series(name="gt")
    for index, row in spa_python.iterrows():
        aoi = pvlib.irradiance.aoi(
            surface_tilt=surface_tilt,
            surface_azimuth=surface_azimuth,
            solar_zenith=row["zenith"],
            solar_azimuth=row["azimuth"],
        )
        # calculate optical losses
        CPVSystem = cpvsystem.CPVSystem()
        gt = CPVSystem.glass_transmission_losses(aoi=aoi)

        aoi_list[index] = aoi
        gt_list[index] = gt

    # add alignement losses and glass transmission losses to DNI and perez_diffuse
    alignement_transmission = 0.95
    weather_loc["dni"] = weather_loc["dni"] * gt_list * alignement_transmission
    weather_loc["perez_diffuse"] = (
        weather_loc["perez_diffuse"] * gt * alignement_transmission
    )

    smallaoi = aoi_list[aoi_list < 60]
    bigaoi = aoi_list[aoi_list > 60]

    weather_smallaoi = weather_loc[weather_loc.index.isin(smallaoi.index)]
    # set DNI= 0 for AOI < 60°
    weather_smallaoi["dni"] = 0
    weather_bigaoi = weather_loc[
        weather_loc.index.isin(bigaoi.index)
    ]                                                                           # todo: adjust DNI for AOI>60°: Durchnitt über Fläche
    calculated_weather = pd.concat([weather_smallaoi, weather_bigaoi])
    calculated_weather.sort_index(inplace=True)

    calculated_weather["aoi"] = aoi_list
    calculated_weather["airmass"] = absulute_airmass

    # calculate ground diffuse                                                  todo: is ground diffuse relevant here?
    ground_diffuse = pvlib.irradiance.get_ground_diffuse(
        surface_tilt=surface_tilt,
        ghi=calculated_weather["ghi"],
        albedo=0.25,
        surface_type=None,
    )
    calculated_weather["dhi"] = calculated_weather["perez_diffuse"] + ground_diffuse
    calculated_weather.fillna(0, inplace=True)

    return calculated_weather[[
            "wind_speed",
            "temp_air",
            "ghi",
            "dhi",
            "dni",
            "aoi",
            "airmass"]]


def create_si_time_series(lat, lon, weather, surface_azimuth, surface_tilt):

    #load example module from sandia library
    sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
    module = sandia_modules['Canadian_Solar_CS5P_220M___2009_']
    # prepare dictionary with poa_global, poa_direct_poa_diffuse, absolute_airmass, aoi
    hybrid_weather = hybrid_weather_data(
        weather,
        lat=lat,
        lon=lon,
        surface_tilt=surface_tilt,
        surface_azimuth=surface_azimuth,
    )

    # todo: adapt function for effective irradiance to increase diffuse fraction
    # todo: adjust sapm_spectral_loss(airmass_absolute, module) and sapm_aoi_loss(aoi, module) that are included in the effective_irradiance

    effective_irradiance = pvsystem.sapm_effective_irradiance(
        poa_direct=hybrid_weather["dni"],
        poa_diffuse=hybrid_weather["dhi"],
        airmass_absolute=hybrid_weather["airmass"],
        aoi=hybrid_weather["aoi"],
        module=module,
    )

    temp_params = TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_glass']

    temp_cell = pvlib.temperature.sapm_cell(
        poa_global=hybrid_weather["ghi"],
        wind_speed=hybrid_weather["wind_speed"],
        temp_air=hybrid_weather["temp_air"],
        **temp_params
    ).fillna(0)

    output = pvsystem.sapm(
        effective_irradiance=effective_irradiance,
        temp_cell=temp_cell,
        module=module,
    )

    return output.p_mp


if __name__ == "__main__":

    #load example weather data
    filename = os.path.abspath(
        "/home/adminlocal/Dokumente/greco_env/pvcompare/pvcompare/data/inputs/weatherdata.csv"
    )
    weather_df = pd.read_csv(
        filename, index_col=0, date_parser=lambda idx: pd.to_datetime(idx, utc=True)
    )
    weather_df.index = pd.to_datetime(weather_df.index).tz_convert("Europe/Berlin")
    weather_df["dni"] = weather_df["ghi"] - weather_df["dhi"]

    ds = create_si_time_series(
            lat=52.11113,
            lon=12.48062,
            weather=weather_df,
            surface_azimuth=180,
            surface_tilt=30,
        )
    plt.plot(ds)

    # dh = create_hybrid_time_series(
    #     lat=40.3, lon=5.4, weather=weather_df, surface_tilt=25,
    #     surface_azimuth=180)
    # plt.plot(dh, "-b", alpha=0.7)


    plt.show()