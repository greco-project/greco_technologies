import matplotlib.pyplot as plt

import pvlib
import pandas as pd
import cpvlib.cpvlib
from greco_technologies.cpv.inputs import mod_params_cpv, mod_params_diffuse



def create_cpv_time_series(lat, lon, weather, surface_azimuth, surface_tilt):
    """

    creates a time series for a type of cpv module

    :param lat: num
        latitude
    :param lon: num
        longitude
    :param weather: pd.DataFrame()
        weather dataframe according to pvlib standards
    :param surface_azimuth: int
        surface azimuth
    :param surface_tilt: int
        surface tilt
    :param cpv_type: str
        possible cpv_types integrated up to this point: "ins", "m300"
    :return: pd.DataFrame()
    """

    location = pvlib.location.Location(
        latitude=lat, longitude=lon, altitude=695, tz='utc')

    weather.index = pd.to_datetime(weather.index)
    spa = pvlib.solarposition.spa_python(
        time=weather.index, latitude=lat, longitude=lon
    )
    weather["dni"] = pvlib.irradiance.dirint(
        weather["ghi"], solar_zenith=spa["zenith"], times=weather.index
    )


    solar_zenith = location.get_solarposition(weather.index).zenith
    solar_azimuth = location.get_solarposition(weather.index).azimuth

    #%%
    # StaticHybridSystem
    static_hybrid_sys = cpvlib.cpvlib.StaticHybridSystem(
        surface_tilt=surface_tilt,
        surface_azimuth=surface_azimuth,
        module_cpv=None,
        module_diffuse=None,
        module_parameters_cpv=mod_params_cpv,
        module_parameters_diffuse=mod_params_diffuse,
        modules_per_string=1,
        strings_per_inverter=1,
        inverter=None,
        inverter_parameters=None,
        racking_model="insulated",
        losses_parameters=None,
        name=None,
    )


    # get_effective_irradiance
    weather['dii_effective'], weather['poa_diffuse_static_effective'] = static_hybrid_sys.get_effective_irradiance(
        solar_zenith,
        solar_azimuth,
        iam_param=0.7,
        aoi_limit=55,
        dii=None,
        ghi=weather['ghi'],
        dhi=weather['dhi'],
        dni=weather['dni']
    )

    # pvsyst_celltemp
    weather['temp_cell_35'], weather['temp_cell_diffuse'] = static_hybrid_sys.pvsyst_celltemp(
        dii=weather['dii_effective'],
        poa_diffuse_static=weather['poa_diffuse_static_effective'],
        temp_air=weather['temp_air'],
        wind_speed=weather['wind_speed']
    )

    # calcparams_pvsyst
    diode_parameters_cpv, diode_parameters_diffuse = static_hybrid_sys.calcparams_pvsyst(
        dii=weather['dii_effective'],
        poa_diffuse_static=weather['poa_diffuse_static_effective'],
        temp_cell_cpv=weather['temp_cell_35'],
        temp_cell_diffuse=weather['temp_cell_diffuse'],
    )

    # singlediode
    dc_cpv, dc_diffuse = static_hybrid_sys.singlediode(
        diode_parameters_cpv, diode_parameters_diffuse)


    # uf_global (uf_am, uf_temp_air)
    weather['am'] = location.get_airmass(weather.index).airmass_absolute

    uf_cpv = static_hybrid_sys.get_global_utilization_factor_cpv(weather['am'], weather['temp_air'])

    # Power
    #dc_cpv.p_mp.plot(label="cpv")
    #dc_diffuse.p_mp.plot(label="diffuse")

    # Energy
    energy_cpv = (dc_cpv['p_mp'] * uf_cpv)
    energy_diffuse = dc_diffuse['p_mp']
    total=energy_cpv +energy_diffuse

    return total
