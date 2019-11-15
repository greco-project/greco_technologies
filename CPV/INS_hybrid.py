import INS_CPV as cpv
import INS_Si as si


def create_hybrid_timeseries(lat, lon, weather, surface_tilt, surface_azimuth):

    """

    This function counts together the output power of the cpv and the backplane si module.
    :param lat:
    :param lon:
    :param weather:
    :param surface_tilt:
    :param surface_azimuth:
    :return:
    """

    power_cpv=cpv.create_cpv_timeseries(lat=lat, lon=lon, weather=weather,
                                        surface_azimuth=surface_azimuth,
                                        surface_tilt=surface_tilt)
    power_si=si.create_si_timeseries(lat=lat, lon=lon, weather=weather,
                                        surface_azimuth=surface_azimuth,
                                        surface_tilt=surface_tilt)

    return power_cpv + power_si

