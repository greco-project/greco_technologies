import greco_technologies.cpv.cpv as cpv
import greco_technologies.cpv.si as si
import pandas as pd
import os

import matplotlib.pyplot as plt


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
                                        surface_tilt=surface_tilt,
                                        type='ins')
    plt.plot(power_cpv, '-r', alpha=0.7)
#    plt.show()
    power_si=si.create_si_timeseries(lat=lat, lon=lon, weather=weather,
                                        surface_azimuth=surface_azimuth,
                                        surface_tilt=surface_tilt)
    plt.plot(power_si, '-y', alpha=0.7)
#    plt.show()
    added=power_cpv + power_si.fillna(0)
    return added

if __name__ == '__main__':

    filename = os.path.abspath('/home/local/RL-INSTITUT/inia.steinbach/Dokumente/greco-project/pvcompare/pvcompare/Data/weatherdata.csv')
    weather_df = pd.read_csv(filename, index_col=0,
                             date_parser=lambda idx: pd.to_datetime(idx,
                                                                    utc=True))
    weather_df.index = pd.to_datetime(weather_df.index).tz_convert(
        'Europe/Berlin')
    weather_df['dni']=weather_df['ghi']-weather_df['dhi']

    ds=create_hybrid_timeseries(lat=40.3, lon=5.4, weather=weather_df,
                             surface_tilt=25, surface_azimuth=180)
    plt.plot(ds, '-b', alpha=0.7)
    plt.show()