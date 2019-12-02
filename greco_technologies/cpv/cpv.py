import pandas as pd
import os
import numpy as np
import pvlib

import cpvtopvlib.cpvsystem as cpv



def create_cpv_timeseries(lat, lon, weather, surface_azimuth, surface_tilt,
                          system, calc_uf_aoi=True):

    panel_location = pvlib.location.Location(latitude=lat, longitude=lon)

    #calculate weather parameters: airmass, aoi
    spa = panel_location.get_solarposition(times=weather.index, pressure=None,
                                           temperature=weather['temp_air'])
    airmass = panel_location.get_airmass(weather.index)
    relative_airmass = airmass['airmass_relative'].fillna(0)

    # calculate AOI
    aoi_list = pd.Series(name='aoi')
    for index, row in spa.iterrows():
        aoi = pvlib.irradiance.aoi(surface_tilt=surface_tilt,
                                   surface_azimuth=surface_azimuth,
                                   solar_zenith=row['zenith'],
                                   solar_azimuth=row['azimuth'])
        aoi_list[index] = aoi

    weather['aoi'] = aoi_list

    #initialize static cpv system class
    cpv_sys=system

    celltemp = cpv_sys.pvsyst_celltemp(weather['ghi'],
                                    weather['temp_air'],
                                    weather['wind_speed'])

    (photocurrent, saturation_current, resistance_series,
     resistance_shunt, nNsVth) = (cpv_sys.calcparams_pvsyst(weather['dni'],
                                                         celltemp))

    cpv_sys.diode_params = (photocurrent, saturation_current, resistance_series,
                         resistance_shunt, nNsVth)

    cpv_sys.dc = cpv_sys.singlediode(photocurrent, saturation_current,
                               resistance_series,
                               resistance_shunt, nNsVth)

    estimation = cpv_sys.dc['p_mp']

    IscDNI_top = 0.96 / 1000
    thld_aoi = 61.978505569631494
    m_low_aoi = -2.716773886925838e-07
    m_high_aoi = -1.781998474992582e-05
    thld_am = 4.574231933073185
    m_low_am = 3.906372068620377e-06
    m_high_am = -3.0335768119184845e-05
    thld_temp = 50
    m_low_temp = 4.6781224141650075e-06
    m_high_temp = 0
    weight_am = 0.2
    weight_temp = 0.8

    uf_am = []
    for i, v in relative_airmass.items():
        uf_am.append(cpv.get_simple_util_factor(v, thld_am,
                                                m_low_am / IscDNI_top,
                                                m_high_am / IscDNI_top))
    uf_temp = []
    for i, v in weather['temp_air'].items():
        uf_temp.append(cpv.get_simple_util_factor(v, thld_temp,
                                                  m_low_temp / IscDNI_top,
                                                  m_high_temp / IscDNI_top))

    uf_am_temp = np.multiply(weight_am, uf_am) + np.multiply(weight_temp,
                                                             uf_temp)
    if calc_uf_aoi==True:
        uf_aoi = []
        for i, v in weather['aoi'].items():
            uf_aoi.append(
                cpv.get_simple_util_factor(v, thld_aoi, m_low_aoi / IscDNI_top,
                                       m_high_aoi / IscDNI_top))

        uf_aoi_top = cpv.get_simple_util_factor(0, thld_aoi,
                                            m_low_aoi / IscDNI_top,
                                            m_high_aoi / IscDNI_top)

        uf_aoi_norm = np.divide(uf_aoi, uf_aoi_top)
        uf_global = np.multiply(uf_am_temp, uf_aoi_norm)
    else:
        uf_global=uf_am_temp

    return estimation * uf_global

if __name__ == '__main__':

    filename = os.path.abspath(
        "/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/163_Open_FRED/03-Projektinhalte/AP2 Wetterdaten/open_FRED_TestWetterdaten_csv/fred_data_test_2016.csv")
    weather_df = pd.read_csv(filename, skiprows=range(1, 50), nrows=(5000),
                             index_col=0,
                             date_parser=lambda idx: pd.to_datetime(idx,
                                                                    utc=True))
    weather_df.index = pd.to_datetime(weather_df.index).tz_convert(
        'Europe/Berlin')


    module_params = {'gamma_ref': 5.524, 'mu_gamma': 0.003, 'I_L_ref': 0.96,
                     'I_o_ref': 0.00000000017, 'R_sh_ref': 5226,
                     'R_sh_0': 21000, 'R_sh_exp': 5.50, 'R_s': 0.01,
                     'alpha_sc': 0.00, 'EgRef': 3.91, 'irrad_ref': 1000,
                     'temp_ref': 25, 'cells_in_series': 12, 'eta_m': 0.32,
                     'alpha_absorption': 0.9, 'Area':0.0001, 'Impo':8.3,
                         'Vmpo':43.9}

    cpv_sys = cpv.StaticCPVSystem(surface_tilt=30, surface_azimuth=180,
                           module=None, module_parameters=module_params,
                           modules_per_string=1, strings_per_inverter=1,
                           inverter=None, inverter_parameters=None,
                           racking_model='insulated',
                           losses_parameters=None, name=None)

    create_cpv_timeseries(lat=40.3, lon=5.4, weather=weather_df,
                          surface_azimuth=180, surface_tilt=30, system=cpv_sys)