import pandas as pd
import os
import numpy as np
import pvlib
import logging


import cpvtopvlib.cpvsystem as cpv
import greco_technologies.cpv.inputs as inputs
import matplotlib.pyplot as plt


def create_cpv_timeseries(
    lat, lon, weather, surface_azimuth, surface_tilt, type
):

    if type == "ins":
        module_params = inputs.create_ins_cpv_dict()
        calc_uf_aoi = True
    elif type == "m300":
        module_params = inputs.create_m300_cpv_dict()
        calc_uf_aoi = False
    else:
        logging.error("The type is not known.")

    cpv_sys = cpv.StaticCPVSystem(
        surface_tilt=surface_tilt,
        surface_azimuth=surface_azimuth,
        module=None,
        module_parameters=module_params,
        modules_per_string=1,
        strings_per_inverter=1,
        inverter=None,
        inverter_parameters=None,
        racking_model="insulated",
        losses_parameters=None,
        name=None,
    )

    panel_location = pvlib.location.Location(latitude=lat, longitude=lon)

    # calculate weather parameters: airmass, aoi
    spa = panel_location.get_solarposition(
        times=weather.index,
        pressure=None,                                                          # todo:if pressure is not None..
        temperature=weather["temp_air"],
    )
    airmass = panel_location.get_airmass(weather.index)
    relative_airmass = airmass["airmass_relative"].fillna(0)

    # calculate AOI
    aoi_list = pd.Series(name="aoi")
    ot_list = pd.Series(name='ot')
    gt_list = pd.Series(name='gt')

    for index, row in spa.iterrows():
        aoi = pvlib.irradiance.aoi(
            surface_tilt=surface_tilt,
            surface_azimuth=surface_azimuth,
            solar_zenith=row["zenith"],
            solar_azimuth=row["azimuth"],
        )
        CPVSystem=cpv.CPVSystem()
        ot_list[index] = CPVSystem.optical_transmission_losses(aoi=aoi)
        gt_list[index] = CPVSystem.glass_transmission_losses(aoi=aoi)
        aoi_list[index] = aoi

    weather["aoi"] = aoi_list
    alignement_transmission = 0.95  # emperical parameter for Insolight module
    weather['glass_transmission'] = gt_list
    #if uf_aoi is calculated, the optical transmission losses are neglected at
    # this point
    if calc_uf_aoi==True:
        weather['dni'] = weather['dni'] * gt_list * alignement_transmission
    else:
        weather['dni'] = weather['dni'] * ot_list * gt_list * alignement_transmission

    celltemp = cpv_sys.pvsyst_celltemp(
        weather["ghi"], weather["temp_air"], weather["wind_speed"]
    )

    (
        photocurrent,
        saturation_current,
        resistance_series,
        resistance_shunt,
        nNsVth,
    ) = cpv_sys.calcparams_pvsyst(weather["dni"], celltemp)

    cpv_sys.diode_params = (
        photocurrent,
        saturation_current,
        resistance_series,
        resistance_shunt,
        nNsVth,
    )

    cpv_sys.dc = cpv_sys.singlediode(
        photocurrent, saturation_current, resistance_series, resistance_shunt, nNsVth
    )

    estimation = cpv_sys.dc["p_mp"]

    #calculate utilization factors
    uf_am = []
    for i, v in relative_airmass.items():
        uf_am.append(
            cpv.get_simple_util_factor(
                v,
                module_params["thld_am"],
                module_params["m_low_am"] / module_params["IscDNI_top"],
                module_params["m_high_am"] / module_params["IscDNI_top"],
            )
        )
    uf_temp = []
    for i, v in weather["temp_air"].items():
        uf_temp.append(
            cpv.get_simple_util_factor(
                v,
                module_params["thld_temp"],
                module_params["m_low_temp"] / module_params["IscDNI_top"],
                module_params["m_high_temp"] / module_params["IscDNI_top"],
            )
        )

    uf_am_temp = np.multiply(module_params["weight_am"], uf_am) + np.multiply(
        module_params["weight_temp"], uf_temp
    )
    if calc_uf_aoi == True:
        uf_aoi = []
        for i, v in weather["aoi"].items():
            uf_aoi.append(
                cpv.get_simple_util_factor(
                    v,
                    module_params["thld_aoi"],
                    module_params["m_low_aoi"] / module_params["IscDNI_top"],
                    module_params["m_high_aoi"] / module_params["IscDNI_top"],
                )
            )

        uf_aoi_top = cpv.get_simple_util_factor(
            0,
            module_params["thld_aoi"],
            module_params["m_low_aoi"] / module_params["IscDNI_top"],
            module_params["m_high_aoi"] / module_params["IscDNI_top"],
        )

        uf_aoi_norm = np.divide(uf_aoi, uf_aoi_top)
        uf_global = np.multiply(uf_am_temp, uf_aoi_norm)
    else:
        uf_global = uf_am_temp

    return estimation * uf_global


if __name__ == "__main__":
    #load example weather_data
    filename = os.path.abspath(
        "/home/local/RL-INSTITUT/inia.steinbach/Dokumente/greco-project/pvcompare/pvcompare/data/inputs/weatherdata.csv"
    )
    weather_df = pd.read_csv(
        filename, index_col=0, date_parser=lambda idx: pd.to_datetime(idx, utc=True)
    )
    weather_df.index = pd.to_datetime(weather_df.index).tz_convert("Europe/Berlin")
    weather_df["dni"] = weather_df["ghi"] - weather_df["dhi"]

    ds = create_cpv_timeseries(
        lat=40.3,
        lon=5.4,
        weather=weather_df,
        surface_azimuth=180,
        surface_tilt=30,
        type="m300",
    )
    plt.plot(ds)
    plt.show()
