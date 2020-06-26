
import pandas as pd
import logging
import sys
import os
import pvlib
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
import math

import pvlib
import pvlib_smarts as smarts
import era5
#import SMARTS.era5

# Reconfiguring the logger here will also affect test running in the PyCharm IDE
log_format = "%(asctime)s %(levelname)s %(filename)s:%(lineno)d %(message)s"
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=log_format)


def calculate_smarts_parameters(year, lat, lon, number_hours, WLMN, WLMX, cell_type, input_directory):

    """

    :param year:
    :return:
    """
    # load atmos data from era5
    atmos_data= era5.load_era5_weatherdata(lat, lon, year, variable="perosi")
    delta = pd.to_timedelta(30, unit='m')
    atmos_data.index = atmos_data.index + delta
    atmos_data['davt'] = atmos_data["temp_air"].resample('D').mean()
    atmos_data=atmos_data.fillna(method='ffill')

    #define constant
    q = 1.602176634 / 10 ** (19)  # in Coulomb = A/s
    # define output data format
    iout = '4 12'
    #define counter for number of hours to be calculated
    c = 0
    #define time interval of one year
    time = pd.date_range(start=f'1/1/{year}', end=f'31/12/{year}', freq='H')
    # calculate Jsc for every timestep
    result=pd.DataFrame()
    for timestep in time:
        if timestep.month in range(3, 8):
            season = 'SUMMER'
        else:
            season = 'WINTER'

        # load spectral data from SMARTS
        spectrum = smarts.SMARTSSpectra(
            IOUT=iout, YEAR=str(year),
            MONTH=str(timestep.month),
            DAY=str(timestep.day), HOUR=str(timestep.hour),
            LATIT=lat,
            LONGIT=str(lon), WLMN=WLMN,
            WLMX=WLMX,
            TAIR=str(atmos_data.at[timestep, 'temp_air']),
            TDAY=str(atmos_data.at[timestep, 'davt']),
            SEASON=season,
            ZONE=1)

        # load EQE data
        for x in cell_type:
            if x == "Korte_pero":
                import data.cell_parameters_korte_pero as param
            elif x == "Korte_si":
                import data.cell_parameters_korte_si as param
            else:
                logging.error("The cell type is not recognized. Please "
                              "choose either 'Korte_si' or 'Korte_pero'.")
            EQE_filename = param.EQE_filename
            EQE = pd.read_csv(os.path.join(input_directory, EQE_filename),
                              index_col=0)
            EQE = EQE / 100

            if spectrum.empty == True:
                result.at[timestep, "Jsc_" + str(x)] = 0
                result.at[timestep, "ghi"] = 0

            else:
                if not spectrum.index.name == "Wvlgth":
                    spectrum.set_index('Wvlgth', inplace=True)
                Jsc_lambda = (spectrum["Global_tilt_photon_irrad"] * EQE[
                    "EQE"]) * q
                Jsc_lambda.fillna(0, inplace=True)
                result.at[timestep, "Jsc_" + str(x)] = Jsc_lambda.sum()  # in A/cm²
                result.at[timestep, "ghi"] = spectrum["Global_horizn_irradiance"].sum()

        result.at[timestep, "temp"] = atmos_data.at[timestep, 'temp_air']
        result.at[timestep, "wind_speed"] = atmos_data.at[timestep, "wind_speed"]

        # check if number of hours is reached
        c = c + 1
        if c == number_hours:
            break
    return result

def create_timeseries(
        lat, lon, surface_azimuth, surface_tilt, year, input_directory=None, cell_type="Korte"
    ):
    """

    creates a time series for a type of pero-si 3T Tandem module

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

    if input_directory == None:
        input_directory = os.path.join(
    os.path.dirname(__file__), "data/"
)
    # define output data format
    iout='4 12'

    q = 1.602176634 / 10 ** (19)  # in Coulomb = A/s
    kB = 1.380649 / 10 ** 23  # J/K

    #calculate spectral parameters from smarts and era5
    smarts_parameters=calculate_smarts_parameters(year=year, lat=lat, lon=lon, WLMN=280, WLMX=1200, number_hours=50, cell_type=cell_type, input_directory=input_directory)

    #calculate cell temperature characteristics
    t_cell=pvlib.temperature.pvsyst_cell(smarts_parameters["ghi"], smarts_parameters['temp'], smarts_parameters['wind_speed'])
    #temperature effect on saturation current
    result=pd.DataFrame()
    for x in cell_type:
        if x == "Korte_pero":
            import data.cell_parameters_korte_pero as param
        elif x == "Korte_si":
            import data.cell_parameters_korte_si as param
        j0=pd.Series()
        for index, value in t_cell.items():
            j0[index]=param.j0_ref * ((value / 25)**3) * math.exp(((q*param.eg)/(param.n*kB))*((1/25)-(1/value)))
        #temperature effect on short circuit current
        smarts_parameters["Jsc_temp"]=smarts_parameters["Jsc_" + str(x)] * param.alpha * (t_cell - param.temp_ref)
        nNsVth= param.n * param.Ns * (kB * (t_cell+273.15)/q)

#       Jsc_in_m = spectral_parameters["Jsc"]*10000
        Isc = smarts_parameters["Jsc_temp"] * param.A

        singlediode=pvlib.pvsystem.singlediode(photocurrent=Isc, saturation_current=param.j0_ref,
                               resistance_series=param.rs, resistance_shunt=param.rsh, nNsVth=nNsVth,
                               ivcurve_pnts=None, method='lambertw')
        result[str(x) + "_p_mp"] = singlediode["p_mp"]

    return result


    analysis = pd.DataFrame()
    analysis["ghi"] = smarts_parameters["ghi"]
    analysis["temp"] = smarts["temp"]
    analysis["p_mp"] = singlediode["p_mp"]

    plt.plot(singlediode['p_mp'], marker='o', label="power in W")
    plt.xlabel("Time")
    plt.ylabel("Power in W")
    plt.legend()
    plt.show()

    fig, ax = plt.subplots()
    im = ax.scatter(analysis["ghi"], analysis["p_mp"], c=analysis["temp"])
    ax.set_xlabel('Irradiance in W/m²')
    ax.set_ylabel('Power in W')
    # Add a colorbar
    cbar=plt.colorbar(im, ax=ax)
    cbar.set_label('Temperatute in °C', rotation=270)
    plt.show()

    fig, ax = plt.subplots()
    im = ax.scatter(analysis["temp"], analysis["p_mp"], c=analysis["ghi"])
    ax.set_xlabel('Temperatute in °C')
    ax.set_ylabel('Power in W')
    # Add a colorbar
    cbar=plt.colorbar(im, ax=ax)
    cbar.set_label('Irradiance in W/m²', rotation=270)
    plt.show()

    return singlediode



def create_pero_si_timeseries(lat="45.", lon=5.87, surface_azimuth=180,
        surface_tilt=30, year=2015,
        input_directory=None):
    """

    :param lat:
    :param lon:
    :param surface_azimuth:
    :param surface_tilt:
    :param year:
    :param input_directory:
    :return:
    """
    pero=create_timeseries(
        lat="45.", lon=5.87, surface_azimuth=180,
        surface_tilt=30, year=2015,
        input_directory=None, cell_type="KortePero"
    )
    si= pero=create_timeseries(
        lat="45.", lon=5.87, surface_azimuth=180,
        surface_tilt=30, year=2015,
        input_directory=None, cell_type="KorteSi"
    )
    output=pero["p_mp"]+ si["p_mp"]

    plt.plot(pero['p_mp'], marker='o', label="power in W")
    plt.plot(si['p_mp'], marker='o', label="power in W")
    plt.xlabel("Time")
    plt.ylabel("Power in W")
    plt.legend()
    plt.show()

    plt.plot(output, marker='o', label="power in W")
    plt.xlabel("Time")
    plt.ylabel("Power in W")
    plt.legend()
    plt.show()

    return output






if __name__ == "__main__":


    output = create_timeseries(
        lat="45.", lon=5.87, surface_azimuth=180,
        surface_tilt=30, year=2015,
        input_directory=None, cell_type=["Korte_si"]
    )
    print(output)