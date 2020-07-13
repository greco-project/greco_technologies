
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
import greco_technologies.perosi.pvlib_smarts as smarts
import greco_technologies.perosi.era5 as era5

# Reconfiguring the logger here will also affect test running in the PyCharm IDE
log_format = "%(asctime)s %(levelname)s %(filename)s:%(lineno)d %(message)s"
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=log_format)


def calculate_smarts_parameters(year, lat, lon, number_hours, cell_type,
                                input_directory,
                                surface_tilt, surface_azimuth,
                                WLMN=280, WLMX=1200):

    """

    :param year: int
        year of interest
    :param lat: str
        latitude ending with a ".", e.g. "45."
    :param lon: int
        longitude
    :param number_hours: int
        number of hours until simulation stops. For one year enter 8760.
    :param WLMN: int
        minimum wavelength of the spectrum. By default this is 280 nm.
    :param WLMX: int
        maximum wavelength of the spectrum. By default this is 1200 nm.
    :param cell_type: list
        list of cells for which the Jsc should be calculated.
    :param input_directory: str
        name of the input directory
    :return: :pd.Dataframe()
           including ghi, temperature, wind_speed, Jsc_"cell_type"
    """
    # load atmos data from era5
    atmos_data= era5.load_era5_weatherdata(lat, lon, year, variable="perosi")
    delta = pd.to_timedelta(30, unit='m')
    atmos_data.index = atmos_data.index + delta
    atmos_data['davt'] = atmos_data["temp_air"].resample('D').mean()
    atmos_data=atmos_data.fillna(method='ffill')

    #define constant
    q = 1.602176634 / 10 ** (19)  # in Coulomb = A*s
    # define output data format
    iout = '8 12'
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
        import decimal
        d = decimal.Decimal(str(lat))
        decimals_lat=d.as_tuple().exponent
        lat_spectrum=str(lat)[:decimals_lat]
        spectrum = smarts.SMARTSSpectra(
            IOUT=iout, YEAR=str(year),
            MONTH=str(timestep.month),
            DAY=str(timestep.day), HOUR=str(timestep.hour),
            LATIT=lat_spectrum,
            LONGIT=str(lon), WLMN=WLMN,
            WLMX=WLMX,
            TAIR=str(atmos_data.at[timestep, 'temp_air']),
            TDAY=str(atmos_data.at[timestep, 'davt']),
            SEASON=season,
            ZONE=1,
            TILT=str(surface_tilt),
            WAZIM=str(surface_azimuth))

        # load EQE data
        for x in cell_type:
            if x == "Korte_pero":
                import greco_technologies.perosi.data.cell_parameters_korte_pero as param
            elif x == "Korte_si":
                import greco_technologies.perosi.data.cell_parameters_korte_si as param
            else:
                logging.error("The cell type is not recognized. Please "
                              "choose either 'Korte_si' or 'Korte_pero'.")
            EQE_filename = param.EQE_filename
            EQE = pd.read_csv(os.path.join(input_directory, EQE_filename),
                              index_col=0)
            EQE = EQE / 100

            # return Jsc and ghi = 0 if the spectrum is empty
            if spectrum.empty == True:
                result.at[timestep, "Jsc_" + str(x)] = 0
                result.at[timestep, "ghi"] = 0

            else:
                if not spectrum.index.name == "Wvlgth":
                    spectrum.set_index('Wvlgth', inplace=True)
                # calculate Jsc
                Jsc_lambda = (spectrum["Global_tilt_photon_irrad"] * EQE[
                    "EQE"]) * q
                Jsc_lambda.fillna(0, inplace=True)
                result.at[timestep, "Jsc_" + str(x)] = Jsc_lambda.sum() # in A/m²
                result.at[timestep, "ghi"] = spectrum["Global_tilted_irradiance"].sum() # in W/m²

        result.at[timestep, "temp"] = atmos_data.at[timestep, 'temp_air']
        result.at[timestep, "wind_speed"] = atmos_data.at[timestep, "wind_speed"]

        # check if number of hours is reached
        c = c + 1
        if c == number_hours:
            break
    return result


def create_timeseries(
        lat, lon, surface_azimuth, surface_tilt, year,
        cell_type, number_hours,
        input_directory=None, plot=True
    ):

    """
    :param year: int
        year of interest
    :param lat: str
        latitude ending with a ".", e.g. "45."
    :param lon: int
        longitude
    :param number_hours: int
        number of hours until simulation stops. For one year enter 8760.
    :param cell_type: list
        list of cells for which the Jsc should be calculated.
    :param input_directory: str
        name of the input directory
    :param surface_azimuth:
    :param surface_tilt:

    :return: :pd.Dataframe()
        maximum power point of each time step for each cell type
    """
    if input_directory == None:
        input_directory = os.path.join(
    os.path.dirname(__file__), "data/"
)
    q = 1.602176634 / 10 ** (19)  # in Coulomb = A/s
    kB = 1.380649 / 10 ** 23  # J/K

    #calculate spectral parameters from smarts and era5
    smarts_parameters=calculate_smarts_parameters(year=year, lat=lat, lon=lon,
                                                  number_hours=number_hours,
                                                  cell_type=cell_type,
                                                  input_directory=input_directory,
                                                  surface_tilt=surface_tilt,
                                                  surface_azimuth=surface_azimuth
    )

    #calculate cell temperature characteristics
    t_cell=pvlib.temperature.pvsyst_cell(smarts_parameters["ghi"], smarts_parameters['temp'], smarts_parameters['wind_speed'])
    #temperature effect on saturation current
    result=pd.DataFrame()
    for x in cell_type:
        if x == "Korte_pero":
            import greco_technologies.perosi.data.cell_parameters_korte_pero as param
        elif x == "Korte_si":
            import greco_technologies.perosi.data.cell_parameters_korte_si as param
#        j0=pd.Series()
#        for index, value in t_cell.items():
 #           j0[index]=param.j0_ref * ((value / 25)**3) * math.exp(((q*param.eg)/(param.n*kB))*((1/25)-(1/value)))
        #temperature effect on short circuit current
#        smarts_parameters["Jsc_temp"]= smarts_parameters["Jsc_" + str(x)] + smarts_parameters["ghi"]/1000 *(param.alpha * (t_cell - param.temp_ref))
        nNsVth= param.n * param.Ns * (kB * (t_cell+273.15)/q)

#       Jsc_in_m = spectral_parameters["Jsc"]*10000
        Isc = smarts_parameters["Jsc_" + str(x)]* param.A
        singlediode=pvlib.pvsystem.singlediode(photocurrent=Isc, saturation_current=param.I_0,
                               resistance_series=param.rs, resistance_shunt=param.rsh, nNsVth=nNsVth,
                               ivcurve_pnts=None, method='lambertw')
        result[str(x) + "_p_mp"] = singlediode["p_mp"]
        result[str(x) + "_p_mp" + "_temp"]=result[str(x) + "_p_mp"] * (1 + (param.alpha * (t_cell - param.temp_ref)))

        if plot == True:
            fig, ax1 = plt.subplots()

            color = 'tab:red'
            ax1.set_xlabel('time')
            ax1.set_ylabel('Power in mW / Temperature in °C', color=color)
            ax1.plot(result[str(x) + "_p_mp" + "_temp"]*1000, color="r", alpha=0.5, label="power")
            ax1.plot(t_cell, color="g", alpha=0.5, label="Cell temperature")
            ax1.tick_params(axis='y', labelcolor=color)

            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

            color = 'tab:blue'
            ax2.set_ylabel('irradiance in W/m²',
                           color=color)  # we already handled the x-label with ax1
            ax2.plot(smarts_parameters["ghi"], color=color, alpha=0.5, label="Irradiance")
            ax2.tick_params(axis='y', labelcolor=color)
            fig.legend()

            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()

            fig, ax1 = plt.subplots()

            ax1.set_xlabel('time')
            ax1.set_ylabel('Power in mW', color="b")
            ax1.plot(result[str(x) + "_p_mp" + "_temp"] * 1000, color="purple",
                     alpha=0.5, label="power_t-corrected")
            ax1.plot((singlediode['p_mp'] * 1000), color=color, alpha=0.5,
                     label="power_uncorrected")
            #ax1.plot(t_cell, color=color, alpha=0.5, label="Temperature")
            ax1.tick_params(axis='y', labelcolor=color)

            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

            color = 'tab:blue'
            ax2.set_ylabel('Temperature in °C',
                           color=color)  # we already handled the x-label with ax1
            ax2.plot(smarts_parameters["temp"], color="orange", alpha=0.5, label="Air Temperature")
            ax2.plot(t_cell, color="red", alpha=0.5, label="Cell temperature")
            ax2.tick_params(axis='y', labelcolor=color)
            fig.legend()

            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            plt.show()


    return result



def create_pero_si_timeseries(year, lat, lon, surface_azimuth,
        surface_tilt, number_hours,
        input_directory=None):

    """
    creates a time series for the output power of a pero-si module

    :param year: int
        year of interest
    :param lat: str
        latitude ending with a ".", e.g. "45."
    :param lon: int
        longitude
    :param number_hours: int
        number of hours until simulation stops. For one year enter 8760.
    :param input_directory: str
        name of the input directory
    :param surface_azimuth:
    :param surface_tilt:

    :return:pd. series()
        time series of the output power
    """
    timeseries=create_timeseries(
        lat=lat, lon=lon, surface_azimuth=surface_azimuth,
        surface_tilt=surface_tilt, year=year,
        cell_type=["Korte_pero", "Korte_si"], number_hours=number_hours,
        input_directory=input_directory, plot=False
    )
    output = timeseries.iloc[:,0] + timeseries.iloc[:,1]

    return output






if __name__ == "__main__":


    # output = create_timeseries(
    #     lat="45.", lon=5.87, surface_azimuth=180,
    #     surface_tilt=30, year=2015,
    #     input_directory=None, number_hours=300, cell_type=["Korte_pero"], plot=True
    # )
    output = create_pero_si_timeseries(
        lat=45.0, lon=5.87, surface_azimuth=180,
        surface_tilt=30, year=2015,
        input_directory=None, number_hours=400
    )
    print(output)