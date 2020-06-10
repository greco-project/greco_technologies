
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
#import SMARTS.era5

# Reconfiguring the logger here will also affect test running in the PyCharm IDE
log_format = "%(asctime)s %(levelname)s %(filename)s:%(lineno)d %(message)s"
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, format=log_format)


def create_perosi_timeseries(
        lat, lon, surface_azimuth, surface_tilt, year, input_directory=None
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
    #load EQE data
    EQE_psi = pd.read_csv(os.path.join(input_directory, "EQE_3T_peroskit_measured_corrected.csv"), index_col=0)
    EQE=EQE_psi/100

    # define output data format
    iout='4 12'


    #Psi parameters from Tockhorne und Korte
    j0_ref=3.5248/10**(15)
    rs=7.6
    rsh=6230
    eg=1.636 *1.602176634 / 10**19#eV
    n=1.5
    kB=1.380649 / 10**23 #J/K
    Ns= 50 #Number of cells in series
    q = 1.602176634 / 10 ** (19)  # in Coulomb = A/s
    A = 0.78 *Ns # Area of the solar cell in cm²

    #calculate spectral parameters from smarts and era5
    spectral_parameters=smarts.calculate_Jsc_from_smarts(year=year, lat=lat, lon=lon,EQE=EQE, stop=96)

    #calculate cell temperature
    t_cell=pvlib.temperature.pvsyst_cell(spectral_parameters["ghi"], spectral_parameters['temp'], spectral_parameters['wind_speed'])
    j0=pd.Series()
    for index, value in t_cell.items():
        j0[index]=j0_ref * ((value / 25)**3) * math.exp(((q*eg)/(n*kB))*((1/25)-(1/value)))


    nNsVth= n * Ns * (kB * t_cell/q)

    Isc = spectral_parameters["Jsc"] * A

    singlediode=pvlib.pvsystem.singlediode(photocurrent=Isc, saturation_current=j0,
                               resistance_series=rs, resistance_shunt=rsh, nNsVth=nNsVth,
                               ivcurve_pnts=None, method='lambertw')
    plt.plot(singlediode['i_sc'], singlediode["v_oc"], marker='o')
    plt.show()
    print(singlediode)
#        def integrand(x):
#            return x ** 2
#        integrand = stat
#        ans, err = quad(integrand, 0, 1)
 #       print
#        ans
  #      I_sc= q * int( b_s(E) QE(E) d(E)







if __name__ == "__main__":

    output=create_perosi_timeseries(
            lat="32.", lon=-110.92, surface_azimuth=180,
            surface_tilt=30, year=2000,
            input_directory=None
    )
    print(output)