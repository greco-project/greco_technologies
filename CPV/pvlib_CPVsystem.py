
import sys

import numpy as np
import pandas as pd
from numpy import nan
import os
from pvlib.modelchain import ModelChain
from pvlib.pvsystem import PVSystem
from pvlib.tracking import SingleAxisTracker
from pvlib.location import Location

from pvlib._deprecation import pvlibDeprecationWarning
import pvlib.solarposition
import pvlib.atmosphere
import pvlib.irradiance as irrad



def UF_corrected_DNI(ufdam, ufdni, uftamb, c1, c2, c3, dni):

    '''
    The approach follows the model of Gerstmaier (Quelle), promoting a Utilization Faktor
    that reflects the spectral (AM), temperature (Tamb) and DNI dependencies of
    multijunction cells. The Utilization Factor is multiplied with the DNI.

    UF = c1 * UF(AM) + c2 * UF(DNI) + c3 * UF(Tamb)

    :param ufdam: Utilization Faktor of AM:
    :param ufdni: Utilization Faktor of DNI
    :param uftamb: Utilization Faktor of Tamb
    :param dni: incoming DNI
    :return: The altered DNI that includes the parametrization of multijunction
            cells
    '''

    UF = c1 * ufdam + c2 * ufdni + c3 * uftamb

    return dni * UF

