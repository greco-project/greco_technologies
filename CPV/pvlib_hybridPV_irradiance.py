from __future__ import division


import numpy as np
import pandas as pd


# the hybrid module accounts for the case, if a solar cell is only exposed to
# DHI irradiance and no DNI. In this case the following functions need to be
# adjusted:
# - sapm_spectral_loss()
# - sapm_aoi_loss()
# - sapm_effective_irradiance()
# - sapm_celltemp() todo: does the cell temperature change?

def hybrid_spectral_loss(airmass_absolute, module):
    """
    Calculates the hybrid spectral loss coefficient, F1.

    Parameters
    ----------
    airmass_absolute : numeric
        Absolute airmass

    module : dict-like
        A dict, Series, or DataFrame defining the SAPM performance
        parameters. See the :py:func:`sapm` notes section for more
        details.

    Returns
    -------
    F1 : numeric
        The SAPM spectral loss coefficient.

    Notes
    -----
    nan airmass values will result in 0 output.
    """

    am_coeff = [module['A4'], module['A3'], module['A2'], module['A1'],
                module['A0']]

    spectral_loss = np.polyval(am_coeff, airmass_absolute)

    spectral_loss = np.where(np.isnan(spectral_loss), 0, spectral_loss)

    spectral_loss = np.maximum(0, spectral_loss)

    if isinstance(airmass_absolute, pd.Series):
        spectral_loss = pd.Series(spectral_loss, airmass_absolute.index)

    return spectral_loss


def hybrid_aoi_loss(aoi, module, upper=None):
    """
    Calculates the hybrid angle of incidence loss coefficient, F2.

    Parameters
    ----------
    aoi : numeric
        Angle of incidence in degrees. Negative input angles will return
        zeros.

    module : dict-like
        A dict, Series, or DataFrame defining the SAPM performance
        parameters. See the :py:func:`sapm` notes section for more
        details.

    upper : None or float, default None
        Upper limit on the results.

    Returns
    -------
    F2 : numeric
        The SAPM angle of incidence loss coefficient.

    Notes
    -----
    The SAPM traditionally does not define an upper limit on the AOI
    loss function and values slightly exceeding 1 may exist for moderate
    angles of incidence (15-40 degrees). However, users may consider
    imposing an upper limit of 1.

    References
    ----------
    [1] King, D. et al, 2004, "Sandia Photovoltaic Array Performance
    Model", SAND Report 3535, Sandia National Laboratories, Albuquerque,
    NM.

    [2] B.H. King et al, "Procedure to Determine Coefficients for the
    Sandia Array Performance Model (SAPM)," SAND2016-5284, Sandia
    National Laboratories (2016).

    [3] B.H. King et al, "Recent Advancements in Outdoor Measurement
    Techniques for Angle of Incidence Effects," 42nd IEEE PVSC (2015).
    DOI: 10.1109/PVSC.2015.7355849
    """

    aoi_coeff = [module['B5'], module['B4'], module['B3'], module['B2'],
                 module['B1'], module['B0']]

    aoi_loss = np.polyval(aoi_coeff, aoi)
    aoi_loss = np.clip(aoi_loss, 0, upper)
    # nan tolerant masking
    aoi_lt_0 = np.full_like(aoi, False, dtype='bool')
    np.less(aoi, 0, where=~np.isnan(aoi), out=aoi_lt_0)
    aoi_loss = np.where(aoi_lt_0, 0, aoi_loss)

    if isinstance(aoi, pd.Series):
        aoi_loss = pd.Series(aoi_loss, aoi.index)

    return aoi_loss


def hybrid_effective_irradiance(poa_direct, poa_diffuse, airmass_absolute, aoi,
                              module, reference_irradiance=1000):
    """
    Calculates the SAPM effective irradiance using the SAPM spectral
    loss and SAPM angle of incidence loss functions.

    Parameters
    ----------
    poa_direct : numeric
        The direct irradiance incident upon the module.

    poa_diffuse : numeric
        The diffuse irradiance incident on module.

    airmass_absolute : numeric
        Absolute airmass.

    aoi : numeric
        Angle of incidence in degrees.

    module : dict-like
        A dict, Series, or DataFrame defining the SAPM performance
        parameters. See the :py:func:`sapm` notes section for more
        details.

    reference_irradiance : numeric, default 1000
        Reference irradiance by which to divide the input irradiance.

    Returns
    -------
    effective_irradiance : numeric
        The SAPM effective irradiance.
    """

    F1 = sapm_spectral_loss(airmass_absolute, module)
    F2 = sapm_aoi_loss(aoi, module)

    E0 = reference_irradiance

    Ee = F1 * (poa_direct*F2 + module['FD']*poa_diffuse) / E0

    return Ee

def sapm_celltemp(poa_global, wind_speed, temp_air,
                  model='open_rack_cell_glassback'):
    '''
    Estimate cell and module temperatures per the Sandia PV Array
    Performance Model (SAPM, SAND2004-3535), from the incident
    irradiance, wind speed, ambient temperature, and SAPM module
    parameters.

    Parameters
    ----------
    poa_global : float or Series
        Total incident irradiance in W/m^2.

    wind_speed : float or Series
        Wind speed in m/s at a height of 10 meters.

    temp_air : float or Series
        Ambient dry bulb temperature in degrees C.

    model : string, list, or dict, default 'open_rack_cell_glassback'
        Model to be used.

        If string, can be:

            * 'open_rack_cell_glassback' (default)
            * 'roof_mount_cell_glassback'
            * 'open_rack_cell_polymerback'
            * 'insulated_back_polymerback'
            * 'open_rack_polymer_thinfilm_steel'
            * '22x_concentrator_tracker'

        If dict, supply the following parameters
        (if list, in the following order):

            * a : float
                SAPM module parameter for establishing the upper
                limit for module temperature at low wind speeds and
                high solar irradiance.

            * b : float
                SAPM module parameter for establishing the rate at
                which the module temperature drops as wind speed increases
                (see SAPM eqn. 11).

            * deltaT : float
                SAPM module parameter giving the temperature difference
                between the cell and module back surface at the
                reference irradiance, E0.

    Returns
    --------
    DataFrame with columns 'temp_cell' and 'temp_module'.
    Values in degrees C.

    References
    ----------
    [1] King, D. et al, 2004, "Sandia Photovoltaic Array Performance
    Model", SAND Report 3535, Sandia National Laboratories, Albuquerque,
    NM.

    See Also
    --------
    sapm
    '''

    temp_models = TEMP_MODEL_PARAMS['sapm']

    if isinstance(model, str):
        model = temp_models[model.lower()]

    elif isinstance(model, (dict, pd.Series)):
        model = [model['a'], model['b'], model['deltaT']]

    a = model[0]
    b = model[1]
    deltaT = model[2]

    E0 = 1000.  # Reference irradiance

    temp_module = pd.Series(poa_global * np.exp(a + b * wind_speed) + temp_air)

    temp_cell = temp_module + (poa_global / E0) * (deltaT)

    return pd.DataFrame({'temp_cell': temp_cell, 'temp_module': temp_module})


def ot_losses(optical_transmission, DNI, DHI):

    '''
    This function accounts for the case that the optical transmission of irradiation
    is changes e.g. by a lens.
    '''
