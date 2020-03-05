import pandas as pd


def create_cpv_dict(cpvtype):

    """
    returns the dict of module parameters and utilization factor for the
    insolight cpv module
    type: str
        'ins' or 'm300'
    :return: dict
    """
    if cpvtype == 'ins':
        module_params = {
            "gamma_ref": 5.524,
            "mu_gamma": 0.003,
            "I_L_ref": 0.96,
            "I_o_ref": 0.00000000017,
            "R_sh_ref": 5226,
            "R_sh_0": 21000,
            "R_sh_exp": 5.50,
            "R_s": 0.01,
            "alpha_sc": 0.00,
            "EgRef": 3.91,
            "irrad_ref": 1000,
            "temp_ref": 25,
            "cells_in_series": 12,
            "cells_in_parallel": 48,
            "eta_m": 0.32,
            "alpha_absorption": 0.9,
            "Area": 0.103,
            "Impo": 8.3,
            "Vmpo": 43.9,
            "v_mp": 33.5,
            "i_mp": 0.893
        }

        UF_parameters = {
            "IscDNI_top": 0.96 / 1000,
            "thld_aoi": 61.978505569631494,
            "m_low_aoi": -2.716773886925838e-07,
            "m_high_aoi": -1.781998474992582e-05,
            "thld_am": 4.574231933073185,
            "m_low_am": 3.906372068620377e-06,
            "m_high_am": -3.0335768119184845e-05,
            "thld_temp": 50,
            "m_low_temp": 4.6781224141650075e-06,
            "m_high_temp": 0,
            "weight_am": 0.2,
            "weight_temp": 0.8,
        }

        module_params.update(UF_parameters)


    elif cpvtype == 'm300':

        module_params = {'gamma_ref': 4.456, 'mu_gamma': 0.0012,
                         'I_L_ref': 3.346,
                         'I_o_ref': 0.000000000004, 'R_sh_ref': 4400,
                         'R_sh_0': 17500, 'R_sh_exp': 5.50, 'R_s': 0.736,
                         'alpha_sc': 0.00, 'irrad_ref': 1000, 'temp_ref': 25,
                         'cells_in_series': 42, 'v_mp': 116.63, 'i_mp': 3.082,
                         'Area': 1.269}

        UF_parameters = {
            "IscDNI_top": 1,
            "thld_am": 2.022411098853249,
            "m_low_am": 0.0423037910485609,
            "m_high_am": -0.0210539236615148,
            "thld_temp": 200,
            "m_low_temp": 0.000923828521724516,
            "m_high_temp": 0.0,
            "weight_am": 0.2,
            "weight_temp": 0.8,
        }

        module_params.update(UF_parameters)

    return module_params
