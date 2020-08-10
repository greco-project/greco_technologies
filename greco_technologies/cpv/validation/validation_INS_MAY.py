import pandas as pd
import pvlib.atmosphere
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
import math
from sklearn.metrics import mean_squared_error

import sys

sys.path.append(
    "/home/local/RL-INSTITUT/inia.steinbach/Dokumente/greco_technologies_to_pvlib/CPV/"
)
from cpvtopvlib import cpvsystem as cpv


df = pd.read_csv("../inputs/InsolightMay2019_filtered.csv", sep=",", index_col=0)


panel_location = pvlib.location.Location(
    latitude=40.453, longitude=-3.727, tz=1, altitude=658
)

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
}

csys = cpv.StaticCPVSystem(
    surface_tilt=30,
    surface_azimuth=180,
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


spa = panel_location.get_solarposition(
    times=df.index, pressure=None, temperature=df["temp_air"]
)


airmass = panel_location.get_airmass(df.index)
relative_airmass = airmass["airmass_relative"].fillna(0)

# calculate AOI
aoi_list = pd.Series(name="aoi")
ot_list = pd.Series(name="ot")
gt_list = pd.Series(name="gt")

for index, row in spa.iterrows():  # todo: correct surface_tilt and surface_azimuth
    aoi = pvlib.irradiance.aoi(
        surface_tilt=30,
        surface_azimuth=180,
        solar_zenith=row["zenith"],
        solar_azimuth=row["azimuth"],
    )
    # calculate optical losses
    aoi_list[index] = aoi
#    ot_list[index]=cpv.optical_transmission_losses(aoi=aoi)
#    gt_list[index]=cpv.glass_transmission_losses(aoi=aoi)

# alignement_transmission = 0.95 #emperical parameter for Insolight module
df["aoi"] = aoi_list
# weather_loc['glass_transmission']=gt_list
# df['DII_new'] = df['dii'] * alignement_transmission * gt_list


celltemp = csys.pvsyst_celltemp(df["gii"], df["temp_air"], df["wind_speed"])

(
    photocurrent,
    saturation_current,
    resistance_series,
    resistance_shunt,
    nNsVth,
) = csys.calcparams_pvsyst(df["dii"], celltemp)

csys.diode_params = (
    photocurrent,
    saturation_current,
    resistance_series,
    resistance_shunt,
    nNsVth,
)

csys.dc = csys.singlediode(
    photocurrent, saturation_current, resistance_series, resistance_shunt, nNsVth
)

real_power = df["Pmp"]
estimation = csys.dc["p_mp"]


# calculate single utilization factors

IscDNI_top = 0.96 / 1000

thld_aoi = 62.33932812381078
m_low_aoi = -3.216303225269502e-07
m_high_aoi = -1.83007860456365e-05


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
    uf_am.append(
        cpv.get_simple_util_factor(
            v, thld_am, m_low_am / IscDNI_top, m_high_am / IscDNI_top
        )
    )


uf_temp = []
for i, v in df["temp_air"].items():
    uf_temp.append(
        cpv.get_simple_util_factor(
            v, thld_temp, m_low_temp / IscDNI_top, m_high_temp / IscDNI_top
        )
    )


uf_aoi = []
for i, v in df["aoi"].items():
    uf_aoi.append(
        cpv.get_simple_util_factor(
            v, thld_aoi, m_low_aoi / IscDNI_top, m_high_aoi / IscDNI_top
        )
    )


uf_aoi_ast = cpv.get_simple_util_factor(
    0, thld_aoi, m_low_aoi / IscDNI_top, m_high_aoi / IscDNI_top
)

uf_aoi_norm = np.divide(uf_aoi, uf_aoi_ast)


uf_am_temp = np.multiply(weight_am, uf_am) + np.multiply(weight_temp, uf_temp)
UF_global = np.multiply(uf_am_temp, uf_aoi_norm)


modeled_power = estimation * UF_global

residualUF = modeled_power - real_power
residualwithoutUF = estimation - real_power

rmsd = math.sqrt(mean_squared_error(real_power, modeled_power))
rmsd1 = math.sqrt(mean_squared_error(real_power, estimation))

print("rmsd real vs modeled power:", rmsd)
print("rmsd real vs estimated power:", rmsd1)

# plt.plot(df['temp_air'], uf_temp, 'b.', label='UF(temp)')
# plt.xlabel("Temperature in C")
# plt.ylabel("Utilization Factor")
# plt.legend()
# plt.show()
#
# plt.plot(relative_airmass, uf_am, 'r.', label='UF(AM)')
# plt.xlabel("Airmass")
# plt.ylabel("Utilization Factor")
# plt.legend()
# plt.show()
#
# plt.plot(df['aoi'], uf_aoi_norm, 'g.', label='UF(aoi)')
# plt.xlabel("AOI in Degrees")
# plt.ylabel("Utilization Factor")
# plt.legend()
# plt.show()
#
#
# plt.plot(modeled_power, 'b', label='modeled power with UF_aoi')
# plt.plot(real_power, 'r', label='real_power')
# plt.plot(estimation, 'g', label='pvlib-calculated power')
# plt.xlabel("Time in Days")
# plt.ylabel("Power in W")
# plt.legend()
# plt.show()

p1 = poly.polyfit(real_power, modeled_power, 1)

plt.plot(
    real_power, real_power.index, "bo", markersize=1, label="with utilization factor"
)
plt.plot(
    modeled_power,
    modeled_power.index,
    "ro",
    markersize=1,
    label="without utilization factor",
)
plt.xlabel("measured power in W")
plt.ylabel("modeled power in W")
plt.legend()
plt.show()

plt.plot(real_power, modeled_power, "bo", markersize=1, label="with utilization factor")
plt.plot(real_power, estimation, "ro", markersize=1, label="without utilization factor")
# plt.plot(real_power, poly.polyval(real_power, p1), 'y-', label='model_power_fit')
plt.plot(real_power, real_power, "g", label="_nolegend_")
plt.xlabel("measured power in W")
plt.ylabel("modeled power in W")
plt.legend()
plt.show()


plt.plot(
    airmass["airmass_relative"],
    residualwithoutUF,
    "go",
    markersize=1,
    label="Airmass residual without UF",
)
plt.plot(
    airmass["airmass_relative"].fillna(0),
    residualUF,
    "ro",
    markersize=1,
    label="Airmass residual with UF",
)
plt.xlabel("Airmass")
plt.ylabel("Residual Pmpp in %")
plt.legend()
plt.show()


plt.plot(
    df["temp_air"], residualUF, "ro", markersize=1, label="Temperature residual with UF"
)
plt.plot(
    df["temp_air"],
    residualwithoutUF,
    "go",
    markersize=1,
    label="Temperature residual without UF",
)
plt.xlabel("Air Temperature in T")
plt.ylabel("Residual Pmpp in %")
plt.legend()
plt.show()

plt.plot(df["aoi"], residualUF, "ro", markersize=1, label="AOI residual with UF")
plt.plot(
    df["aoi"], residualwithoutUF, "go", markersize=1, label="AOI residual without UF"
)
plt.xlabel("Air Temperature in T")
plt.ylabel("Residual Pmpp in %")
plt.legend()
plt.show()
