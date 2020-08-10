import pandas as pd
from pvlib.pvsystem import PVSystem
import matplotlib.pyplot as plt
import pvlib.atmosphere
import pvlib.pvsystem as pvsystem
import greco_technologies.cpv.hybrid as hybrid

"""
note that Ins-Si cannot be validated yet, because there is not sufficient 
measurement data available
"""

df = pd.read_csv("../inputs/InsolightMay2019_filtered.csv", sep=",", index_col=0)
df.index = pd.to_datetime(df.index)

df["dhi"] = df["ghi"] - df["dni"]

sandia_modules = pvlib.pvsystem.retrieve_sam("SandiaMod")
module_ref = sandia_modules["Canadian_Solar_CS5P_220M___2009_"]
module_ref["Cells_in_Series"] = 4
module_ref["Cells_in_Parallel"] = 1
module_ref["Isco"] = 1
module_ref["Voco"] = 4
module_ref["Impo"] = 6.52
module_ref["Vmpo"] = 5

system_ref = PVSystem(
    surface_tilt=30,
    surface_azimuth=180,
    module_parameters=module_ref,
    inverter_parameters=None,
)

# prepare dictionary with poa_global, poa_direct_poa_diffuse, absolute_airmass, aoi
hybrid_weather = hybrid.hybrid_weather_data(
    df, lat=45.641603, lon=-3.727, surface_tilt=30, surface_azimuth=180
)


# todo: adapt function for effective irradiance to increase diffuse fraction
# todo: adjust sapm_spectral_loss(airmass_absolute, module) and sapm_aoi_loss(aoi, module) that are included in the effective_irradiance


effective_irradiance = pvsystem.sapm_effective_irradiance(
    poa_direct=hybrid_weather["dni"],
    poa_diffuse=hybrid_weather["dhi"],
    airmass_absolute=hybrid_weather["airmass"],
    aoi=hybrid_weather["aoi"],
    module=module_ref,
)

df["effective_irrandiance"] = effective_irradiance
from pvlib.temperature import sapm_cell, TEMPERATURE_MODEL_PARAMETERS

temp_params = TEMPERATURE_MODEL_PARAMETERS["sapm"]["open_rack_glass_glass"]

temp_cell = pvsystem.temperature.sapm_cell(
    poa_global=hybrid_weather["ghi"],
    wind_speed=hybrid_weather["wind_speed"],
    temp_air=hybrid_weather["temp_air"],
    **temp_params
)  # todo: adjust model
temp_cell = temp_cell.fillna(0)
output = pvsystem.sapm(
    effective_irradiance=effective_irradiance, temp_cell=temp_cell, module=module_ref
)

modeled_power = output.p_mp
modeled_isc = output.i_sc

real_power = df["Pmp_Si"]
real_isc = df["Isc_Si"]

plt.plot(
    hybrid_weather["aoi"], effective_irradiance, "b.", label="effective irradiance"
)
# plt.plot(aoi_list, weather_loc['perez_diffuse'], 'g.', label='DHI')
plt.xlabel("AOI")
plt.ylabel("Irradiance")
plt.legend()
plt.show()


plt.plot(modeled_power, "b", label="modeled power")
plt.plot(real_power, "r", label="measured power")
plt.xlabel("time")
plt.ylabel("power in W")
plt.legend()
plt.show()

plt.plot(hybrid_weather["aoi"], modeled_power, "b.", label="modeled_power")
plt.plot(hybrid_weather["aoi"], real_power, "r.", label="measured power")
plt.xlabel("AOI")
plt.ylabel("power in W")
plt.legend()
plt.show()
