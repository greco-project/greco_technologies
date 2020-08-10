import cpvtopvlib.uf_preprocessing as pre
from cpvtopvlib import cpvsystem as cpv

import numpy as np
import pvlib
import pandas as pd

import matplotlib.pyplot as plt

data = pd.read_csv(
    "../inputs/filtered_dataset_November_Marcos.csv", sep=",", index_col=0
)
df_fixairmass = pd.read_csv(
    "../inputs/filtered_dataset_November_fixairmass.csv", sep=",", index_col=0
)
df_fixtemp = pd.read_csv(
    "../inputs/filtered_dataset_November_fixtemperature.csv", sep=",", index_col=0
)


panel_location = pvlib.location.Location(
    latitude=40.453, longitude=-3.727, tz=1, altitude=658
)


fixtemp_IscDNI = df_fixtemp["Isc/DNI"]
fixtemp_airmass = df_fixtemp["relative_airmass"]


median_df = pd.Series()
for j in np.arange(2, 5, 0.1):
    am_data = df_fixtemp[
        (df_fixtemp["relative_airmass"] > j - 0.05)
        & (df_fixtemp["relative_airmass"] > j + 0.05)
    ]
    median_df[j] = am_data["Isc/DNI"].median()
    median_Isc = median_df.tolist()

m_low_am, n_low_am, m_high_am, n_high_am, thld_am = pre.calc_two_regression_lines(
    median_df.index, median_Isc, limit=4.0
)

x1 = np.arange(2, 5.1, 0.1)
y1 = m_low_am * x1 + n_low_am
x2 = np.arange(4, 8, 0.1)
y2 = m_high_am * x2 + n_high_am

plt.plot(
    df_fixtemp["relative_airmass"],
    df_fixtemp["Isc/DNI"],
    "b+",
    median_df.index,
    median_Isc,
    "r.",
    x1,
    y1,
    "g",
    x2,
    y2,
    "r",
)
plt.xlabel("Airmass")
plt.ylabel("Isc/DNI")
plt.show()

IscDNI_top = 0.96 / 1000
print(
    "thld_am = ", thld_am, "\n" "m_low_am = ", m_low_am, "\n" "m_high_am = ", m_high_am
)

uf_am = pd.Series()
for i, row in data.iterrows():
    uf_am[i] = cpv.get_simple_util_factor(
        row["relative_airmass"], thld_am, m_low_am / IscDNI_top, m_high_am / IscDNI_top
    )


m_low_temp, n_low_temp, m_high_temp, n_high_temp, thld_temp = pre.calc_uf_lines(
    df_fixairmass["temp"], df_fixairmass["Isc/DNI"], datatype="temp_air"
)
x = np.arange(10, 25, 1)
y1 = m_low_temp * x + n_low_temp

# calculate median for each temperature
# Isc_median_temp = data.groupby([data['temp']]).median()


plt.plot(df_fixairmass["temp"], df_fixairmass["Isc/DNI"], "g+", x, y1, "b")
plt.xlabel("Temperature")
plt.ylabel("Isc/DNI")
plt.show()

print(
    "thld_temp = ",
    thld_temp,
    "\n" "m_low_temp = ",
    m_low_temp,
    "\n" "m_high_temp = ",
    m_high_temp,
)

uf_temp = pd.Series()
for i, row in data.iterrows():
    uf_temp[i] = cpv.get_simple_util_factor(
        row["temp"], thld_temp, m_low_temp / IscDNI_top, m_high_temp / IscDNI_top
    )
