import numpy as np
import cpvtopvlib.uf_preprocessing as pre
from cpvtopvlib import cpvsystem as cpv

import numpy as np
import pvlib
import datetime

import matplotlib.pyplot as plt

data = np.loadtxt("../inputs/m300_data_filtered.txt", delimiter=",")

datetimestring = np.genfromtxt(
    "../inputs/m300_datetime.txt", dtype="str", delimiter="\n"
)

datetimeobject = []
for i in range(len(datetimestring)):
    datetimeobject.append(
        datetime.datetime.strptime(datetimestring[i], "%d-%b-%Y %H:%M:%S")
    )

panel_location = pvlib.location.Location(
    latitude=45.641603, longitude=5.875387, tz=1, altitude=234
)

Airmass = panel_location.get_airmass(times=datetimeobject)

airmass_array = np.array(Airmass["airmass_relative"])
relative_airmass = np.zeros((len(airmass_array), 1))
for i in range(len(airmass_array)):
    relative_airmass[i, 0] = airmass_array[i]

data = np.append(data, relative_airmass, 1)

np.savetxt(
    fname="../inputs/m300_data_filtered_complete.txt",
    X=data,
    delimiter=",",
    fmt="%.10f",
)


nontemp_data = np.loadtxt("../inputs/m300_nontemp_measurements.txt", delimiter=",")

nontemp_IscDNI = nontemp_data[:, 25]
nontemp_airmass = nontemp_data[:, 33]

import statistics as stats

IscDNI_medians = []
Airmass_aux = []
for i in np.arange(1, 2.8, 0.1):
    array_aux1 = []
    for j in range(len(nontemp_IscDNI)):
        if nontemp_airmass[j] > i - 0.05 and nontemp_airmass[j] < i + 0.05:
            array_aux1.append(nontemp_IscDNI[j])
    if len(array_aux1) > 0:
        IscDNI_medians.append(stats.median(array_aux1))
        Airmass_aux.append(i)

m_low, n_low, m_high, n_high, thld = pre.calc_two_regression_lines(
    Airmass_aux, IscDNI_medians, limit=2.1
)
# m_low, n_low, error1 = calc_regression_line(Airmass_aux[:10],
# IscDNI_medians[:10])
# m_high, n_high, error2 = calc_regression_line(Airmass_aux[10:],
# IscDNI_medians[10:])
# thld = (n_high - n_low) / (m_low - m_high)

x = np.arange(1, 5, 0.1)
y1 = m_low * x + n_low
y2 = m_high * x + n_high

plt.plot(
    nontemp_airmass,
    nontemp_IscDNI,
    "b.",
    Airmass_aux,
    IscDNI_medians,
    "g.",
    x,
    y1,
    "g",
    x,
    y2,
    "r",
)
plt.show()

IscDNI_ast = 3.346 / 1000
uf_am = []

print(
    "thld_am = ",
    thld,
    "\n" "m_low_am = ",
    m_low / IscDNI_ast,
    "\n" "m_high_am = ",
    m_high / IscDNI_ast,
)

for i in range(len(airmass_array)):
    uf_am.append(
        cpv.get_simple_util_factor(
            airmass_array[i], thld, m_low / IscDNI_ast, m_high / IscDNI_ast
        )
    )

plt.plot(uf_am)
plt.show()

# Carga y Procesado de Datos sin influencia de la Masa de Aire

nonairmass_data = np.loadtxt(
    "../inputs/m300_nonairmass_measurements.txt", delimiter=","
)

nonairmass_IscDNI = nonairmass_data[:, 25]
nonairmass_temp = nonairmass_data[:, 10]
m_low, n_low, m_high, n_high, thld = pre.calc_uf_lines(
    nonairmass_temp, nonairmass_IscDNI, "temp_air"
)

AmbientTemp = data[:, 10]

uf_at = []
for i in range(len(airmass_array)):
    uf_at.append(
        cpv.get_simple_util_factor(
            AmbientTemp[i], thld, m_low / IscDNI_ast, m_high / IscDNI_ast
        )
    )
print(
    "thld_temp = ",
    thld,
    "\n" "m_low_temp = ",
    m_low / IscDNI_ast,
    "\n" "m_high_temp = ",
    m_high / IscDNI_ast,
)
