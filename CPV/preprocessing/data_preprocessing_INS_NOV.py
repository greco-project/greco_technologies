# Carga y Procesado de Datos sin influencia de la Temperatura Ambiente
import numpy as np
import regression_analysis as reg

# Cálculos de la Masa de Aire
import numpy as np
import pvlib
import pandas as pd
import datetime

data= pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/INS/INS_Nov_twoaxistracking/filtered_dataset_November_Marcos.csv', sep=',', index_col=0)
df_fixairmass=pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/INS/INS_Nov_twoaxistracking/filtered_dataset_November_fixairmass.csv', sep=',', index_col=0)
df_fixtemp=pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/INS/INS_Nov_twoaxistracking/filtered_dataset_November_fixtemperature.csv', sep=',', index_col=0)


# Converting the index as date
panel_location = pvlib.location.Location(latitude=40.453,longitude=-3.727,
                                         tz=1, altitude=658)


fixtemp_IscDNI = df_fixtemp['Isc/DNI']
fixtemp_airmass = df_fixtemp['relative_airmass']


# power = data['Isc/DNI']
# airmass = data['relative_airmass']
# temperature = data['temp']
#
# new_data = pd.DataFrame()
# new_data['airmass'] = airmass
# new_data['temp'] = temperature
# new_data['Isc/DNI'] = power

median_df = pd.Series()
for j in np.arange(1, 5, 0.1):
    am_data = df_fixtemp[
        (df_fixtemp['relative_airmass'] > j - 0.05) & (df_fixtemp['relative_airmass'] > j + 0.05)]
    median_df[j] = am_data['Isc/DNI'].median()
    median_Isc = median_df.tolist()

m_low_am, n_low_am, m_high_am, n_high_am, thld_am = reg.calc_two_regression_lines(
    median_df.index,
    median_Isc, limit=4.0)

x1 = np.arange(1, 5, 0.1)
y1 = m_low_am * x1 + n_low_am
x2= x = np.arange(4, 5, 0.1)
y2 = m_high_am * x2 + n_high_am


import matplotlib.pyplot as plt

plt.plot(df_fixtemp['relative_airmass'], df_fixtemp['Isc/DNI'], 'b+', median_df.index, median_df, 'r.', x1, y1,
         'g', x2, y2, 'r')
plt.xlabel("Airmass")
plt.ylabel("Isc/DNI")
plt.show()

IscDNI_top = 0.96/1000
print("thld_am = ", thld_am, '\n'
    'm_low_am = ', m_low_am / IscDNI_top, '\n'
    'm_high_am = ', m_high_am / IscDNI_top)

uf_am = pd.Series()
for i, row in data.iterrows():
    uf_am[i] = reg.get_single_util_factor(row['relative_airmass'], thld_am,
                                      m_low_am / IscDNI_top,
                                      m_high_am / IscDNI_top)

# Carga y Procesado de Datos sin influencia de la Masa de Aire

m_low_temp, n_low_temp, m_high_temp, n_high_temp, thld_temp = reg.calc_uf_lines(
    df_fixairmass['temp'],
    df_fixairmass['Isc/DNI'],
    datatype='temp_air')
x = np.arange(10, 25, 1)
y1 = m_low_temp * x + n_low_temp

# calculate median for each temperature
#Isc_median_temp = data.groupby([data['temp']]).median()


plt.plot(df_fixairmass['temp'], df_fixairmass['Isc/DNI'], 'g+', x, y1, 'b')
plt.xlabel("Temperature")
plt.ylabel("Isc/DNI")
plt.show()

print("thld_temp = ", thld_temp, '\n'
    'm_low_temp = ', m_low_temp / IscDNI_top, '\n'
    'm_high_temp = ', m_high_temp / IscDNI_top)

uf_temp = pd.Series()
for i, row in data.iterrows():
    uf_temp[i] = reg.get_single_util_factor(row['temp'], thld_temp,
                                        m_low_temp / IscDNI_top,
                                        m_high_temp / IscDNI_top)



