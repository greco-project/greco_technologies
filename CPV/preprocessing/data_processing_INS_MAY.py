# Carga y Procesado de Datos sin influencia de la Temperatura Ambiente
import numpy as np
import regression_analysis as reg

# Cálculos de la Masa de Aire
import numpy as np
import pvlib
import pandas as pd
import datetime

df= pd.read_csv('/home/local/RL-INSTITUT/inia.steinbach/rl-institut/04_Projekte/220_GRECO/03-Projektinhalte/AP4_High_Penetration_of_Photovoltaics/T4_3_CPV/INS/MAY/InsolightMay2019_filtered.csv', sep=',', index_col=0)

# Converting the index as date
panel_location = pvlib.location.Location(latitude=40.453,longitude=-3.727,
                                         tz=1, altitude=658)


IscDNI_top=0.96/1000
aoi=df['aoi']
#df['Isc/dii'] = df['Isc']/df['dii']
df=df.replace([np.inf, -np.inf], np.nan)
df=df.fillna(0)
m_low, n_low, m_high, n_high, thld = reg.calc_uf_lines(df['aoi'] ,df['Isc/dii'],datatype = 'airmass', limit=60)


x1 = np.arange(10,70,0.1)
y1 = m_low * x1 + n_low
x2 = np.arange(55,85,0.1)
y2 = m_high * x2 + n_high


print("thld_aoi = ", thld, '\n'
       'm_low_aoi = ',  m_low/IscDNI_top, '\n'
      'm_high_aoi = ', m_high/IscDNI_top)

import matplotlib.pyplot as plt
plt.plot(aoi, df['Isc/dii'], 'b.', x1, y1, 'g', x2, y2, 'r')
plt.xlabel('Ángulo de Incidencia (º)')
plt.ylabel('Isc/DII (A/(W/m2))')
plt.title('Análisis de Isc/DII en función del Ángulo de Incidencia')
#plt.savefig("grafica1.png", dpi=300)
plt.show()