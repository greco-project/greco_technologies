import cpvtopvlib.uf_preprocessing as pre
from cpvtopvlib import cpvsystem as cpv

import numpy as np
import pvlib
import pandas as pd

import matplotlib.pyplot as plt

df= pd.read_csv('../inputs/InsolightMay2019_filtered.csv', sep=',', index_col=0)

# Converting the index as date
panel_location = pvlib.location.Location(latitude=40.453,longitude=-3.727,
                                         tz=1, altitude=658)


IscDNI_top=1
aoi=df['aoi']
#df['Isc/dii'] = df['Isc']/df['dii']
df=df.replace([np.inf, -np.inf], np.nan)
df=df.fillna(0)
m_low, n_low, m_high, n_high, thld = pre.calc_uf_lines(df['aoi'] ,df['Isc/dii'],datatype = 'airmass', limit=60)


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