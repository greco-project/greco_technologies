import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math


def optical_transmission_losses(aoi):

    '''
    reduces the dni due to optical transmission losses by the lens

    :param numeric: dni at a certain timestep
    :param numeric: aoi at a certain timestep
    :return: float, new dni minus optical tranmission losses
    '''



    c0=4.54545455e-09
    c1=-2.21212121e-06
    c2=8.09090909e-05
    c3=-2.68463203e-03
    c4=8.62948052e-01

    if aoi >=60:
        return 0
    else:
        ot= c4 + c3* aoi + c2*aoi**2+ c1*aoi**3 + c0 *aoi**4
        return ot


def glass_transmission_losses(aoi): #todo: find mistake for aoi=59°

    # this is an insolight specific anti-reflection coating
    glass_ar_offset = 0.015
    n1 = 1.0
    n2 = 1.5

    if aoi <= 90:
        theta = aoi
    else:
        return 0

    Rs = abs((n1 * math.cos(math.radians(theta)) - n2 * np.sqrt(1-((n1/n2) * np.sin(math.radians(theta)))**2.0))/
             (n1 * np.cos(math.radians(theta)) + n2 * np.sqrt(1-((n1/n2) * np.sin(math.radians(theta)))**2.0)))**2

    Rp= abs((n1 * np.sqrt(1 - ((n1/n2) * np.sin(math.radians(theta)))**2.0) - n2 * np.cos(math.radians(theta)))/
            (n1 * np.sqrt(1 - ((n1/n2) * np.sin(math.radians(theta)))**2.0) + n2 * np.cos(math.radians(theta))))**2

    Reff= 0.5 *(Rs + Rp)
    glass_transmission = (1.0 - Reff + glass_ar_offset)**2
    return glass_transmission

def ufam(am):

    a0 =1.0015288
    a1=-0.00878182
    a2=-0.00146694

    return a0+a1*am + a2*am**2

def ufdni(dni):

    d0=-2.22040018e-01
    d1=3.08399124e-03
    d2=-1.95350121e-06

    return d0 + d1 * dni + d2 * dni ** 2

def uftamb(t_ambient):

    t0=8.43234532e-01
    t1=1.04128138e-02
    t2=-1.73841969e-04

    return t0 + t1 * t_ambient + t2 * t_ambient ** 2

def UF_corrected_DNI(am, t_ambient, dni):

    '''
    The approach follows the model of Gerstmaier (Quelle), promoting a Utilization Faktor
    that reflects the spectral (AM), temperature (Tamb) and DNI dependencies of
    multijunction cells. The Utilization Factor is multiplied with the DNI.

    UF = c1 * UF(AM) + c2 * UF(DNI) + c3 * UF(Tamb)

    :param dni: incoming DNI
    :return: The altered DNI that includes the parametrization of multijunction
            cells
    '''

    uf_am=ufam(am)
    uf_dni=ufdni(dni)
    uf_tamb=uftamb(t_ambient)
    c1=1
    c2=1
    c3=1

    df_new = pd.concat([uf_am, uf_dni, uf_tamb], axis=1)
    df=df_new.sum(axis=1)

#    UF = c1 * uf_am + c2 * uf_dni + c3 * uf_tamb

    return dni * df

if __name__ == "__main__":

    aoi_list={}

    aoi=1
    while aoi < 200:
        gt=glass_transmission_losses(aoi)
        aoi_list[aoi]=gt
        aoi=aoi+1


    #gt=glass_transmission_losses(aoi)
    #print(gt)
    get=pd.DataFrame.from_dict(aoi_list, orient='index')
    print(get)
    plt.plot(get)
    plt.show()