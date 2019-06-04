


def optical_transmission_losses(dni, aoi):

    '''
    reduces the dni due to optical transmission losses by the lense

    :param numeric: dni at a certain time step
    :param numeric: aoi at a certain time step
    :return: float, new dni minus optical tranmission losses
    '''
    c0=4.54545455e-09
    c1=-2.21212121e-06
    c2=8.09090909e-05
    c3=-2.68463203e-03
    c4=8.62948052e-01


    dni_ot= dni * (c4 + c3* aoi + c2*aoi**2+ c1*aoi**3 + c0 *aoi**4)
    return dni_ot


def UF_corrected_DNI(c1, c2, c3, am, t_ambient, dni):

    '''
    The approach follows the model of Gerstmaier (Quelle), promoting a Utilization Faktor
    that reflects the spectral (AM), temperature (Tamb) and DNI dependencies of
    multijunction cells. The Utilization Factor is multiplied with the DNI.

    UF = c1 * UF(AM) + c2 * UF(DNI) + c3 * UF(Tamb)

    :param dni: incoming DNI
    :return: The altered DNI that includes the parametrization of multijunction
            cells
    '''

    ufam=ufam(am)
    ufdni(dni)
    uftamb(t_ambient)

    UF = c1 * ufam + c2 * ufdni + c3 * uftamb

    return dni * UF

def ufam(am):

    a0 = 2.54991088e-03
    a1 = 6.08723215e-04
    a2 = -1.27238047e-04
    a3 = 8.42391410e-06

    return a0+a1*am + a2*am**2 + a3*am**3

def ufdni(dni):

    d0=-7.05120964e-03
    d1=4.02286998e-05
    d2=-5.01069512e-08
    d3=1.98629580e-11

    return d0 + d1 * am + d2 * am ** 2 + d3 * am ** 3

def uftamb(t_ambient):

    t0=3.51270357e-03
    t1=1.25173402e-05
    t2=-2.12813771e-06
    t3=4.37897948e-08

    return t0 + t1 * am + t2 * am ** 2 + t3 * am ** 3