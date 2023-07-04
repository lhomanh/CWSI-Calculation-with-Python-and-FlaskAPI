
from numba import jit_module
from numpy import empty, full

CACHE = False

cp = 1.011e-3 # MJ * kg^-1 * C^-1 (Campbell & Norman)

def Rnc(Rn):
    '''net radiation at crop surface (MJ*m^-2)


    Arguments
    ---------
    Rn : float
        net radiation at crop surface (MJ*m^-2)


    Returns
    -------
    float
    '''
    return 0.9 * Rn

def cwsi(Ta, Tc, lower, upper):
    '''crop water stress index (unitless)


    Arguments
    ---------
    Ta : float
        air temperature (C)
    Tc : float
        canopy temperature (C)
    lower : float
        lower limit of Tc - Ta (C)
    upper : float
        upper limit of Tc - Ta (C)


    Returns
    -------
    float
    '''
    res = (Tc - Ta - lower) / (upper - lower)
    if res < 0:
        return 0
    elif res > 1:
        return 1
    else:
        return res

def ll(DELTA, ea, es, gamma, ra, rc, upper):
    '''theoretical lower limit of Tc - Ta (C)


    Arguments
    ---------
    DELTA : float
        slope of the saturation vapor pressure curve (kPa*C^-1)
    ea : float
        actual vapor pressure (kPa)
    es : float
        saturation vapor pressure (kPa)
    gamma : float
        psychrometric constant (kPa*C^-1)
    ra : float
        aerodynamic resistance (s*m^-1)
    rc : float
        canopy resistance (s*m^-1)
    upper : float
        upper limit of Tc - Ta (C)


    Returns
    -------
    float
    '''
    n = gamma * (1 + rc / ra) # numerator
    d = DELTA + n # denominator
    return upper * n / d - (es - ea) / d

def ra(DELTA, Rnc, a, b, rho):
    '''semi-empirical aerodynamic resistance (s*m^-1)


    Arguments
    ----------
    DELTA : float
        slope of the saturation vapor pressure to temperature curve (kPa*C^-1)
    Rnc : float
        net radiation at crop surface (MJ*m^-2*h-1)
    a : float
        intercept of non-water-stressed baseline (C)
    b : float
        slope of non-water-stressed baseline (C*kPa^-1)
    rho : float
        density of air (kg*m^-3)


    Returns
    -------
    float
    '''
    # 3600 converts from h*m^-1 to s*m^-1
    return 3600 * rho * cp * a / Rnc / b / (DELTA + 1 / b)

def rc(DELTA, b, gamma, ra):
    '''semi-empirical canopy resistance (s*m^-1)


    Arguments
    ---------
    DELTA : float
        slope of the saturation vapor pressure to temperature curve (kPa*C^-1)
    b : float
        slope of non-water-stressed baseline (C*kPa^-1)
    gamma : float
        psychrometric constant (kPa*C^-1)
    ra : float
        aerodynamic resistance (s*m^-1)


    Returns
    -------
    float
    '''
    return -ra * ((DELTA + 1 / b) / gamma + 1)

def rho(P, Ta):
    '''density of air (kg*m^-3)


    Arguments
    ---------
    P : float
        atmospheric pressure (kPa)
    Ta : float
        air temperature (C)


    Returns
    -------
    float
    '''
    return 3.484 * P / (Ta + 273.15)

def ul(Rnc, ra, rho):
    '''theoretical upper limit of Tc - Ta (C)


    Arguments
    ---------
    Rnc  : float
        net radiation at crop surface (MJ*m^-2*h-1)
    ra : float
        aerodynamic resistance (s*m^-1)
    rho : float
        density of air (kg*m^-3)


    Returns
    -------
    float
    '''
    return Rnc * ra / rho / cp / 3600

jit_module(nopython=True, inline='always', cache=CACHE)
