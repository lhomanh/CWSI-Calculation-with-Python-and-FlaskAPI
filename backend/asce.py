
from numba import jit_module
from numpy import abs, arccos, arcsin, cos, empty, exp, isfinite, pi, \
    radians, sin, sqrt, tan, where

import numpy as np
CACHE = False

def Cd(Rn):
    '''tall reference denominator constant (s*m^-1)


    Arguments
    ---------
    Rn : float
        net radiation at crop surface (MJ*m^-2)


    Returns
    -------
    float
    '''
    if Rn > 0:
        return 0.25
    else:
        return 1.7

def DELTA(Ta):
    '''slope of the saturation vapor pressure to temperature curve (kPa*C^-1)


    Arguments
    ---------
    Ta : float
        mean air temperature (C)


    Returns
    -------
    float
    '''
    return 2503 * exp(17.27 * Ta / (Ta + 237.3)) / (Ta + 237.3) ** 2

def DOY(D, M, Y):
    '''day of year


    Arguments
    ---------
    D : int
        day [1, 31]
    M : int
        month [1, 12]
    Y : int
        year, i.e. 1999


    Returns
    -------
    int
    '''
    return D - 32 + int(275 * M / 9) + 2 * int(3 / (M + 1)) \
        + int(M / 100 - Y % 4 / 4 + 0.975)

def ETsz(Cd, Cn, DELTA, G, Rn, Ta, ea, es, gamma, u2):
    '''reference evapotranspiration (mm)


    Arguments
    ---------
    Cd : float
        denominator constant (s*m^-1)
    Cn : int
        numerator constant (K*mm*s^3*Mg^-1)
    DELTA : float
        slope of the saturation vapor pressure to temperature curve (kPa*C^-1)
    G : float
        soil heat flux density (MJ*m^-2)
    Rn : float
        net radiation at crop surface (MJ*m^-2)
    Ta : float
        air temperature (C)
    ea : float
        actual vapor pressure (kPa)
    es : float
        saturation vapor pressure (kPa)
    gamma : float
        psychrometric constant (kPa*C^-1)
    u2 : float
        wind speed at 2 m (m*s^-1)


    Returns
    -------
    float
    '''
    return (0.408 * DELTA * (Rn - G) + gamma * Cn / (Ta + 273) * u2 * \
        (es - ea)) / (DELTA + gamma * (1 + Cd * u2))

def G(Rn):
    '''tall reference soil heat flux (MJ*m^-2)


    Arguments
    ---------
    Rn : float
        net radiation at crop surface (MJ*m^-2)


    Returns
    -------
    float
    '''
    if Rn > 0:
        return Rn * 0.04
    else:
        return Rn * 0.2

def Lm(lon):
    ''' longitude of measurement site (positive degrees west of Greenwich England)


    Arguments
    ---------
    lon : float
        longitude of measurement site (degrees)


    Returns
    -------
    float
    '''
    if lon < 0:
        return abs(lon)
    else:
        return 360 - lon

def P(z):
    '''atmospheric pressure (kPa)


    Arguments
    ---------
    z : int
        station elevation (m)


    Returns
    -------
    float
    '''
    return 101.3 * ((293 - 0.0065 * z) / 293) ** 5.26

def Ra(dr, delta, omega1, omega2, phi):
    '''extraterrestrial radiation (MJ*m^-2)


    Arguments
    ---------
    dr : float
        inverse relative earth-sun distance
    delta : float
        solar declination (rad)
    omega1 : float
        solar time angle at beginning of measurement period (rad)
    omega2 : float
        solar time angle at end of measurement period (rad)
    phi : float
        latitude (rad)


    Returns
    -------
    float
    '''
    return 12 / pi * 4.92 * dr * ((omega2 - omega1) * sin(phi) * sin(delta) \
        + cos(phi) * cos(delta) * (sin(omega2) - sin(omega1)))

def Rn(Rnl, Rns):
    '''net radiation at crop surface (MJ*m^-2)


    Arguments
    ---------
    Rnl : float
        net long-wave radiation, positive skyward (MJ*m^-2)
    Rns : float
        net short-wave radiation, positive skyward (MJ*m^-2)


    Returns
    -------
    float
    '''
    return Rns - Rnl

def Rnl(Ta, ea, fcd):
    '''net long-wave radiation, positive skyward (MJ*m^-2)


    Arguments
    ---------
    Ta : float
        mean air temperature (C)
    ea : float
        actual vapor pressure (kPa)
    fcd : float
        cloudiness factor


    Returns
    -------
    float
    '''
    return 2.042e-10 * fcd * (0.34 - 0.14 * sqrt(ea)) * (Ta + 273.16) ** 4

def Rns(Rs):
    '''net short-wave radiation, positive skyward (MJ*m^-2)


    Arguments
    ---------
    Rs : float
        measured solar radiation (MJ*m^-2)


    Returns
    -------
    float
    '''
    return 0.77 * Rs

def Rso(Ra, z):
    '''clear-sky radiation (MJ*m^-2)


    Arguments
    ---------
    Ra : float
        extraterrestrial radiation (MJ*m^-2)
    z : int
        station elevation (m)


    Returns
    -------
    float
    '''
    return (0.75 + 2e-5 * z) * Ra

def Sc(DOY):
    '''seasonal correction for solar time (h)


    Arguments
    ---------
    DOY : int
        day of year


    Returns
    -------
    float
    '''
    b = 2 * pi * (DOY - 81) / 364
    return 0.1645 * sin(2 * b) - 0.1255 * cos(b) - 0.025 * sin(b)

def beta(delta, omega, phi):
    '''angle of the sun above the horizon (rad)


    Arguments
    ---------
    delta : float
        solar declination (rad)
    omega : float
        solar time angle (rad)
    phi : float
        latitude (rad)


    Returns
    -------
    float
    '''
    return arcsin(sin(phi) * sin(delta) + cos(phi) * cos(delta) * cos(omega))

def dr(DOY):
    '''inverse relative earth-sun distance


    Arguments
    ---------
    DOY : int
        day of year


    Returns
    -------
    float
    '''
    return 1 + 0.033 * cos(DOY * 2 * pi / 365)

def delta(DOY):
    '''solar declination (rad)


    Arguments
    ---------
    DOY : int
        day of year


    Returns
    -------
    float
    '''
    return 0.409 * sin(DOY * 2 * pi / 365 - 1.39)

def ea(RH, es):
    '''actual vapor pressure (kPa)


    Arguments
    ---------
    RH : float
        mean fractional relative humidity
    es : float
        mean saturation vapor pressure (kPa)


    Returns
    -------
    float
    '''
    return RH * es

def es(Ta):
    '''saturation vapor pressure (kPa)


    Arguments
    ---------
    Ta : float
        mean air temperature (C)


    Returns
    -------
    float
    '''
    return 0.6108 * exp(17.27 * Ta / (Ta + 237.3))

def fcd(Rs, Rso, beta, prev):
    '''cloudiness factor


    Arguments
    ---------
    Rs : float
        solar radiation (MJ*m^-2)
    Rso : float
        clear-sky radiation (MJ*m^-2)
    beta : float
        angle of the sun above the horizon (rad)
    prev : float
        previous value of fcd


    Returns
    -------
    float


    Notes
    -----
    To initialize for a clear night, use prev=0.055.
    '''
    if beta < 0.3:
        return prev
    elif Rso > 0: # Will divide by zero at night (shouldn't happen)
        res = Rs / Rso
        if res < 0.3: res = 0.3
        elif res > 1: res = 1
        return 1.35 * res - 0.35
    else: # any ammount of radiation over zero should is the full ammount
        return 1

def gamma(P):
    '''psychrometric constant (kPa*C^-1)


    Arguments
    ---------
    P : float
        atmospheric pressure (kPa)


    Returns
    -------
    float
    '''
    return 0.000665 * P

def omega(Lm, Lz, Sc, h, m, t):
    '''solar time angle at midpoint of measurement period (rad)


    Arguments
    ---------
    Lm : float
        longitude of measurement site (positive degrees west of Greenwich, England)
    Lz : float
        longitude of center of local time zone (positive degrees west of Greenwich,
        England)
    Sc : float
        seasonal correction for solar time (h)
    h : int
        hour of measurement (1, 24)
    m : int
        minute of measurement (1, 60)
    t : float
        length of measurement period (h) (0, 1]
        i.e.: 0.5 for half-hourly, 1 for hourly


    Returns
    -------
    float
    '''
    c = h + m % 60 / 60 - t / 2 # fractional midpoint hour
    return pi / 12 * ((c + 0.06667 * (Lz - Lm) + Sc) - 12)

def omega1(omega, omega2, omegas, t):
    '''solar time angle at beginning of measurement period (rad)


    Arguments
    ----------
    omega : float
        solar time angle (rad)
    omega2 : float
        solar time angle at end of measurement period (rad)
    omegas : float
        sunset hour angle (rad)
    t : float
        length of measurement period (h) (0, 1]
        i.e.: 0.5 for half-hourly, 1 for hourly


    Returns
    -------
    float
    '''
    omega1 = omega - t * np.pi / 24
    if omega1 < -omegas:
        omega1 = -omegas
    elif omega1 > omegas:
        omega1 = omegas
    if omega1 > omega2:
        return omega2
    else:
        return omega1

def omega2(omega, omegas, t):
    '''solar time angle at end of measurement period (rad)


    Arguments
    ---------
    omega : float
        solar time angle (rad)
    omegas : float
        sunset hour angle (rad)
    t : float
        length of measurement period (h) (0, 1]
        i.e.: 0.5 for half-hourly, 1 for hourly


    Returns
    -------
    float
    '''
    omega2 = omega + t * np.pi / 24
    if omega2 < -omegas:
        return -omegas
    elif omega2 > omegas:
        return omegas
    else:
        return omega2

def omegas(delta, phi):
    '''sunset hour angle (rad)


    Arguments
    ---------
    delta :
        solar declination (rad)
    phi :
        latitude (rad)


    Returns
    -------
    float
    '''
    return arccos(-tan(phi) * tan(delta))

def phi(lat):
    '''latitude of measurement site (rad)


    Arguments
    ---------
    lat : float
        latitude of measurement site (degrees)


    Returns
    -------
    float
    '''
    return radians(lat)

jit_module(nopython=True, inline='always', cache=CACHE)
