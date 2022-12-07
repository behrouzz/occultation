import spiceypy as sp
import numpy as np
from datetime import datetime, timedelta


d2r = np.pi/180
r2d = 180/np.pi

RADIUS_MARS = 3389.5
RADIUS_MOON = 1737.4


def angular_separation(r1, d1, r2, d2):
    """
    Calculate angular separation between two point
    Arguments
    ---------
        r1 (float): right ascension of the first point in degrees
        d1 (float): declination of the first point in degrees
        r2 (float): right ascension of the second point in degrees
        d2 (float): declination of the second point in degrees
    Returns
    -------
        angular sepration in degrees
    """
    r1 = (np.pi/180) * r1
    d1 = (np.pi/180) * d1
    r2 = (np.pi/180) * r2
    d2 = (np.pi/180) * d2

    radical = np.sqrt(
        (np.cos(d2)**2)*(np.sin(r2-r1)**2) + ((np.cos(d1)*np.sin(d2) - np.sin(d1)*np.cos(d2)*np.cos(r2-r1))**2)
        )
    
    kasr = radical / ( (np.sin(d1)*np.sin(d2)) + (np.cos(d1)*np.cos(d2)*np.cos(r2-r1)) )

    sep = (180/np.pi) * np.arctan(kasr)
    
    return sep


def lonlat2cart(obs_loc):
    if len(obs_loc)==2:
        obs_loc = (obs_loc[0], obs_loc[1], 0)
    lon, lat, alt = obs_loc
    re = 6378.1366
    rp = 6356.7519
    f = (re-rp)/re
    obspos = sp.pgrrec(body='earth', lon=lon*d2r, lat=lat*d2r, alt=alt/1000, re=re, f=f)
    return obspos


def get_apparent_bodies(bodies, t, obs_loc, kernels, abcorr='LT+S'):

    bodies = [str(i) for i in bodies]

    for k in kernels:
        sp.furnsh(k)

    et = sp.str2et(str(t))
    obspos = lonlat2cart(obs_loc)

    r_az_alt = np.zeros((len(bodies),3))

    for i in range(len(bodies)):
        state, lt  = sp.azlcpo(
            method='ELLIPSOID',
            target=bodies[i],
            et=et,
            abcorr=abcorr,
            azccw=False,
            elplsz=True,
            obspos=obspos,
            obsctr='earth',
            obsref='ITRF93')
        r, az, alt = state[:3]
        r_az_alt[i,:] = r, az*r2d, alt*r2d

    sp.kclear()

    return r_az_alt


def create_range(t0, steps, dt):
    t1 = t0 - timedelta(seconds=dt)
    t2 = t0 + timedelta(seconds=dt)
    rng = t2 - t1
    dt = rng / steps
    return [t1 + dt*i for i in range(steps+1)]


def find_exact_t0(t0, bodies, obs_loc, kernels):
    """
    Finds the best moment of occulation based on an inial guess (t0)
    """
    
    def limit_range(t0, rng, n):
        t_win = [t0 + timedelta(seconds=i) for i in np.linspace(-rng, rng, n)]
        new_t_win = []
        dist = []
        for i in range(len(t_win)):
            [_,az1,alt1], [_,az2,alt2] = get_apparent_bodies(bodies, t_win[i], obs_loc, kernels)

            if (alt1>0) and (alt2>0):
                new_t_win.append(t_win[i])
                d = angular_separation(az1, alt1, az2, alt2)
                dist.append(d)
        dist = np.array(dist)
        if len(dist)>0:
            ind = np.argmin(np.abs(dist))
            return new_t_win[ind]
        else:
            return None
    
    t = limit_range(t0, 86400*2, 5)
    if t is None:
        return None
    else:
        t = limit_range(t, 3600*24, 24)
        t = limit_range(t, 60*60, 60)
        t = limit_range(t, 1*60, 60)
        return t


def get_occultation(bodies, obs_loc, t0, kernels, n_inter=100, dt=4000):

    time_window = create_range(t0, n_inter, dt)
    
    pos_mars = []
    pos_moon = []
    ls_size_moon = []
    ls_size_mars = []
    delta = []

    for t in time_window:
        
        moon, mars = get_apparent_bodies(bodies, t, obs_loc, kernels, abcorr='LT+S')
        pos_mars.append(mars)
        pos_moon.append(moon)

        dist_mars, az_mars, alt_mars = mars
        dist_moon, az_moon, alt_moon = moon

        size_moon = (2 * np.arctan(RADIUS_MOON/dist_moon)*r2d)
        size_mars = (2 * np.arctan(RADIUS_MARS/dist_mars)*r2d)

        ls_size_moon.append(size_moon)
        ls_size_mars.append(size_mars)

        dist_tang = (size_moon/2) + (size_mars/2)

        d = angular_separation(az_moon, alt_moon, az_mars, alt_mars)

        delta.append(d - dist_tang)

    immersion = min([i for i in delta if i>0])
    emersion = max([i for i in delta if i<0])

    ind_im = delta.index(immersion)
    ind_em = delta.index(emersion)

    t_im = time_window[ind_im]
    t_em = time_window[ind_em]
    
    return time_window, ind_im, ind_em, pos_moon, pos_mars, np.array(ls_size_moon), np.array(ls_size_mars)
