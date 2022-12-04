import spiceypy as sp
import numpy as np
from datetime import datetime, timedelta


d2r = np.pi/180
r2d = 180/np.pi

RADIUS_MARS = 3389.5
RADIUS_MOON = 1737.4



def lonlat2cart(obs_loc):
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
    # inner function
    def limit_range(t0, rng, n):
        t_win = [t0 + timedelta(seconds=i) for i in np.linspace(-rng, rng, n)]
        d = np.zeros(len(t_win))
        for i in range(len(t_win)):
            [_,x1,y1], [_,x2,y2] = get_apparent_bodies(bodies, t_win[i], obs_loc, kernels)
            d[i] = ((x1-x2)**2 + (y1-y2)**2) ** 0.5
        ind = np.argmin(np.abs(d))
        return t_win[ind]
    
    t = limit_range(t0, 86400*2, 5)
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

        size_moon = (2 * np.arctan(RADIUS_MOON/dist_moon)*r2d) * 3600
        size_mars = (2 * np.arctan(RADIUS_MARS/dist_mars)*r2d) * 3600

        ls_size_moon.append(size_moon)
        ls_size_mars.append(size_mars)

        dist_tang = (size_moon/2) + (size_mars/2)

        d = ((az_mars-az_moon)**2 + (alt_mars-alt_moon)**2)**0.5
        d = d*3600

        delta.append(d - dist_tang)

    immersion = min([i for i in delta if i>0])
    emersion = max([i for i in delta if i<0])

    ind_im = delta.index(immersion)
    ind_em = delta.index(emersion)

    t_im = time_window[ind_im]
    t_em = time_window[ind_em]
    
    return time_window, ind_im, ind_em, pos_moon, pos_mars, np.array(ls_size_moon), np.array(ls_size_mars)
