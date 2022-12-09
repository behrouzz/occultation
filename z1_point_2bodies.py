from bspice import gfsep
import numpy as np
from datetime import timedelta

d2r = np.pi/180
r2d = 180/np.pi

t1 = '2023-01-01'
t2 = '2024-01-01'


adr = 'C:/Moi/_py/Astronomy/Solar System/kernels/'

kernels = [
    adr + 'naif0012.tls',
    adr + 'pck00010.tpc',
    adr + 'de440s.bsp',
    ]


def time_intersection(t1start, t1end, t2start, t2end):
    return (t1start <= t2start <= t1end) or (t2start <= t1start <= t2end)
    

def conj2sso(body1, body2, t1, t2, maxval, kernels):
    """
    Conjucation times of two Solar System Objects (SSO) wrt Earth

    Arguments
    ---------
        body1 (str)    : first SSO
        body2 (str)    : second SSO
        t1 (str)       : start time
        t2 (str)       : end time
        maxval (float) : maximum angular distant (degree)
        kernels (str)  : address of kernel files
        
    Returns
    -------
        times  : list of times distance of the SSOs is less than maxval
                 Each member of the list is a tuple containing:
                   - time that conjunction starts
                   - time that conjunction finishes
                   - duration of conjunction (in hours)
                   - median time (between start and end)
    """

    t = gfsep(t1=t1, t2=t2,
              targ1=body1, targ2=body2,
              shape1='POINT', shape2='POINT',
              abcorr='LT+S',
              relate='<',
              refval=maxval*d2r,
              step=3600,
              kernels=kernels,
              return_utc=True)
    times = []
    for i in t:
        dt = (i[1]-i[0]).total_seconds()/3600
        t_mean = i[0] + timedelta(hours=dt/2)
        times.append((i[0], i[1], dt, t_mean))
    return times



times = conj2sso(body1='moon', body2='4', t1=t1, t2=t2, maxval=10, kernels=kernels)

for i in times:
    print(str(i[0])[:19], '|', str(i[1])[:19], '|', round(i[2]))

print()

