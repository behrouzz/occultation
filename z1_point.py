from bspice import gfsep
import numpy as np
from datetime import timedelta

d2r = np.pi/180
r2d = 180/np.pi

t1 = '2023-01-01'
t2 = '2024-01-01'
targ1 = 'moon'
targ2 = '4'

adr = 'C:/Moi/_py/Astronomy/Solar System/kernels/'

kernels = [
    adr + 'naif0012.tls',
    adr + 'pck00010.tpc',
    adr + 'de440s.bsp',
    ]


def get_time(targ1, targ2, t1, t2, refval, kernels):

    times = gfsep(t1=t1, t2=t2,
                  targ1=targ1, targ2=targ2,
                  shape1='POINT', shape2='POINT',
                  abcorr='LT+S',
                  relate='<',
                  refval=refval*d2r,
                  step=3600,
                  kernels=kernels,
                  return_utc=True)
    return times


times = get_time(targ1, targ2, t1, t2, refval=1, kernels=kernels)

for i in times:
    dt = (i[1]-i[0]).total_seconds()/3600
    t_mean = i[0] + timedelta(hours=dt)
    print(i[0], '|', i[1], '|', t_mean, '|', round(dt))

