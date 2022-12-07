from bspice import gfsep
import numpy as np

d2r = np.pi/180
r2d = 180/np.pi

t1 = '2022-01-01'
t2 = '2023-01-01'

adr = 'C:/Moi/_py/Astronomy/Solar System/kernels/'

kernels = [
    adr + 'naif0012.tls',
    adr + 'pck00010.tpc',
    adr + 'de440s.bsp',
    ]


for obj in [2, 4, 5, 6, 7]:
    print('NAME:', obj)
    targ1 = 'moon'
    targ2 = str(obj)
    refval = 0.6*d2r

    times = gfsep(t1=t1, t2=t2,
                  targ1=targ1, targ2=targ2,
                  shape1='POINT', shape2='POINT',
                  abcorr='LT+S',
                  relate='<',
                  refval=refval,
                  step=3600,
                  kernels=kernels,
                  return_utc=True)

    times = [i[0] for i in times]

    for i in times:
        print(i)
    print('-'*70)
