import spiceypy as sp
import numpy as np
from datetime import datetime


d2r = np.pi/180
r2d = 180/np.pi

t1 = '2023-01-01'
t2 = '2024-01-01'

adr = 'C:/Moi/_py/Astronomy/Solar System/kernels/'

kernels = [
    adr + 'naif0012.tls',
    adr + 'pck00010.tpc',
    adr + 'de421.bsp',
    ]


targ1 = 'moon'
targ2 = '499'
abcorr='LT+S'
relate='ABSMIN'
refval=0.0
step=3600
return_utc = True


#def gfsep(t1, t2, targ1, targ2, shape1, shape2, abcorr, relate, refval, step, kernels, return_utc):

for k in kernels:
    sp.furnsh(k)

MAXWIN = 1000
result = sp.Cell_Double(2*MAXWIN)
cnfine = sp.Cell_Double(2)

et1 = sp.str2et(t1)
et2 = sp.str2et(t2)

sp.wninsd(et1, et2, cnfine)

sp.gfsep(
    targ1=targ1,
    shape1='SPHERE',
    inframe1='NULL',
    targ2=targ2,
    shape2='SPHERE',
    inframe2='NULL',
    abcorr=abcorr,
    obsrvr='EARTH',
    relate=relate,
    refval=refval,
    adjust=0.0,
    step=step,
    nintvls=1000,
    cnfine=cnfine,
    result=result
    )

count = sp.wncard(result)

times = []

for i in range(count):
    t1_tdb, t2_tdb = sp.wnfetd(result, i)
    
    if return_utc:
        t1_utc = sp.et2utc(t1_tdb, 'ISOC', 0, 20)
        t2_utc = sp.et2utc(t2_tdb, 'ISOC', 0, 20)
        t1_utc = datetime.strptime(t1_utc, '%Y-%m-%dT%H:%M:%S')
        t2_utc = datetime.strptime(t2_utc, '%Y-%m-%dT%H:%M:%S')
        times.append((t1_utc, t2_utc))
    else:
        times.append((t1_tdb, t2_tdb))
    
sp.kclear()

for i in times:
    print(i[0])
