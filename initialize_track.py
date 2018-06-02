## Import modules
import math; import datetime; import numpy as np; import matlab.engine; import os;

def solve_kepler(M,e):
    ## Purpose: Given mean anomaly and eccentricity, return eccentric anomaly
    ## Inputs:
    ##     M: mean anomaly
    ##     e: eccentricity

    ## Define equations
    f = lambda x: x - e*math.sin(x) - M
    f_prime = lambda x: 1 - e*math.cos(x)

    ## Pick guess
    if M < math.pi:
        E = M + e/2
    else:
        E = M - e/2

    ## Loop until we are close to the root
    ratio = f(E)/f_prime(E)
    while abs(ratio) > 1e-8:
        E -= ratio
        ratio = f(E)/f_prime(E)

    if E > 2.0*math.pi:
        two_pi = 2.0*math.pi
        rem = E % two_pi
        E = rem

    return E

def calc_w(E,e):
    ## Purpose: Given eccentric anomaly and eccentricity, return true anomaly

    ## Inputs:
    ##   E: eccentric anomaly
    ##   e: eccentricity

    e_coef = math.sqrt((1-e)/(1+e))
    w = 2.0*math.atan(math.tan(E/2.0)/e_coef)

    if w < 0:
        w += 2.0*math.pi

    return w

def peri2geo(O,i,wp):
    ## Purpose: Given RAAN, inclination, and argument of perigee, return rotation matrix from perifocal to geocentric frame

    ## Inputs:
    ##   O: RAAN (rad)
    ##   i: inclination (rad)
    ##   wp: argument of perigee (rad)

    q11 = -math.sin(O)*math.cos(i)*math.sin(wp) + math.cos(O)*math.cos(wp)
    q12 = -math.sin(O)*math.cos(i)*math.cos(wp) - math.cos(O)*math.sin(wp)
    q13 = math.sin(O)*math.sin(i)
    q21 = math.cos(O)*math.cos(i)*math.sin(wp) + math.sin(O)*math.cos(wp)
    q22 = math.cos(O)*math.cos(i)*math.cos(wp) - math.sin(O)*math.sin(wp)
    q23 = -math.cos(O)*math.sin(i)
    q31 = math.sin(i)*math.sin(wp)
    q32 = math.sin(i)*math.cos(wp)
    q33 = math.cos(i)

    Q = np.array([[q11,q12,q13],[q21,q22,q23],[q31,q32,q33]])
    return Q

def geo2ecef(time):
    ## Purpose: Return the rotation matrix from geocentrix to ecef at the input times

    y = time.year
    m = time.month
    d = time.day
    hr = time.hour
    min = time.minute
    sec = time.second

    j0 = 367.0*y - math.trunc((7.0*(y+math.trunc((m+9)/12.0)))/4.0) + math.trunc(275.0*m/9.0) + d + 1721013.5
    T0 = (j0 - 2451545.0)/36525.0
    thet_g0 = 100.4606184 + 36000.77004*T0 + 0.000387933*math.pow(T0,2) - 2.583e-8*math.pow(T0,3)

    thet_g0 = thet_g0 % 360.0
    if thet_g0 < 0:
        thet_g0 += 360

    UT = hr + min/60.0 + sec/3600.0
    thet_g = thet_g0 + 360.98564724*UT/24
    thet_g = thet_g*math.pi/180.0

    Q = [[math.cos(thet_g),math.sin(thet_g),0],[-math.sin(thet_g),math.cos(thet_g),0],[0,0,1]]
    return Q

def get_latlon(r):
    ## Purpose: Given radius vector in ECEF, return latiude and longitude

    ## Inputs:
    ##   r: radius vector in ECEF frame (m)

    x = r[0]
    y = r[1]
    z = r[2]

    r_mag = math.sqrt(math.pow(x,2)+math.pow(y,2)+math.pow(z,2))

    l = x/r_mag
    m = y/r_mag
    n = z/r_mag

    dec = math.asin(n)
    if m > 0:
        ra = math.acos(l/math.cos(dec))
    else:
        ra = 2.0*math.pi - math.acos(l/math.cos(dec))

    return dec,ra

## Define constants
J2 = 1.08263e-3;
G = 6.67e-11;
Me = 5.972e24;
Mew = G*Me;
Re = 6378.37e3; #m
eng = matlab.engine.start_matlab()
eng.addpath(os.getcwd()+"/")

## Read text file with orbit initial conditions (ic.kep)
with open("ic.kep","r") as f:
    ## Format is:
    ##   epoch (UTC)
    ##   semi-major axis (m)
    ##   eccentricity
    ##   inclination (deg)
    ##   true anomaly (deg)
    ##   RAAN (deg)
    ##   argument of perigee (deg)
    ic = f.read()
    ic = ic.split("\n")
    epoch = ic[0]
    a0 = float(ic[1])
    e0 = float(ic[2])
    i0 = float(ic[3])
    w0 = float(ic[4])
    O0 = float(ic[5])
    wp0 = float(ic[6])

## Determine rate of change of O and wp
i_rad = i0*math.pi/180.0
coef = -((3*math.sqrt(Mew)*J2*math.pow(Re,2))/(2*pow((1-pow(e0,2)),2)*pow(a0,3.5)))
#dO = coef*math.cos(i_rad)
dO = 0
#dwp = coef*(2.5*pow(math.sin(i_rad),2)-2)
dwp = 0

## Calculate time since perigee passage
w0_rad = w0*math.pi/180.0
E_coef = math.sqrt((1-e0)/(1+e0))*math.tan(w0_rad/2.0)
E0 = 2.0*math.atan(E_coef)
if E0 < 0:
    E0 += 2*math.pi
M0 = E0-e0*math.sin(E0)
T = 2.0*math.pi*pow(a0,1.5)/math.sqrt(Mew)
dt = (M0*T)/2.0/math.pi ## units of seconds
print(str(dt)+" seconds since perigee passage.")

## Define times
t0 = datetime.datetime.strptime(epoch,"%Y-%m-%d %H:%M:%S")
print("Epoch: " + str(t0))
tp = t0 - datetime.timedelta(seconds=dt)
print("Time of perigee passage: " + tp.strftime("%Y-%m-%d %H:%M:%S"))
h = math.sqrt(Mew*a0*(1-pow(e0,2)))
P = pow(h,2)/Mew

## Main loop -> generate matrix of rs and vs
lats = []
longs = []
N = int(24*3600/60)
times = np.linspace(0,24*3600,N);
for time in times:
    ## Define the time
    l_time = tp + datetime.timedelta(seconds=time)

    ## Define M, E, and w (true anomaly)
    M = 2.0*math.pi*time/T
    E = solve_kepler(M,e0)
    w = calc_w(E,e0)

    ## Define change in RAAN and wp
    O = O0 + dO*time
    wp = wp0 + dwp*time

    ## Define coordinates in perifocal and geocentric frames
    r = P/(1+e0*math.cos(w))
    r_p = np.array([r*math.cos(w),r*math.sin(w),0])
    qp_g = peri2geo(O,i0,wp)
    r_geo = np.matmul(qp_g,r_p)

    ## Define coordinates in ecef frame
    bar = [l_time.year,l_time.month,l_time.day,l_time.hour,l_time.minute,l_time.second]
    foo = matlab.double(bar)

    #qg_ecef = eng.dcmeci2ecef('IAU-2000/2006',foo)
    qg_ecef = geo2ecef(l_time)
    qg_ecef = np.asarray(qg_ecef)
    r_ecef = np.matmul(qg_ecef,r_geo)

    ## Calculate instanteous latitude and longitude
    lat,lon = get_latlon(r_ecef)

    ## Append
    lats.append(lat*180/math.pi)
    longs.append(lon*180/math.pi)

## Plot
print(max(lats))
mat_lats = matlab.double(lats)
mat_lons = matlab.double(longs)
eng.plot_track(mat_lats,mat_lons,nargout=0)
