## Import modules
import math; import datetime; import numpy as np;

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

    return E

def calc_w(E,e0):
    ## Purpose: Given eccentric anomaly and eccentricity, return true anomaly

    e_coef = math.sqrt((1-e)/(1+e))

## Define constants
J2 = 1.08263e-3;
G = 6.67e-11;
Me = 5.972e24;
Mew = G*Me;
Re = 6378.37e3; #m

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
dO = coef*math.cos(i_rad)
dwp = coef*(2.5*pow(math.sin(i_rad),2)-2)

## Calculate time since perigee passage
w0_rad = w0*math.pi/180.0
E_coef = math.sqrt((1-e0)/(1+e0))*math.tan(w0_rad/2.0)
E0 = 2.0*math.atan(E_coef)
if E0 < 0:
    E0 += 2*math.pi
M0 = E0-e0*math.sin(E0)
T = 2*math.pi*pow(a0,1.5)/math.sqrt(Mew)
dt = (M0*T)/2/math.pi ## units of seconds
print(str(dt)+" seconds since perigee passage.")

## Define times
t0 = datetime.datetime.strptime(epoch,"%Y-%m-%d %H:%M:%S")
print("Epoch: " + str(t0))
tp = t0 - datetime.timedelta(seconds=dt)
print("Time of perigee passage: " + tp.strftime("%Y-%m-%d %H:%M:%S"))

## Main loop -> generate matrix of rs and vs
rs = []
vs = []
times = np.linspace(0,86400,60);
for time in times:
    M = 2*math.pi*time/T
    E = solve_kepler(M,e0)
    w = calc_w(E,e0)
