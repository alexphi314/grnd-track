## Read text file with orbit initial conditions (ic.kep)
with open("ic.kep","r") as f:
    ## Format is:
    ##   epoch
    ##   semi-major axis
    ##   eccentricity
    ##   inclination
    ##   true anomaly
    ##   RAAN
    ##   argument of perigee
    ic = f.read()
    ic = ic.split("\n")
    epoch = ic[0]
    a = float(ic[1])
    e = float(ic[2])
    i = float(ic[3])
    w = float(ic[4])
    O = float(ic[5])
    wp = float(ic[6])
