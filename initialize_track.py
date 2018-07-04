## Import modules
import math
import datetime
import numpy as np
import os
import argparse

import matlab.engine


#################
### Functions ###
#################
def parse_rcfile(line):
    coords = line.split(',')[0:2]
    ref_a = float(line.split(',')[-1].strip())

    ref_coords = []
    for coord in coords:
        coord = coord.strip()

        sign = 1
        if coord[0] == "-":
            sign = -1
            coord = coord[1:]

        comps = coord.split('-')
        foo = float(comps[0]) + float(comps[1]) / 60 + float(comps[2]) / 3600
        ref_coords.append(math.radians(foo * sign))

    ref_coords.append(ref_a)
    return ref_coords

def ic_calc_time_since_perigee(e,w,a,epoch):
    E_coef = math.sqrt((1 - e) / (1 + e)) * math.tan(w / 2.0)
    E0 = 2.0 * math.atan(E_coef)
    if E0 < 0:
        E0 += 2 * math.pi
    M0 = E0 - e * math.sin(E0)
    T = 2.0 * math.pi * pow(a, 1.5) / math.sqrt(Mew)
    dt = (M0 * T) / 2.0 / math.pi  ## units of seconds
    print(str(dt) + " seconds since perigee passage.")
    h = math.sqrt(Mew * a * (1 - pow(e, 2)))
    P = pow(h, 2) / Mew

    ## Define times
    t0 = datetime.datetime.strptime(epoch, "%Y-%m-%d %H:%M:%S")
    print("Epoch: " + str(t0))
    tp = t0 - datetime.timedelta(seconds=dt)

    return tp, h, P, T

def parse_tle(line1,line2):
    ## Read tle
    sat_num = line1[2:7]
    epoch_year = line1[18:20]
    epoch_day = line1[20:32]
    i0 = math.radians(float(line2[8:16]))
    O0 = math.radians(float(line2[17:25]))
    e0 = float("." + line2[26:33])
    wp0 = math.radians(float(line2[34:42]))
    M0 = math.radians(float(line2[43:51]))
    n = float(line2[52:63]) * 2 * math.pi / 86400  # rad/s

    ## Define times
    if int(epoch_year) > datetime.datetime.now().year:
        year = "19" + str(epoch_year)
    else:
        year = "20" + str(epoch_year)

    frac, doy = math.modf(float(epoch_day))
    frac, hour = math.modf(frac * 24)
    frac, min = math.modf(frac * 60)
    frac, sec = math.modf(frac * 60)

    if doy < 10:
        doy = "00" + str(doy)
    elif doy < 100:
        doy = "0" + str(doy)

    epoch = str(year) + "-" + str(int(doy)) + " " + str(int(hour)) + ":" + str(int(min)) + ":" + str(
        int(sec)) + "." + str(frac)[2:6]

    return sat_num, epoch, i0, O0, e0, wp0, M0, n

def tle_calc_time_since_perigee(M, e, n, i, epoch):
    ## Calculate values
    E0 = solve_kepler(M, e)
    w0 = calc_w(E0, e)
    a0 = math.pow(Mew / math.pow(n, 2), float(1 / 3))
    T = 2.0 * math.pi * pow(a0, 1.5) / math.sqrt(Mew)
    dt = M / n
    dO, dwp = get_dO_dwp(a0, e, i)
    print(str(dt) + " seconds since perigee passage.")
    h = math.sqrt(Mew * a0 * (1 - math.pow(e, 2)))
    P = math.pow(h, 2) / Mew

    t0 = datetime.datetime.strptime(epoch, "%Y-%j %H:%M:%S.%f")
    print("Epoch: " + str(t0))
    tp = t0 - datetime.timedelta(seconds=dt)

    return tp, h, P, dO, dwp, T

def solve_kepler(M, e):
    ## Purpose: Given mean anomaly and eccentricity, return eccentric anomaly
    ## Inputs:
    ##     M: mean anomaly
    ##     e: eccentricity

    ## Define equations
    f = lambda x: x - e * math.sin(x) - M
    f_prime = lambda x: 1 - e * math.cos(x)

    ## Pick guess
    if M < math.pi:
        E = M + e / 2
    else:
        E = M - e / 2

    ## Loop until we are close to the root
    ratio = f(E) / f_prime(E)
    while abs(ratio) > 1e-8:
        E -= ratio
        ratio = f(E) / f_prime(E)

    if E > 2.0 * math.pi:
        two_pi = 2.0 * math.pi
        rem = E % two_pi
        E = rem

    return E

def calc_w(E, e):
    ## Purpose: Given eccentric anomaly and eccentricity, return true anomaly

    ## Inputs:
    ##   E: eccentric anomaly
    ##   e: eccentricity

    e_coef = math.sqrt((1 - e) / (1 + e))
    w = 2.0 * math.atan(math.tan(E / 2.0) / e_coef)

    if w < 0:
        w += 2.0 * math.pi

    return w

def peri2geo(O, i, wp):
    ## Purpose: Given RAAN, inclination, and argument of perigee, return rotation matrix from perifocal to geocentric frame

    ## Inputs:
    ##   O: RAAN (rad)
    ##   i: inclination (rad)
    ##   wp: argument of perigee (rad)

    q11 = -math.sin(O) * math.cos(i) * math.sin(wp) + math.cos(O) * math.cos(wp)
    q12 = -math.sin(O) * math.cos(i) * math.cos(wp) - math.cos(O) * math.sin(wp)
    q13 = math.sin(O) * math.sin(i)
    q21 = math.cos(O) * math.cos(i) * math.sin(wp) + math.sin(O) * math.cos(wp)
    q22 = math.cos(O) * math.cos(i) * math.cos(wp) - math.sin(O) * math.sin(wp)
    q23 = -math.cos(O) * math.sin(i)
    q31 = math.sin(i) * math.sin(wp)
    q32 = math.sin(i) * math.cos(wp)
    q33 = math.cos(i)

    Q = np.array([[q11, q12, q13], [q21, q22, q23], [q31, q32, q33]])
    return Q

def calc_greenwich_sidereal(time):
    ## Purpose: Calculate the greenwich sidereal time at a given time

    y = time.year
    m = time.month
    d = time.day
    hr = time.hour
    min = time.minute
    sec = time.second

    j0 = 367.0 * y - math.trunc((7.0 * (y + math.trunc((m + 9) / 12.0))) / 4.0) + math.trunc(
        275.0 * m / 9.0) + d + 1721013.5
    T0 = (j0 - 2451545.0) / 36525.0
    thet_g0 = 100.4606184 + 36000.77004 * T0 + 0.000387933 * math.pow(T0, 2) - 2.583e-8 * math.pow(T0, 3)

    thet_g0 = thet_g0 % 360.0
    if thet_g0 < 0:
        thet_g0 += 360

    UT = hr + min / 60.0 + sec / 3600.0
    thet_g = thet_g0 + 360.98564724 * UT / 24
    thet_g = thet_g * math.pi / 180.0

    return thet_g

def geo2ecef(time):
    ## Purpose: Return the rotation matrix from geocentric to ecef at the input times

    thet_g = calc_greenwich_sidereal(time)

    Q = np.array([[math.cos(thet_g), math.sin(thet_g), 0], [-math.sin(thet_g), math.cos(thet_g), 0], [0, 0, 1]])
    return Q

def get_latlon(r):
    ## Purpose: Given radius vector in ECEF, return latitude and longitude in radians

    ## Inputs:
    ##   r: radius vector in ECEF frame (m)

    x = r[0]
    y = r[1]
    z = r[2]

    r_mag = math.sqrt(math.pow(x, 2) + math.pow(y, 2) + math.pow(z, 2))

    l = x / r_mag
    m = y / r_mag
    n = z / r_mag

    dec = math.asin(n)
    if m > 0:
        ra = math.acos(l / math.cos(dec))
    else:
        ra = 2.0 * math.pi - math.acos(l / math.cos(dec))

    return dec, ra

def find_rgeo(P, e, w, i, wp, O):
    ## Purpose: Calculate radius in ECI
    ##
    ## Inputs:
    ##   P: semi-latus rectum
    ##   e: eccentricity
    ##   w: true anomaly (rad)
    ##   i: inclination (rad)
    ##   wp: argument of perigee (rad)

    r = P / (1 + e * math.cos(w))
    r_p = np.array([[r * math.cos(w)], [r * math.sin(w)], [0]])
    qp_g = peri2geo(O, i, wp)
    r_geo = np.matmul(qp_g, r_p)

    return r_geo

def get_dO_dwp(a0, e0, i0):
    ## Purpose: Calculate rate of change of RAAN and wp

    ## Inputs:
    ##   a0: sma (m)
    ##   e0: eccentricity
    ##   i0: inclination (rad)

    coef = -((3 * math.sqrt(Mew) * J2 * math.pow(Re, 2)) / (2 * pow((1 - pow(e0, 2)), 2) * pow(a0, 3.5)))
    dO = coef * math.cos(i0)
    dwp = coef * (2.5 * pow(math.sin(i0), 2) - 2)

    return dO, dwp


def haversine(ref_coords, coords):
    ## Purpose: Haversine formulat to calculate distance between two lat lon pairs

    ## Inputs:
    ##   ref_coords: lat lon pair one
    ##   coords: lat lon pair two

    s1 = math.radians(ref_coords[0])
    s2 = math.radians(coords[0])
    l1 = math.radians(ref_coords[1])
    l2 = math.radians(coords[1])
    l12 = l2 - l1

    y = math.sqrt(math.pow(math.cos(s1) * math.sin(s2) - math.sin(s1) * math.cos(s2) * math.cos(l12), 2) + math.pow(
        math.cos(s2) * math.sin(l12), 2))
    x = math.sin(s1) * math.sin(s2) + math.cos(s1) * math.cos(s2) * math.cos(l12)

    s12 = math.atan2(y, x)
    d = Re * s12

    return d

def lla2geo(ref_coords, time):
    ## Purpose: Given coordinates in lla, convert to geocentric/eci coordinate frame
    ##
    ## Inputs:
    ##   ref_coords: lat, lon, alt coordinates
    ##   time: current time

    lat = ref_coords[0]
    lon = ref_coords[1]
    h = ref_coords[2]
    local_theta = calc_greenwich_sidereal(time) + lon

    C = 1/math.sqrt(1 + WGS84_f*(WGS84_f-2)*math.pow(math.sin(lat),2))
    S = math.pow(1-WGS84_f,2)*C

    x_geo = WGS84_a*C*math.cos(lat)*math.cos(local_theta)
    y_geo = WGS84_a*C*math.cos(lat)*math.sin(local_theta)
    z_geo = WGS84_a*S*math.sin(lat)

    r_geo = np.array([[x_geo], [y_geo], [z_geo]])

    return r_geo

def geo2topo(theta, lat):
    ## Purpose: Calculate transformation matrix from geocentric (ECI) to topocentric horizon

    ## Inputs:
    ##   theta: local sidereal time
    ##   lat: latitude

    Q1 = [math.sin(lat) * math.cos(theta), math.sin(lat) * math.sin(theta), -math.cos(lat)]
    Q2 = [-math.sin(theta), math.cos(theta), 0]
    Q3 = [math.cos(lat) * math.cos(theta), math.cos(lat) * math.sin(theta), math.sin(lat)]

    Q = np.array([Q1, Q2, Q3])
    return Q

def calc_sat_subpoint(lat,lon,r_geo):
    ## Purpose: Given geodetic satellite latitude and longitude, calculate geocentric latitude, longitude, and altitude
    ##
    ## Inputs:
    ##   lat: geodetic latitude
    ##   lon: geodetic longitude
    ##   r_geo: ECI coordinates

    e2 = 2*WGS84_f-math.pow(WGS84_f,2)
    lat_i = lat
    tol = 1e-8
    R = math.sqrt(math.pow(r_geo[0],2) + math.pow(r_geo[1],2))

    f_C = lambda x: 1/math.sqrt(1-e2*math.pow(math.sin(x),2))
    f_lat = lambda x, C: math.atan2(r_geo[2] + WGS84_a*C*e2*math.sin(x),R)
    C = f_C(lat_i)
    lat_n = f_lat(lat_i,C)

    while abs(lat_i - lat_n) > tol:
        lat_i = lat_n
        C = f_C(lat_i)
        lat_n = f_lat(lat_i,C)

    geoc_lat = lat_n
    geoc_lon = lon
    h = R/math.cos(geoc_lat) - WGS84_a*C

    r_lla = [geoc_lat,geoc_lon,h]
    return r_lla

#####################################################
### 01. Initializing variables and getting inputs ###
#####################################################

## Define constants
J2 = 1.08263e-3
Mew = 3986004.418e8 #m3/s2
WGS84_a = 6378137.0 #m
WGS84_f = 1/298.257223563
WGS84_b = WGS84_a*(1-WGS84_f)
WGS84_w = 7292115e-11 #rad/s
Re = 6378.37e3 # m
eng = matlab.engine.start_matlab()
eng.addpath(os.getcwd() + "/")

if __name__ == "__main__":

    ## Define argument parsing
    parser = argparse.ArgumentParser(description="Process initial conditions")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--kep", '-k', dest="tle", action="store_false",
                   default="true",
                   help="Initial conditions are in the ic.kep format (default: tle):\n  epoch (UTC)\n  semi-major axis (m)\n  eccentricity\n  inclination (deg)\n  true anomaly (deg)\n  RAAN (deg)\n  argument of perigee (deg)\n")
    group.add_argument("--tle", '-t', dest="tle", action="store_true",
                   default="true", help="Initial conditions are in the tle.txt format (default: tle)\n")
    parser.add_argument('--init_cond', '-i', metavar="ic_file",
                    help="file holding initial conditions", required=True)
    parser.add_argument("--end_time", '-e', dest="end_time", default="T",
                    help="Number of minutes to simulate the ground track. Default: one revolution")
    parser.add_argument("--ref_coord", '-r',
                    help="file holding reference coordinates. Format: LAT (N), LON (E), ALT (m). For example, Denver would be: 39-45-43, -104-52-52 (DMS)",
                    default="None")
    args = vars(parser.parse_args())

    ## Define input args
    ic_file = args["init_cond"]
    ic_tle = args["tle"]
    end_time = args["end_time"]
    rc_file = args["ref_coord"]

    # print(ic_file)
    # print(ic_tle)
    # print(end_time)
    print(rc_file)

    if rc_file is not None:
        with open(rc_file, "r") as f:
            line = f.readline()
        ref_coords = parse_rcfile(line)

    if not ic_tle:
        print("Now parsing input Keplerian orbit elements")
        ## Read text file with orbit initial conditions (ic.kep)
        with open(ic_file, "r") as f:
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
            i0 = math.radians(float(ic[3]))
            w0 = math.radians(float(ic[4]))
            O0 = math.radians(float(ic[5]))
            wp0 = math.radians(float(ic[6]))

        #######################################
        ### 02. Calculate starting position ###
        #######################################

        ## Determine rate of change of O and wp
        dO, dwp = get_dO_dwp(a0, e0, i0)

        ## Calculate time since perigee passage
        tp, h, P, T = ic_calc_time_since_perigee(e0,w0,a0,epoch)
        print("Time of perigee passage: " + tp.strftime("%Y-%m-%d %H:%M:%S"))

    else:
        print("Now parsing input tle")
        with open(ic_file, "r") as f:
            line1 = f.readline()
            line2 = f.readline()

        sat_num, epoch, i0, O0, e0, wp0, M0, n = parse_tle(line1,line2)
        tp, h, P, dO, dwp, T = tle_calc_time_since_perigee(M0, e0, n, i0, epoch)
        print("Time of perigee passage: " + tp.strftime("%Y-%m-%d %H:%M:%S"))

    #################################################################
    ### 03. Main loop: step through times and get a lat/lon point ###
    #################################################################

    ## Main loop -> generate matrix of lats and longs
    lats = []
    longs = []

    if end_time == "T":
        end_time = T
        print("Generating ground track for one revolution")
    else:
        print("Generating ground track for " + end_time + " minutes.")
        end_time = float(end_time) * 60

    loop_time = 86400 * 1
    N = int(loop_time / 60)
    times = np.linspace(0, loop_time, N)
    vis_t = []
    vis_a = []
    vis_lats = []
    vis_lons = []
    vis_az = []
    elevs = []
    prev_elev = 999
    lim = math.radians(10)
    for time in times:
        ## Define the time
        l_time = tp + datetime.timedelta(seconds=time)

        ## Define M, E, and w (true anomaly)
        M = 2.0 * math.pi * time / T
        E = solve_kepler(M, e0)
        w = calc_w(E, e0)

        ## Define change in RAAN and wp
        O = O0 + dO * time
        wp = wp0 + dwp * time

        ## Define coordinates in geocentric frames
        r_geo = find_rgeo(P, e0, w, i0, wp, O)

        ## Define coordinates in ecef frame
        # bar = [l_time.year,l_time.month,l_time.day,l_time.hour,l_time.minute,l_time.second]
        # foo = matlab.double(bar)

        # qg_ecef = eng.dcmeci2ecef('IAU-2000/2006',foo)

        qg_ecef = geo2ecef(l_time)
        qg_ecef = np.asarray(qg_ecef)
        r_ecef = np.matmul(qg_ecef, r_geo)

        ## Calculate instanteous latitude and longitude
        lat, lon = get_latlon(r_ecef)

        ## Append
        if time < end_time:
            lats.append(lat * 180 / math.pi)
            longs.append(lon * 180 / math.pi)

        ## See if satellite is visible from ref coords (if given)
        if rc_file != None:
            o_geo = lla2geo(ref_coords, l_time)
            theta_g = calc_greenwich_sidereal(l_time)
            theta_l = theta_g + ref_coords[1]

            rsat_o = np.subtract(r_geo, o_geo)
            Qgeo_topo = geo2topo(theta_l, lat)
            r_topo = np.matmul(Qgeo_topo, rsat_o)

            rs = float(r_topo[0])
            re = float(r_topo[1])
            rz = float(r_topo[2])

            range = math.sqrt(math.pow(rs, 2) + math.pow(re, 2) + math.pow(rz, 2))
            elev = math.asin(rz / range)
            az = math.atan2(rs, re) + math.radians(90)

            if elev > lim:
                vis_t.append(l_time)
                vis_a.append(elev)
                vis_lats.append(math.degrees(lat))
                vis_lons.append(math.degrees(lon))
                vis_az.append(math.degrees(az))

            if elev < lim and prev_elev > lim and prev_elev != 999:
                caption = "Visible at reference location from " + vis_t[0].strftime("%Y-%m-%d %H:%M:%S") + " to " + vis_t[
                    -1].strftime("%Y-%m-%d %H:%M:%S") + " with maximum elevation of " + str(
                    round(math.degrees(max(vis_a)), 0)) + " and azimuth " + str(vis_az[0]) + " to " + str(vis_az[-1])
                print(caption)
                line1 = caption[:78]
                line2 = caption[78:]
                name = 'Plots/' + vis_t[0].strftime('%Y-%m-%dT%H:%M:%S') + '_' + vis_t[-1].strftime(
                    '%Y-%m-%dT%H:%M:%S') + '.jpg'
                eng.plot_track(matlab.double(vis_lats), matlab.double(vis_lons),
                               matlab.double([math.degrees(ref_coords[0]), math.degrees(ref_coords[1])]), name, line1,
                               line2, nargout=0)
                vis_t = []
                vis_a = []
                vis_lats = []
                vis_lons = []
                vis_az = []

            prev_elev = elev

    ## Plot
    mat_lats = matlab.double(lats)
    mat_lons = matlab.double(longs)
    eng.plot_track(mat_lats, mat_lons, matlab.double([]), 'ground_track.jpg', '', '', nargout=0)
