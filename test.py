import unittest
import math
import datetime as dt
import numpy as np

import initialize_track

def assert_within_bounds(testCase,test_num,num,tol):
    ## Purpose: Assert that the test number is within a tolerance of the true number
    ##
    ## Inputs:
    ##  test_num: test number from script
    ##  num: true number
    ##  tol: tolerance test number must be within

    testCase.assertGreaterEqual(test_num,num-tol)
    testCase.assertLessEqual(test_num,num+tol)

def assert_matrix_almost_equal(testCase,true,test):
    ## Purpose: Assert that all elements of input matrices are equal
    ##
    ## Inputs:
    ##   True array
    ##   Test array

    r,c = true.shape
    r_test,c_test = test.shape

    testCase.assertEqual(r,r_test)
    testCase.assertEqual(c,c_test)

    for i in range(0,r):
        for j in range(0,c):
            testCase.assertAlmostEqual(test[i][j],true[i][j],2)

class Test_parse_rcfile(unittest.TestCase):
    def test_function(self):
        line = '53-46-5, -269-52-32'
        coords = np.array([math.radians(53.7181),math.radians(-269.88)])
        coords_out = np.array(initialize_track.parse_rcfile(line))

        self.assertAlmostEqual(coords[0],coords_out[0],2)
        self.assertAlmostEqual(coords[1],coords_out[1],2)

class Test_calc_time_since_perigee(unittest.TestCase):
    e = 0.37255
    w = math.radians(120)
    a = 15300e3
    epoch = '2018-06-07 23:54:00'
    i = 0
    tp = 4077
    n = math.sqrt(initialize_track.Mew/math.pow(a,3))
    E = 2 * math.atan2(math.sqrt((1 - e) / (1 + e)) * math.tan(w / 2), 1)
    M = E - e*math.sin(E)

    def test_ic(self):
        tp_out, h, P, T = initialize_track.ic_calc_time_since_perigee(Test_calc_time_since_perigee.e,
                                                                      Test_calc_time_since_perigee.w,
                                                                      Test_calc_time_since_perigee.a,
                                                                      Test_calc_time_since_perigee.epoch)
        tp_out = (dt.datetime.strptime(Test_calc_time_since_perigee.epoch,'%Y-%m-%d %H:%M:%S')
                                       - tp_out).total_seconds()
        self.assertAlmostEqual(Test_calc_time_since_perigee.tp,tp_out,0)

    def test_tle(self):
        epoch = '2018-180 23:54:00.56'
        tp_out, h, P, dO, dwp, T = initialize_track.tle_calc_time_since_perigee(Test_calc_time_since_perigee.M,
                                                                                Test_calc_time_since_perigee.e,
                                                                                Test_calc_time_since_perigee.n,
                                                                                Test_calc_time_since_perigee.i,
                                                                                epoch)
        tp_out = (dt.datetime.strptime(epoch,'%Y-%j %H:%M:%S.%f')
                                       - tp_out).total_seconds()
        self.assertAlmostEqual(Test_calc_time_since_perigee.tp, tp_out, 0)

class Test_parse_tle(unittest.TestCase):
    def test_function(self):
        line1 = '1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927'
        line2 = '2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537'

        sat_num_out, epoch_out, i_out, O_out, e_out, wp_out, M_out, n_out = initialize_track.parse_tle(line1,line2)
        sat_num = '25544'
        epoch = '2008-264 12:25:40.1041'
        i = math.radians(51.6416)
        O = math.radians(247.4627)
        e = 0.0006703
        wp = math.radians(130.5360)
        M = math.radians(325.0288)
        n = 15.72125391*2*math.pi/86400

        self.assertEqual(sat_num,sat_num_out)
        self.assertEqual(epoch,epoch_out)
        self.assertEqual(i,i_out)
        self.assertEqual(O,O_out)
        self.assertEqual(e,e_out)
        self.assertEqual(wp,wp_out)
        self.assertEqual(M,M_out)
        self.assertEqual(n,n_out)

class Test_solve_kepler(unittest.TestCase):
    def test_function(self):
        E = math.radians(85)
        e = 0.8
        M = 0.68657

        E_out = initialize_track.solve_kepler(M,e)
        self.assertAlmostEqual(E,E_out,2)

class Test_calc_w(unittest.TestCase):
    def test_function(self):
        w = math.radians(85)
        e = 0.8
        E = 2*math.atan2(math.sqrt((1-e)/(1+e))*math.tan(w/2),1)

        w_out = initialize_track.calc_w(E,e)
        self.assertAlmostEqual(w,w_out,2)

class Test_peri2geo(unittest.TestCase):
    def test_function(self):
        i = math.radians(30)
        O = math.radians(40)
        wp = math.radians(60)

        Q = np.array([[-0.0991,0.89593,0.43301],[-0.94175,-0.22496,0.25],[0.32139,-0.38302,0.86603]])
        Q = np.transpose(Q)
        Q_out = initialize_track.peri2geo(O,i,wp)

        assert_matrix_almost_equal(self,Q,Q_out)

class Test_calc_greenwich_sidereal_time(unittest.TestCase):
    def test_function(self):
        time = dt.datetime(year=2004,month=3,day=3,hour=4,minute=30)
        lon = math.radians(139.8)

        theta = initialize_track.calc_greenwich_sidereal(time)
        self.assertAlmostEqual(theta,math.radians(228.79354),2)

class Test_geo2ecef(unittest.TestCase):
    def test_function(self):
        time = dt.datetime(year=2018, month=6, day=15, hour=19, minute=39, second=26)
        Q = np.array([[-0.9474,-0.3200,0.0017],[0.3200,-0.9474,-0.0006],[0,0,1.0000]])
        Q_out = initialize_track.geo2ecef(time)

        assert_matrix_almost_equal(self,Q,Q_out)

class Test_get_lat_lon(unittest.TestCase):
    w = math.radians(30)
    O = math.radians(40)
    i = math.radians(30)
    wp = math.radians(60)
    h = 80e3*1e6
    e = 1.4

    P = math.pow(h, 2) / initialize_track.Mew

    output_rgeo = initialize_track.find_rgeo(P, e, w, i, wp, O)
    r_geo = [-4040e3,4815e3,3629e3]

    def test_comparison(self):
        r_eci = Test_get_lat_lon.r_geo
        l_time = dt.datetime(year=2018, month=6, day=15, hour=19, minute=39, second=26)
        qg_ecef = initialize_track.geo2ecef(l_time)
        qg_ecef = np.asarray(qg_ecef)
        r_ecef = np.matmul(qg_ecef, r_eci)

        lat, lon = initialize_track.get_latlon(r_ecef)
        lat2 = math.atan2(r_eci[2],math.sqrt(math.pow(r_eci[0],2)+math.pow(r_eci[1],2)))
        theta_g = initialize_track.calc_greenwich_sidereal(l_time)
        lon2 = math.atan2(r_eci[1],r_eci[0]) - theta_g

        while lon2 < math.radians(0):
            lon2 += 2*math.pi

        while lon2 > math.radians(360):
            lon2 -= 2*math.pi

        while lat2 < math.radians(-90):
            lat2 += 2*math.pi

        while lat2 > math.radians(90):
            lat2 -= 2*math.pi

        self.assertAlmostEqual(lat,lat2,2)
        self.assertAlmostEqual(lon,lon2,2)

    def test_function_find_r_geo(self):
        r_geo = Test_get_lat_lon.r_geo
        output_rgeo = Test_get_lat_lon.output_rgeo
        tol = 1

        assert_within_bounds(self,output_rgeo[0][0]/1e3,r_geo[0]/1e3,tol)
        assert_within_bounds(self,output_rgeo[1][0]/1e3,r_geo[1]/1e3,tol)
        assert_within_bounds(self,output_rgeo[2][0]/1e3,r_geo[2]/1e3,tol)

    def test_function_get_latlon(self):
        r_eci = [-5368e3,-1784e3,3691e3]
        lat, lon = initialize_track.get_latlon(r_eci)

        self.assertAlmostEqual(lat,math.radians(33.12),2)
        self.assertAlmostEqual(lon,math.radians(198.4),2)

class Test_get_dO_dwp(unittest.TestCase):
    def test_function(self):
        e = 0.008931
        a = 6718e3
        i = math.radians(51.43)

        dO = -1.0465e-6
        dwp = 7.9193e-7

        dO_out, dwp_out = initialize_track.get_dO_dwp(a,e,i)

        self.assertAlmostEqual(dO,dO_out,2)
        self.assertAlmostEqual(dwp,dwp_out,2)

# class Test_haversine(unittest.TestCase):
#     def test_function(self):
#         ref_coords = [50.76,-83.4]
#         coords = [-60.2,89]
#
#         dist = 18860e3
#         dist_out = initialize_track.haversine(ref_coords,coords)
#
#         self.assertAlmostEqual(dist,dist_out,0)

class Test_lla2geo(unittest.TestCase):
    def test_function(self):
        lat = math.radians(40)
        lon = math.radians(-75)
        time = dt.datetime(year=1995,month=10,day=1,hour=9)

        r_geo = initialize_track.lla2geo([lat,lon],time)
        x = r_geo[0][0]
        y = r_geo[1][0]
        z = r_geo[2][0]

        true_x = int(1703.295e3)
        true_y = int(4586.650e3)
        true_z = int(4077.984e3)
        tol = 1.6

        # self.assertLessEqual(x,true_x+tol) # Check that values align within 2 m
        # self.assertGreaterEqual(x,true_x-tol)
        # self.assertLessEqual(y,true_y+tol)
        # self.assertGreaterEqual(y,true_y-tol)
        # self.assertLessEqual(z,true_z+tol)
        # self.assertGreaterEqual(z,true_z-tol)
        assert_within_bounds(self,x,true_x,tol)
        assert_within_bounds(self,y,true_y,tol)
        assert_within_bounds(self,z,true_z,tol)

class Test_geo2topo(unittest.TestCase):
    def test_function(self):
        theta = math.radians(60)
        lat = math.radians(30)

        Q = np.array([[0.25,0.43301,-0.866],[-0.866,0.5,0],[0.43301,0.75,0.5]])
        Q_out = initialize_track.geo2topo(theta,lat)

        assert_matrix_almost_equal(self,Q,Q_out)

class Test_calc_sat_subpoint(unittest.TestCase):
    def test_function(self):
        lat = math.radians(45)
        lon = math.radians(-93)
        r_geo = [-4400.594e3, 1932.87e3, 4760.712e3]

        r_lla_out = initialize_track.calc_sat_subpoint(lat, lon, r_geo)
        r_lla = [math.radians(44.91),math.radians(-92.31),397.507e3]
        tol = 2

        self.assertAlmostEqual(r_lla[0],r_lla_out[0],1)
        self.assertAlmostEqual(r_lla[1],r_lla_out[1],1)
        self.assertGreaterEqual(r_lla_out[2],r_lla[2]-tol)
        self.assertLessEqual(r_lla_out[2],r_lla[2]+tol)

class Test_get_look_angles(unittest.TestCase):
    def test_function(self):
        ref_coords = [math.radians(45),math.radians(-93)]
        time = dt.datetime(year=1995, month=11, day=18, hour=12, minute=46)
        sat_geo = np.array([[-4400.594e3], [1932.87e3], [4760.712e3]])

        elev_out, az_out = initialize_track.get_look_angles(ref_coords, time, sat_geo)
        elev = 81.52
        az = 100.36

        self.assertAlmostEqual(elev,math.degrees(elev_out),2)
        self.assertAlmostEqual(az,math.degrees(az_out),2)

if __name__ == '__main__':
    unittest.main()

