import unittest
import math
import datetime as dt
import initialize_track

class Test_lla2geo(unittest.TestCase):

    def test_convert(self):
        lat = math.radians(40)
        lon = math.radians(-75)
        h = 0
        time = dt.datetime(year=1995,month=10,day=1,hour=9)

        r_geo = initialize_track.lla2geo([lat,lon,h],time)
        x = r_geo[0][0]
        y = r_geo[1][0]
        z = r_geo[2][0]

        true_x = int(1703.295e3)
        true_y = int(4586.650e3)
        true_z = int(4077.984e3)
        tol = 1.6

        self.assertLessEqual(x,true_x+tol) # Check that values align within 2 m
        self.assertGreaterEqual(x,true_x-tol)
        self.assertLessEqual(y,true_y+tol)
        self.assertGreaterEqual(y,true_y-tol)
        self.assertLessEqual(z,true_z+tol)
        self.assertGreaterEqual(z,true_z-tol)

class Test_calc_greenwich_sidereal_time(unittest.TestCase):

    def test_time(self):
        time = dt.datetime(year=2004,month=3,day=3,hour=4,minute=30)
        lon = math.radians(139.8)

        theta = initialize_track.calc_greenwich_sidereal(time)
        self.assertAlmostEqual(theta,math.radians(228.79354),2)

class Test_get_lat_lon(unittest.TestCase):
#
#     def test_conversion(self):
#         r_eci = [-5778543.80350671,-3450942.43214222,-822041.74691368]
#         r_ecef = [6584735.38392856,1393496.35185912,-822041.74691368]
#         l_time = dt.datetime(year=2018,month=6,day=15,hour=19,minute=39,second=26,microsecond=836033)
#
#         lat, lon = initialize_track.get_latlon(r_ecef)
#         lat2 = math.atan2(r_eci[2],math.sqrt(math.pow(r_eci[0],2)+math.pow(r_eci[1],2)))
#         theta_g = initialize_track.calc_greenwich_sidereal(l_time)
#         lon2 = math.atan2(r_eci[1],r_eci[0]) - theta_g
#
#         self.assertAlmostEqual(lat,lat2,2)
#         self.assertAlmostEqual(lon,lon2,2)

    def test_r_geo(self):
        w = math.radians(30)
        O = math.radians(40)
        i = math.radians(30)
        wp = math.radians(60)
        h = 80e3*1e6
        e = 1.4

        P = math.pow(h, 2) / initialize_track.Mew

        output_rgeo = initialize_track.find_rgeo(P, e, w, i, wp, O)
        r_geo = [-4040e3,4815e3,3629e3]

        self.assertAlmostEqual(output_rgeo[0][0]/1e3,r_geo[0]/1e3,2)
        self.assertAlmostEqual(output_rgeo[1][0]/1e3,r_geo[1]/1e3,2)
        self.assertAlmostEqual(output_rgeo[2][0]/1e3,r_geo[2]/1e3,2)

class Test_calc_sat_subpoint(unittest.TestCase):

    def test_conversion(self):
        lat = math.radians(45)
        lon = math.radians(-93)
        time = dt.datetime(year=1995,month=11,day=18,hour=12,minute=46)


if __name__ == '__main__':
    unittest.main()

