import unittest
import math
import datetime as dt
import initialize_track

class TestLla2Geo(unittest.TestCase):

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

        self.assertLessEqual(int(x),true_x+1) # Check that values align within 1 m
        self.assertGreaterEqual(int(x),true_x-1)
        self.assertLessEqual(int(y),true_y+1)
        self.assertGreaterEqual(int(y),true_y-1)
        self.assertLessEqual(int(z),true_z+1)
        self.assertGreaterEqual(int(z),true_z-1)

if __name__ == '__main__':
    unittest.main()

