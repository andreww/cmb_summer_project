#!/usr/bin/env python

import geographical as gt
import unittest
import numpy.testing as npt # This has support for array-like test comparisons
import random

class TestGeographicalFunctions(unittest.TestCase):

    def setUp(self):
        # These are all "known points" as (x, y, z) or (r, lat, lon)
        self.north_pole_geog = (6378.0, 90.0, 0.0)
        self.north_pole_cart = (0.0, 0.0, 6378.0)
        self.south_pole_geog = (6378.0, -90.0, 0.0)
        self.south_pole_cart = (0.0, 0.0, -6378.0)
        self.x0_geog = (6378.0, 0.0, 0.0)
        self.x0_cart = (6378.0, 0.0, 0.0)
        self.x1_geog = (6378.0, 0.0, 180.0)
        self.x1_cart = (-6378.0, 0.0, 0.0)
        self.y0_geog = (6378.0, 0.0, 90.0)
        self.y0_cart = (0.0, 6378.0, 0.0)
        self.y1_geog = (6378.0, 0.0, -90.0)
        self.y1_cart = (0.0, -6378.0, 0.0)

    def test_cart2geog_elements(self):
        self.assertEqual(gt.cart2geog(self.north_pole_cart[0],self.north_pole_cart[1],
                           self.north_pole_cart[2]), self.north_pole_geog)
        self.assertEqual(gt.cart2geog(self.south_pole_cart[0],self.south_pole_cart[1],
                           self.south_pole_cart[2]), self.south_pole_geog)
        self.assertEqual(gt.cart2geog(self.x0_cart[0],self.x0_cart[1],self.x0_cart[2]), self.x0_geog)
        self.assertEqual(gt.cart2geog(self.x1_cart[0],self.x1_cart[1],self.x1_cart[2]), self.x1_geog)
        self.assertEqual(gt.cart2geog(self.y0_cart[0],self.y0_cart[1],self.y0_cart[2]), self.y0_geog)
        self.assertEqual(gt.cart2geog(self.y1_cart[0],self.y1_cart[1],self.y1_cart[2]), self.y1_geog)

    def test_geog2cart_elements(self):
        npt.assert_almost_equal(gt.geog2cart(self.north_pole_geog[0],self.north_pole_geog[1],self.north_pole_geog[2]), self.north_pole_cart)
        npt.assert_almost_equal(gt.geog2cart(self.south_pole_geog[0],self.south_pole_geog[1],self.south_pole_geog[2]), self.south_pole_cart)
        npt.assert_almost_equal(gt.geog2cart(self.x0_geog[0],self.x0_geog[1],self.x0_geog[2]), self.x0_cart)
        npt.assert_almost_equal(gt.geog2cart(self.x1_geog[0],self.x1_geog[1],self.x1_geog[2]), self.x1_cart)
        npt.assert_almost_equal(gt.geog2cart(self.y0_geog[0],self.y0_geog[1],self.y0_geog[2]), self.y0_cart)
        npt.assert_almost_equal(gt.geog2cart(self.y1_geog[0],self.y1_geog[1],self.y1_geog[2]), self.y1_cart)

    def test_geog2cart2geog_fuse(self):
        random.seed(987.987)
        for i in range(1000):
            r = random.uniform(0.0, 1000) 
            lat = random.uniform(-90,90)
            lon = random.uniform(-180, 180)
            (x, y, z) = gt.geog2cart(r, lat, lon)
            npt.assert_almost_equal(gt.cart2geog(x, y, z), (r, lat, lon))
        
    def test_cart2geog2cart_fuse(self):
        random.seed(987.987)
        for i in range(1000):
            x = random.uniform(-1000, 1000) 
            y = random.uniform(-1000,1000)
            z = random.uniform(-1000, 1000)
            (r, lat, lon) = gt.cart2geog(x, y, z)
            npt.assert_almost_equal(gt.geog2cart(r, lat, lon), (x, y, z))

class TestVincentyFunctions(unittest.TestCase):
   
    # Table II of Vincenty's paper (T. Vincenty 1975, "Direct
    # and inverse solutions of geodesics on the ellipsoid with
    # application of nested equations" Survey Review XXII 
    # pp.88-93) has five test examples for the forward and inverse
    # problem (with results rounded to 0.00001 seconds of arc and 
    # 1 mm). The inverse versions of these are implemented here. Note
    # the non-standard (old) ellipsoid usage.

    #def test_vincenty_inverse_tabulated_a(self):
    #    correct_result = (14110.526170, gt.dms2dec(96.0,36.0,8.79960), 
    #                      gt.dms2dec(137.0,52.0,22.01454))
    #    lat1 = gt.dms2dec(55.0,45.0,0.0)
    #    lat2 = gt.dms2dec(-33.0,26.0,0.0)
    #    lon2 = gt.dms2dec(108.0,13.0,0.0)
    #    npt.assert_almost_equal(gt.vincenty(lat1,0.0,lat2,lon2,
    #                  r_major=6377.397155,r_minor=6356.078963),correct_result)

    def test_vincenty_inverse_tabulated_b(self):
        correct_result = (4085.966703, gt.dms2dec(95.0,27.0,59.63089), 
                          gt.dms2dec(118,5.0,58.96161))
        lat1 = gt.dms2dec(37.0,19.0,54.95367)
        lat2 = gt.dms2dec(26.0,7.0,42.83946)
        lon2 = gt.dms2dec(41.0,28.0,35.50729)
        npt.assert_almost_equal(gt.vincenty(lat1,0.0,lat2,lon2,
                  r_major=6378.388000,r_minor=6356.911946),correct_result,6)


    def test_vincenty_direct_tabulated_b(self):
        lat1 = gt.dms2dec(37.0,19.0,54.95367)
        lon1 = 0.0
        lat2 = gt.dms2dec(26.0,7.0,42.83946)
        lon2 = gt.dms2dec(41.0,28.0,35.50729)
        dist = 4085.966703
        azi = gt.dms2dec(95.0,27.0,59.63089)
        correct_result = (lat2, lon2)
        npt.assert_almost_equal(gt.vincenty_direct(lat1, lon1, azi, dist,
             r_major=6378.388000,r_minor=6356.911946),correct_result,6)

    def test_vincenty_inverse_tabulated_c(self):
        correct_result = (8084.823839, gt.dms2dec(15.0,44.0,23.74850), 
                          gt.dms2dec(144.0,55.0,39.92147))
        lat1 = gt.dms2dec(35.0,16.0,11.24862)
        lat2 = gt.dms2dec(67.0,22.0,14.77638)
        lon2 = gt.dms2dec(137.0,47.0,28.31435)
        npt.assert_almost_equal(gt.vincenty(lat1,0.0,lat2,lon2,
                  r_major=6378.388000,r_minor=6356.911946),correct_result,6)

    def test_vincenty_direct_tabulated_c(self):
        lat1 = gt.dms2dec(35.0,16.0,11.24862)
        lon1 = 0.0
        lat2 = gt.dms2dec(67.0,22.0,14.77638)
        lon2 = gt.dms2dec(137.0,47.0,28.31435)
        dist = 8084.823839
        azi = gt.dms2dec(15.0,44.0,23.74850)
        correct_result = (lat2, lon2)
        npt.assert_almost_equal(gt.vincenty_direct(lat1, lon1, azi, dist,
             r_major=6378.388000,r_minor=6356.911946),correct_result,6)

    # FIXME
    # Skip tabulated data d - very sensive to lattitude...

    def test_vincenty_inverse_tabulated_e(self):
        correct_result = (19780.006558, gt.dms2dec(4.0,59.0,59.99995), 
                          gt.dms2dec(174.0,59.0,59.88481))
        lat1 = gt.dms2dec(1.0,0.0,0.0)
        lat2 = gt.dms2dec(1.0,1.0,15.18952)
        lon2 = gt.dms2dec(179.0,46.0,17.84244)
        npt.assert_almost_equal(gt.vincenty(lat1,0.0,lat2,lon2,
                  r_major=6378.388000,r_minor=6356.911946),correct_result,6)

    def test_vincenty_direct_tabulated_c(self):
        lat1 = gt.dms2dec(1.0,0.0,0.0)
        lon1 = 0.0
        lat2 = gt.dms2dec(1.0,1.0,15.18952)
        lon2 = gt.dms2dec(179.0,46.0,17.84244)
        dist = 19780.006558
        azi = gt.dms2dec(4.0,59.0,59.99995)
        correct_result = (lat2, lon2)
        npt.assert_almost_equal(gt.vincenty_direct(lat1, lon1, azi, dist,
             r_major=6378.388000,r_minor=6356.911946),correct_result,6)

if __name__ == '__main__':
    unittest.main()

