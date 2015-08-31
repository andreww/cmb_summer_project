"""Test cases from ttimel fortran code

   Each test case corresponds to a different
   combination of event location, delta and
   azimuth. All for ak135. 'basic' phases 
   are tested.
"""
import numpy.testing as npt

import ellippy

def parse_ttimel_data(data):
    """ttimeel data into phase and correction"""

    phases = []
    corrections = []
    for line in data.splitlines():
        words = line.split()
        phases.append(words[1])
        corrections.append(float(words[3]) - float(words[2]))
        
    return (phases, corrections)

def test_ellippy_nophase():
    npt.assert_raises(ValueError, ellippy.ellip_correct,
                       10.0, 10.0, 10.0, 10.0, "Q")

def test_ellippy_ttimel_lat_0_depth_0_azi_0_delta_45():

    src_lat = 0.0
    src_depth = 0.0
    depth = 0.0
    bazim = 0.0
    delta = 45.0

    phases, corrections = parse_ttimel_data(
        """1  P          497.1021   497.1409     7.9602  -1.57E-01  -4.24E-03
           2  PcP        598.2408   598.3131     3.4483  -1.70E-01   6.53E-03
           3  PP         602.0785   602.2909    10.6359  -1.43E-01  -4.85E-03
           4  PP         605.1327   605.3451     9.1824  -1.51E-01  -2.47E-02
           5  ScP        831.9159   832.2742     4.1221  -2.87E-01   7.86E-03
           6  S          896.6008   896.6761    14.4881  -2.58E-01  -3.01E-03
           7  PKiKP     1016.9644  1017.2147     0.9614  -1.72E-01   1.57E-02
           8  PKKPdf    1891.0547  1895.3062    -0.9828  -1.72E-01  -1.46E-02
           9  SKKPdf    2104.2104  2108.3516    -0.9346  -2.89E-01  -1.51E-02
          10  P'P'df    2396.2705  2394.2739    -1.2510  -1.72E-01  -1.24E-02
          11  P'P'ab    2459.6294  2457.6328    -4.3179  -1.68E-01   3.53E-02""")

    for phase, correction in zip(phases, corrections):
        correction_calc = ellippy.ellip_correct(src_lat, src_depth, 
                                            bazim, delta, phase)
        npt.assert_almost_equal(correction_calc, correction, decimal=3,
                                err_msg="(phase is " + phase + ")")


if __name__ == "__main__":
    test_ellippy_ttimel_lat_0_depth_0_azi_0_delta_45()
