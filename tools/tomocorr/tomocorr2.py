#!/usr/bin/env python

from future.utils import native_str

import numpy as np
import obspy.taup as taup
import geographiclib.geodesic as geod


# Geographical raypaths...
# ========================

# We need to be able to find ray paths with the 
# latitude, longitude and depth of each point
# (not just angular distance and depth) to do
# this we can create a subclass of TauPyModel
# which makes use of geographiclib to calculate
# the location of each point. TauPyModelGeo 
# works exactly the same as TauPyModel but has
# an additional optional argument when created
# (to specify the ellipsoid) and a new method
# (get_ray_paths_geo) which takes the latitude and
# logitude for each raypath, which is calculated from
# the location of the source and recever, not the 
# distance.


"""
Holds the ray parameter, time and distance increments, and optionally a
depth, latitude and longitude for a ray passing through some layer.
"""
TimeDistLoc = np.dtype([
    (native_str('p'), np.float_),
    (native_str('time'), np.float_),
    (native_str('dist'), np.float_),
    (native_str('depth'), np.float_),
    (native_str('lat'), np.float_),
    (native_str('lon'), np.float_)
])

# Subclass TauPyModel - if merged this could be 
# a real subclass, or could be added to TauPyModel
# directly (with some more functionality in
# obspy's geodesics stuff.

class TauPyModelGeo(taup.TauPyModel):
    """TauPyModel with geographical coordinates of a path

       Make use of the geographiclib module to
       solve the inverse and direct problem for each point
       along a path.
    """

    def __init__(self, ellipsoid=geod.Geodesic.WGS84, **kwargs):

        taup.TauPyModel.__init__(self, **kwargs)
        self.ellipsoid = ellipsoid

    def get_ray_paths_geo(self, source_depth_in_km, source_latitude_in_deg,
                          source_longitude_in_deg, station_latitude_in_deg,
                          station_longitude_in_deg, phase_list=("ttall",)):
        """
        Return ray paths of every given phase with geographical info.

        :param source_depth_in_km: Source depth in km
        :type source_depth_in_km: float
        :param source_latitude_in_deg: Source location latitude in degrees
        :type source_latitude_in_deg: float
        :param source_longitude_in_deg: Source location longitue in degrees
        :type source_longitude_in_deg: float
        :param station_latitude_in_deg: Station location latitude in degrees
        :type station_latitude_in_deg: float
        :param station_longitude_in_deg: Station location longitude in degrees
        :type station_longitude_in_deg: float
        :param phase_list: List of phases for which travel times should be
            calculated. If this is empty, all phases will be used.
        :type phase_list: list of str
        :return: List of ``Arrival`` objects, each of which has the time,
            corresponding phase name, ray parameter, takeoff angle, etc. as
            attributes.
        :rtype: :class:`Arrivals`
        """

        # Find the geodesic, azimuth and angular distance between epicenter
        # and station 
        g = self.ellipsoid.Inverse(source_latitude_in_deg, source_longitude_in_deg,
                              station_latitude_in_deg, station_longitude_in_deg)
        distance_in_deg = g['a12']
        azimuth = g['azi1']

        # Find the path(s) and decorate with lat and lon found from the direct
        # geodesic problem. We just change the path attribute of the
        # arrival object - everything else stays the same.
        arrivals = self.get_ray_paths(source_depth_in_km, distance_in_deg, 
                                      phase_list)
        line = self.ellipsoid.Line(source_latitude_in_deg, 
                                   source_longitude_in_deg, azimuth)
       
        pathList = []
        for arrival in arrivals:
            for path_point in arrival.path:
                pos = line.ArcPosition(np.degrees(path_point['dist']))
                diffTDG = np.array([(
                    path_point['p'],
                    path_point['time'],
                    path_point['dist'],
                    path_point['depth'],
                    pos['lat2'],
                    pos['lon2'])], dtype=TimeDistLoc)
                pathList.append(diffTDG)
            arrival.path = np.concatenate(pathList)

        return arrivals


    def get_pierce_points_geo(self, source_depth_in_km, source_latitude_in_deg,
                              source_longitude_in_deg, station_latitude_in_deg,
                              station_longitude_in_deg, phase_list=("ttall",)):
        """
        Return ray paths of every given phase with geographical info.

        :param source_depth_in_km: Source depth in km
        :type source_depth_in_km: float
        :param source_latitude_in_deg: Source location latitude in degrees
        :type source_latitude_in_deg: float
        :param source_longitude_in_deg: Source location longitue in degrees
        :type source_longitude_in_deg: float
        :param station_latitude_in_deg: Station location latitude in degrees
        :type station_latitude_in_deg: float
        :param station_longitude_in_deg: Station location longitude in degrees
        :type station_longitude_in_deg: float
        :param phase_list: List of phases for which travel times should be
            calculated. If this is empty, all phases will be used.
        :type phase_list: list of str
        :return: List of ``Arrival`` objects, each of which has the time,
            corresponding phase name, ray parameter, takeoff angle, etc. as
            attributes.
        :rtype: :class:`Arrivals`
        """

        # Find the geodesic, azimuth and angular distance between epicenter
        # and station 
        g = self.ellipsoid.Inverse(source_latitude_in_deg, source_longitude_in_deg,
                              station_latitude_in_deg, station_longitude_in_deg)
        distance_in_deg = g['a12']
        azimuth = g['azi1']

        # Find the path(s) and decorate with lat and lon found from the direct
        # geodesic problem. We just change the path attribute of the
        # arrival object - everything else stays the same.
        arrivals = self.get_pierce_points(source_depth_in_km, distance_in_deg, 
                                      phase_list)
        line = self.ellipsoid.Line(source_latitude_in_deg, 
                                   source_longitude_in_deg, azimuth)
       
        pathList = []
        for arrival in arrivals:
            for path_point in arrival.pierce:
                pos = line.ArcPosition(np.degrees(path_point['dist']))
                diffTDG = np.array([(
                    path_point['p'],
                    path_point['time'],
                    path_point['dist'],
                    path_point['depth'],
                    pos['lat2'],
                    pos['lon2'])], dtype=TimeDistLoc)
                pathList.append(diffTDG)
            arrival.pierce = np.concatenate(pathList)

        return arrivals



# Tomographic correction...
# =========================
#
# Given this has some one-time setup, it seems sensible
# to also give this an OO interface. However, we can only have
# one instance (see line 25 of tomo_predict.f90. What fun. 

import tomo_predict

class TomographicCorrection(object):

    def __init__(self, file_1d, file_3d, ellipsoid=geod.Geodesic.WGS84, 
                    taup_model="iasp91"):

        # Something to calculate the path
        self.earth_model = TauPyModelGeo(ellipsoid=ellipsoid,model=taup_model)

        # Setup the Fortran...
        assert not tomo_predict.tomo_predict.setup_done, \
            "Only one instance is permitted by the Fortran"
        tomo_predict.tomo_predict.setup(file_1d, file_3d)

    def calculate(self, evtlat, evtlon, evtdep, stalat, stalon, phase_list):

        arrivals = self.earth_model.get_ray_paths_geo(evtdep, evtlat, evtlon,
                    stalat, stalon, phase_list)

        dts = []
        for arrival in arrivals:
            dts.append(tomo_predict.tomo_predict.tomo_delay(arrival.path['lat'], 
                                  arrival.path['lon'], arrival.path['depth']))

        return dts




if __name__ == "__main__":      

        # Example useage

        # Data - two events...
        slat = 0.0
        slon = 90.0
        elat = 0.0
        elon = 0.0
        edep = 100.0

        corrector = TomographicCorrection('ak135.1D_vp', 'vdh3D_1999', ellipsoid=geod.Geodesic(a=6371000.0, f=0))

        dt = corrector.calculate(elat, elon, edep, slat, slon, ['PcP'])
        print 'PcP: ', dt




