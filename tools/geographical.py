#!/usr/bin/env python
#
# Key primitives for cartesian to geographical 
#  coordinate transforms
#
# Copyright (C) Andrew Walker, 2010, 2011
#               <andrew.walker@bristol.ac.uk>
import math as m

def geog2cart(r, lat, lon):
    """Converts from geograpical to 
       cartesian coordinates. Input lat
       and lon in degrees (lat: 90 -> 0 -> -90 
       north to south, lon: 0 -> 360 or 0 -> 
       90 -> 180 -> -90 -> 0 E -> W). Output in 
       units of r."""

    theta = m.radians(90.0 - lat)
    if (lon < 0.0): 
        lon = 360.0 + lon # -ve lon is degrees W of grenwich
    phi = m.radians(lon)

    return (sper2cart(r, phi, theta))

def cart2geog(x, y, z):
    """Converts from our cartesian system (X3 to N pole,
       X1 to 0 E 0 N, X2 to 90 E 0 N) to geographical 
       coordinates (r, lat, lon). Radius units is the same
       as the cartesian units system (typically km)."""
    (r, theta, phi) = cart2sper(x, y, z)
    lon = m.degrees(phi)
    lat = (90.0 - m.degrees(theta))
    return(r, lat, lon)

def cart2sper(x, y, z):
    """Converts from cartesian to spherical coordiates,
       angular results in radians, radial in units of 
       cart system. Theta is the angle between the point
       and the +ve Z-axis, phi is the angle between the 
       XZ plane and the point measured paralell to the XY
       plane."""
    r = m.sqrt(x**2 + y**2 + z**2)
    # NB arctan2(a, b) is arctan(a/b) sorting out the quadrent
    phi = m.atan2(y, x)
    theta = m.acos(z/r)
    return(r, theta, phi)

def sper2cart (r, phi, theta):
    """Converts from sperical coordinates 
       to cartesian coordinates. Input theta
       and phi are in radians"""

    x = r * m.cos(phi) * m.sin(theta)
    y = r * m.sin(phi) * m.sin(theta)
    z = r * m.cos(theta)
    return (x, y, z)

def geog2str(r, lat, lon):
    return "r = " + str(r) + " lat = " + str(lat) + " lon = " + str(lon)

def cart2str(x, y, z):
    return "x = " + str(x) + " y = " + str(y) + " z = " + str(z)

def vincenty_direct(lat, lon, azimuth, distance,
             r_major=6378.1370, r_minor=6356.752314, r_sphere=None):
    """
    Calculate the location of a point on an ellipsoid given the azimuth
    and distance from a starting point using Vincenty's direct approach
    (T. Vincenty 1975, "Direct and inverse solutions of geodesics on the 
    ellipsoid with application of nested equations" Survey Review XXII 
    pp.88-93).

    Arguments are the lattitude and longitude of the starting point (in
    degrees) and the azimuth (in degrees, clockwise from north) and 
    distance (normally in km) of the geodesic. Optiionally the semi-major
    and semi-minor radius of the ellipsoid or the radius of a sphere 
    (if polar flattening is ignored) can be supplied. Units of distance
    must be consistant between these and the distance along the geodesic.
    These default to the WGS 84 ellipsoid in km.

    A length two tuple is returned with the first element containing the
    latititude and the second the longitude of the point (in degrees).
    """
   
    lat = m.radians(lat)
    lon = m.radians(lon)
    azi = m.radians(azimuth)
    

    if (r_sphere is not None):
        r_major = r_sphere
        r_minor = r_sphere
        f = 0.0
    else:
        f = (r_major-r_minor)/r_major

    
    U1 = m.atan((1.0-f) * m.tan(lat))

    sigma1 = m.atan2(m.tan(U1), m.cos(azi))
    sin_apha = m.cos(U1)*m.sin(azi)
    cos_2_alpha = 1 - sin_apha**2
    u_2 = cos_2_alpha * ((r_major**2 - r_minor**2) / r_minor**2) 
    A = 1 + (u_2 / 16384)*(4096 + u_2 * (-768 + u_2 * (320 - 175 * u_2)))
    B = (u_2/1024)*(256 + u_2*(-128+u_2*(74-47*u_2)))
    sig_0 = distance/(r_minor*A)
    sig = sig_0

    epsilon = 1E-12 # Accuracy 
    max_iter = 500
    for i in range(max_iter):
        sig_old = sig
        sig_2_m = 2.0*sigma1 + sig
        d_sig = B * m.sin(sig) * (
                  m.cos(sig_2_m) + 0.25*B*(
                    m.cos(sig)*(-1.0 + 2.0*m.cos(sig_2_m)**2.0) 
                  - (1.0/6.0)*B*m.cos(sig_2_m)*(-3.0+4.0*m.sin(sig)**2.0)
                  * (-3.0+4.0*m.cos(sig_2_m)**2)))
        sig = sig_0 + d_sig

        if (m.fabs(sig-sig_old) <= epsilon):
            # Found a solution in i iters...
            break

    lat2 = m.atan2(m.sin(U1)*m.cos(sig)+m.cos(U1)*m.sin(sig)*m.cos(azi),
               (1.0-f)*m.sqrt(sin_apha**2.0+(m.sin(U1)*m.sin(sig)-m.cos(U1)
                          *m.cos(sig)*m.cos(azi))**2.0))

    lam = m.atan2(m.sin(sig)*m.sin(azi),
                  m.cos(U1)*m.cos(sig) - m.sin(U1)*m.sin(sig)*m.cos(azi))

    C = (f/16.0) * cos_2_alpha * (4.0 + f*(4.0-3.0*cos_2_alpha))

    L = lam - (1.0-C)*f*sin_apha*(
        sig+C*m.sin(sig)*(m.cos(sig_2_m)
         +C*m.cos(sig)*(-1.0+2.0*m.cos(sig_2_m)**2.0)))

    lon2 = L + lon

    return(m.degrees(lat2),m.degrees(lon2))


def vincenty(lat1, lon1, lat2, lon2, 
             r_major=6378.1370, r_minor=6356.752314, r_sphere=None):
    """
    Calulates the distance and direction between two points on an 
    ellipsoid using Vincenty's inverse approach (T. Vincenty 1975,
    "Direct and inverse solutions of geodesics on the ellipsoid with
    application of nested equations" Survey Review XXII pp.88-93).

    Arguments are the lattitude and longitude of the two points (in 
    degrees) and, optionally, the semi-major and semi-minor radius of the
    ellipsoid or the radius of a sphere (if polar flattening is ignored). 
    These default to the WGS 84 ellipsoid. Changing the units of the 
    radii changes the units of the resultant distance, but note that the
    other arguments shouls always be in degrees.

    A length three tuple is returned. The first element is the distance 
    between the two points, the second and third are the azimuths of the 
    geodesic in degrees clockwise from north. 
    """

    lat1 = m.radians(lat1)
    lat2 = m.radians(lat2)
    lon1 = m.radians(lon1)
    lon2 = m.radians(lon2)

    if (r_sphere is not None):
        r_major = r_sphere
        r_minor = r_sphere
        f = 0.0
    else:
        f = (r_major-r_minor)/r_major

    U1 = m.atan((1.0-f) * m.tan(lat1))
    U2 = m.atan((1.0-f) * m.tan(lat2))
    L = lon2 - lon1

    epsilon = 1E-12 # Accuracy (10E-12 -> ~ 0.06mm)
    max_iter = 500 
    lam = L

    cU1 = m.cos(U1)
    cU2 = m.cos(U2)
    sU1 = m.sin(U1)
    sU2 = m.sin(U2)

    for i in range(max_iter):
        lam_old = lam
        sLam = m.sin(lam)
        cLam = m.cos(lam)
        sin_sig = m.sqrt((cU2*sLam)**2 + (cU1*sU2 - sU1*cU2*cLam)**2)
        cos_sig = sU1*sU2 + cU1*cU2*cLam
        sig = m.atan2(sin_sig,cos_sig)
        sin_alp = (cU1*cU2*sLam) / sin_sig
        cos2_alp = 1.0 - sin_alp**2
        if (cos2_alp == 0.0): 
            # equitorial line
            cos_2sigm = 100
            C = 0.0
        else:
            cos_2sigm = cos_sig - (2.0*sU1*sU2)/cos2_alp
            C = f/16.0 * cos2_alp * (4.0 + f*(4.0-3.0*cos2_alp))
        lam = L + (1.0 - C) * f * sin_alp * \
            (sig + C * sin_sig * (cos_2sigm + C * cos_sig * \
            (-1.0 + 2.0 * cos_2sigm**2)))
        if ((m.fabs(lam - lam_old)) <= epsilon):
            # Found a solution in i iters...
            break

    usq = cos2_alp * ((r_major**2 - r_minor**2) / r_minor**2)
    A = 1 + usq/16384 * (4096 + usq*(-768 + usq*(320 - 175*usq)))
    B = usq/1024 * (256 + usq*(-128 + usq*(74 - 47*usq)))
    del_sig = B * sin_sig * (cos_2sigm + 0.25*B*(cos_sig*( \
        -1 + 2*cos_2sigm**2) - (1.0/6.0)*B*cos_2sigm * ( \
        -3 + 4*sin_sig**2) * (-3 + 4 * cos_2sigm**2)))
    s = r_minor * A * (sig - del_sig)
    alp1 = m.atan2(cU2*m.sin(lam),(cU1*sU2-sU1*cU2*m.cos(lam)))
    alp2 = m.atan2(cU1*m.sin(lam),(cU1*sU2*m.cos(lam)-sU1*cU2))

    return (s, m.degrees(alp1), m.degrees(alp2))

def dms2dec(degs,mins,secs):
        """Converts angle given in degrees, mins and secs to
           angle in decimal degrees"""
        return (degs + mins/60.0 + secs/3600.0)


if __name__ == "__main__":
    print "Testing co-ordinate conversion"

    print "Geographical point:"
    print geog2str(1000, 0.0, 0.0)
    print "Cartesian point:"
    (x, y, z) = geog2cart(1000, 0.0, 0.0)
    print cart2str(x, y, z)
