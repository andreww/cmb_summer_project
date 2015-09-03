# python script tomocorr.py
# Python script to call c-shell script and fortran code to generate travel time
# corrections of a raypath through the 3D velocity model Van der Hilst et al (1999).
# Travel time corrections are produced for phases P and PcP.

# Written by H Bentham Aug 2015

# Requires:
# station lat and lon (slat, slon) earthquake lat, lon, depth (elat, elon, edep)
# and station (or path) identifier
# Output parameters:
# dtP and dtPcP (travel time variations for dtP and dtPcP)

# Other scripts and files required:
# c.tomo_predict_vdh	(c-shell script)
# tomo_predict_vdh.f	(fortran code)
# ak135.1D_vp		(1D velocity model ak135)
# vdh3D_1999		(3D velocity model Van der Hilst et al (1999))


# function run_tomocorr
def run_tomocorr(slat, slon, elat, elon, edep, statn):
	import csv
	import numpy
	import subprocess
	
	# produce a 2D array containing station and earthquake coordinate & identifier
	# information. Write array out as input file for fortran code
	TCinfo = numpy.column_stack([slat,slon,elat,elon,edep,statn])
	numpy.savetxt('input.latlon_tomo_predict', TCinfo, delimiter=" ", fmt='%s' )
	
	## Run shell scripts to call fortran code
	# For P
	subprocess.call(["./c.tomo_predict_vdh", "P"])
	f = open('vdh1999.P.dt','r')
	dtPf = f.read()
	f.close()
	
	# For PcP
	subprocess.call(["./c.tomo_predict_vdh", "PcP"])
	f = open('vdh1999.PcP.dt','r')
	dtPcPf = f.read()
	f.close()	

	return (dtPf, dtPcPf)
	
	
################################################
# Example of usage in main python script

# define arrays
dtP = []
dtPcP = []
slat1 = [ ]
slon1 = [ ]
elat1 = [ ]
elon1 = [ ]
edep1 = [ ]
statn1 = [ ]

# example input values
slat1 = [13.60, 11.60]
slon1 = [77.44, 80.44]
elat1 = [93.0, 91.0]
elon1 = [8.0, 7.0]
edep1 = [100, 120]
statn1 = ['GBA', 'HLMB']

# extract 3D travel time variations dtP & dtPcP from function run_tomocorr
dtP, dtPcP = run_tomocorr(slat1,slon1,elat1,elon1,edep1,statn1)
