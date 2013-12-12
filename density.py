#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from optparse import OptionParser
import numpy as np

np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=120)

def define_plane_by_planar_stereonet_coordinates(x, y):
	dir = np.degrees(np.atan2(x,y))
	l = sqrt(x**2+y**2)
	dip = 90 - np.degrees (2 * np.atan(l))
	return [dir, dip]

def return_coss(dir, dip):
	cos3 = np.sin(np.radians(dir)) * np.cos(np.radians(dip))
	cos2 =                          -np.sin(np.radians(dip))
	cos1 = np.cos(np.radians(dir)) * np.cos(np.radians(dip))
	return cos1, cos2, cos3

def return_angle_between((dir_a, dip_a, cos1_a, cos2_a, cos3_a), \
                         (dir_b, dip_b, cos1_b, cos2_b, cos3_b)):
	angle = (cos1_a * cos1_b + cos2_a * cos2_b + cos3_a * cos3_b)
	angle = np.degrees(np.arccos(angle))

	nulls = np.ones(dir_b.shape) * 0
	half_pi = np.ones(dir_b.shape)*180

	angle = np.where(np.all(dir_a - dir_b) and np.all(dip_a - dip_b), \
	                    angle, nulls)
	angle = np.where(np.all(np.abs(dir_a - dir_b)-180) and np.all(dip_a+dip_b), 
	                    angle, half_pi)
	return angle

def gauss(m, sigma):
	return np.exp(-1 * m**2/(2*sigma**2))

#--------------------------------------------------
# Parsing of an options
optparser = OptionParser()
optparser.add_option('-s', '--sigma', dest='sigma',
                     help='size of a net cell',
                     default=5)
optparser.add_option('-r', '--resolution', dest='resolution',
                     help='density (resolution) of the grid',
                     default=50)
optparser.add_option('-i', '--input', dest='ifile',
                     help='input data file name')
optparser.add_option('-o', '--output', dest='ofile',
                     help='output data density file name')
options, args = optparser.parse_args(sys.argv[1:])

resolution = float(options.resolution) # amount of rasters
sigma = float(options.sigma) # sigma in degrees;
                             # defines the size of average winidow
ifile = options.ifile
ofile = options.ofile

# parsing of the data
data     = np.loadtxt(ifile)
data_dbl = (data + [180,0]) % 360 * [1, -1]
data     = np.append(data, data_dbl, axis=0)

#---------------------------------------------------
# prepare an empty grid of rasters with its attitudes

x     = np.linspace(-1,1,resolution)
y     = np.linspace(-1,1,resolution)
xx,yy = np.meshgrid(x,y)

# calculate dir and dip of a plane based on their
# x, y coords on a stereonet diagramm
dir   = np.degrees(np.arctan2(xx,yy))
l     = np.sqrt(xx**2+yy**2)
dip   = 90 - np.degrees(2*np.arctan(l))
nulls = np.empty(dip.shape)*np.NaN
dip   = np.where(dip>=0,dip,nulls)
dir   = np.where(dip>=0,dir,nulls) % 360

# calculate directional cosinuse of a planes of grid
cos1_b, cos2_b, cos3_b = return_coss(dir, dip)

density = np.ones(dir.shape)*0
sigma   = np.ones(dir.shape)*sigma

#----------------------------------------------------
# computing density for each plane of data

for line in data:
	data_dir = (line[0] -90) % 360 # ugly string but it works
	data_dip = line[1]
	if data_dip>=0:
		data_dip = 90 - data_dip
	else:
		data_dip = -1 * (90 + data_dip)
	print data_dir, data_dip
	cos1_a, cos2_a, cos3_a = return_coss(data_dir, data_dip)
	angle = return_angle_between((data_dir, data_dip, cos1_a, cos2_a, cos3_a),
								 (dir,      dip,      cos1_b, cos2_b, cos3_b))
	density += gauss(angle, sigma)

density  = np.nan_to_num(density)
density /= density.max()

#---------------------------------------------------
# save density to file
np.savetxt(ofile, density[::-1], fmt='%0.4f')
