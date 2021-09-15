#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import scipy
import numpy as np
import netCDF4 as netcdf
from crocotools_param import param_class
import matplotlib.pyplot as plt
ct = param_class()

taux_file = ct.coads_dir + 'taux.cdf'
taux_name = 'taux'
tauy_file = ct.coads_dir + 'tauy.cdf'
tauy_name = 'tauy'
#
#  Heat fluxes w3
#
shf_file = ct.coads_dir + 'netheat.cdf'
shf_name = 'netheat'
#
#  Fresh water fluxes (evaporation - precipitation)
#
swf_file = ct.coads_dir + 'emp.cdf'
swf_name = 'emp'
#
#  Sea surface temperature and heat flux sensitivity to the
#  sea surface temperature (dQdSST).
#  To compute dQdSST we need:
#    sat     : Surface atmospheric temperature
#    airdens : Surface atmospheric density
#    w3      : Wind speed at 10 meters
#    qsea    : Sea level specific humidity
#
sst_file = ct.coads_dir + 'sst.cdf'
sst_name = 'sst'
sat_file = ct.coads_dir + 'sat.cdf'
sat_name = 'sat'
airdens_file = ct.coads_dir + 'airdens.cdf'
airdens_name = 'airdens'
w3_file = ct.coads_dir + 'w3.cdf'
w3_name = 'w3'
qsea_file = ct.coads_dir + 'qsea.cdf'
qsea_name = 'qsea'

#  Sea surface salinity
#
sss_file = ct.coads_dir + 'sss.cdf'
sss_name = 'salinity'
#
#  Short wave radiation
#
srf_file = ct.coads_dir + 'shortrad.cdf'
srf_name = 'shortrad'
#
#
################### END USERS DEFINED VARIABLES #######################
#
# Title
#
print('\n {} \n'.format(ct.CROCO_title))
#
# Read in the grid
#
print(' Read in the grid...')
nc = netcdf.Dataset(ct.grdname, 'r+', 'NC_')
Lp = len(nc['xi_rho'][:])
Mp = len(nc['eta_rho'][:])
lon = nc['lon_rho'][:]
lat = nc['lat_rho'][:]
angle = nc['angle'][:]
nc.close()
cosa = np.cos(angle)
sina = np.sin(angle)


#
