#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import netCDF4 as netcdf
import time


def create_grid(L, M, grdname, title):
    """
    Create an empty netcdf gridfile
    L: total number of psi points in x
    M: total number of psi points in y direction
    grdname: name of the grid file
    title: title in the netcdf file

    Further Information:
    http://www.croco-ocean.org

    Adapted from CROCOTOOLS
    """

    Lp = L + 1
    Mp = M + 1

    # File Creation
    print('Creating netCDF grid file: {}'.format(grdname))
    nw = netcdf.Dataset(grdname, 'w', format='NETCDF4')
    print('Initialization')
    # Set global attributes

    nw.title = title
    nw.date = time.asctime(time.localtime())
    nw.type = 'CROCO grid file'
    nw.set_fill_off()

    # Create dimensions
    xi_u = nw.createDimension('xi_u', L)
    eta_u = nw.createDimension('eta_u', Mp)

    xi_v = nw.createDimension('xi_v', Lp)
    eta_v = nw.createDimension('eta_v', M)

    xi_rho = nw.createDimension('xi_rho', Lp)
    eta_rho = nw.createDimension('eta_rho', Mp)

    xi_psi = nw.createDimension('xi_psi', L)
    eta_psi = nw.createDimension('eta_psi', M)

    one = nw.createDimension('one', 1)
    two = nw.createDimension('two', 2)
    four = nw.createDimension('four', 4)
    bath = nw.createDimension('bath', 1)

    # Create variables and attributes

    xl = nw.createVariable('xl', 'f8', ('one', ))
    xl.units = 'meter'
    xl.long_name = 'domain length in the XI-direction'

    el = nw.createVariable('el', 'f8', ('one', ))
    el.units = 'meter'
    el.long_name = 'domain length in the ETA-direction'

    depthmin = nw.createVariable('depthmin', 'f8', ('one', ))
    depthmin.units = 'meter'
    depthmin.long_name = 'Shallow bathymetry clipping depth'

    depthmax = nw.createVariable('depthmax', 'f8', ('one', ))
    depthmax.units = 'meter'
    depthmax.long_name = 'Deep bathymetry clipping depth'

    # ---------------

    spherical = nw.createVariable('spherical', 'S1', ('one', ))
    spherical.option_T = 'spherical'
    spherical.long_name = 'Grid type logical switch'

    angle = nw.createVariable('angle', 'f8', ('eta_rho', 'xi_rho'))
    angle.units = 'radian'
    angle.long_name = 'angle between xi axis and east'

    h = nw.createVariable('h', 'f8', ('eta_rho', 'xi_rho'))
    h.units = 'meter'
    h.long_name = 'Final bathymetry at RHO-points'

    hraw = nw.createVariable('hraw', 'f8', ('bath', 'eta_rho', 'xi_rho'))
    hraw.units = 'meter'
    hraw.long_name = 'Working bathymetry at RHO-points'

    alpha = nw.createVariable('alpha', 'f8', ('eta_rho', 'xi_rho'))
    alpha.long_name = 'Weights between coarse and fine grids at RHO-points'

    f = nw.createVariable('f', 'f8', ('eta_rho', 'xi_rho'))
    f.long_name = 'Coriolis parameter at RHO-points'
    f.units = 'second-1'

    # ---------------

    pm = nw.createVariable('pm', 'f8', ('eta_rho', 'xi_rho'))
    pm.long_name = 'curvilinear coordinate metric in XI'
    pm.units = 'meter-1'

    pn = nw.createVariable('pn', 'f8', ('eta_rho', 'xi_rho'))
    pn.long_name = 'curvilinear coordinate metric in ETA'
    pn.units = 'meter-1'

    dndx = nw.createVariable('dndx', 'f8', ('eta_rho', 'xi_rho'))
    dndx.long_name = 'xi derivative of inverse metric factor pn'
    dndx.units = 'meter'

    dmde = nw.createVariable('dmde', 'f8', ('eta_rho', 'xi_rho'))
    dmde.long_name = 'eta derivative of inverse metric factor pm'
    dmde.units = 'meter'

    # ---------------

    x_rho = nw.createVariable('x_rho', 'f8', ('eta_rho', 'xi_rho'))
    x_rho.long_name = 'x location of RHO-points'
    x_rho.units = 'meter'

    x_u = nw.createVariable('x_u', 'f8', ('eta_u', 'xi_u'))
    x_u.long_name = 'x location of U-points'
    x_u.units = 'meter'

    x_v = nw.createVariable('x_v', 'f8', ('eta_v', 'xi_v'))
    x_v.long_name = 'x location of V-points'
    x_v.units = 'meter'

    x_psi = nw.createVariable('x_psi', 'f8', ('eta_psi', 'xi_psi'))
    x_psi.long_name = 'x location of PSI-points'
    x_psi.units = 'meter'

    # ---------------

    y_rho = nw.createVariable('y_rho', 'f8', ('eta_rho', 'xi_rho'))
    y_rho.long_name = 'y location of RHO-points'
    y_rho.units = 'meter'

    y_u = nw.createVariable('y_u', 'f8', ('eta_u', 'xi_u'))
    y_u.long_name = 'y location of U-points'
    y_u.units = 'meter'

    y_v = nw.createVariable('y_v', 'f8', ('eta_v', 'xi_v'))
    y_v.long_name = 'y location of V-points'
    y_v.units = 'meter'

    y_psi = nw.createVariable('y_psi', 'f8', ('eta_psi', 'xi_psi'))
    y_psi.long_name = 'y location of PSI-points'
    y_psi.units = 'meter'

    # ---------------

    lon_rho = nw.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    lon_rho.long_name = 'longitude of RHO-points'
    lon_rho.units = 'degree_east'

    lon_u = nw.createVariable('lon_u', 'f8', ('eta_u', 'xi_u'))
    lon_u.long_name = 'longitude of U-points'
    lon_u.units = 'degree_east'

    lon_v = nw.createVariable('lon_v', 'f8', ('eta_v', 'xi_v'))
    lon_v.long_name = 'longitude of V-points'
    lon_v.units = 'degree_east'

    lon_psi = nw.createVariable('lon_psi', 'f8', ('eta_psi', 'xi_psi'))
    lon_psi.long_name = 'longitude of PSI-points'
    lon_psi.units = 'degree_east'

    # ---------------

    lat_rho = nw.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    lat_rho.long_name = 'latitude of RHO-points'
    lat_rho.units = 'degree_north'

    lat_u = nw.createVariable('lat_u', 'f8', ('eta_u', 'xi_u'))
    lat_u.long_name = 'latitude of U-points'
    lat_u.units = 'degree_north'

    lat_v = nw.createVariable('lat_v', 'f8', ('eta_v', 'xi_v'))
    lat_v.long_name = 'latitude of V-points'
    lat_v.units = 'degree_north'

    lat_psi = nw.createVariable('lat_psi', 'f8', ('eta_psi', 'xi_psi'))
    lat_psi.long_name = 'latitude of PSI-points'
    lat_psi.units = 'degree_north'

    # ---------------

    mask_rho = nw.createVariable('mask_rho', 'f8', ('eta_rho', 'xi_rho'))
    mask_rho.long_name = 'mask on RHO-points'
    mask_rho.option_0 = 'land'
    mask_rho.option_1 = 'water'

    mask_u = nw.createVariable('mask_u', 'f8', ('eta_u', 'xi_u'))
    mask_u.long_name = 'mask on U-points'
    mask_u.option_0 = 'land'
    mask_u.option_1 = 'water'

    mask_v = nw.createVariable('mask_v', 'f8', ('eta_v', 'xi_v'))
    mask_v.long_name = 'mask on V-points'
    mask_v.option_0 = 'land'
    mask_v.option_1 = 'water'

    mask_psi = nw.createVariable('mask_psi', 'f8', ('eta_psi', 'xi_psi'))
    mask_psi.long_name = 'mask on PSI-points'
    mask_psi.option_0 = 'land'
    mask_psi.option_1 = 'water'

    # ---------------

    z0b = nw.createVariable('z0b', 'f8', ('eta_rho', 'xi_rho'))
    z0b.long_name = 'Bottom friction on RHO-points (length)'
    z0b.units = 'meter'


    # Close file
    print('Closing netCDF grid file')
    nw.close()


if __name__ == '__main__':
    import sys
    create_grid(int(sys.argv[3]), int(sys.argv[2]), sys.argv[1], 'grid test py')
