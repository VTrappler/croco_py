#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import netCDF4 as netcdf
from crocotools_param import param_class
from create_grid import create_grid
from get_angle import get_angle
import sys

ct = param_class()


def rho2uvp(rfield):
    vfield = 0.5 * (rfield[:-1, :] + rfield[1:, :])
    ufield = 0.5 * (rfield[:, :-1] + rfield[:, 1:])
    pfield = 0.5 * (ufield[:-1, :] + ufield[1:, :])
    return ufield, vfield, pfield


def spheric_dist(lat1, lat2, lon1, lon2):
    """compute distances for a simple spheric earth"""
    R = 6367442.76
    lshift = np.abs(lon2 - lon1)
    lshift[lshift >= 180] = 360 - lshift[lshift >= 180]
    deg2rad = np.pi / 180.0
    lat1 = lat1 * deg2rad
    lat2 = lat2 * deg2rad
    lshift = lshift * deg2rad

    dist = R
    A = ((np.sin(lshift) * np.cos(lat2)) ** 2) + (
        ((np.sin(lat2) * np.cos(lat1)) - (np.sin(lat1) * np.cos(lat2) * np.cos(lshift)))
        ** 2
    )
    dist = R * np.arcsin(np.sqrt(A))
    return dist


def get_metrics(inputfile, close=True):
    if isinstance(inputfile, str):
        nc = netcdf.Dataset(inputfile, "r+", "NC_")
    else:
        nc = inputfile
    latu = nc["lat_u"][:]
    lonu = nc["lon_u"][:]
    latv = nc["lat_v"][:]
    lonv = nc["lon_v"][:]

    if close:
        nc.close()

    Mp, L = latu.shape
    M, Lp = latv.shape
    Lm = L - 1
    Mm = M - 1
    dx = np.zeros((Mp, Lp))
    dy = np.zeros((Mp, Lp))

    dx[:, 1:-1] = spheric_dist(latu[:, :-1], latu[:, 1:], lonu[:, :-1], lonu[:, 1:])
    dx[:, 0] = dx[:, 1]
    dx[:, -1] = dx[:, -2]

    dy[1:-1, :] = spheric_dist(latv[:-1, :], latv[1:, :], lonv[:-1, :], lonv[1:, :])
    dy[0, :] = dy[1, :]
    dy[-1, :] = dy[-2, :]
    pm = 1.0 / dx
    pn = 1.0 / dy

    dndx = np.zeros((Mp, Lp))
    dmde = np.zeros((Mp, Lp))

    dndx[1:-1, 1:-1] = 0.5 * (1 / pn[1:-1, 2:] - 1 / pn[1:-1, :-2])
    dmde[1:-1, 1:-1] = 0.5 * (1 / pm[2:, 1:-1] - 1 / pm[:-2, 1:-1])
    return pm, pn, dndx, dmde


def make_grid(z0b=None, grdname=ct.grdname):
    """
    Build a CROCO grid file with value of z0b of shape (166, 141)
    """
    print("\n Making the grid: {}".format(grdname))
    print(" Title: {}".format(ct.CROCO_title))
    print(" Resolution: 1/{}".format(str(1.0 / ct.dl)))
    if z0b is None:
        print("z0b set to {} by default".format(z0b))
        z0b = 0.008

    lonr = np.arange(ct.lonmin, ct.lonmax + ct.dl, ct.dl)

    i = 0
    latr = []
    latr.append(ct.latmin)
    while latr[i] <= ct.latmax:
        to_append = latr[i] + ct.dl * np.cos(latr[i] * np.pi / 180.0)
        latr.append(to_append)
        i += 1

    Lonr, Latr = np.meshgrid(lonr, latr)
    Lonu, Lonv, Lonp = rho2uvp(Lonr)
    Latu, Latv, Latp = rho2uvp(Latr)

    print("\n Create the grid file...")
    M, L = Latp.shape
    print(" LLm = {}\n MMm = {}".format(L - 1, M - 1))
    create_grid(L, M, grdname, ct.CROCO_title)

    print("\n Fill the grid file...")
    nc = netcdf.Dataset(grdname, "r+", format="NETCDF4")
    nc["lat_u"][:] = Latu
    nc["lon_u"][:] = Lonu
    nc["lat_v"][:] = Latv
    nc["lon_v"][:] = Lonv
    nc["lat_rho"][:] = Latr
    nc["lon_rho"][:] = Lonr
    nc["lat_psi"][:] = Latp
    nc["lon_psi"][:] = Lonp
    nc["z0b"][:] = z0b
    nc.close()
    if isinstance(z0b, (list, tuple, np.ndarray)):
        print("z0b ranges from {} to {}".format(z0b.min(), z0b.max()))
    else:
        print(" z0b set to {}".format(z0b))

    # ---------------------------------------------------------------------------------
    print("\n Compute the metrics...")

    pm, pn, dndx, dmde = get_metrics(grdname)
    xr = np.zeros_like(pm)
    yr = np.zeros_like(pm)

    for i in range(L):
        xr[:, i + 1] = xr[:, i] + 2 / (pm[:, i + 1] + pm[:, i])

    for j in range(M):
        yr[j + 1, :] = yr[j, :] + 2 / (pn[j + 1, :] + pn[j, :])

    xu, xv, xp = rho2uvp(xr)
    yu, yv, yp = rho2uvp(yr)
    dx = 1.0 / pm
    dy = 1.0 / pn
    dxmax = np.max(dx / 1000.0)
    dxmin = np.min(dx / 1000.0)

    dymax = np.max(dy / 1000.0)
    dymin = np.min(dy / 1000.0)
    print("\n Min dx = {} km - Max dx = {} km".format(dxmin, dxmax))
    print(" Min dy = {} km - Max dy = {} km".format(dymin, dymax))

    # Angle between XI and EAST at RHO points
    angle = get_angle(Latu, Lonu)

    # Coriolis
    f = 4 * np.pi * np.sin(np.pi * Latr / 180.0) * 366.25 / (24 * 3600 * 365.25)

    print(" Fill the grid file...")
    nc = netcdf.Dataset(grdname, "r+", format="NETCDF4")

    nctarget = netcdf.Dataset(
        "/home/victor/croco_dahu2/Run/CROCO_FILES/croco_grd_template.nc",
        "r",
        format="NETCDF4",
    )
    nc["pm"][:] = pm
    nc["pn"][:] = pn
    nc["dndx"][:] = dndx
    nc["dmde"][:] = dmde
    nc["x_u"][:] = xu
    nc["y_u"][:] = yu
    nc["x_v"][:] = xv
    nc["y_v"][:] = yv
    nc["x_rho"][:] = xr
    nc["y_rho"][:] = yr
    nc["x_psi"][:] = xp
    nc["y_psi"][:] = yp
    nc["angle"][:] = angle
    nc["f"][:] = f
    nc["spherical"][:] = "T"

    nc["h"][:] = nctarget["h"][:]
    nc["hraw"][:] = nctarget["hraw"][:]
    nc["mask_u"][:] = nctarget["mask_u"][:]
    nc["mask_v"][:] = nctarget["mask_v"][:]
    nc["mask_psi"][:] = nctarget["mask_psi"][:]
    nc["mask_rho"][:] = nctarget["mask_rho"][:]

    nc.close()
    nctarget.close()


def change_M2_amplitude(u, scale, in_file, out_file, idx=[0]):
    """Take the product of the initial amplitude of the first tide component in the in_file
    u should be normalized between 0 and 1:
    u=0.5 => no change of amplitude
    u=0.0 => amp *= (1-scale)
    u=1.0 => amp *= (1+scale)
    """
    nc_in = netcdf.Dataset(in_file, "r", format="NETCDF4")
    nc_out = netcdf.Dataset(out_file, "w", format="NETCDF4")

    toexclude = ["tide_Eamp"]
    # copy global attributes all at once via dictionary
    nc_out.setncatts(nc_in.__dict__)
    # copy dimensions
    for name, dimension in nc_in.dimensions.items():
        nc_out.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None)
        )
    # copy all file data except for the excluded
    for name, variable in nc_in.variables.items():
        if name not in toexclude:
            x = nc_out.createVariable(name, variable.datatype, variable.dimensions)
            nc_out[name][:] = nc_in[name][:]
            # copy variable attributes all at once via dictionary
            nc_out[name].setncatts(nc_in[name].__dict__)
        else:
            x = nc_out.createVariable(name, variable.datatype, variable.dimensions)
            uvec = np.ones(5)
            for j, index in enumerate(idx):
                uvec[index] = uvec[index] + scale * (2 * u[j] - 1)
            print(uvec)
            nc_out[name][:] = nc_in[name][:] * uvec[:, np.newaxis, np.newaxis]
            nc_out[name].setncatts(nc_in[name].__dict__)
    nc_in.close()
    nc_out.close()


# in_file = '/home/victor/croco_dahu/croco-tap/croco/RunC/CROCO_FILES/croco_frc_5.nc'
# out_file = '/home/victor/croco_dahu/croco-tap/croco/RunC/CROCO_FILES/croco_frc_test.nc'
# change_M2_amplitude(1.0, 0.05, in_file, out_file)


if __name__ == "__main__":
    N = len(sys.argv)
    if N == 2:
        z0b = float(sys.argv[1])
    make_grid(z0b)
