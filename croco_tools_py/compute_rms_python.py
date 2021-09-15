#!/usr/bin/python3


import netCDF4 as netcdf
import os
import numpy as np
import sys
import argparse

RUN_NAME = '/bettik/trapplev/croco/croco_SA/Run/CROCO_FILES/croco_SA_his.nc'
OBS_NAME = '/bettik/trapplev/croco/croco_SA/Run/obs.nc'
out_file = '/bettik/trapplev/croco/croco_SA/SA_rms.txt'


def main(ind, additional, RUN_NAME, OBS_NAME, outfile):
    index_start = 388
    index_end = index_start + 49
    zeta_obs = np.asarray(netcdf.Dataset(OBS_NAME, 'r', format="NETCDF4")['zeta'])
    zeta_run = np.asarray(netcdf.Dataset(RUN_NAME, 'r', format="NETCDF4")['zeta'])
    chi2 = (zeta_obs[index_start:index_end] - zeta_run)**2  # Get last 49 written obs in the obs file
    rss = chi2.sum()

    print('rss => ', rss)
    with open(outfile, 'a+') as fi:
        fi.write('{}, {}, {}\n'.format(ind, additional, rss))


def compute_rms(RUN_NAME, OBS_NAME, outfile, i, additional, start=388):
    index_start = start
    index_end = index_start + 49
    zeta_obs = np.asarray(netcdf.Dataset(OBS_NAME, 'r', format="NETCDF4")['zeta'])
    zeta_run = np.asarray(netcdf.Dataset(RUN_NAME, 'r', format="NETCDF4")['zeta'])
    chi2 = (zeta_obs[index_start:index_end] - zeta_run)**2  # Get last 49 written obs in the obs file
    rss = chi2.sum()
    if outfile is not None:
        with open(outfile, 'a+') as fi:
            fi.write('{}, {}, {}\n'.format(i, additional, rss))
    return rss


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', default=0, type=int)
    parser.add_argument('--add', default=0, type=str)
    parser.add_argument('--nc_his', default=RUN_NAME, type=str)
    parser.add_argument('--nc_obs', default=OBS_NAME, type=str)
    parser.add_argument('--output', default=out_file, type=str)
    args = parser.parse_args()
    main(ind=args.index,
         additional=args.add,
         RUN_NAME=args.nc_his,
         OBS_NAME=args.nc_obs,
         outfile=args.output)
