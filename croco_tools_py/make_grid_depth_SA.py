#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import netCDF4 as netcdf
import csv
import sys

# sys.path.append('/bettik/trapplev/croco/croco_SA/croco_tools_py/')
sys.path.append("/home/victor/croco_dahu2/croco_tools_py/")
import make_grid
import pickle
import argparse

target_path = "/bettik/trapplev/croco/croco_SA/Run/CROCO_FILES/croco_grd_template.nc"
target_path = "/home/victor/croco_dahu2/Run/CROCO_FILES/croco_grd_template.nc"
nctarget = netcdf.Dataset(target_path, "r", format="NETCDF4")
bathy = nctarget["h"][:]
nctarget.close()

# depth_bins = [10., 20, 30, 50, 100, 200, 5000]
depth_bins = [0.0, 30, 60, 100, 5000]
with open("/home/victor/indexes_sediments.pk", "rb") as handle:
    idx = pickle.load(handle)
# depth_bins = [0., 5000]
# depth_bins = [0., 60., 100., 200., 5000.]
# depth_bins = [0., 60., 5000.]
idx_reduced = dict()
segm = {
    "Roches": ["Roche"],
    "Cailloutis": ["C", "CG", "CS", "CV"],
    "Graviers": ["G", "GC", "GV", "GS"],
    "Sables": ["S", "SC", "SG", "SGV", "SSi", "SV"],
    "Sables fins": ["SF", "SFC", "SFV"],
    "Silts": ["Si", "SiA"],
    "Argiles": ["ASi", "A"],
    "Vases": ["V", "VC", "VG", "VS", "VSF"],
}

segm_2 = {
    "Roches": ["Roche"],
    "Cailloutis": ["C", "CG", "CS", "CV"],
    "Graviers": ["G", "GC", "GV", "GS"],
    "Sables": ["S", "SC", "SG", "SGV", "SSi", "SV"],
    "Sables fins": ["SF", "SFC", "SFV"],
    "Silts et Vases": ["Si", "SiA", "V", "VC", "VG", "VS", "VSF"],
}

for key in segm.keys():
    idx_reduced[key] = np.full_like(idx["Roche"], fill_value=False, dtype=bool)
    for v in segm[key]:
        try:
            idx_reduced[key] = np.logical_or(idx_reduced[key], idx[v])
        except (AttributeError, KeyError):
            pass
    print("{}: {}".format(key, np.sum(idx_reduced[key])))


idx_reduced_2 = dict()
for key in segm_2.keys():
    idx_reduced_2[key] = np.full_like(idx["Roche"], fill_value=False, dtype=bool)
    for v in segm_2[key]:
        try:
            idx_reduced_2[key] = np.logical_or(idx_reduced_2[key], idx[v])
        except (AttributeError, KeyError):
            pass
    print("{}: {}".format(key, np.sum(idx_reduced_2[key])))


def map_grid_by_depth(z0bvalues, depth_idx):
    """
    Create a np array of values, from a list of values for z0b,
    based on the indexes of depth_idx
    """
    z0barray = np.zeros_like(bathy)
    # print('depth_idx={}'.format(len((depth_idx))))
    # print('z0bvalues={}'.format(len(z0bvalues)))
    if len(depth_idx) != len(z0bvalues):
        print("TODO: add error message")
    else:
        for i, z in enumerate(z0bvalues):
            # print(i)
            # print(depth_idx[i])
            z0barray[depth_idx[i]] = z
    return z0barray


def make_grid_by_depth(z0bvalues, depth_bins, grd_file_modif=False, grdname=None):
    """Create and np array and nc grid file, which corresponds to the
    list of z0b mapped in the different bins of the depths_bins list
    """
    nctarget = netcdf.Dataset(target_path, "r", format="NETCDF4")
    bathy = nctarget["h"][:]
    nctarget.close()
    idx = np.full(
        (len(depth_bins) - 1, bathy.shape[0], bathy.shape[1]), False, dtype=bool
    )
    print(len(depth_bins))
    for j in range(len(depth_bins) - 1):
        idx[j] = np.logical_and((bathy >= depth_bins[j]), (bathy < depth_bins[j + 1]))
    z0barray = map_grid_by_depth(z0bvalues, depth_idx=idx)
    if grd_file_modif:
        make_grid.make_grid(z0barray, grdname)
    else:
        return z0barray


def map_grid_sediment_dict(z0bdict, idxdict):
    """Create a numpy array based on the dictionary of z0bvalues and the
    dictionary of indices

    """
    if set(z0bdict.keys()) != set(idxdict.keys()):
        raise NameError
    z0barray = np.zeros_like(bathy)
    for ke in z0bdict.keys():
        z0barray[idxdict[ke]] = z0bdict[ke]
    return z0barray


def map_grid_sediment_array(z0b_list, idxdict, order):
    """Create numpy array based on list of values, order given and
    dictionary of indices

    """
    z0bdict = dict(zip(order, z0b_list))
    return map_grid_sediment_dict(z0bdict, idxdict)


# z0barray = make_grid_by_depth(range(len(depth_bins) - 1), depth_bins, plot=False)
in_file = "/bettik/trapplev/croco/croco_SA/Run/CROCO_FILES/croco_frc_5.nc"
out_file = "/bettik/trapplev/croco/croco_SA/Run/CROCO_FILES/croco_frc_SA.nc"
CSV_FILE = "/bettik/trapplev/croco/croco_SA/SA_DoE.csv"


in_file = "/home/victor/croco_dahu2/Run/CROCO_FILES/croco_frc_5.nc"
out_file = "/home/victor/croco_dahu2/Run/CROCO_FILES/croco_frc_SA.nc"
CSV_FILE = "/home/victor/croco_dahu2/Sobol_sediments.csv"
grdname = "/home/victor/croco_dahu2/croco_grd_AD.nc"


def main(
    line_number=1,
    DoE=CSV_FILE,
    grdname=grdname,
    ncfrcin=in_file,
    ncfrcout=out_file,
    dimK=4,
    dimU=2,
):
    """
    Keyword Arguments:
    line_number  -- linenumber of the line to pass to the model
    DoE          -- csv file of design of experiment: K, U
    grdname      -- grid file to create
    in_frc_file  -- forcing file to take and modify
    out_frc_file -- out forcing file modified
    dimK         -- dimension of K
    dimU         -- dimension of U
    """
    reader = csv.reader(open(DoE, "r"))
    order = [
        "Roches",
        "Cailloutis",
        "Graviers",
        "Sables",
        "Sables fins",
        "Silts et Vases",
        # 'Vases',
        # 'Argiles'
    ]

    for i, row in enumerate(reader):
        if (i + 1) != line_number:
            pass
        else:
            print("SA: input:", row)
            z0b = np.asarray([float(r) for r in row[:dimK]])
            uu = np.asarray([float(r) for r in row[dimK:]])
            print("SA z0b: ", z0b)
            print("SA u: ", uu)
            # make_grid_by_depth(z0b, depth_bins, grd_file_modif=True, grdname=grdname)

            z0barray = map_grid_sediment_array(z0b, idx_reduced_2, order)
            make_grid.make_grid(z0barray, grdname)

            print("SA: Modify forcing")
            if dimU > 0:
                make_grid.change_M2_amplitude(
                    uu, scale=0.01, in_file=ncfrcin, out_file=ncfrcout, idx=[0, 1]
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--DoE", help="Design of experiments", type=str, default=CSV_FILE
    )
    parser.add_argument("--line", help="Corresponding line of the input file", type=int)
    parser.add_argument("--ncgrd", help="grd file", type=str)
    parser.add_argument(
        "--ncfrcin", help="forcing file to modify", default=in_file, type=str
    )
    parser.add_argument(
        "--ncfrcout", help="forcing file modified", type=str, default=out_file
    )
    parser.add_argument("--dimK", help="dimK", type=int, default=4)
    parser.add_argument("--dimU", help="dimU", type=int, default=2)
    args = parser.parse_args()
    print("DoE file: {}".format(args.DoE))
    main(
        line_number=args.line,
        DoE=args.DoE,
        grdname=args.ncgrd,
        ncfrcin=args.ncfrcin,
        ncfrcout=args.ncfrcout,
        dimK=args.dimK,
        dimU=args.dimU,
    )


# EOF ----------------------------------------------------------------------
