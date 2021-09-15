#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Usual imports
import os
import numpy as np
import sys

## Plot settings if needed
# exec(open("/home/victor/acadwriting/Manuscrit/plots_settings.py").read())


## croco tools utility rewritten in python
from croco_tools_py import make_grid
from croco_tools_py import compute_rms_python

## Sediment segmentation
import croco_tools_py.make_grid_depth_SA as sediments

## Truth value definition
order = ["Roches", "Cailloutis", "Graviers", "Sables", "Sables fins", "Silts et Vases"]
truth = np.asarray([50, 25, 7, 1, 1.5e-1, 2e-2]) * 1e-3
z0btruth = dict(zip(order, truth))
grdname_truth = "/home/victor/croco_dahu2/Run/CROCO_FILES/croco_grd_truth.nc"

CROCO_DIR = "/home/victor/croco_dahu2"

for i, sed in enumerate(order):
    print("{}, {}mm".format(sed, truth[i] * 1e3))


def generate_truth_sediments():
    z0barray = sediments.map_grid_sediment_array(truth, sediments.idx_reduced_2, order)
    make_grid.make_grid(z0barray, grdname_truth)


## Generate the observations
def generate_observations():
    os.chdir(CROCO_DIR + "/Run")

    os.system("./croco ../AD/CONFIGS/ATLN2/croco-obs-py.in | grep WRT_HIS")
    os.system("mv CROCO_FILES/croco_his.nc CROCO_FILES/obs_truth.nc")


prefix_dir = "/home/victor/croco_dahu2/"
prefix_dir = "/home/victor/croco_dahu/croco-tap/croco/"


## Run the model and compute the rms
def croco_rms(
    z0b=truth, u=[0.5, 0.5], OBS_NAME="/home/victor/croco_dahu2/obs_sediments_true.nc"
):
    """
    Compute the RMS between croco evaluated with z0b and u, and the obs

    Parameters
    ----------
    z0b : Array of size 6, which corresponds to z0b of each sediments, optional

    u : Array of size 2, which indicates error on amplitude , optional

    OBS_NAME : Filename of the observations


    Returns
    -------
    out : the rms

    """
    grdname = prefix_dir + "Run/CROCO_FILES/croco_grd_py.nc"
    z0barray = sediments.map_grid_sediment_array(z0b, sediments.idx_reduced_2, order)
    make_grid.make_grid(z0barray, grdname)
    in_file = "/home/victor/croco_dahu2/Run/CROCO_FILES/croco_frc_5.nc"
    out_file = prefix_dir + "Run/CROCO_FILES/croco_frc_py.nc"
    make_grid.change_M2_amplitude(
        u, scale=0.01, in_file=in_file, out_file=out_file, idx=[0, 1]
    )
    os.chdir(CROCO_DIR + "/Run")
    os.chdir(prefix_dir)
    # os.system('ln -sf ../AD/CONFIGS/ATLN2/croco-py.in')
    # os.system('LD_LIBRARY_PATH=../.install/lib mpirun -np 1 /home/victor/croco_dahu/croco-tap/croco/RunC/croco ../AD/CONFIGS/ATLN2/croco-py.in > log.txt')
    # os.system('LD_LIBRARY_PATH=../.install/lib  mpirun -np 1 ./croco ../AD/CONFIGS/ATLN2/croco-py.in > log.txt')
    os.system("./croco-run-py > log.txt")
    RUN_NAME = "/home/victor/croco_dahu2/Run/CROCO_FILES/croco_py_his.nc"
    RUN_NAME = prefix_dir + "Run/CROCO_FILES/croco_py_his.nc"
    return compute_rms_python.compute_rms(RUN_NAME, OBS_NAME, None, None, None, 388)
