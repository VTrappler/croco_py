#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import geopandas
import matplotlib.ticker as mticker
import netCDF4 as netcdf

filepathMondiale = "/home/victor/shape_files/SEDIM_MONDIALE/nagene_V9_total.dbf"
df_mond = geopandas.read_file(filepathMondiale)

# Get bathy data and coordinates
nctarget = netcdf.Dataset("data/croco_grd_template.nc", "r", format="NETCDF4")
lon_rho = nctarget["lon_rho"][:]
lat_rho = nctarget["lat_rho"][:]
bathy = nctarget["h"][:]
nctarget.close()
points = np.stack([lon_rho, lat_rho]).T.reshape(-1, 2)


# inde = []
# for j, (xpt, ypt) in enumerate(points):
#     print(j)
#     tmp = []
#     for i, pol in enumerate(df_mond['geometry']):
#         if pol.contains(Point(xpt, ypt)):
#             tmp.append(i)
#     if len(tmp) == 0:
#         tmp.append(-1)
#     inde.append(tmp)

inde_1st = np.genfromtxt("data/df_indices.txt")

type_latlon = []
for ind in inde_1st:
    if ind == -1:
        type_latlon.append("Land")
    else:
        type_latlon.append(df_mond.iloc[int(ind)]["TYPE_VALEU"][2:])


mapping = {}
for typ in np.unique(type_latlon):
    mapping[typ] = np.where(np.asarray(type_latlon) == typ)
# np.save('/home/victor/mapping_type_point.npy', mapping)

mapping = np.load("data/mapping_type_point.npy", allow_pickle=True)
mapping = mapping[()]

bathy_ma = np.ma.masked_array(bathy)
bathy_ma[idx["Land"]] = np.ma.masked
plt.contourf(
    bathy_ma,
    levels=np.logspace(np.log(1), np.log(6000), 50),
    locator=mticker.LogLocator(),
)
plt.contour(bathy_ma, levels=[20, 30, 50, 100, 200, 500, 1000, 5000])
plt.colorbar()
plt.show()


for k, v in mapping.items():
    plt.scatter(points[v][:, 0], points[v][:, 1], s=10, alpha=0.7)
plt.show()

idx = np.full(
    (len(mapping.keys()) - 1, bathy.shape[0], bathy.shape[1]), False, dtype=bool
)
idx = dict()
for sed_typ in mapping.keys():
    print(sed_typ)
    idx[sed_typ] = np.full((bathy.shape[0], bathy.shape[1]), False, dtype=bool)
    for i in mapping[sed_typ][0]:
        idx[sed_typ][i % 166, i / 166] = True


def map_grid_sediment_dict(z0bdict, idxdict):
    if set(z0bdict.keys()) != set(idxdict.keys()):
        raise NameError
    z0barray = np.zeros_like(bathy)
    for ke in z0bdict.keys():
        z0barray[idxdict[ke]] = z0bdict[ke]
    return z0barray


def map_grid_sediment_array(z0b_list, idxdict, order):
    z0bdict = dict(zip(order, z0b_list))
    return map_grid_sediment_dict(z0bdict, idxdict)


z0bdict = {}
for i, ke in enumerate(mapping.keys()):
    z0bdict[ke] = i * 1e-3


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
order = [
    "Roches",
    "Cailloutis",
    "Graviers",
    "Sables",
    "Sables fins",
    "Silts",
    "Vases",
    "Argiles",
]
idx_reduced = dict()
for key in segm.keys():
    idx_reduced[key] = np.full_like(idx["Roche"], fill_value=False, dtype=bool)
    for v in segm[key]:
        try:
            idx_reduced[key] = np.logical_or(idx_reduced[key], idx[v])
        except (AttributeError, KeyError):
            pass
    print("{}: {}".format(key, np.sum(idx_reduced[key])))


z0b_dict = {
    "Roches": 5e-2,
    "Cailloutis": 1.5e-2,
    "Graviers": 7e-3,
    "Sables": 1e-3,
    "Sables fins": 1.5e-4,
    "Silts": 2e-5,
    "Vases": 2e-5,
    "Argiles": 2e-5,
}


z0barray = map_grid_sediment_dict(z0b_dict, idx_reduced)
plt.contourf(np.log(z0barray))
plt.show()
# make_grid.make_grid(z0barray, "/home/victor/croco_dahu2/croco_grd_sediments.nc")
