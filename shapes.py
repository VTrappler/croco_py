#!/usr/bin/env python
# -*- coding: utf-8 -*-

## Imports and plot config
import matplotlib.pyplot as plt
import numpy as np
import geopandas
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as netcdf
import pickle
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from collections import OrderedDict
exec(open('/home/victor/acadwriting/Manuscrit/plots_settings.py').read())
plt.rcParams['font.family'] = 'sans-serif'


# nctarget = netcdf.Dataset('/home/victor/sandbox/CROCO_FILES_test/croco_grd_template.nc',
#                           'r', format='NETCDF4')
## paths of sediments map
filepathMondiale = '/home/victor/shape_files/SEDIM_MONDIALE/nagene_V9_total.dbf'
filepathGDG = '/home/victor/shape_files/NATURES_FOND_500/GDG/SR-Golfe-Gascogne.dbf'
filepathMC = '/home/victor/shape_files/NATURES_FOND_500/MC/SR-Mers-celtiques.dbf'
filepathMMN = '/home/victor/shape_files/NATURES_FOND_500/MMN/SR-Mer-du-nord.dbf'

with open('indexes_sediments.pk', 'rb') as handle:
    idx = pickle.load(handle)


## 
def add_all_decorations(ax):
    lonmin = -9.
    lonmax = 1.
    latmin = 43.
    latmax = 51.

    ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND.with_scale('10m'), zorder=100, edgecolor='k')
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), zorder=101)
    ax.add_feature(cfeature.RIVERS.with_scale('10m'), zorder=101)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', zorder=102)
    ax.coastlines(resolution='50m')
    ax.set_ylim([43, 51])
    gl = ax.gridlines(alpha=0.1, draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlocator = mticker.FixedLocator([-10, -8, -6, -4, -2, 0, 2])
    gl.ylocator = mticker.FixedLocator([42, 44, 46, 48, 50, 52])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 6}
    gl.ylabel_style = {'size': 6}
    ax.set_aspect(1)


# dff = []
# for fp in [filepathGDG, filepathMC, filepathMMN]:
#     dff.append(geopandas.read_file(fp))
# df = pd.concat(dff)

df_mond = geopandas.read_file(filepathMondiale)
# df.TYPE_VALEU[:] = df.TYPE_VALEU[:][2:]

co = plt.get_cmap('tab20c')
# print(len(df_mond.TYPE_VALEU.unique()))
sed_type = df_mond.TYPE_VALEU.unique()
# fig = plt.figure(figsize = [5.748031686730317, 3.552478950810724])
# ax = plt.subplot(1, 2, 1, projection=ccrs.PlateCarree())
# df.plot(ax=ax, column='TYPE_VALEU', edgecolor=None)
# add_all_decorations(ax)
# ax.title(r'{} different types of sediments'.format(len(df['TYPE_VALEU'].unique())))

type_of_sediments = {'C': 'cailloutis',
                     'CG': 'cailloutis et graviers',
                     'GC': 'graviers et cailloutis',
                     'CV': 'cailloutis envasés',
                     'CS': 'cailloutis sables',
                     'G': 'graviers',
                     'GS': 'graviers sableux',
                     'SG': 'graviers sableux',
                     'GV': 'graviers envasés',
                     'VG': 'Vases et graviers',
                     'VSF': 'Vases et sable fins',
                     'Sg': 'sables grossiers',
                     'SGV': 'Sables graviers envasés',
                     'SC': 'Sables et cailloutis',
                     'SFC': 'Sables fins et cailloutis',
                     'S': 'sables',
                     'SF': 'sables fins',
                     'SV': 'sables vaseux',
                     'SFV': 'sables fins vaseux',
                     'VS': 'vases sableuses',
                     'V': 'vases',
                     'Roche': 'roches',
                     'VSi': 'vases silteuses',
                     'Si': 'Silts',
                     'A': 'argiles',
                     'ASi': 'argiles silteuses',
                     'SiA': 'silts argileux',
                     'SSi': 'sables et silts',
                     'SFSi': 'sables fins et silts',
                     'M': 'maerl'}

types_of_sediments_EN = {'C': 'Grits',
                         'CG': 'Grits and gravel',
                         'GC': 'Gravel and grits',
                         'CV': 'Muddy grits',
                         'CS': 'Sandy grits',
                         'G': 'Gravel',
                         'GS': 'Sandy gravels',
                         'SG': 'Gravely sand',
                         'GV': 'Muddy gravel',
                         'VG': 'Mud and gravel',
                         'VSF': 'Mud and fine sand',
                         'Sg': 'Coarse sand',
                         'SGV': 'Gravely muddy sand',
                         'SC': 'Sands and grits',
                         'SFC': 'Fine sands and grits',
                         'S': 'Sands',
                         'SF': 'Fine sands',
                         'SV': 'Muddy sands',
                         'SFV': 'Fine and muddy sands',
                         'VS': 'Sandy mud',
                         'V': 'Muds',
                         'Roche': 'Rocks',
                         'VSi': 'Silty muds',
                         'Si': 'Silts',
                         'A': 'Clays',
                         'ASi': 'Silty clay',
                         'SiA': 'Argillaceous silts',
                         'SSi': 'Sands and silts',
                         'SFSi': 'Fine sands and silts',
                         'M': 'Maerl'}

segm = {'Roches': ['Roche'],
        'Cailloutis': ['C', 'CG', 'CS', 'CV'],
        'Graviers': ['G', 'GC', 'GV', 'GS'],
        'Sables': ['S', 'SC', 'SG', 'SGV', 'SSi', 'SV'],
        'Sables fins': ['SF', 'SFC', 'SFV'],
        'Silts': ['Si', 'SiA'],
        'Argiles': ['ASi', 'A'],
        'Vases': ['V', 'VC', 'VG', 'VS', 'VSF']}

segm_SA = {r'$\theta^{(1)}$': ['C', 'CG', 'CS', 'CV'],
           r'$\theta^{(2)}$': ['G', 'GC', 'GV', 'GS'],
           r'$\theta^{(3)}$': ['Roche',
                                 'S', 'SC', 'SG', 'SGV', 'SSi', 'SV',
                                 'SF', 'SFC', 'SFV',
                                 'Si', 'SiA', 'ASi', 'A',
                                 'V', 'VC', 'VG', 'VS', 'VSF']
           }

segm_EN_FR = {'Rocks': ['Roches'],
              'Pebbles': ['Cailloutis'],
              'Gravels': ['Graviers'],
              'Sands': ['Sables'],
              'Fine Sands': ['Sables fins'],
              'Silts': ['Silts'],
              'Clays': ['Argiles'],
              'Muds': ['Vases']}


def invert_dict(d):
    inverse = dict()
    for key in d:
        # Go through the list that is saved in the dict:
        for item in d[key]:
            # Check if in the inverted dict the key exists
            if item not in inverse:
                # If not create a new list
                inverse[item] = [key]
            else:
                inverse[item].append(key)
    return inverse

total = dict()
for key in segm.keys():
    total[key] = np.full_like(idx['Roche'], fill_value=False, dtype=bool)
    for v in segm[key]:
        try:
            total[key] = np.logical_or(total[key], idx[v])
        except (AttributeError, KeyError):
            pass
    print('{}: {}'.format(key, np.sum(total[key])))
 

inv_segm = invert_dict(segm)
inv_segm_SA = invert_dict(segm_SA)

for i, key in enumerate(total.keys()):
    plt.subplot(3, 3, i+1)
    plt.contourf(total[key])
    plt.title(key)
plt.close()

not_present = np.array([9,   # NFCV
                        10,  # NFA
                        12,  # NFGV
                        13,  # NFASi
                        14,  # NFSiA
                        16,  # NFSSi
                        17,  # NFVG
                        18,  # NFCS
                        25]) # NFSFC


to_keep = []
for i, val in enumerate(df_mond.TYPE_VALEU.values):
    if val not in sed_type[not_present]:
        to_keep.append(i)


df_mond = df_mond.iloc[to_keep]

df_mond2 = df_mond.copy()

df_SA = df_mond2.copy()

df_mond2['TYPE_VALEU'] = df_mond2['TYPE_VALEU'].map(lambda x: inv_segm[x[2:]][0])
df_SA['TYPE_VALEU'] = df_SA['TYPE_VALEU'].map(lambda x: inv_segm_SA[x[2:]][0])

df_mond['TYPE_VALEU'] = df_mond['TYPE_VALEU'].map(lambda x: types_of_sediments_EN[x[2:]])

print('Full segmentation: {}'.format(len(df_mond.TYPE_VALEU.unique())))
print('Reduced segmentation: {}'.format(len(df_mond2.TYPE_VALEU.unique())))
print('SA segmentation: {}'.format(len(df_SA.TYPE_VALEU.unique())))



plt.figure(figsize=[5.748031686730316, 3.552478950810724])
ax2 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
df_mond.plot(ax=ax2, categorical=True,
             column='TYPE_VALEU', edgecolor=None, legend=True,
             legend_kwds={'bbox_to_anchor': (1.0, 1.05),
                          'fontsize': 6,
                          'frameon': False,
                          'loc': 'upper left'},
             cmap=co)
add_all_decorations(ax2)
ax2.set_title(r'Repartition of the different types of sediments')
plt.savefig('/home/victor/acadwriting/Manuscrit/Text/Chapter5/img/sediments_full.png',
            dpi=400, bbox_inches='tight')
plt.close()


plt.figure(figsize=[5.748031686730317, 3.552478950810724])
ax2 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
# for key in np.unique(df_mond2['TYPE_VALEU']):
#     df_iter = df_mond2.copy()
#     df_iter[df_iter['TYPE_VALEU'] != key] = np.nan
#     df_iter.plot(ax=ax2, categorical=True,
#                  column='TYPE_VALEU', edgecolor='none', legend=True,
#                  legend_kwds={'bbox_to_anchor': (1.0, 1.05),
#                               'fontsize': 8,
#                               'frameon': False,
#                               'loc': 'upper left'},
#                  cmap=co)
co = plt.get_cmap('tab20c')

color_mapping = {'Argiles': 'black',
                 'Cailloutis': co(0),
                 'Graviers': co(4),
                 'Roches': co(8),
                 'Sables': co(12),
                 'Sables fins': co(16),
                 'Silts': co(20),
                 'Vases': co(24)}
# cbrewer_map = ['#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd', 'black']
# cbrewer_map = ['#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#8c2d04']
# cbrewer_map = ['#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd', 'black']
# cbrewer_map = ['#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016450']
# cbrewer_map = ['#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#0c2c84']
# cbrewer_map = ['#f0f9e8','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#08589e']
cbrewer_map = ['#d73027', '#fc8d59', '#fee090', '#ffffbf', '#e0f3f8', '#91bfdb', '#4575b4']
names_leg = ['Rocks', 'Pebbles', 'Gravel', 'Sand', 'Fine Sand', 'Silts', 'Muds']#, 'Clays']
names_ = ['Roches', 'Cailloutis', 'Graviers', 'Sables', 'Sables fins', 'Silts', 'Vases']#, 'Argiles']
color_mapping = OrderedDict(zip(names_, cbrewer_map))

# color_mapping = {'Cailloutis': '#8DD3C7',
#                  'Graviers': '#FFFFB3',
#                  'Roches': '#BEBADA',
#                  'Sables': '#FB8072',
#                  'Sables fins': '#80B1D3',
#                  'Silts': '#FDB462',
#                  'Vases': '#B3DE69',
#                  'Argiles':'black'}
df_mond2.plot(ax=ax2, categorical=True,
              # column='TYPE_VALEU',
              edgecolor=None, legend=True,
             # legend_kwds={'bbox_to_anchor': (1.0, 1.05),
             #              'fontsize': 8,
             #              'frameon': False,
             #              'loc': 'upper left'},
              color=df_mond2['TYPE_VALEU'].map(color_mapping))
add_all_decorations(ax2)

handles = [mpatches.Patch(color=val, label=invert_dict(segm_EN_FR)[key][0])
           for key, val in color_mapping.items()]
plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
ax2.set_title('Repartition of the sediments')
plt.savefig('/home/victor/acadwriting/Manuscrit/Text/Chapter5/img/sediments_reduced_sserif.png', dpi=400)
# plt.show()
plt.close()

## grouped after SA
plt.figure(figsize=[5.748031686730317, 3.552478950810724])
ax2 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
co = plt.get_cmap('tab20c')
plt.rcParams['font.family'] = 'sans-serif'
cbrewer_map = ['#fc8d59','#ffffbf','#91bfdb']
names_ = [r'$\theta^{{({})}}$'.format(i) for i in [1, 2, 3]]
color_mapping = OrderedDict(zip(names_, cbrewer_map))
print(color_mapping)
print(df_SA['TYPE_VALEU'])
df_mond2.plot(ax=ax2, categorical=True,
              edgecolor=None, legend=True,
              color=df_SA['TYPE_VALEU'].map(color_mapping))
add_all_decorations(ax2)

handles = [mpatches.Patch(color=val, label=key)
           for key, val in color_mapping.items()]
plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=8)
ax2.set_title('Repartition of the sediments')
plt.savefig('/home/victor/acadwriting/Manuscrit/Text/Chapter5/img/sediments_reduced_SA_sserif.png', dpi=400)
plt.show()
plt.close()


# EOF ----------------------------------------------------------------------
