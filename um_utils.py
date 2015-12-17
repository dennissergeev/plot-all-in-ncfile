# -*- coding: utf-8 -*-
from __future__ import print_function, division

import cf_units
import iris
import numpy as np
import os
import pandas as pd

default_interface = 'iris'

def replace_unknown_names(dataset, use=default_interface, default_name='unknown'):
    for ivar in dataset:
        if use.lower() == 'iris':
            if default_name in ivar.name():
                try:
                    stash_id = ivar.attributes['STASH'].__str__()
                    ivar.rename(find_by_stash_id(stash_id))
                except AttributeError:
                    print('Unable to rename, STASH attribute is missing')

        elif use.lower() == 'xray':
            import xray
            newnames_dict = dict()
            if default_name in dataset[ivar].name:
                try:
                    stash_id = dataset[ivar].attrs['um_stash_source']
                    newname = find_by_stash_id(stash_id)
                    newnames_dict[dataset[ivar].name] = newname
                except KeyError:
                    print('Unable to rename '+dataset[ivar].name+': STASH attribute is missing; skipping')

            dataset.rename(newnames_dict, inplace=True)

def find_by_stash_id(stash_id, path_to_stash=None, default_name='unknown'):
    if path_to_stash is None:
        path_to_stash = os.path.join(os.curdir,'stash.csv')
        try:
            df = pd.read_csv(path_to_stash)
        except:
            print('File stash.csv not found')
            print('Trying to download it from http://reference.metoffice.gov.uk ...')
            try:
                import urllib
                f = urllib.URLopener()
                f.retrieve('http://reference.metoffice.gov.uk/um/stash?_format=csv&_view=with_metadata',
                           path_to_stash)
            except:
                print('Download failed')
                print('Default name returned instead')

                return default_name

    df = pd.read_csv(path_to_stash)

    stash_label = df['rdfs:label'][df['@notation']==stash_id]
    if len(stash_label) > 0:
        return stash_label.values[0]
    else:
        print('Match not found, default name returned instead')
        return default_name

def get_cube(cubelist, cube_name, lazy_search=True):
    for icube in cubelist:
        if lazy_search:
            if cube_name.lower() in icube.name().lower():
                return icube
        else:
            if cube_name == icube.name():
                return icube

def unrotate_wind(cubelist,
                  uwind_name = 'x_wind', vwind_name = 'y_wind',
                  newcs=iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS),
                  replace=False, verbose=0):

        u = get_cube(cubelist, uwind_name)
        v = get_cube(cubelist, vwind_name)

        if u is not None or v is not None:
            oldcs = u.coord_system()
            if verbose > 1:
                print('Rotating winds from {}'.format(oldcs) + ' to {}'.format(newcs))
                print()
            u_rot, v_rot = iris.analysis.cartography.rotate_winds(u, v, newcs)
            if replace:
                cubelist[cubelist.index(u)] = u_rot
                cubelist[cubelist.index(v)] = v_rot
            else:
                cubelist.append(u_rot)
                cubelist.append(v_rot)
        else:
            print('u-wind or v-wind cubes not found. No winds rotating.')
