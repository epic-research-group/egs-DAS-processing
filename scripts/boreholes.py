#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 15:52:37 2022

@author: spri902
"""
import os

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

from itertools import cycle
from glob import glob

def parse_surf_boreholes(file_path):
    """
    Parse the surf 4850 xyz file to dict of hole: {[(x, y, z), (x1, y1, z1)]}
    :param file_path: Path to borehole text file
    :return: dict
    """
    well_dict = {}
    with open(file_path, 'r') as f:
        next(f)
        for ln in f:
            ln = ln.rstrip('\n')
            line = ln.split(',')
            xm = float(line[3]) / 3.28084
            ym = float(line[4]) / 3.28084
            zm = float(line[5]) / 3.28084
            dp = float(line[2])
            name = line[0].split('-')[1]
            if name in well_dict:
                well_dict[name] = np.concatenate(
                    (well_dict[name],
                     np.array([xm, ym, zm, dp]).reshape(1, 4)))
            else:
                well_dict[name] = np.array([xm, ym, zm, dp]).reshape(1, 4)
    return well_dict

def depth_to_xyz(well_dict, well, depth):
    """
    Return xyz coords for depth in a given borehole
    :param well: Well string
    :param depth: Depth float
    :return:
    """
    easts, norths, zs, deps = np.hsplit(well_dict[well], 4)
    # Get closest depth point
    dists = np.squeeze(np.abs(depth - deps))
    x = easts[np.argmin(dists)][0]
    y = norths[np.argmin(dists)][0]
    z = zs[np.argmin(dists)][0]
    return (x, y, z)

def make_4100_boreholes(path, plot=False):
    """
    Take as-planned for 4100L boreholes and return array of xyz pts
    :param path: Path to excel file from Paul's drilling plan
    :return:
    """
    df = pd.read_excel(path)
    # Do the iterrows cause who cares
    well_dict = {}
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
    for i, row in df.iterrows():
        strt = row[['Easting (ft)', 'Northing (ft)',
                    'Height (ft)']].values * 0.3048
        bearing = np.deg2rad(row['Bearing (deg)'])
        tilt = np.deg2rad(row['Tilt (deg)'])
        unit_v = np.array([np.sin(bearing) * np.cos(tilt),
                           np.cos(bearing) * np.cos(tilt),
                           np.sin(tilt)])
        unit_v /= np.linalg.norm(unit_v)
        vect = unit_v * row['Length (ft)'] * 0.3048
        end = strt + vect
        pts = np.vstack(
            [np.linspace(strt[0], end[0], 200),
             np.linspace(strt[1], end[1], 200),
             np.linspace(strt[2], end[2], 200)]).T
        well_dict[row['New Name']] = pts
        if plot:
            ax.plot(xs=pts[:, 0], ys=pts[:, 1], zs=pts[:, 2], label=row['New Name'])
    if plot:
        fig.legend()
        ax.set_xlabel('Easting')
        ax.set_ylabel('Northing')
        ax.set_zlabel('Elevation')
        plt.show()
    return well_dict