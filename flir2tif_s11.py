#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez, Michele Cosi, Holly Ellingson, Jeffrey Demieville
Note   : Parts of this code was initially developed by the AgPipeline and TERRA-REF teams.
Date   : 2020-07-09
Purpose: Convert FLIR .bin files to .tif (Season 11)
"""

import argparse
import os
import sys
import logging
import json
import numpy as np
import glob
from terrautils.spatial import scanalyzer_to_utm, scanalyzer_to_latlon, geojson_to_tuples
from terrautils.formats import create_geotiff
import matplotlib.pyplot as plt
from osgeo import gdal, osr


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Season 11 flir2tif',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bin',
                        metavar='str',
                        help='Bin file to be converted to TIF')

    parser.add_argument('-m',
                        '--metadata',
                        help='Cleaned metadata file',
                        metavar='metadata',
                        type=str,
                        required=True)
                        #default='cleanmetadata_out')

    parser.add_argument('-z',
                        '--zoffset',
                        help='Z-axis offset',
                        metavar='z-offset',
                        type=int,
                        default=0.76)#0.76

    parser.add_argument('-o',
                        '--outdir',
                        help='Output directory where .tif files will be saved',
                        metavar='str',
                        type=str,
                        default='flir2tif_out')

    args = parser.parse_args()

    if '/' not in args.outdir:
        args.outdir = args.outdir + '/'

    return args

# --------------------------------------------------
def get_boundingbox(metadata, z_offset):

    with open(metadata) as f:
        meta = json.load(f)['lemnatec_measurement_metadata']

    loc_gantry_x = float(meta['sensor_fixed_metadata']['location in camera box x [m]'])
    loc_gantry_y = float(meta['sensor_fixed_metadata']['location in camera box y [m]'])
    loc_gantry_z = float(meta['sensor_fixed_metadata']['location in camera box z [m]'])

    gantry_x = float(meta['gantry_system_variable_metadata']['position x [m]']) + loc_gantry_x
    gantry_y = float(meta['gantry_system_variable_metadata']['position y [m]']) + loc_gantry_y
    gantry_z = float(meta['gantry_system_variable_metadata']['position z [m]']) + loc_gantry_z + z_offset#offset in m

    fov_x, fov_y = float(meta['sensor_fixed_metadata']['field of view x [m]']), float(meta['sensor_fixed_metadata']['field of view y [m]'])
    
    img_height, img_width = 640, 480
    
    B = gantry_z
    A_x = np.arctan((0.5*float(fov_x))/2)
    A_y = np.arctan((0.5*float(fov_y))/2)
    L_x = 2*B*np.tan(A_x)
    L_y = 2*B*np.tan(A_y)

    x_n = gantry_x + (L_x/2)
    x_s = gantry_x - (L_x/2)
    y_w = gantry_y + (L_y/2)
    y_e = gantry_y - (L_y/2)

    bbox_nw_latlon = scanalyzer_to_latlon(x_n, y_w)
    bbox_se_latlon = scanalyzer_to_latlon(x_s, y_e)

    # TERRA-REF
    lon_shift = 0.000020308287

    # Drone
    lat_shift = 0.000018292 #0.000015258894
    b_box =  ( bbox_se_latlon[0] - lat_shift,
                bbox_nw_latlon[0] - lat_shift,
                bbox_nw_latlon[1] + lon_shift,
                bbox_se_latlon[1] + lon_shift)

    return b_box, img_height, img_width


# --------------------------------------------------
def flirRawToTemperature(rawData, calibP):

    shutter_temp = calibP['sensor_variable_metadata']['shutter temperature [K]']
    T = float(shutter_temp) - 273.15

    P_5_outmean = [1.137440642331793e-11,
                   -7.151963918140453e-07,
                   2.040023288027391e-02,
                   -1.480567234537099e+02]
    P_15_outmean = [1.081311914979629e-11,
                    -7.016010881023338e-07,
                    2.054630019627413e-02,
                    -1.521561215301546e+02]
    P_20_outmean = [7.884866004076222e-12,
                    -5.627752964123624e-07,
                    1.841833557270094e-02,
                    -1.424489740528044e+02]
    P_25_outmean = [9.583147873422692e-12,
                    -6.411047671547955e-07,
                    1.957403307722059e-02,
                    -1.488744387542483e+02]
    P_30_outmean = [7.731929583673130e-12,
                    -5.450000399690083e-07,
                    1.788280850465480e-02,
                    -1.397155089900219e+02]
    P_35_outmean = [9.979352154351443e-12,
                    -6.638673059086900e-07,
                    2.015587753410061e-02,
                    -1.556220395053390e+02]
    P_40_outmean = [1.113388420010232e-11,
                    -7.376131006851630e-07,
                    2.162806444290634e-02,
                    -1.657425341330783e+02]
    P_45_outmean = [8.689237696307418e-12,
                    -6.008401296566917e-07,
                    1.914217995514052e-02,
                    -1.514361986681356e+02]
    T_list = [5, 15, 20, 25, 30, 35, 40, 45]
    a = [P_5_outmean[0], P_15_outmean[0], P_20_outmean[0],
         P_25_outmean[0], P_30_outmean[0], P_35_outmean[0],
         P_40_outmean[0], P_45_outmean[0]]
    b = [P_5_outmean[1], P_15_outmean[1], P_20_outmean[1],
         P_25_outmean[1], P_30_outmean[1], P_35_outmean[1],
         P_40_outmean[1], P_45_outmean[1]]
    c = [P_5_outmean[2], P_15_outmean[2], P_20_outmean[2],
         P_25_outmean[2], P_30_outmean[2], P_35_outmean[2],
         P_40_outmean[2], P_45_outmean[2]]
    d = [P_5_outmean[3], P_15_outmean[3], P_20_outmean[3],
         P_25_outmean[3], P_30_outmean[3], P_35_outmean[3],
         P_40_outmean[3], P_45_outmean[3]]

    # use numpy linear interpolation function to generate calibration coefficients for actual sensor temperature
    P_val = [np.interp(T, T_list, a), np.interp(T, T_list, b),
                       np.interp(T, T_list, c), np.interp(T, T_list, d)]
    im = rawData
    pxl_temp = P_val[0]*im**3 + P_val[1]*im**2 + P_val[2]*im + P_val[3]
    #pxl_temp = pxl_temp.astype(int)

    return pxl_temp


# --------------------------------------------------
def main():
    """Create TIF here"""

    args = get_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    bin_file = args.bin

    if bin_file is not None:
        with open(args.metadata, 'r') as mdf:

            full_md = json.load(mdf)['lemnatec_measurement_metadata']
            extractor_info = None

            if full_md:
                if bin_file is not None:
                    out_file = os.path.join(args.outdir, bin_file.split('/')[-1].replace(".bin", ".tif"))
                    gps_bounds_bin, img_height, img_width = get_boundingbox(args.metadata, args.zoffset)

                    raw_data = np.fromfile(bin_file, np.dtype('<u2')).reshape(
                        [480, 640]).astype('float')
                    raw_data = np.rot90(raw_data, 3)

                    tc = flirRawToTemperature(raw_data, full_md)

                    create_geotiff(tc, gps_bounds_bin, out_file, None,
                                True, extractor_info, None, compress=True)

                    print(f'Done. See output in {args.outdir}')


# --------------------------------------------------
if __name__ == '__main__':
    main()
