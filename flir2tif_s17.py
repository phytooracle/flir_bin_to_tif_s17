#!/usr/bin/env python3
"""
Author : Emmanuel Gonzalez, Michele Cosi, Holly Ellingson, Jeffrey Demieville
Note   : Parts of this code was initially developed by the AgPipeline and TERRA-REF teams. 
         Several comments were added based on those found in TinoDornbusch's original code.
Date   : 2020-07-09
Purpose: Convert FLIR .bin files to .tif (Season 17+)
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
import math


# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='Season 17 flir2tif',
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
                        type=float,
                        required=True)# Check with gantry operator for correct height for specified season

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

    fov_x, fov_y = float(meta['sensor_fixed_metadata']['field of view X [m]']), float(meta['sensor_fixed_metadata']['field of view y [m]'])
    
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
    # Camera-Specific constants output by FLIR camera 
    R = calibP['sensor_fixed_metadata']['calibration R'] # Function of integration time and wavelength; Planck Constant
    B = calibP['sensor_fixed_metadata']['calibration B'] # Function of wavelength; Planck Constant
    F = calibP['sensor_fixed_metadata']['calibration F'] # Positive value (0-1); Planck Constant
    J1 = calibP['sensor_fixed_metadata']['calibration J1'] # Global Gain
    J0 = calibP['sensor_fixed_metadata']['calibration J0'] # Global Offset

    # Constant Atmospheric transmission parameter by Flir
    a1 = calibP['sensor_fixed_metadata']['calibration alpha1']
    a2 = calibP['sensor_fixed_metadata']['calibration alpha2']
    X = calibP['sensor_fixed_metadata']['calibration X']
    b1 = calibP['sensor_fixed_metadata']['calibration beta1']
    b2 = calibP['sensor_fixed_metadata']['calibration beta2']

    # Constant for VPD computation (sqtrH2O)
    H2O_K1 = 1.56
    H2O_K2 = 0.0694
    H2O_K3 = -0.000278
    H2O_K4 = 0.000000685
    
    # Environmental factors
    # According to FLIR, atmospheric absorption under 10m object distance can be neglected, expecially under dry desert climate
	# Assumption: Ambient Temperature ~= Shutter Temperature
    shutter_temp = calibP['sensor_variable_metadata']['shutter temperature [K]']
    T = float(shutter_temp) - 273.15 # Proxy for ambient temperature from the gantry
    H = 0.1 # Gantry Relative Humidity; Try to pull dynamic value if possible
    D = 3.778 # Distance from sensor to target; scans in S17 onwards set camera_ref to 3.5 m from target, resulting in thermal_camera_ref at 3.778 m from target 
    E = 0.98 # Approximate emissivity of vegetation

    # Atmospheric Transmission

    # VPD
    H2OInGperM2 = H*math.exp(H2O_K1 + H2O_K2*T + H2O_K3*math.pow(T, 2) + H2O_K4*math.pow(T, 3))
    
    # Atmospheric Transmission Correction Tau
    tau = (X * math.exp(-math.sqrt(D / 2) * (a1 + b1 * math.sqrt(H2OInGperM2))) + (1 - X) * math.exp(-math.sqrt(D / 2) * (a2 + b2 * math.sqrt(H2OInGperM2))))

    # Object Radiation obj_rad = Theoretical object radiation * emissivity * atmospheric transmission
    obj_rad = im * E * tau

    # Atmospheric Radiation atm_rad= (1 - atmospheric transmission) * Theoretical atmospheric radiation
    theo_atm_rad = (R * J1 / (math.exp(B / shutter_temp) - F)) + J0
    atm_rad = np.tile((1 - tau) * theo_atm_rad, (640, 480))

    # Ambient Reflection Radiation: amb_refl_rad = (1 - emissivity) * atmospheric transmission * Theoretical Ambient Reflection Radiation
    theo_amb_refl_rad = (R * J1 / (math.exp(B / shutter_temp) - F)) + J0
    amb_refl_rad = np.tile((1 - E) * tau * theo_amb_refl_rad, (640, 480))

    # Total Radiation
    corr_pxl_val = obj_rad + atm_rad + amb_refl_rad
    
    im = rawData
    pxl_temp = B / np.log(R /(corr_pxl_val - J0) * J1 + F) - 273.15 # Radial Basis Function (RBF)

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
