# FLIR BIN2TIF - Season 17+
This script converts BIN files to GeoTIFF images for thermal (FLIR) sensors.

This module includes a conversion for digital number to temperature by using camera internal calibration parameters and atmospheric values. Prior seasons have been calibrated using experimentally-derived polynomials relating DN, ambient temperature, and blackbody target temperature. This calibration process was not performed when the camera was replaced at the beginning of Season 17. Camera operation has changed as well, now focusing on full-field scanning. 

Note that some values are collected from the provided metadata, including:

sensor_fixed_metadata
```
    location in camera box x [m]
    location in camera box y [m]
    location in camera box z [m]
    field of view x [m]
    field of view y [m]
```
gantry_system_variable_metadata
```
    position x [m]
    position y [m]
    position z [m]
```
sensor_variable_metadata
```
    shutter temperature [K]
```

Also, a experimentally derived offset is applied:

```
    # TERRA-REF
    lon_shift = 0.000020308287

    # Drone
    lat_shift = 0.000018292 #0.000015258894
```

## Inputs

Raw metadata and BIN files.

## Outputs

Converted GeoTIFs in a single directory.

## Arguments and Flags
- **Required Arguments:** 
    - **Directory of flir files to decompress:** 'bin'
    - **raw metadata files:** '-m', '--metadata'
    - **Z-axis offset:** '-z', '--zoffset'

- **Optional Arguments**
    - **Output directory:** '-o', '--outdir', default='flir2tif_out/'
                                        
## Executing example (using singularity)
`singularity run -B $(pwd):/mnt --pwd /mnt/ docker://phytooracle/flir_bin_to_tif_s10 -m <metadata.json> <bin_dir> -z <zoffset>`

