# Flir bin2tif for s11 (Sorghum)

This script converts BIN files to GeoTIFF images for thermal (FLIR) sensors.

Re-calibrated transformer: converts bin to tif for gantry files (s11).
Differences between s10 and s11:
- S11 only contains center of agricultural plots, therefore no GCP lids are detectable. This is taken into consideration.
- Apart from removing the S10 correction, the code used is very similar to the one for s10.

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
`singularity run -B $(pwd):/mnt --pwd /mnt/ docker://phytooracle/flir_bin_to_tif_s10 -m <metadata.json> <bin_dir>`

