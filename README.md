[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Compensate cross-talk between neighboring cells

Based on the original code found at https://github.com/nolanlab/REDSEA

## Installation

```
pip install git+https://github.com/clemenshug/redseapy.git#egg=redseapy
```

## Usage

```
usage: redsea [-h] [--element-shape {star,square}]
              [--element-size ELEMENT_SIZE]
              [--markers-of-interest [MARKERS_OF_INTEREST ...]]
              IMAGE SEGMENTATION_MASK MARKERS DESTINATION

Correct cross-talk using the REDSEA method developed by Bai et.al.
https://doi.org/10.3389/fimmu.2021.652631

positional arguments:
  IMAGE                 Path to registered image in TIFF format
  SEGMENTATION_MASK     Path to segmentation mask image in TIFF format
  MARKERS               Path to a csv file containing a single column
                        with the marker names called 'marker_name'
  DESTINATION           Path to a new directory where uncorrected and
                        corrected single cell quantifications will be
                        saved

optional arguments:
  -h, --help            show this help message and exit
  --element-shape {star,square}
                        Shape of the element to be used for cross-talk
                        correction
  --element-size ELEMENT_SIZE
                        Size of the element to be used for cross-talk
                        correction
  --markers-of-interest [MARKERS_OF_INTEREST ...]
                        Only correct cross-talk for the specified
                        markers. Default: all markers
```

## Example

Example data are available for [download at mcmicro.org](https://mcmicro.org/datasets.html).

```
redsea exemplar-001/registration/exemplar-001.ome.tif \
  exemplar-001/segmentation/unmicst-exemplar-001/cellMask.tif \
  exemplar-001/markers.csv \
  exemplar-001/corrected-quantification
```
