[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Compensate cross-talk between neighboring cells

Based on the original code found at https://github.com/nolanlab/REDSEA

To take full advantage of the basic science and clinical potential of
multiplexed imaging technologies, various challenges, such as cell segmentation
and cellular feature extraction, must first be addressed. However, the
uncertainty in cell boundaries and technical noise across protein markers
in the image can cause inaccurate cell segmentation in dense tissue and create
conditions in which signals from adjacent cells spatially “bleed” into each
other. This leads to nonsensical cell states as determined by unsupervised
clustering methods. Recent efforts have led to the development of a novel
spatial cross-talk correction method called REinforcement Dynamic Spillover
EliminAtion (REDSEA, PMID: 34295327).

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

Correct cross-talk between neighboring cells in immuno-fluorescence images
using the REDSEA method developed by Bai et.al.
https://doi.org/10.3389/fimmu.2021.652631

positional arguments:
  IMAGE                 Path to registered image in TIFF format
  SEGMENTATION_MASK     Path to segmentation mask image in TIFF format
  MARKERS               Path to a csv file containing a single column with
                        the marker names called 'marker_name'
  DESTINATION           Path to a new directory where uncorrected and
                        corrected single cell quantifications will be saved

optional arguments:
  -h, --help            show this help message and exit
  --element-shape {star,square}
                        Shape of the element to be used for cross-talk
                        correction. See Fig S1C in the READSEA manuscript.
                        The authors recommend the star shape.
  --element-size ELEMENT_SIZE
                        Size in pixels of the element to be used for cross-
                        talk correction. The REDSEA authors recommend 2
  --markers-of-interest [MARKERS_OF_INTEREST ...]
                        Only correct cross-talk for the specified markers.
                        Default: all markers
```

## Example

Example data are available for [download at mcmicro.org](https://mcmicro.org/datasets.html).

```
redsea exemplar-001/registration/exemplar-001.ome.tif \
  exemplar-001/segmentation/unmicst-exemplar-001/cellMask.tif \
  exemplar-001/markers.csv \
  exemplar-001/corrected-quantification
```
