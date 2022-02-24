[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Compensate cross-talk between neighboring cells

To take full advantage of the basic science and clinical potential of
multiplexed imaging technologies, various challenges, such as cell segmentation
and cellular feature extraction, must first be addressed. However, the
uncertainty in cell boundaries and technical noise across protein markers
in the image can cause inaccurate cell segmentation in dense tissue and create
conditions in which signals from adjacent cells spatially “bleed” into each
other. This leads to nonsensical cell states as determined by unsupervised
clustering methods. Recent efforts have led to the development of a novel
spatial cross-talk correction method called REinforcement Dynamic Spillover
EliminAtion (REDSEA, [PMID: 34295327](https://pubmed.ncbi.nlm.nih.gov/34295327/)).

Here, we present a command-line interface and a Docker container that allow REDSEA to be deployed in any execution environment.

## References

* When using the method, please cite [Bai and Zhu, et. al](https://www.frontiersin.org/articles/10.3389/fimmu.2021.652631/full)
* Original implementation: https://github.com/nolanlab/REDSEA

## Installation

* When running as a container, please ensure that you have [Docker installed](https://docs.docker.com/get-docker/). Note that when using Docker Desktop on Mac OS X, you may need to increase the amount of RAM that Docker is allowed to use, [which is set to 2GB by default](https://docs.docker.com/desktop/mac/#resources).
* When running in a local environment,
```
pip install git+https://github.com/labsyspharm/redseapy.git#egg=redseapy
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

### Example

Example data are available for [download at mcmicro.org](https://mcmicro.org/datasets.html).

```
redsea exemplar-001/registration/exemplar-001.ome.tif \
  exemplar-001/segmentation/unmicst-exemplar-001/cellMask.tif \
  exemplar-001/markers.csv \
  exemplar-001/corrected-quantification
```

### Running as a Docker container

From the directory that contains your data files (e.g., inside `exemplar-001`),
```
docker run -v "$PWD":/data --rm labsyspharm/redsea:0.1 redsea data/registration/exemplar-001.ome.tif \
  data/segmentation/unmicst-exemplar-001/cellMask.tif \
  data/markers.csv \
  data/corrected-quantification
```
where
* `-v "$PWD":/data` maps the current working directory (containing the data files) to be visible as `/data` inside the container
* `--rm` cleans up the container after it's done executing
* `0.1` specifies the container [version](https://github.com/labsyspharm/redseapy/releases) to execute.
