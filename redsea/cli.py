import argparse
import logging
import os

import numpy as np
import pandas as pd

from . import redsea as redsea_mod


def redsea():
    parser = argparse.ArgumentParser(
        description="""Correct cross-talk between neighboring cells
        in immuno-fluorescence images using the REDSEA method
        developed by Bai et.al.

        https://doi.org/10.3389/fimmu.2021.652631
        """,
    )
    parser.add_argument(
        "IMAGE",
        help="Path to registered image in TIFF format",
    )
    parser.add_argument(
        "SEGMENTATION_MASK",
        help="Path to segmentation mask image in TIFF format",
    )
    parser.add_argument(
        "MARKERS",
        help="Path to a csv file containing a single column with the marker names called 'marker_name'",
    )
    parser.add_argument(
        "DESTINATION",
        help="Path to a new directory where uncorrected and corrected single cell quantifications will be saved",
    )
    parser.add_argument(
        "--element-shape",
        choices=["star", "square"],
        default="star",
        help="Shape of the element to be used for cross-talk correction. See Fig S1C in the READSEA "
        "manuscript. The authors recommend the star shape.",
    )
    parser.add_argument(
        "--element-size",
        type=int,
        default=2,
        help="Size in pixels of the element to be used for cross-talk correction. The REDSEA authors recommend 2",
    )
    parser.add_argument(
        "--markers-of-interest",
        nargs="*",
        default=None,
        help="Only correct cross-talk for the specified markers. Default: all markers",
    )
    args = parser.parse_intermixed_args()
    logging.basicConfig(
        format="%(processName)s %(asctime)s %(levelname)s: %(message)s",
        level=os.environ.get("LOGLEVEL", "INFO").upper(),
    )
    logging.debug(args)
    # parameters for compensation
    # boundary_mode = boundaryMod = 2 # 2 means boundary, 1 whole cell
    # compensation_mode = REDSEAChecker = 1 # 1 means subtract+ reinforce
    # element_shape = elementShape = 2 # star, 1 == square size
    # element_size = elementSize = 2 # star or square extension size
    redsea_mod.run_redsea(
        args.IMAGE,
        args.SEGMENTATION_MASK,
        args.MARKERS,
        args.DESTINATION,
        element_shape=2 if args.element_shape == "star" else 1,
        element_size=args.element_size,
        markers_of_interest=args.markers_of_interest,
    )
    logging.info("Done")
