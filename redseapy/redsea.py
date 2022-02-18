import logging
import pathlib
import numpy as np
import pandas as pd
import skimage
from skimage.io import imread
import skimage.measure
import skimage.morphology


# helper function 1


def ismember(a, b):
    bind = {}
    for i, elt in enumerate(b):
        if elt not in bind:
            bind[elt] = i
    return [
        bind.get(itm, None) for itm in a
    ]  # None can be replaced by any other "not in b" value


# helper function 2


def printProgressBar(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=100,
    fill="█",
    printEnd="\r",
):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print()


def run_redsea(
    tiff,
    seg_mask,
    markers_csv,
    output_dir,
    markers_of_interest=None,
    boundary_mode=2,  # Ignored atm
    compensation_mode=1,
    element_shape=2,
    element_size=2,
):
    # parameters for compensation
    # boundary_mode = boundaryMod = 2 # 2 means boundary, 1 whole cell
    # compensation_mode = REDSEAChecker = 1 # 1 means subtract+ reinforce
    # element_shape = elementShape = 2 # star, 1 == square size
    # element_size = elementSize = 2 # star or square extension size
    massDS = pd.read_csv(markers_csv)  # read the mass csv
    if "marker_name" not in massDS.columns:
        # Assume that csv only has a single column with the marker names
        # and no header
        massDS = pd.read_csv(markers_csv, header=None)
        massDS.columns = ["marker_name"]

    # remove the wavelenth
    markers = []
    for m in massDS.marker_name:
        for s in ["_488", "_555", "_570", "_647", "_660"]:
            m = m.replace(s, "")
        markers.append(m)
    massDS["marker_name"] = markers

    if markers_of_interest is None:
        markers_of_interest = markers
    normChannels = markers_of_interest
    #### should be inside the function
    normChannelsInds = ismember(normChannels, massDS["marker_name"])
    channelNormIdentity = np.zeros((len(massDS["marker_name"]), 1))
    # make a flag for compensation
    for i in range(len(normChannelsInds)):
        channelNormIdentity[normChannelsInds[i]] = 1

    clusterChannels = massDS["marker_name"]  # only get the label column
    clusterChannelsInds = np.where(np.isin(clusterChannels, massDS["marker_name"]))[
        0
    ]  # channel indexes

    logging.info("Reading image")
    tiff_img = imread(tiff)

    countsNoNoise = np.swapaxes(np.swapaxes(tiff_img, 0, 1), 1, 2)

    logging.info("Reading segmentation mask")
    segMat = imread(seg_mask)

    # rename the labels as this is a subset from the full sample
    lb = 0
    for i in sorted(np.unique(segMat)):
        segMat = np.where(segMat == i, lb, segMat)
        lb += 1

    logging.info("Quantifying markers before correction")
    labelNum = np.max(segMat)
    stats = skimage.measure.regionprops(segMat)
    newLmod = segMat

    ##### stuff related to mat finisehd
    channelNum = len(clusterChannels)  # how many channels

    data = np.zeros((labelNum, channelNum))
    dataScaleSize = np.zeros((labelNum, channelNum))
    cellSizes = np.zeros((labelNum, 1))

    centroid_xs = []
    centroid_ys = []

    for i in range(labelNum):  # for each cell (label)
        label_counts = [
            countsNoNoise[coord[0], coord[1], :] for coord in stats[i].coords
        ]  # all channel count for this cell
        data[i, 0:channelNum] = np.sum(
            label_counts, axis=0
        )  #  sum the counts for this cell
        dataScaleSize[i, 0:channelNum] = (
            np.sum(label_counts, axis=0) / stats[i].area
        )  # scaled by size
        cellSizes[i] = stats[i].area  # cell sizes
        centroid_xs.append(stats[i].centroid[0])
        centroid_ys.append(stats[i].centroid[1])

    [rowNum, colNum] = newLmod.shape
    cellNum = labelNum
    cellPairMap = np.zeros(
        (cellNum, cellNum)
    )  # cell-cell shared perimeter matrix container

    ## need to add border to the segmentation mask (newLmod in this case)
    newLmod_border = np.pad(newLmod, pad_width=1, mode="constant", constant_values=0)

    logging.info("Creating cell-cell contact matrix")
    # start looping the mask and produce the cell-cell contact matrix
    for i in range(rowNum):
        for j in range(colNum):
            if newLmod[i, j] == 0:
                tempMatrix = newLmod_border[
                    i : i + 3, j : j + 3
                ]  # the 3x3 window, xy shifted +1 due to border
                tempFactors = np.unique(tempMatrix)  # unique
                tempFactors = tempFactors - 1  # minus one for python index
                if len(tempFactors) == 3:  # means only two cells
                    cellPairMap[tempFactors[1], tempFactors[2]] = (
                        cellPairMap[tempFactors[1], tempFactors[2]] + 1
                    )  # count zero
                elif len(tempFactors) == 4:  # means three cells, three pairs
                    cellPairMap[tempFactors[1], tempFactors[2]] = (
                        cellPairMap[tempFactors[1], tempFactors[2]] + 1
                    )  # count zero
                    cellPairMap[tempFactors[1], tempFactors[3]] = (
                        cellPairMap[tempFactors[1], tempFactors[3]] + 1
                    )  # count zero
                    cellPairMap[tempFactors[2], tempFactors[3]] = (
                        cellPairMap[tempFactors[2], tempFactors[3]] + 1
                    )  # count zero
                elif len(tempFactors) == 5:  # means four cells, 6 pairs
                    cellPairMap[tempFactors[1], tempFactors[2]] = (
                        cellPairMap[tempFactors[1], tempFactors[2]] + 1
                    )  # count zero
                    cellPairMap[tempFactors[1], tempFactors[3]] = (
                        cellPairMap[tempFactors[1], tempFactors[3]] + 1
                    )  # count zero
                    cellPairMap[tempFactors[1], tempFactors[4]] = (
                        cellPairMap[tempFactors[1], tempFactors[4]] + 1
                    )  # count zero

                    cellPairMap[tempFactors[2], tempFactors[3]] = (
                        cellPairMap[tempFactors[2], tempFactors[3]] + 1
                    )  # count zero
                    cellPairMap[tempFactors[2], tempFactors[4]] = (
                        cellPairMap[tempFactors[2], tempFactors[4]] + 1
                    )  # count zero

                    cellPairMap[tempFactors[3], tempFactors[4]] = (
                        cellPairMap[tempFactors[3], tempFactors[4]] + 1
                    )  # count zero

    # double direction
    cellPairMap = cellPairMap + np.transpose(cellPairMap)

    ###############
    cellBoundaryTotal = np.sum(cellPairMap, axis=0)  # count the boundary

    # Cells without neighbors cause division by zero problems later one
    # removing them for now
    no_neighbor_cells = np.where(cellBoundaryTotal == 0)[0]
    if len(no_neighbor_cells) > 0:
        cellPairMap = np.delete(cellPairMap, no_neighbor_cells, axis=0)
        cellPairMap = np.delete(cellPairMap, no_neighbor_cells, axis=1)
        cellNum = cellNum - len(no_neighbor_cells)
        labelNum = cellNum
        cellBoundaryTotal = np.delete(cellBoundaryTotal, no_neighbor_cells, axis=0)
        data = np.delete(data, no_neighbor_cells, axis=0)
        dataScaleSize = np.delete(dataScaleSize, no_neighbor_cells, axis=0)
        cellSizes = np.delete(cellSizes, no_neighbor_cells, axis=0)
        centroid_xs = np.delete(np.array(centroid_xs), no_neighbor_cells, axis=0)
        centroid_ys = np.delete(np.array(centroid_ys), no_neighbor_cells, axis=0)

    ############### this step might cause error in ark version, double check with YH

    # devide to get fraction
    cellBoundaryTotalMatrix = np.tile(cellBoundaryTotal, (cellNum, 1))
    # cellBoundaryTotalMatrix = repmat(cellBoundaryTotal',[1 cellNum]);
    cellPairNorm = (
        compensation_mode * np.identity(cellNum) - cellPairMap / cellBoundaryTotalMatrix
    )
    cellPairNorm = np.transpose(
        cellPairNorm
    )  ### this is a werid bug in python, need to transpose
    # now starts the calculation of signals from pixels along the boudnary of cells
    MIBIdataNearEdge1 = np.zeros((cellNum, channelNum))

    logging.info("Performing correction")
    ##### A List of Items
    items = list(range(cellNum))
    l = len(items)
    printProgressBar(
        0, l, prefix="Progress:", suffix="Complete", length=50
    )  # progress bar
    #####

    ######pre-calculated shape
    if element_shape == 1:  # square
        square = skimage.morphology.square(2 * element_size + 1)
        square_loc = np.where(square == 1)
    elif element_shape == 2:  # diamond
        diam = skimage.morphology.diamond(
            element_size
        )  # create diamond shapte based on elementSize
        diam_loc = np.where(diam == 1)
    else:
        print("Error elementShape Value not recognized.")
    ############

    for i in range(cellNum):
        label = i + 1  # python problem
        [tempRow, tempCol] = np.where(newLmod == label)
        # sequence in row not col, should not affect the code
        for j in range(len(tempRow)):
            label_in_shape = []  # empy list in case
            # make sure not expand outside
            if (
                (element_size - 1 < tempRow[j])
                and (tempRow[j] < rowNum - element_size - 2)
                and (element_size - 1 < tempCol[j])
                and (tempCol[j] < colNum - element_size - 2)
            ):
                ini_point = [
                    tempRow[j] - element_size,
                    tempCol[j] - element_size,
                ]  # corrected top-left point

                if element_shape == 1:  # square
                    square_loc_ini_x = [item + ini_point[0] for item in square_loc[0]]
                    square_loc_ini_y = [item + ini_point[1] for item in square_loc[1]]

                    label_in_shape = [
                        newLmod[square_loc_ini_x[k], square_loc_ini_y[k]]
                        for k in range(len(square_loc_ini_x))
                    ]

                elif element_shape == 2:  # diamond
                    diam_loc_ini_x = [item + ini_point[0] for item in diam_loc[0]]
                    diam_loc_ini_y = [item + ini_point[1] for item in diam_loc[1]]
                    # finish add to ini point

                    label_in_shape = [
                        newLmod[diam_loc_ini_x[k], diam_loc_ini_y[k]]
                        for k in range(len(diam_loc_ini_x))
                    ]

            if 0 in label_in_shape:
                MIBIdataNearEdge1[i, :] = (
                    MIBIdataNearEdge1[i, :] + countsNoNoise[tempRow[j], tempCol[j], :]
                )

        # Update Progress Bar
        printProgressBar(i + 1, l, prefix="Progress:", suffix="Complete", length=50)

    ## fome final formatting

    MIBIdataNorm2 = np.transpose(np.dot(np.transpose(MIBIdataNearEdge1), cellPairNorm))
    # this is boundary signal subtracted by cell neighboor boundary
    MIBIdataNorm2 = (
        MIBIdataNorm2 + data
    )  # reinforce onto the whole cell signal (original signal)
    MIBIdataNorm2[MIBIdataNorm2 < 0] = 0  # clear out the negative ones
    # flip the channelNormIdentity for calculation
    rev_channelNormIdentity = np.ones_like(channelNormIdentity) - channelNormIdentity
    # composite the normalized channels with non-normalized channels
    # MIBIdataNorm2 is the matrix to return
    MIBIdataNorm2 = data * np.transpose(
        np.tile(rev_channelNormIdentity, (1, cellNum))
    ) + MIBIdataNorm2 * np.transpose(np.tile(channelNormIdentity, (1, cellNum)))
    # scale by size
    dataCompenScaleSize = MIBIdataNorm2 / cellSizes
    # some last steps
    ############ SKIP THE POSITIVE NUCLEAR IDENTITY FILTER
    ############ SHOULD ADD by user's choice

    labelIdentityNew2 = np.ones(cellNum)  ####### this part is the skipped line
    sumDataScaleSizeInClusterChannels = np.sum(
        dataScaleSize[:, clusterChannelsInds], axis=1
    )  # add all the cluster channels
    labelIdentityNew2[
        sumDataScaleSizeInClusterChannels < 0.1
    ] = 2  # remove the cells that does not have info in cluster channels
    # the function should return 4 varaibles
    # matrix
    dataCells = data[labelIdentityNew2 == 1, :]
    dataScaleSizeCells = dataScaleSize[labelIdentityNew2 == 1, :]
    dataCompenCells = MIBIdataNorm2[labelIdentityNew2 == 1, :]
    dataCompenScaleSizeCells = dataCompenScaleSize[labelIdentityNew2 == 1, :]

    # create the final matrixs ( 4 types of them)

    labelVec = np.where(labelIdentityNew2 == 1)
    labelVec = [
        item + 1 for item in labelVec
    ]  # python indexing difference need to add 1

    # get cell sizes
    cellSizesVec = cellSizes[labelIdentityNew2 == 1].flatten()

    # produce the matrices

    ## first dataframe
    dataL = pd.DataFrame({"CellID": labelVec[0].tolist(), "cell_size": cellSizesVec})
    dataCells_df = pd.DataFrame(dataCells)
    dataCells_df.columns = clusterChannels
    dataL_full = pd.concat((dataL, dataCells_df), axis=1)
    ### second
    dataScaleSizeL_df = pd.DataFrame(dataScaleSizeCells)
    dataScaleSizeL_df.columns = clusterChannels
    dataScaleSizeL_full = pd.concat((dataL, dataScaleSizeL_df), axis=1)
    ### third
    dataCompenL_df = pd.DataFrame(dataCompenCells)
    dataCompenL_df.columns = clusterChannels
    dataCompenL_full = pd.concat((dataL, dataCompenL_df), axis=1)
    ### forth
    dataCompenScaleSizeL_df = pd.DataFrame(dataCompenScaleSizeCells)
    dataCompenScaleSizeL_df.columns = clusterChannels
    dataCompenScaleSizeL_full = pd.concat((dataL, dataCompenScaleSizeL_df), axis=1)

    for d in [
        dataL_full,
        dataScaleSizeL_full,
        dataCompenL_full,
        dataCompenScaleSizeL_full,
    ]:
        d["x_centroid"] = centroid_xs
        d["y_centroid"] = centroid_ys

    output_dir = pathlib.Path(output_dir).resolve()
    logging.info(f"Writing output files to {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    dataScaleSizeL_full.to_csv(
        f"{output_dir}/single_cell_before_redsea.csv"
    )
    dataCompenScaleSizeL_full.to_csv(
        f"{output_dir}/single_cell_after_redsea.csv"
    )
