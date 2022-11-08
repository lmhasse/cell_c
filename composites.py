import argparse

import numpy as np
import os
import matplotlib.pyplot as plt
from PIL import Image
import tifffile
import csv
import PIL
from tqdm import tqdm


def generate_map(image_file, x, y, downsample):
    """ Generates np array from binary image

    :param image_file: string, path to binary image
    :param x: int, x of upper left corner of first tile in QuPath
    :param y: int, y of upper left corner of first tile in QuPath
    :param downsample: int, downsample factor of QuPath export
    :return: nparray
    """
    # turns off decompression bomb error
    Image.MAX_IMAGE_PIXELS = None

    # binary image needed
    image = Image.open(image_file)

    nparray = np.asarray(image)

    # scaling up np array in x and y direction
    np2 = np.repeat(nparray, downsample, axis=0)
    np3 = np.repeat(np2, downsample, axis=1)

    # cut off upper and left borders that are not included in tiles
    np4 = np3[x:, y:]

    # returns map where erythrocytes areas are 1, normal areas 0
    return np4 


def composite_prediction(pred_path, outfile_name, map_file, x, y, offset_x, offset_y, downsample):
    """Generates composite file out of several prediciton.txt generated from immunet

    :param pred_path: string, path to folder containing predictions
    :param outfile_name: string, composite files name
    :param map_file: string, path to binary image file containing erythrocyte map
    :param x: int, x coord of first tile in QuPath
    :param y: int, y coord of first tile in QuPath
    :param offset_x: int, name of first tile x,y
    :param offset_y: int, name of first tile x,y
    :param downsample: int, factor of downsampling in QuPath map export
    :return: void
    """

    # composite file of all prediction.txt, "a" for appending
    outfile = open(outfile_name, "a")

    # generates np array map out of binary image
    erythrocyte_map = generate_map(map_file, x, y, downsample)

    # iterate over all tiles in project (e.g. "454,213")
    for tile in tqdm(os.listdir(pred_path)):

        # iterate over files in tile folder
        for file in os.listdir(pred_path + tile):

            if file == "prediction.txt":
                # for identifying cells to tiles write tile
                outfile.write("tile:" + tile + "\n")

                with open(pred_path + tile + "/" + file) as f:

                    lines = f.readlines()

                    tile = tile.split(",")

                    for line in lines:
                        line = line.split("\t")

                        line[0] = int(line[0]) - offset_y + (int(tile[1]) - offset_y)
                        line[1] = int(line[1]) - offset_x + (int(tile[0]) - offset_x)

                        # check if cell is in erythrocyte area, write only if no ery. area
                        if not (erythrocyte_map[line[0], line[1]]):
                            outfile.write(str(line[0]) + "\t" + str(line[1]) + "\t" + line[2] + "\t" + line[3])

                outfile.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--function', type=str, default="plot", required=False,
                        help="type 'comp_txt' for generating a composite text from ImmuNets prediction.txt\n")
    parser.add_argument('--x', type=int, default=1908, required=False,
                        help="type x coord of upper left tile in qptiff\n")
    parser.add_argument('--y', type=int, default=2006, required=False,
                        help="type y coord of upper left tile in qptiff\n")
    parser.add_argument('--offset_x', type=int, default=5012, required=False,
                        help="upper left tiles name (x)\n")
    parser.add_argument('--offset_y', type=int, default=26886, required=False,
                        help="upper left tiles name (y)\n")
    parser.add_argument('--path', type=str,
                        default="immunet_modified_pipeline/demo-output/output_tilecache_data2_cluster/output_data2/HALO-N2125_11/",
                        required=False,
                        help="type path for projects output data\n")
    parser.add_argument('--outfile', type=str, default="composite_pred_11_wo_eries.txt", required=False,
                        help="Composite output txt \n")
    parser.add_argument('--binary', type=str, default="binary_N2125_eries.png", required=False,
                        help="Input binary image \n")
    parser.add_argument('--downsample', type=int, default=2, required=False,
                        help="Downsample factor of QuPath binary image export\n")

    args = parser.parse_args()

    function_type = args.function
    path = args.path
    x = args.x
    y = args.y
    offset_x = args.offset_x
    offset_y = args.offset_y
    outfile = args.outfile
    binary = args.binary
    downsample = args.downsample

    if function_type == "comp_txt":
        composite_prediction(path, outfile, binary, x, y, offset_x, offset_y, downsample)
