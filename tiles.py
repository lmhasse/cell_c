import argparse
import random
import numpy as np
import os.path
import tifffile
import csv
import json
import re
from pathlib import Path


def extract_tiff(filename):
    print("Starting to extract page 1 out of qptiff", filename, "...")

    image = tifffile.imread(filename, squeeze=False)
    print(image.shape)

    pattern = image.shape[1:3]
    print(pattern)

    newfile = filename.replace(".qptiff", "_extracted_cropped.tiff")
    # tifffile.imwrite(newfile, tifffile.TiffFile(imagepath).pages[0].asarray(), bigtiff=True, photometric="minisblack")

    # Works with three colored qptiff files
    series = [tifffile.TiffFile(filename).pages[0].asarray(), tifffile.TiffFile(filename).pages[1].asarray(),
              tifffile.TiffFile(filename).pages[2].asarray(), tifffile.TiffFile(filename).pages[3].asarray()]

    with tifffile.TiffWriter(directory + newfile) as tif:
        for s in series:
            tif.write(s, photometric="minisblack")

    image = tifffile.imread(newfile, squeeze=False)
    print(image.shape)

    print("Extraction successful:", newfile, "created.")
    print("-----------------------------------------------")
    return newfile


def make_tiles(filename, crop_side_px=2000, size_px=500, amt_tiles=15, autofl_layer=3):
    """ Generates given amount of crops of given size out of qptiff file and saves them as tiffs without
    autofluorescence

    :param filename: string, name of the file to make tiles of
    :param crop_side_px: int, cuts off the sides of the qptiff to reduce tiles consisting only of background
    :param size_px: int, determines size of the tiles (height and width)
    :param amt_tiles: int, amount of tiles created
    :param autofl_layer: int, index (starting from 1) of autofluorescence layer
    :return: void
    """
    print("Starting to generate", amt_tiles, "tiles from", filename, "...")
    file_path = path_qptiff + filename

    image = tifffile.imread(file_path, squeeze=False)
    print(image.shape)

    if autofl_layer == -1: autofl_layer = image.shape[0]

    step_y = crop_side_px
    step_x = crop_side_px

    tiles = []  # list containing tile coordinates

    while step_y < image.shape[1] - crop_side_px:
        while step_x < image.shape[2] - crop_side_px:
            tiles.append([step_x, step_y])
            step_x += size_px
        step_x = crop_side_px
        step_y += size_px

    rand_tiles = random.choices(tiles, k=amt_tiles)

    series = [] # list containing all the layers of the tiff

    for i in range(image.shape[0]):
        if i == autofl_layer-1:
            continue
        series.append(tifffile.TiffFile(file_path).pages[i].asarray())

    filename = filename.replace(".qptiff", "")

    # Creating tiff file in folder ~/tilecache/test/<scan name>/<xxxx,yyyy>/components.tiff
    for tile in rand_tiles:
        new_folder = directory + "/tilecache/test/" + filename + "/" + ",".join(str(s) for s in tile)
        tile_cache = Path(new_folder)
        tile_cache.mkdir(parents=True, exist_ok=True)

        data = []

        for layer in series:
            data.append(layer[int(tile[1]):int(tile[1] + size_px), int(tile[0]):int(tile[0] + size_px)])

        tifffile.imwrite(new_folder + "/components.tiff", np.stack(tuple(data)), photometric="minisblack")

        print("The tile", new_folder + "/" + "components.tiff", "is saved")

    print("-----------------------------------------------")


def csv_to_json():
    """ Creates json file from several csv files containing info about cells position and type.

    :return: void
    """

    print("Starting to create annotations_train.json from tilecache...")
    slide_names = os.listdir(path_tilecache)

    slide_list = []

    for slide_name in slide_names:  # iterate over all tiles in the ds id folder
        print(slide_name)
        tile_names = os.listdir(path_tilecache + slide_name)
        tile_list = []

        for tile_name in tile_names:  # iterate over all tiles in the scan folder
            csv_list = os.listdir(path_tilecache + slide_name + "/" + tile_name)
            annotations_list = []
            if len(csv_list) <= 1:
                continue

            for csv_name in csv_list:  # iterate over specific csv for each tile

                if csv_name.endswith(".csv"):  # exclude components.tiff

                    with open(path_tilecache + slide_name + "/" + tile_name + "/" + csv_name) as csv_file:
                        csv_reader = csv.reader(csv_file)

                        for row in csv_reader:  # iterate over each row in csv file
                            if row == [' ', 'Label', 'Area', 'Mean', 'Min', 'Max', 'X', 'Y', 'Ch'] or \
                                    row == [' ', 'Label', 'Area', 'Mean', 'Min', 'Max', 'X', 'Y', 'Slice'] or \
                                    row == [' ', 'Area', 'Mean', 'Min', 'Max', 'X', 'Y', 'Ch'] or \
                                    row == [' ', 'Area', 'Mean', 'Min', 'Max', 'X', 'Y', 'Slice']:
                                continue

                            if len(row) == 8:
                                x, y = float(row[5]), float(row[6])
                            else:
                                x, y = float(row[6]), float(row[7])

                            if csv_name.endswith("_prob_pos.csv"):
                                if csv_name.startswith("DAPI"):
                                    annotations_list.append({"type": "Other cell", "x": x, "y": y,
                                                             "positivity": [4, 1, 1]})
                                if csv_name.startswith("CD3"):
                                    annotations_list.append({"type": "T cell", "x": x, "y": y,
                                                             "positivity": [3, 4, 1]})
                                if csv_name.startswith("CD20"):
                                    annotations_list.append({"type": "B cell", "x": x, "y": y,
                                                             "positivity": [3, 1, 4]})

                            elif csv_name.endswith("_pos.csv"):
                                if csv_name.startswith("DAPI"):
                                    annotations_list.append({"type": "Other cell", "x": x, "y": y,
                                                             "positivity": [5, 1, 1]})
                                if csv_name.startswith("CD3"):
                                    annotations_list.append({"type": "T cell", "x": x, "y": y,
                                                             "positivity": [3, 5, 1]})
                                if csv_name.startswith("CD20"):
                                    annotations_list.append({"type": "B cell", "x": x, "y": y,
                                                             "positivity": [3, 1, 5]})

                            elif csv_name.endswith("unclear.csv"):
                                if csv_name.startswith("DAPI"):
                                    annotations_list.append({"type": "Other cell", "x": x, "y": y,
                                                             "positivity": [3, 1, 1]})
                                if csv_name.startswith("CD3"):
                                    annotations_list.append({"type": "T cell", "x": x, "y": y,
                                                             "positivity": [3, 3, 1]})
                                if csv_name.startswith("CD20"):
                                    annotations_list.append({"type": "B cell", "x": x, "y": y,
                                                             "positivity": [3, 1, 3]})

                            elif csv_name.endswith("_prob_neg.csv"):
                                if csv_name.startswith("DAPI"):
                                    annotations_list.append({"type": "Other cell", "x": x, "y": y,
                                                             "positivity": [2, 1, 1]})
                                if csv_name.startswith("CD3"):
                                    annotations_list.append({"type": "T cell", "x": x, "y": y,
                                                             "positivity": [3, 2, 1]})
                                if csv_name.startswith("CD20"):
                                    annotations_list.append({"type": "B cell", "x": x, "y": y,
                                                             "positivity": [3, 1, 2]})

                            elif csv_name.endswith("_neg.csv"):
                                annotations_list.append({"type": "No cell", "x": x, "y": y,
                                                         "positivity": [1, 1, 1]})

            tile_dict = {"tile_id": tile_name, "annotations": annotations_list}
            tile_list.append(tile_dict)

        slide_dict = {"slide_id": slide_name, "tiles": tile_list}
        slide_list.append(slide_dict)

    json_dict = {"ds_id": ds_id_name.replace(".qptiff", ""), "slides": slide_list}

    with open('annotations_train.json', 'w') as outfile:
        json.dump(json_dict, outfile)

    print("annotations_train.json is saved at", directory)
    print("-----------------------------------------------")


def preprocess_qptiff(filepath, filename, autofl_layer):
    print("Preprocessing qptiff", filename, "...")

    image = tifffile.imread(filepath + filename)
    if autofl_layer == -1: autofl_layer = image.shape[0]

    series = [] # list containing all the layers of the tiff

    for i in range(image.shape[0]):
        if i == autofl_layer-1:
            continue
        series.append(tifffile.TiffFile(filepath + "/" + filename).pages[i].asarray())

    data = []

    for layer in series:
        data.append(layer[:, :])

    filename.replace(".qptiff", ".tiff")

    tifffile.imwrite(directory + filename, np.stack(tuple(data)), photometric="minisblack")

    print("The tiff file", filename, "is saved and autofluorescence is excluded")
    print("-----------------------------------------------")


def preprocess_tif(filepath, filename, autofl_layer):
    print("Preprocessing tif", filename, "...")

    image = tifffile.imread(filepath + filename, squeeze=False)
    if autofl_layer == -1: autofl_layer = image.shape[0]

    series = [] # list containing all the layers of the tiff
    print(image.shape)
    for i in range(image.shape[0]):
        if i == autofl_layer-1:
            continue
        series.append(tifffile.TiffFile(filepath + filename).pages[i].asarray())

    data = []

    for layer in series:
        data.append(layer[:, :])

    filename.replace(".tif", ".tiff")

    matchObj = re.search(r'\d+,\d+', filename)

    new_folder = path_tilecache + matchObj.group()
    tile_cache = Path(new_folder)
    tile_cache.mkdir(parents=True, exist_ok=True)

    if matchObj:
        print(matchObj.group())
    else:
        print("oops")

    tifffile.imwrite(path_tilecache + matchObj.group() + "/" + filename, np.stack(tuple(data)), photometric="minisblack")

    print("The tiff file", filename, "is saved and autofluorescence is excluded")
    print("-----------------------------------------------")

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--function', type=str, required=True,
                        help="'tiles' for creating tiles, 'csv' for creating a csv file containing cells from "
                             "tiles or 'qptiff' for preprocessing qptiff for following analysis ")
    parser.add_argument('--path_qptiff', type=str, default=os.getcwd() + "/", required=False,
                        help="directory where the qptiff files are located (default: current working directory)")
    parser.add_argument('--path_tilecache', type=str, default=os.getcwd() + "/tilecache/", required=False,
                        help="directory where the tilecache folder is located (default: current working directory "
                             "+ '/tilecache/test/')")
    parser.add_argument('--amount_tiles', type=int, default=10, required=False,
                        help="amount of tiles created")
    parser.add_argument('--tile_size', type=int, default=500, required=False,
                        help="height and length of tiles in pixel")
    parser.add_argument('--project_name', type=str, default="test", required=False,
                        help="name of the project containing all qptiff files")
    parser.add_argument('--autofl_layer', type=int, default=-1, required=False,
                        help="layer on which autofluorescence is located")
    parser.add_argument('--qptiff', type=str, default="image.qptiff", required=False,
                        help="name of the qptiff file to be preprocessed")

    args = parser.parse_args()

    function = args.function
    path_qptiff = args.path_qptiff
    ds_id_name = args.project_name
    path_tilecache = args.path_tilecache
    amount_tiles = args.amount_tiles
    tile_size = args.tile_size
    autofl_layer = args.autofl_layer
    qptiff = args.qptiff

    directory = os.getcwd()
    filenames = os.listdir(path_qptiff)

    if function == "tiles":
        for file in filenames:
            if file.endswith(".qptiff"):
                make_tiles(file, amt_tiles=amount_tiles, size_px=tile_size, autofl_layer=autofl_layer)

    elif function == "csv":
        csv_to_json()

    elif function == "qptiff":
        preprocess_qptiff(path_qptiff, qptiff, autofl_layer)

    elif function == "tif":
        print(directory)
        for file in filenames:
            if file.endswith("component_data.tif"):
                preprocess_tif(path_qptiff, file, autofl_layer)

    else:
        print("No valid function!")
        print("You entered", function)
        print("Type 'python3 tiles.py --help' for help")
