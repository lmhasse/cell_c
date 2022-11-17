import argparse
import math
import random
import numpy as np
import os.path
import tifffile
import csv
import json
import re
from pathlib import Path
from tqdm import tqdm


def merge_cells(a_list, allowed_distance):
    i = 0
    while i < len(a_list)-1:
        if math.sqrt((a_list[i]["x"] - a_list[i + 1]["x"]) ** 2
                     + (a_list[i]["y"] - a_list[i + 1]["y"]) ** 2) <= allowed_distance:
            print("----------")
            print(math.sqrt((a_list[i]["x"] - a_list[i + 1]["x"]) ** 2 + (a_list[i]["y"] - a_list[i + 1]["y"]) ** 2))
            print(a_list[i], a_list[i + 1])
            print("----------")

            if (a_list[i]["type"] in ["T cell", "B cell", "Macrophage", "Dendritic cell"] and
                a_list[i + 1]["type"] in ["Other cell"]) or \
                    (a_list[i + 1]["type"] in ["T cell", "B cell", "Macrophage", "Dendritic cell"] and
                     a_list[i]["type"] in ["Other cell"]):
                # Cell type with PD-1

                a_list[i] = {"type": (a_list[i]["type"] if a_list[i]["type"] != "Other cell" else a_list[i + 1]),
                             "x": (a_list[i]["x"] + a_list[i + 1]["x"]) // 2,
                             "y": (a_list[i]["y"] + a_list[i + 1]["y"]) // 2,
                             "positiviy": list((a_list[i]["positivity"][j]
                                               if a_list[i]["positivity"][j] >= a_list[i + 1]["positivity"][j]
                                               else a_list[i + 1]["positivity"][j])
                                               for j in range(7))}

                print(a_list[i], a_list[i+1])

                a_list.remove(a_list[i + 1])

                print(a_list[i], "case 1")

                continue

            elif (a_list[i]["type"] == "Macrophage" and a_list[i + 1]["type"] == "Dendritic cell") or \
                    (a_list[i + 1]["type"] == "Macrophage" and a_list[i]["type"] == "Dendritic cell"):
                # Macrophage

                a_list[i] = {"type": (a_list[i]["type"] if a_list[i]["type"] == "Macrophage" else a_list[i + 1]),
                             "x": a_list[i]["x"] + a_list[i + 1]["x"] // 2,
                             "y": a_list[i]["y"] + a_list[i + 1]["y"] // 2,
                             "positiviy": list(a_list[i]["positivity"][j]
                                               if a_list[i]["positivity"][j] > a_list[i + 1]["positivity"][j]
                                               else a_list[i + 1]["positivity"][j]
                                               for j in range(6))}

                a_list.remove(a_list[i + 1])

                print(a_list[i], "case 2")

                continue

        i += 1

    return a_list


def expressions_to_list(row):
    list_expressions = [4,
                        5 if row[3] == "vWF" else (2 if row[3] == "vWF_prob_neg" else 1),
                        5 if row[3] == "PD-1" else (2 if row[3] == "PD-1_prob_neg" else 1),
                        5 if row[3] == "CD20" else (2 if row[3] == "CD20_prob_neg" else 1),
                        5 if row[3] == "CD163" else (2 if row[3] == "CD163_prob_neg" else 1),
                        5 if row[3] == "Lamp3" else (2 if row[3] == "Lamp3_prob_neg" else 1),
                        5 if row[3] == "CD3" else (2 if row[3] == "CD3_prob_neg" else 1)]

    type_cell = "Vessel cell" if list_expressions[1] != 1 else (
        "B cell" if list_expressions[3] != 1 else (
            "Macrophage" if list_expressions[4] != 1 else (
                "Dendritic cell" if list_expressions[5] != 1 else (
                    "T cell" if list_expressions[6] != 1 else
                        "Other cell"  # cells with only DAPI or DAPI + PD-1
                ))))

    return {"type": type_cell, "x": round(float(row[0])), "y": round(float(row[1])), "positivity": list_expressions}


"""
                    "type": "Other cell", "x": x, "y": y,
                    "positivity": [4, 1, 1]})
"""


def tsv_to_json():
    """ Creates json file from several tsv files containing info about cells position and type.
    Folder structure needed:
    annotations_qupath folder
    └── slide 1 folder
        ├── tile_[1,1].tsv - [x, y, class, name, color]
        ├── tile_[2,2].tsv - [x, y, class, name, color]
        ├── tile_[3,3].tsv - [x, y, class, name, color]
        ...

    :return: void
    """

    print("Starting to create annotations_train.json from tilecache...")

    slide_names = os.listdir(path_annotations)
    slide_list = []

    for slide_name in tqdm(slide_names):  # iterate over all slide folders in path_annotations
        tile_names = os.listdir(path_annotations + slide_name)
        tile_list = []

        for tile_name in tqdm(tile_names):  # iterate over all tiles tsv in the scan's folder

            tile_coords = tile_name.split("]")[0].split("[")[1]  # parsing tile coordinates out of file name

            annotations_list = []

            with open(path_annotations + slide_name + "/" + tile_name) as tsv_file:
                print(tile_name)
                tsv_reader = csv.reader(tsv_file, delimiter="\t")

                for row in tsv_reader:  # iterate over each row in tsv file
                    if row[0] == "x":
                        continue

                    annotations_list.append(expressions_to_list(row))

            tile_dict = {"tile_id": tile_coords, "annotations": merge_cells(annotations_list, 5)}
            tile_list.append(tile_dict)

        slide_dict = {"slide_id": slide_name, "tiles": tile_list}
        slide_list.append(slide_dict)

    json_dict = {"ds_id": ds_id_name, "slides": slide_list}

    with open('annotations_train.json', 'w') as outfile:
        json.dump(json_dict, outfile)

    print("annotations_train.json is saved at", directory)
    print("-----------------------------------------------")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--function', type=str, required=False, default="tsv",
                        help="'tiles' for creating tiles, 'csv' for creating a csv file containing cells from "
                             "tiles or 'qptiff' for preprocessing qptiff for following analysis ")
    parser.add_argument('--path_annotations', type=str, default="/home/agimkeller/mnt_AG_Imkeller/group_members/hasse/immunet_data/annotations_qupath/", required=False,
                        help="directory where the annotations folder is located (default: current working directory "
                             "+ '/tilecache/test/')")
    parser.add_argument('--project_name', type=str, default="test", required=False,
                        help="name of the project containing all qptiff files")

    args = parser.parse_args()

    function = args.function
    path_annotations = args.path_annotations
    ds_id_name = args.project_name
    directory = os.getcwd()

    if function == "tsv":
        tsv_to_json()

    else:
        print("No valid function!")
        print("You entered", function)
        print("Type 'python3 tiles.py --help' for help")
