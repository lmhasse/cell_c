import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import tifffile
import csv


def plot_predictions(txt_path, img_name, threshold):
    coordinates = []

    with open(txt_path) as f:
        lines = f.readlines()
        for row in lines:
            if row == "\n" or row.startswith("tile:"):
                continue
            row = row.split("\t")
            coordinates.append([int(row[1]), int(row[0]), float(row[2]), float(row[3].strip())])

    x_CD3, y_CD3 = [], []
    x_CD20, y_CD20 = [], []

    print("All cells:", len(coordinates))

    for cord in coordinates:
        if cord[2] > 0 and cord[3] < cord[2] / threshold:
            # T Cell
            x_CD3.append(cord[0])
            y_CD3.append(cord[1])
        elif cord[2] < cord[3] / threshold and cord[3] > 0:
            # B Cell
            x_CD20.append(cord[0])
            y_CD20.append(cord[1])

    print("CD3", len(x_CD3), "CD20", len(x_CD20))

    x_CD3 = np.array(x_CD3)
    y_CD3 = np.array(y_CD3)
    x_CD20 = np.array(x_CD20)
    y_CD20 = np.array(y_CD20)

    fig = plt.figure(figsize=(10, 12), dpi=1200)
    ax = fig.add_subplot()

    ax.text(8, 8, 'All cells: ' + str(len(coordinates)) + ', CD3: ' + str(len(x_CD3)) + ', CD20: ' + str(len(x_CD20)),
            bbox={'facecolor': 'lightblue', 'alpha': 0.5, 'pad': 10})

    plt.scatter(x_CD3, y_CD3, s=0.01, marker=",", color="limegreen", edgecolors=None)
    plt.scatter(x_CD20, y_CD20, s=0.01, marker=",", color="r", edgecolors=None)

    # plt.imshow(image, cmap=plt.get_cmap('binary'))
    # plt.show()
    plt.savefig(img_name, dpi=fig.dpi)


def plot_csv(csv_path, img_name):
    coordinates = {"CD20": [], "CD3": []}

    with open(csv_path) as csv_file:
        csv_reader = csv.reader(csv_file)
        line_count = 0

        for row in csv_reader:
            if line_count == 0:
                # exclude header
                line_count += 1
            else:
                line_count += 1

                cell_coordinates = [((int(row[1]) + int(row[2])) / 2), ((int(row[3]) + int(row[4])) / 2)]

                if row[9] == "1":
                    # CD20
                    coordinates["CD20"].append(cell_coordinates)
                elif row[13] == "1":
                    # CD3
                    coordinates["CD3"].append(cell_coordinates)

    fig = plt.figure(figsize=(10, 12), dpi=1200)

    x_CD3 = [i[0] for i in coordinates["CD3"]]
    y_CD3 = [i[1] for i in coordinates["CD3"]]
    x_CD20 = [i[0] for i in coordinates["CD20"]]
    y_CD20 = [i[1] for i in coordinates["CD20"]]

    plt.scatter(x_CD3, y_CD3, s=0.01, marker=",", color="limegreen")
    plt.scatter(x_CD20, y_CD20, s=0.01, marker=",", color="r")

    # plt.imshow(image, cmap=plt.get_cmap('binary'))
    # plt.show()
    plt.savefig(img_name, dpi=fig.dpi)


if __name__ == '__main__':

    plot_predictions("compare_output/composites/composite_pred_rad11_wo_eries.txt", "plot_output/radius_11/immunet_rad11_thr60_3_20_wo_eries.png", 60)
    # immunet_radxx_thrxx_3_20.png --- rad - radius size --- thr - threshold --- 20_3 - CD20 dann CD 3 oder umgekehrt
    # plot_csv("N2125-14-1.csv", "csv_img_less_resolution_resized.png")
    