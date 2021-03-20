from sys import argv
from scipy.io import FortranFile
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import sys
import os
import pandas as pd


class Planet:

    def __init__(self, df):
        self.x = [float(df["x"])]
        self.y = [float(df["y"])]
        self.z = [float(df["z"])]
        self.m = float(df["m"])
        self.color = self.get_color()
        self.color = "yellow"

    def add(self, df):
        self.x.append(float(df["x"]))
        self.y.append(float(df["y"]))
        self.z.append(float(df["z"]))

    def get_color(self):
        if self.y[0] > 0:
            return "red"
        else:
            return "yellow"


def read_binary(path):
    # path = "out_10.bin"
    x = np.fromfile(path, dtype=np.float64)

    x = x.reshape((x.shape[0] // 9, 9))

    d = {"m": x[:, 0],
         "x": x[:, 1],
         "y": x[:, 2],
         "z": x[:, 3],
         "vx": x[:, 4],
         "vy": x[:, 5],
         "vz": x[:, 6],
         }
    df = pd.DataFrame(data=d)
    return df


def T(a, i, e):
    au = 1.5e11
    a_j = 5.204 * au
    return a_j / a + 2 * np.cos(i) * np.sqrt(a / a_j * (1 - e ^ 2))


planets = []
j = 0

while True:
    if not os.path.exists(f'output/out_{j:05d}.dat'):
        break
    print(f'reading out_{j:05d}.bin')
    master_file = f'output/out_{j:05d}.dat'

    if j % 100 == 0:
        df = read_binary(master_file)
        if j == 0:
            for i in range(len(df)):
                planets.append(Planet(df.loc[i]))
        else:
            for i in range(len(df)):
                planets[i].add(df.loc[i])
    j += 1
if j == 0:
    print("no files found, exiting...")
    exit()
else:
    print("found {} files".format(j))

body_i = 30
plt.plot(planets[body_i].x, planets[body_i].y)
plt.axis("equal")
plt.show()
