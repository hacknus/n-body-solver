import numpy as np
import matplotlib.pyplot as plt
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
    return a_j / a + 2 * np.cos(i) * np.sqrt(a / a_j * (1 - e ** 2))


def ellipse(x, a, b):
    return b * np.sqrt(1 - x / a ** 2)


planets = []
j = 0

while True:
    if j % 10 == 0:
        if not os.path.exists(f'../output/out_{j:07d}.dat'):
            break
        print(f'reading out_{j:07d}.dat')
        master_file = f'../output/out_{j:07d}.dat'
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

comets = range(10, len(planets))
num_orbits = 2
j_orbit = 365 * 12 // 10 * num_orbits
orbits = range(0, len(planets[0].x), j_orbit)
tisserand = np.zeros((len(comets), len(orbits)))
ind = 0
for body_i in comets:
    x = np.array(planets[body_i].x[:j_orbit])
    y = np.array(planets[body_i].y[:j_orbit])
    plt.plot(x, y)
plt.show()

for body_i in comets:
    ind_orb = 0
    for i in orbits:
        x = np.array(planets[body_i].x[i:i + j_orbit])
        y = np.array(planets[body_i].y[i:i + j_orbit])
        r = np.sqrt(x ** 2 + y ** 2)
        a = (np.min(r) + np.max(r)) / 2
        b = np.sqrt(np.min(r) * np.max(r))
        e = np.sqrt(1 - (b / a) ** 2)
        print(f"a: {a}, b: {b}, e: {e:.4f}")
        tisserand[ind, ind_orb] = T(a, 0, e)
        ind_orb += 1
    ind += 1
plt.plot(np.arange(len(orbits)) * num_orbits, tisserand.T)
plt.plot(np.arange(len(orbits)) * num_orbits, np.mean(tisserand, axis=0), color="black", ls="--")
print(f"T_mean = {np.mean(tisserand[:, -1], axis=0):.2f} +/- {np.std(tisserand[:, -1], axis=0):.2f}")
plt.axhline(2, color="red", ls="--")
plt.axhline(3, color="red", ls="--")
plt.ylim(0, 5)
plt.xlabel(r"$n$ Jupiter Orbits")
plt.ylabel(r"$T_{Jupit  er}$")
plt.savefig("comets.pdf")
plt.clf()
plt.plot(np.arange(len(orbits)) * num_orbits, np.mean(tisserand, axis=0), color="black", ls="--")
plt.axhline(2, color="red", ls="--")
plt.axhline(3, color="red", ls="--")
plt.ylim(0, 5)
plt.xlabel(r"$n$ Jupiter Orbits")
plt.ylabel(r"$T_{Jupiter}$")
plt.savefig("comets_mean.pdf")
