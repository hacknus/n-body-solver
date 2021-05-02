import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from scipy.optimize import least_squares, curve_fit, leastsq
from python.kepler import xyzele

au = 1.5e11

mode = "_rev"


# mode = ""


class Planet:

    def __init__(self, df):
        self.x = [float(df["x"])]
        self.y = [float(df["y"])]
        self.z = [float(df["z"])]
        self.vx = [float(df["vx"])]
        self.vy = [float(df["vy"])]
        self.vz = [float(df["vz"])]
        self.ax = [float(df["ax"])]
        self.ay = [float(df["ay"])]
        self.az = [float(df["az"])]
        self.dt = [float(df["dt"])]
        self.t = [float(df["t"])]
        self.ekin = [float(df["ekin"])]
        self.epot = [float(df["epot"])]
        self.m = float(df["m"])
        self.color = self.get_color()
        self.color = "yellow"

    def add(self, df):
        self.x.append(float(df["x"]))
        self.y.append(float(df["y"]))
        self.z.append(float(df["z"]))
        self.vx.append(float(df["vx"]))
        self.vy.append(float(df["vy"]))
        self.vz.append(float(df["vz"]))
        self.ax.append(float(df["ax"]))
        self.ay.append(float(df["ay"]))
        self.az.append(float(df["az"]))
        self.dt.append(float(df["dt"]))
        self.t.append(float(df["t"]))
        self.ekin.append(float(df["ekin"]))
        self.epot.append(float(df["epot"]))

    def get_color(self):
        if self.y[0] > 0:
            return "red"
        else:
            return "yellow"


def read_binary(path):
    # path = "out_10.bin"
    x = np.fromfile(path, dtype=np.float64)

    x = x.reshape((x.shape[0] // 14, 14))

    d = {"m": x[:, 0],
         "x": x[:, 1],
         "y": x[:, 2],
         "z": x[:, 3],
         "vx": x[:, 4],
         "vy": x[:, 5],
         "vz": x[:, 6],
         "ax": x[:, 7],
         "ay": x[:, 8],
         "az": x[:, 9],
         "dt": x[:, 10],
         "t": x[:, 11],
         "epot": x[:, 12],
         "ekin": x[:, 13],
         }
    df = pd.DataFrame(data=d)
    return df


def T(a, i, e):
    a_j = 5.204 * au
    return a_j / a + 2 * np.cos(i) * np.sqrt(a / a_j * (1 - e ** 2))


def ellipse(x, a, b):
    return b * np.sqrt(1 - x / a ** 2)


planets = []
j = 0
save_interval = 1000000
increment = 1 * save_interval

while True:
    try:
        if j % increment == 0:
            if not os.path.exists(f'../output{mode}/out_{j:09d}.dat'):
                break
            print(f'reading out_{j:09d}.dat')
            master_file = f'../output{mode}/out_{j:09d}.dat'
            df = read_binary(master_file)
            if j == 0:
                for i in range(len(df)):
                    planets.append(Planet(df.loc[i]))
            else:
                for i in range(len(df)):
                    planets[i].add(df.loc[i])
        j += 1
        # if j >= 365 * 12 * save_interval:
        #     break
    except KeyboardInterrupt:
        break

if j == 0:
    print("no files found, exiting...")
    exit()
else:
    print("found {} files".format(j))

a = np.zeros((len(planets), len(planets[0].t)))
for i in range(len(planets)):
    a[i, :] = np.sqrt(np.array(planets[i].ax) ** 2 + np.array(planets[i].ay) ** 2 + np.array(planets[i].az) ** 2)
a = np.max(a, axis=0)

plt.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), planets[0].dt)
plt.show()
print(f"t_end: {planets[0].t[-1] / (3600 * 365 * 12 * 24)}")

E = np.zeros((len(planets), len(planets[0].t)))
ekin = np.zeros((len(planets), len(planets[0].t)))
epot = np.zeros((len(planets), len(planets[0].t)))
for i in range(len(planets)):
    ekin[i, :] = np.array(planets[i].ekin)
    epot[i, :] = np.array(planets[i].epot)

E = epot + ekin
E = E[9:, :]
print(E.shape, epot.shape, ekin.shape)

plt.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(E, axis=0))
#plt.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(ekin, axis=0), ls="-")
#plt.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(epot, axis=0), ls="-")
plt.ylabel(r"$E$")
plt.xlabel(r"$n$ Jupiter Orbits")
plt.savefig(f"full_comets_energy_{mode}.png")
plt.show()

E = np.mean(E, axis=0)
dE = np.diff(E)

plt.plot(np.array(planets[0].t[1:]) / (3600 * 365 * 12 * 24), dE/E[1:])
plt.ylabel(r"$dE/E$")
plt.xlabel(r"$n$ Jupiter Orbits")
plt.savefig(f"full_comets_delta_energy_{mode}.png")
plt.show()
print(f"t_end: {planets[0].t[-1] / (3600 * 365 * 12 * 24)}")

comets = np.array(range(10, len(planets)))
n_per_orbit = 365 * 12 * 24 // increment
ind = 0
n_data = len(planets[0].x)
tisserand = np.zeros((len(comets), n_data))
semi_major_a = np.zeros((len(comets), n_data))
semi_minor_a = np.zeros((len(comets), n_data))
es = np.zeros((len(comets), n_data))
inclinations = np.zeros((len(comets), n_data))
rs = np.zeros((len(comets), n_data))
vs = np.zeros((len(comets), n_data))
disabled = np.zeros((len(comets)))

for body_i in comets:
    for j in range(n_data):
        x = np.array(planets[body_i].x[j]) - planets[0].x[j]
        y = np.array(planets[body_i].y[j]) - planets[0].y[j]
        z = np.array(planets[body_i].z[j]) - planets[0].z[j]
        vx = np.array(planets[body_i].vx[j]) - planets[0].vx[j]
        vy = np.array(planets[body_i].vy[j]) - planets[0].vy[j]
        vz = np.array(planets[body_i].vz[j]) - planets[0].vz[j]
        r = np.array([x, y, z])
        v = np.array([vx, vy, vz])
        GM = 1.32712440018e20  # sun
        # ele = xyzele(GM, r, v)

        h = np.cross(r, v)
        eps = np.dot(v, v) / 2 - GM / np.linalg.norm(r)
        e = np.sqrt(1 + 2 * eps * np.dot(h, h) / GM ** 2)
        a = np.dot(h, h) / (GM * (1 - e ** 2))
        b = a * np.sqrt(1 - e ** 2)
        i = np.arctan2(np.sqrt(h[0] ** 2 + h[1] ** 2), h[2])
        if e >= 1 or np.linalg.norm(v) > np.sqrt(2 * GM / np.linalg.norm(r)) or disabled[ind] == 1:
            disabled[ind] = 1
            tisserand[ind, j] = tisserand[ind, j - 1]
            semi_major_a[ind, j] = semi_major_a[ind, j - 1]
            semi_minor_a[ind, j] = semi_minor_a[ind, j - 1]
            es[ind, j] = es[ind, j - 1]
            inclinations[ind, j] = inclinations[ind, j - 1]
            rs[ind, j] = rs[ind, j - 1]
            vs[ind, j] = vs[ind, j - 1]
        else:
            tisserand[ind, j] = T(a, i, e)
            semi_major_a[ind, j] = a
            semi_minor_a[ind, j] = b
            es[ind, j] = e
            inclinations[ind, j] = i
            rs[ind, j] = np.sqrt(x ** 2 + y ** 2 + z ** 2)
            vs[ind, j] = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    ind += 1

n_per_orbit_plot = 365 * 12 * 24 // increment * 100

comets = comets[disabled == 0]

for body_i in comets:
    x = np.array(planets[body_i].x[:n_per_orbit_plot]) - planets[0].x[:n_per_orbit_plot]
    y = np.array(planets[body_i].y[:n_per_orbit_plot]) - planets[0].y[:n_per_orbit_plot]
    plt.plot(x / au, y / au, color="yellow")
    break
for body_i in comets:
    x = np.array(planets[body_i].x[-n_per_orbit_plot:]) - planets[0].x[-n_per_orbit_plot:]
    y = np.array(planets[body_i].y[-n_per_orbit_plot:]) - planets[0].y[-n_per_orbit_plot:]
    plt.plot(x / au, y / au, color="red")
    break

x = np.array(planets[5].x[-n_per_orbit_plot:]) - planets[0].x[-n_per_orbit_plot:]
y = np.array(planets[5].y[-n_per_orbit_plot:]) - planets[0].y[-n_per_orbit_plot:]
plt.plot(x / au, y / au, color="blue")
x = np.array(planets[5].x[:n_per_orbit_plot]) - planets[0].x[:n_per_orbit_plot]
y = np.array(planets[5].y[:n_per_orbit_plot]) - planets[0].y[:n_per_orbit_plot]
plt.plot(x / au, y / au, color="brown")

plt.savefig(f"full_comets_orbits{mode}.png")
plt.show()
plt.clf()

tisserand = tisserand[disabled == 0, :]
es = es[disabled == 0, :]
semi_major_a = semi_major_a[disabled == 0, :]
semi_minor_a = semi_minor_a[disabled == 0, :]

semi_minor_a /= au
semi_major_a /= au

fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True)

ax0.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), tisserand.T)
ax0.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(tisserand, axis=0), color="black", ls="--")
print(f"T_mean = {np.mean(tisserand[:, -1], axis=0):.2f} +/- {np.std(tisserand[:, -1], axis=0):.2f}")
ax0.axhline(2, color="red", ls="--")
ax0.axhline(3, color="red", ls="--")
ax0.set_ylim(1, 4)
ax0.set_ylabel(r"$T_{Jupiter}$")

ax1.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), semi_major_a.T)
ax1.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(semi_major_a, axis=0), color="black", ls="--")
ax1.set_ylabel(r"$a$ [a.u.]")

ax2.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), semi_minor_a.T)
ax2.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(semi_minor_a, axis=0), color="black", ls="--")
ax2.set_ylabel(r"$b$ [a.u.]")

ax3.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), es.T)
ax3.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(es, axis=0), color="black", ls="--")
ax3.set_ylabel(r"$e$")

ax3.set_xlabel(r"$n$ Jupiter Orbits")
plt.savefig(f"full_comets{mode}.png")
plt.show()
plt.clf()

fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=4, sharex=True)

ax0.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(tisserand, axis=0), color="black", ls="--")
ax0.set_ylabel(r"$T_{Jupiter}$")
ax0.axhline(2, color="red", ls="--")
ax0.axhline(3, color="red", ls="--")
ax0.set_ylim(1, 4)

ax1.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(semi_major_a, axis=0), color="black", ls="--")
ax1.set_ylabel(r"$a$ [a.u.]")

ax2.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(semi_minor_a, axis=0), color="black", ls="--")
ax2.set_ylabel(r"$b$ [a.u.]")

ax3.plot(np.array(planets[0].t) / (3600 * 365 * 12 * 24), np.mean(es, axis=0), color="black", ls="--")
ax3.set_ylabel(r"$e$")
ax3.set_xlabel(r"$n$ Jupiter Orbits")
plt.savefig(f"full_comets_mean{mode}.png")
plt.show()
