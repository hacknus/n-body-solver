import numpy as np
import pandas as pd


def read_collision_data(filename="dubinski.tab"):
    # read milkyway andromeda collision dataset
    header = ["m", "x", "y", "z", "vx", "vy", "vz"]
    df = pd.read_csv(filename, sep="\s", names=header, engine="python")

    fac = 1 / 128

    subset = 81920 * fac
    subset /= 2
    f2 = int(subset / 5)
    f1 = 2 * f2

    print(f1, f2)

    gal_disk = list(range(0, f1))
    and_disk = list(range(16384, 16384 + f1))
    gal_bulge = list(range(32768, 32768 + f2))
    and_bulge = list(range(40960, 40960 + f2))
    gal_halo = list(range(49152, 49152 + f1))
    and_halo = list(range(65536, 65536 + f1))

    mask = gal_disk + and_disk + gal_bulge + and_bulge + gal_halo + and_halo
    df = df.loc[mask]
    df["m"] = df.m * 1 / fac
    #df.reset_index(inplace=True)
    print(df.head())
    print(f"saving {len(df)} objects to cdata.csv")
    df.to_csv("../input/milky_way_andromeda.csv")
    return


if __name__ == "__main__":
    read_collision_data()
