# This file is part of DXMClib.

# DXMClib is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# DXMClib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

# Copyright 2024 Erlend Andersen


import seaborn as sns
import pandas as pd
from matplotlib import pylab as plt
import numpy as np
from pathlib import Path
import os
import sys
import re

HUE_ORDER = list(["NoneLC", "Livermore", "IA", "TG195"])


def readData(filename="validationTable.txt"):
    converters = {"Result": float, "Stddev": float, "SimulationTime": float}
    dt = pd.read_csv(filename, sep=", ", engine="python", converters=converters)
    nNaN = dt.isnull().sum().sum()
    if nNaN > 0:
        print(
            "Warning: Found {} NaN values or missing data, dropping rows".format(nNaN)
        )
        dt = dt.dropna()

    # Making TG195 results an own model
    m_series = dt["Model"].value_counts().sort_values(ascending=False)
    model = m_series.index[0]
    dm = dt[dt["Model"] == model].copy()
    dm["Result"] = dm["TG195Result"]
    dm["Model"] = ["TG195" for _ in dm["TG195Result"]]
    for key in ["Stddev", "nEvents", "SimulationTime"]:
        dm[key] = [0 for _ in dm[key]]

    dm_min = dt[dt["Model"] == model].copy()
    dm_min["Result"] = dm["TG195Result_min"]
    dm_min["Model"] = ["TG195" for _ in dm_min["TG195Result_min"]]
    for key in ["Stddev", "nEvents", "SimulationTime"]:
        dm_min[key] = [0 for _ in dm_min[key]]

    dm_max = dt[dt["Model"] == model].copy()
    dm_max["Result"] = dm_max["TG195Result"]
    dm_max["Model"] = ["TG195" for _ in dm_max["TG195Result"]]
    for key in ["Stddev", "nEvents", "SimulationTime"]:
        dm_max[key] = [0 for _ in dm_max[key]]

    df = pd.concat([dt, dm, dm_min, dm_max], ignore_index=True)

    present_models = list([m for m in df["Model"].unique()])

    cmodels = list([m for m in HUE_ORDER])
    for m in cmodels:
        if m not in present_models:
            HUE_ORDER.remove(m)

    # adding errors for seaborn
    df_max = df[df["Model"] != "TG195"][df["Case"] != "Case 4.2"].copy()
    df_max["Result"] += 1.96 * df_max["Stddev"]
    df_min = df[df["Model"] != "TG195"].copy()
    df_min["Result"] -= 1.96 * df_max["Stddev"]
    return pd.concat([df, df_min, df_max], ignore_index=True)


def readAllData(filelist: list):
    dfs = []
    for file in filelist:
        dfs.append(readData(file))
    return pd.concat(dfs, ignore_index=True)


def fix_axis(fg, rotate_labels=True, ylabel="Energy [eV/history]", set_y0=False):
    if set_y0:
        fg.set(ylim=(0, None))
    fg.set_axis_labels(y_var=ylabel)
    if rotate_labels:
        for ax in fg.axes_dict.values():
            for tick in ax.get_xticklabels():
                tick.set_rotation(70)
    fg.tight_layout()


def plotCase1(dt_full, kind="strip", show=False):
    dt = dt_full[dt_full["Case"] == "Case 1"]
    if dt.size == 0:
        return

    labels = ["NoneFilter", "HVLFilter", "QVLFilter"]

    for m in ["monoenergetic", "polyenergetic"]:
        dtt = dt[dt["Mode"] == m]
        g = sns.catplot(
            x="Volume",
            y="Result",
            hue="Model",
            col="Specter",
            order=labels,
            hue_order=HUE_ORDER,
            data=dtt,
            kind=kind,
            errorbar=lambda x: (x.min(), x.max()),
        )

        fix_axis(g, ylabel="KERMA per history")

        plt.savefig("plots/Case1_{}.png".format(m), dpi=300)
        if show:
            plt.show()
        plt.clf()
        plt.close()


def plotCase2(dt_full, kind="strip", show=False):
    dt = dt_full[dt_full["Case"] == "Case 2"]
    if dt.size == 0:
        return
    dt_vol = dt[dt["Volume"] != "Total body"]

    labels = list(set([v for v in dt_vol["Volume"]]))
    labels.sort()

    g = sns.catplot(
        x="Volume",
        y="Result",
        hue="Model",
        row="Specter",
        col="Mode",
        order=labels,
        hue_order=HUE_ORDER,
        data=dt_vol,
        kind=kind,
        errorbar=lambda x: (x.min(), x.max()),
    )

    fix_axis(g)

    plt.savefig("plots/Case2_volume.png", dpi=300)
    if show:
        plt.show()
    plt.clf()
    plt.close()

    dt_tot = dt[dt["Volume"] == "Total body"]
    g = sns.catplot(
        x="Specter",
        y="Result",
        hue="Model",
        col="Mode",
        hue_order=HUE_ORDER,
        data=dt_tot,
        kind=kind,
        errorbar=lambda x: (x.min(), x.max()),
    )
    fix_axis(g, False)
    plt.savefig("plots/Case2_total_body.png", dpi=300)
    if show:
        plt.show()
    plt.clf()
    plt.close()


def plotCase3(dt_full, kind="strip", show=False):
    dt = dt_full[dt_full["Case"] == "Case 3"]
    if dt.size == 0:
        return
    dt_vol = dt[dt["Volume"] != "Total body"]

    labels = list(set([v for v in dt_vol["Volume"]]))
    labels.sort()

    g = sns.catplot(
        x="Volume",
        y="Result",
        hue="Model",
        row="Specter",
        col="Mode",
        order=labels,
        hue_order=HUE_ORDER,
        data=dt_vol,
        kind=kind,
        errorbar=lambda x: (x.min(), x.max()),
    )

    fix_axis(g)

    plt.savefig("plots/Case3_volume.png", dpi=300)
    if show:
        plt.show()
    plt.clf()
    plt.close()

    dt_tot = dt[dt["Volume"] == "Total body"]
    labels = list(set([v for v in dt_tot["Specter"]]))
    labels.sort()

    g = sns.catplot(
        x="Specter",
        y="Result",
        hue="Model",
        col="Mode",
        order=labels,
        hue_order=HUE_ORDER,
        data=dt_tot,
        kind=kind,
        errorbar=lambda x: (x.min(), x.max()),
    )
    fix_axis(g, False)
    plt.savefig("plots/Case3_total_body.png", dpi=300)
    if show:
        plt.show()
    plt.clf()
    plt.close()


def plotCase41(dt_full, kind="strip", show=False):
    dt = dt_full[dt_full["Case"] == "Case 4.1"]
    if dt.size == 0:
        return

    labels = list(set([v for v in dt["Volume"]]))
    labels.sort()

    g = sns.catplot(
        x="Volume",
        y="Result",
        hue="Model",
        row="Specter",
        col="Mode",
        order=labels,
        hue_order=HUE_ORDER,
        data=dt,
        kind=kind,
        errorbar=lambda x: (x.min(), x.max()),
    )
    fix_axis(g, False)
    plt.savefig("plots/Case41.png", dpi=300)
    if show:
        plt.show()
    plt.close()


def plotCase42(dt_full, show=False):
    dt = dt_full[dt_full["Case"] == "Case 4.2"]
    if dt.size == 0:
        return
    dt["Volume [angle]"] = [float(d) for d in dt["Volume"]]
    dtp_ind = ["Cent" in e for e in dt["Mode"]]
    dtp = dt[dtp_ind]
    g = sns.relplot(
        x="Volume [angle]",
        y="Result",
        hue="Model",
        row="Specter",
        col="Mode",
        hue_order=HUE_ORDER,
        data=dtp,
    )
    fix_axis(g)
    plt.savefig("plots/Case42_Cent.png", dpi=300)
    if show:
        plt.show()

    dtp_ind = ["Pher" in e for e in dt["Mode"]]
    dtp = dt[dtp_ind]
    g = sns.relplot(
        x="Volume [angle]",
        y="Result",
        hue="Model",
        row="Specter",
        col="Mode",
        hue_order=HUE_ORDER,
        data=dtp,
    )
    fix_axis(g)
    plt.savefig("plots/Case42_Pher.png", dpi=300)
    if show:
        plt.show()
    plt.close()


def plotCase5(dt_full, kind="strip", show=False):
    dt = dt_full[dt_full["Case"] == "Case 5"]
    if dt.size == 0:
        return

    dt["Mode"] = [int(v) for v in dt["Mode"]]
    labels = list(set(dt["Mode"]))
    labels.sort()

    vols = set([m for m in dt["Volume"]])

    for m in vols:
        dtm = dt[dt["Volume"] == m]
        g = sns.catplot(
            x="Mode",
            y="Result",
            hue="Model",
            row=None,
            col="Specter",
            order=labels,
            hue_order=HUE_ORDER,
            data=dtm,
            kind=kind,
            errorbar=lambda x: (x.min(), x.max()),
        )
        fix_axis(g)
        plt.savefig("plots/Case5_{}.png".format(m), dpi=300)
        if show:
            plt.show()
        plt.close()


def plotRuntimes(dt, kind="strip", show=False):
    df = dt[dt["Model"] != "TG195"]
    for cas in set(df["Case"]):
        dff = df[df["Case"] == cas]
        labels = list(set(dff["Mode"]))
        if all([m.isdigit() for m in labels]):
            labels.sort(key=lambda x: int(x))
        else:
            labels.sort()
        g = sns.catplot(
            x="Mode",
            y="SimulationTime",
            row="Specter",
            hue="Model",
            data=dff,
            kind=kind,
            order=labels,
            hue_order=HUE_ORDER[:-1],
        )
        fix_axis(g, ylabel="Simulation time [ms]", set_y0=False)
        plt.savefig("plots/Runtimes_{}.png".format(cas), dpi=300)
        if show:
            plt.show()
        plt.clf()
        plt.close()


def merge_files(filenames: list, output="validationTable.txt", aggregate=True):
    have_header = os.path.isfile(output)
    if not aggregate:
        with open(output, "a") as out:
            for fname in filenames:
                if os.path.isfile(fname) and fname != output:
                    with open(fname, "r") as f:
                        if have_header:
                            lines = f.readlines()
                            if len(lines) > 1:
                                out.writelines(lines[1:])
                        else:
                            out.write(f.read())
                            have_header = True
                else:
                    print("No file named {}, skipping".format(fname))
    else:
        tables = list()
        for fname in filenames:
            if os.path.isfile(fname):
                converters = {"Result": float, "Stddev": float, "SimulationTime": float}
                dt = pd.read_csv(
                    fname, sep=", ", engine="python", converters=converters
                )
                tables.append(dt)
            else:
                print("No file named {}, skipping".format(fname))

        if len(tables) == 0:
            print("No files found, exiting.")
            return

        df = tables[0]
        for i in range(1, len(tables)):
            df["Result"] += tables[i]["Result"]
        df["Result"] /= len(tables)
        df.to_csv(output, index=False, sep=",", mode="w")

        ##fixing seps
        with open(output, "r") as inn:
            txt = inn.read()
        txt.replace(",", ", ")
        with open(output, "w") as out:
            out.write(txt)


def expand_filenames(filenames: list):
    expanded = list()
    files_in_path = os.listdir(os.path.dirname(os.path.abspath(__file__)))
    for fname in filenames:
        if "*" in fname:
            for p in Path(".").glob(fname):
                expanded.append(str(p))
        else:
            expanded.append(fname)

    return list(set(expanded))


if __name__ == "__main__":
    # Setting current path to this file folder
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    infile = "validationTable.txt"
    if len(sys.argv) > 2:
        infile = sys.argv[1]
        merge_files(expand_filenames(sys.argv[2:]), output=infile, aggregate=False)
        # dt = readAllData(expand_filenames(sys.argv[1:]))
    elif len(sys.argv) == 2:
        merge_files(expand_filenames(sys.argv[1:]), output=infile, aggregate=False)

    dt = readData(infile)

    try:
        os.mkdir("plots")
    except FileExistsError:
        pass

    kind = "bar"
    plotCase1(dt, kind=kind)
    plotCase2(dt, kind=kind)
    plotCase3(dt, kind=kind)
    plotCase41(dt, kind=kind)
    plotCase42(dt)
    plotCase5(dt, kind=kind)
    plotRuntimes(dt, kind=kind)
