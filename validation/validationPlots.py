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

HUE_ORDER = list(["EGSnrc", "Geant4", "MCNP", "Penelope", "IA", "Livermore", "NoneLC"])


plt.rcParams["font.family"] = "Times New Roman"
sns.set_theme(palette="husl", font="Times New Roman", style="ticks")


def readTG195Data(path):
    converters = {
        "Case": str,
        "Volume": str,
        "Specter": str,
        "Mode": str,
        "EGSnrc_Result": float,
        "EGSnrc_Uncertainty": float,
        "Geant4_Result": float,
        "Geant4_Uncertainty": float,
        "MCNP_Result": float,
        "MCNP_Uncertainty": float,
        "Penelope_Result": float,
        "Penelope_Uncertainty": float,
    }
    dt_long = pd.read_csv(path, sep=",", engine="python", converters=converters)

    ##adding mean results
    agg_names = ["EGSnrc_Result", "Geant4_Result", "MCNP_Result", "Penelope_Result"]
    agg = pd.DataFrame({n: dt_long[n] for n in agg_names})
    dt_long["TG195"] = agg.mean(axis=1, skipna=True)

    return dt_long


def readdxmcData(path):
    converters = {
        "Case": str,
        "Volume": str,
        "Specter": str,
        "Mode": str,
        "Model": str,
        "Result": float,
        "Uncertainty": float,
    }
    cols = list(converters.keys())
    dx = pd.read_csv(
        path, sep=",", engine="python", converters=converters, usecols=cols
    )

    # Transform from long to wide form
    models = list(set(dx["Model"]))
    longs = list([(m, dx[dx["Model"] == m]) for m in models])
    for _, l in longs:
        del l["Model"]
    if len(longs) > 1:
        m, dxl = longs.pop(0)
        dxl.rename(
            columns={
                "Result": "{}_Result".format(m),
                "Uncertainty": "{}_Uncertainty".format(m),
            },
            inplace=True,
        )

        for m, dxr in longs:
            dxl = dxl.merge(
                dxr,
                how="outer",
            )
            dxl.rename(
                columns={
                    "Result": "{}_Result".format(m),
                    "Uncertainty": "{}_Uncertainty".format(m),
                },
                inplace=True,
            )
    else:
        m = models[0]
        dxl = dx.rename(
            columns={
                "Result": "{}_Result".format(m),
                "Uncertainty": "{}_Uncertainty".format(m),
            }
        )
        del dxl["Model"]
    return dxl


def readData(dxmc_path, TG195_path, relative_percent=True):

    dx = readdxmcData(dxmc_path)
    dtg = readTG195Data(TG195_path)
    dt_long = dtg.merge(dx, how="outer")

    # dropping non existing data
    dt_long.dropna(
        how="any",
        subset=[k + "_Result" for k in ["IA", "Livermore", "NoneLC"]],
        inplace=True,
    )

    col = dt_long.columns
    models = ["EGSnrc", "Geant4", "MCNP", "Penelope", "IA", "Livermore", "NoneLC"]
    models = [m for m in models if "{}_Result".format(m) in col]

    ## making uncertainty absolute values
    for m in models:
        u = "{}_Uncertainty".format(m)
        r = "{}_Result".format(m)
        dt_long[u] = dt_long[r] * dt_long[u]

    ## making relative data
    for m in models:
        r = "{}_Result".format(m)
        dt_long["{}_relative".format(r)] = dt_long[r] / dt_long["TG195"]
        u = "{}_Uncertainty".format(m)
        dt_long["{}_relative".format(u)] = dt_long[u] / dt_long["TG195"]

    ## making long form
    melted = list()
    for c in ["Result", "Uncertainty", "Result_relative", "Uncertainty_relative"]:
        melt = dt_long.melt(
            id_vars=["Case", "Volume", "Specter", "Mode"],
            value_vars=["{}_{}".format(m, c) for m in models],
            var_name="Model",
            value_name=c,
        )
        melt["Model"] = melt["Model"].replace(
            to_replace={"{}_{}".format(m, c): m for m in models}, inplace=False
        )
        melted.append(melt)

    dt = melted.pop(0)
    for melt in melted:
        dt = dt.merge(melt, how="outer")

    ##dropping missing TG195 results
    dt.dropna(how="any", inplace=True)

    ##Remove non existing models
    for h in HUE_ORDER.copy():
        if h not in models:
            HUE_ORDER.remove(h)

    if relative_percent:
        dt["Result_relative"] = (dt["Result_relative"] - 1) * 100
        dt["Uncertainty_relative"] = (dt["Uncertainty_relative"]) * 100

    ## adding conf intervals
    dt_max = dt.copy()
    dt_min = dt.copy()
    dt_max["Result"] = dt["Result"] + dt["Uncertainty"]
    dt_max["Result_relative"] = dt["Result_relative"] + dt["Uncertainty_relative"]
    dt_min["Result"] = dt["Result"] - dt["Uncertainty"]
    dt_min["Result_relative"] = dt["Result_relative"] - dt["Uncertainty_relative"]
    dt = pd.concat([dt, dt_min, dt_max])
    return dt


def fix_axis(
    fg,
    rotate_labels=True,
    ylabel=r"Energy imparted $\left[ \frac{\mathdefault{eV}}{\mathdefault{history}} \right]$",
    set_y0=False,
):
    if set_y0:
        fg.set(ylim=(0, None))
    fg.set_axis_labels(y_var=ylabel)
    if rotate_labels:
        for ax in fg.axes_dict.values():
            for tick in ax.get_xticklabels():
                tick.set_rotation(70)
    fg.tight_layout()


def plotCase1(dt_full, kind="strip", show=False, relative=False):
    dt = dt_full[dt_full["Case"] == "Case 1"]
    if dt.size == 0:
        return

    labels = ["NoneFilter", "HVLFilter", "QVLFilter"]
    g = sns.catplot(
        x="Volume",
        y="Result_relative" if relative else "Result",
        hue="Model",
        col="Specter",
        col_wrap=2,
        order=labels,
        hue_order=HUE_ORDER,
        data=dt,
        kind=kind,
        errorbar=lambda x: (x.min(), x.max()),
    )
    if relative:
        ylabel = "Difference from mean of TG195 [%]"
    else:
        ylabel = r"$\mathdefault{KERMA}_{\mathdefault{Air}} $ per history"

    fix_axis(g, ylabel=ylabel)

    plt.savefig(
        "plots/Case1_{}.png".format("relative" if relative else "value"), dpi=300
    )
    if show:
        plt.show()
    plt.clf()
    plt.close()


def plotCase2and3(case, dt_full, kind="strip", show=False, relative=False):
    casename = case.replace(" ", "")
    dt = dt_full[dt_full["Case"] == case]
    if dt.size == 0:
        return
    if relative:
        dt_vol = dt
    else:
        dt_vol = dt[dt["Volume"] != "Total body"]

    labels = list(set([v for v in dt_vol["Volume"]]))
    labels.sort()

    g = sns.catplot(
        x="Volume",
        y="Result_relative" if relative else "Result",
        hue="Model",
        row="Specter",
        col="Mode",
        order=labels,
        hue_order=HUE_ORDER,
        data=dt_vol,
        kind=kind,
        errorbar=lambda x: (x.min(), x.max()),
    )

    if relative:
        ylabel = "Difference from mean of TG195 [%]"
    else:
        ylabel = r"Energy imparted $\left[ \frac{\mathdefault{eV}}{\mathdefault{history}} \right]$"

    fix_axis(g, ylabel=ylabel)
    plt.savefig(
        "plots/{}_{}.png".format(casename, "relative" if relative else "value"), dpi=300
    )
    if show:
        plt.show()
    plt.clf()
    plt.close()

    if not relative:
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
        plt.savefig("plots/{}_value_TotalBody.png".format(casename), dpi=300)
        if show:
            plt.show()
        plt.clf()
        plt.close()


def plotCase41(dt_full, kind="strip", show=False, relative=False):
    dt = dt_full[dt_full["Case"] == "Case 4.1"]
    if dt.size == 0:
        return

    labels = list(set([v for v in dt["Volume"]]))
    labels.sort()

    g = sns.catplot(
        x="Volume",
        y="Result_relative" if relative else "Result",
        hue="Model",
        row="Specter",
        col="Mode",
        order=labels,
        hue_order=HUE_ORDER,
        data=dt,
        kind=kind,
        errorbar=lambda x: (x.min(), x.max()),
    )
    if relative:
        ylabel = "Difference from mean of TG195 [%]"
    else:
        ylabel = r"Energy imparted $\left[ \frac{\mathdefault{eV}}{\mathdefault{history}} \right]$"

    fix_axis(g, rotate_labels=False, ylabel=ylabel)
    plt.savefig(
        "plots/Case41_{}.png".format("relative" if relative else "value"), dpi=300
    )
    if show:
        plt.show()
    plt.close()


def plotCase42(dt_full, show=False, kind="strip", relative=False):
    dt = dt_full[dt_full["Case"] == "Case 4.2"]
    if dt.size == 0:
        return
    dt["Volume [angle]"] = [int(d) for d in dt["Volume"]]
    dtp_ind = ["Cent" in e for e in dt["Mode"]]
    dtp = dt[dtp_ind]
    for mode in set(dtp["Mode"]):
        g = sns.catplot(
            x="Volume [angle]",
            y="Result_relative" if relative else "Result",
            hue="Model",
            row="Specter",
            # col="Mode",
            col=None,
            kind=kind,
            hue_order=HUE_ORDER,
            data=dtp[dtp["Mode"] == mode],
            errorbar=lambda x: (x.min(), x.max()),
            aspect=2,
        )
        if relative:
            ylabel = "Difference from mean of TG195 [%]"
        else:
            ylabel = r"Energy imparted $\left[ \frac{\mathdefault{eV}}{\mathdefault{history}} \right]$"

        fix_axis(g, ylabel=ylabel)
        plt.savefig(
            "plots/Case42_Cent_{}_{}.png".format(
                mode, "relative" if relative else "value"
            ),
            dpi=300,
        )
        if show:
            plt.show()
        plt.close()

    dtp_ind = ["Pher" in e for e in dt["Mode"]]
    dtp = dt[dtp_ind]
    for mode in set(dtp["Mode"]):
        g = sns.catplot(
            x="Volume [angle]",
            y="Result_relative" if relative else "Result",
            hue="Model",
            row="Specter",
            col=None,
            kind=kind,
            hue_order=HUE_ORDER,
            data=dtp[dtp["Mode"] == mode],
            errorbar=lambda x: (x.min(), x.max()),
            aspect=2,
        )
        fix_axis(g, ylabel=ylabel)
        plt.savefig(
            "plots/Case42_Pher_{}_{}.png".format(
                mode, "relative" if relative else "value"
            ),
            dpi=300,
        )
        if show:
            plt.show()
        plt.close()


def plotCase5(dt_full, kind="strip", show=False, relative=False):
    dt = dt_full[dt_full["Case"] == "Case 5"]
    if dt.size == 0:
        return

    # dt["Volume"] = [int(v) for v in dt["Volume"]]
    labels = list(set(dt["Volume"]))
    labels.sort()

    vols = set([m for m in dt["Mode"]])

    for m in vols:
        dtm = dt[dt["Mode"] == m]
        g = sns.catplot(
            x="Volume",
            y="Result_relative" if relative else "Result",
            hue="Model",
            col=None,
            row="Specter",
            order=labels,
            hue_order=HUE_ORDER,
            data=dtm,
            kind=kind,
            errorbar=lambda x: (x.min(), x.max()),
            aspect=2,
        )
        if relative:
            ylabel = "Difference from mean of TG195 [%]"
        else:
            ylabel = r"Energy imparted $\left[ \frac{\mathdefault{eV}}{\mathdefault{history}} \right]$"
        fix_axis(g, ylabel=ylabel)
        plt.savefig(
            "plots/Case5_{}_{}.png".format(m, "relative" if relative else "value"),
            dpi=300,
        )
        if show:
            plt.show()
        plt.close()


if __name__ == "__main__":
    # Setting current path to this file folder
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    data = readData("ValidationTable.txt", "TG195Results.txt")

    try:
        os.mkdir("plots")
    except FileExistsError:
        pass

    # kind = "bar"
    kind = "point"
    plotCase1(data, kind=kind, relative=False)
    plotCase1(data, kind=kind, relative=True)
    plotCase2and3("Case 2", data, kind=kind, relative=False)
    plotCase2and3("Case 2", data, kind=kind, relative=True)
    plotCase2and3("Case 3", data, kind=kind, relative=False)
    plotCase2and3("Case 3", data, kind=kind, relative=True)
    plotCase41(data, kind=kind, relative=False)
    plotCase41(data, kind=kind, relative=True)
    plotCase42(data, kind=kind, relative=False)
    plotCase42(data, kind=kind, relative=True)
    plotCase5(data, kind=kind, relative=True)
    plotCase5(data, kind=kind, relative=False)
