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
import os


ROW = "Material"
HUE = "Model"
COL = "Energy"


def readData():
    converters = {
        "Model": str,
        "InteractionType": str,
        "x": float,
        "y": float,
        "Energy": int,
        "Material": str,
    }
    dt = pd.read_csv(
        "validationScatterTable.txt", sep=";", engine="python", converters=converters
    )
    return dt


def plotComptonScatter(dtt, show=True):

    dt = dtt[dtt["InteractionType"] == "ComptonAngle"]

    g = sns.relplot(
        x="x",
        y="y",
        hue=HUE,
        col=COL,
        row=ROW,
        kind="line",
        data=dt,
    )
    g.set_axis_labels(y_var="Intensity [A.U]", x_var="Scatter angle [deg]")
    plt.tight_layout()
    plt.savefig("ComptonScatterAngle.png", dpi=300)
    if show:
        plt.show()


def plotComptonEnergy(dtt, show=True):

    dt = dtt[dtt["InteractionType"] == "ComptonEnergy"]

    g = sns.relplot(
        x="x",
        y="y",
        hue=HUE,
        col=COL,
        row=ROW,
        kind="line",
        data=dt,
    )
    g.set_axis_labels(y_var="Intensity [A.U]", x_var="Energy fraction")
    plt.tight_layout()
    plt.savefig("ComptonScatterEnergy.png", dpi=300)
    if show:
        plt.show()


def plotRayleightScatter(dtt, show=True):

    dt = dtt[dtt["InteractionType"] == "RayleighAngle"]

    g = sns.relplot(
        x="x",
        y="y",
        hue=HUE,
        col=COL,
        row=ROW,
        kind="line",
        data=dt,
    )
    g.set_axis_labels(y_var="Intensity [A.U]", x_var="Scatter angle [deg]")
    plt.tight_layout()
    plt.savefig("RayleighScatterAngle.png", dpi=300)
    if show:
        plt.show()


if __name__ == "__main__":
    # Setting current path to this file folder
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    dt = readData()

    plotComptonScatter(dt, show=False)
    plotComptonEnergy(dt, show=False)
    plotRayleightScatter(dt, show=False)
