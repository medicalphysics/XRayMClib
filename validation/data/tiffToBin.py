"""
This file is part of XRayMClib.

XRayMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

XRayMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with XRayMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2020 Erlend Andersen
"""

import numpy as np
from skimage import io

# Path to tiff file from report TG195 case 5 describing the phantom
TIFF_FILE = r"./TG 195 Case 5 Voxelized Volume.tif"
TIFF_FILE = r"./TG 195 Case 5 Voxelized Volume.tif"


# Procedure to convert the tiff file to binary file that can be easily read by c++
def convert_tiff(path):
    im = io.imread(path)
    flat = im.flatten().astype(np.uint8)
    flat.tofile("case5world.bin")
    return


def readBinaryArray(path, dim, dtype=np.float32, reshape=True):
    size = np.prod(dim)
    vec = np.fromfile(path, dtype=dtype)
    if not reshape:
        return vec
    arr = vec.reshape(dim)
    return arr


def showArray():
    from itertools import permutations
    from matplotlib import pylab as plt

    dir_dose = r"."
    path_mat = r"./case5world.bin"
    dim = np.array([500, 320, 260], dtype=np.int)

    mat_flat = readBinaryArray(path_mat, dim, np.uint8, False)

    angles = [0, 45, 90, 135, 180, 225, 270, 315]
    for angle in angles:
        path_dose = dir_dose + "/doseAngle_{}Mono.bin".format(angle)
        dose_flat = readBinaryArray(path_dose, dim, np.float32, False)
        shape = (260, 320, 500)
        dose = dose_flat.reshape(shape)
        mat = mat_flat.reshape(shape)
        plt.subplot(3, 2, 1)
        plt.imshow(dose[shape[0] // 2, :, :], vmin=0, vmax=dose.max() * 0.01)
        plt.subplot(3, 2, 2)
        plt.imshow(mat[shape[0] // 2, :, :], vmin=0, vmax=19)
        plt.subplot(3, 2, 3)
        plt.imshow(dose[:, shape[1] // 2, :], vmin=0, vmax=dose.max() * 0.01)
        plt.subplot(3, 2, 4)
        plt.imshow(mat[:, shape[1] // 2, :], vmin=0, vmax=19)
        plt.subplot(3, 2, 5)
        plt.imshow(dose[:, :, shape[2] // 2], vmin=0, vmax=dose.max() * 0.01)
        plt.subplot(3, 2, 6)
        plt.imshow(mat[:, :, shape[2] // 2], vmin=0, vmax=19)
        plt.show()


if __name__ == "__main__":
    convert_tiff(TIFF_FILE)
    # showArray()
