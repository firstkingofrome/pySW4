"""
This file contains the binary routines for flushing data to the rFile

I want to make it clear that everything in this file comes from the pySW4 project (I made a few very small bug fixes/changes)
The reason that I have made this is because this was the only PYSW4 dependency that my code had so I just made it independent of PYSW4 (for now)
Also you really should use pySW4 to see what youre rFiles actually look like
Thank You to Shahar and everyone else whom has contributed to pySW4
"""



from __future__ import absolute_import, print_function, division
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from warnings import warn

flush = sys.stdout.flush()


def write_hdr(f, magic=1, precision=4, attenuation=1,
              az=0., lon0=33.5, lat0=28.0,
              proj_str=('+proj=utm +zone=36 +datum=WGS84 '
                        '+units=m +no_defs'),
              nb=1):
    """
    Write rfile header.

    Parameters
    ----------
    f : file
        Open file handle in ``'wb'`` mode

    magic : int
        Determine byte ordering in file. Defaults to 1.

    precision : int
        The number of bytes per entry in the data section.

        4 - single precision (default)

        8 - double precision

    attenuation : int
        Indicates whether the visco-elastic attenuation parameters QP
        and QS are included in the data section.
        0 - no visco-elastic attenuation parameters included
        1 - visco-elastic attenuation parameters included (default)

    az : float
        Angle in degrees between North and the positive x-axis.
        Defaults to 0. See the
        `SW4 User Guide <https://geodynamics.org/cig/software/sw4/>`_.

    lon0, lat0 : float
        Longitude and Latitude of the origin of the data.
        Defaults to 33.5, 28.0 .

    proj_str : str
        Projection string which is read by the Proj4 library if SW4 was
        built with Proj4 support. See the
        `SW4 User Guide <https://geodynamics.org/cig/software/sw4/>`_
        and the `Proj4 <https://trac.osgeo.org/proj/wiki/GenParms>`_
        documentation. Defaults to
        '+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs'.

    nb : int
        The number of blocks in the data section. Must be > 0.
        Defaults to 1.
    """

    magic        = np.int32(magic)
    precision    = np.int32(precision)
    attenuation  = np.int32(attenuation)
    az           = np.float64(az)
    lon0         = np.float64(lon0)
    lat0         = np.float64(lat0)
    mlen         = np.int32(len(proj_str))
    nb           = np.int32(nb)

    hdr = [magic, precision, attenuation, az, lon0, lat0, mlen, proj_str, nb]
    for val in hdr:
        f.write(val)
    return

def write_block_hdr(f, hh, hv, z0, nc, ni, nj, nk):
    """
    Write rfile block header

    Block headers are appended after the rfile header has been written.
    All block headers are written one after the other.

    Parameters
    ----------
    f : file
        Open file handle in ``'wa'`` mode.

    hh, hv : numpy.float64
        Grid size in the horizontal (x and y) and vertical (z)
        directions  in meters.

    z0 : numpy.float64
        The base z-level of the block. Not used for the
        first block which holds the elevation of the topography/
        bathymetry.

    nc : int
        The number of components:

        The first block holds the elevation of the topography/
        bathymetry, so ``nc=1``.
        The following blocks must have either 3 if only rho, vp, and vs
        are present (``attenuation=0``) or 5 if qp and qs are pressent
        (``attenuation=1``).

    ni, nj, nk : int
        Number of grid points in the i, j, k directions.

        Because the topography/bathymetry is (only) a function of the
        horizontal coordinates, the first block must have ``nk=1``.
    """

    f.write(np.float64(hh))
    f.write(np.float64(hv))
    f.write(np.float64(z0))
    f.write(np.int32(nc))
    f.write(np.int32(ni))
    f.write(np.int32(nj))
    f.write(np.int32(nk))

def write_properties(f,vp, nc, vs=None, rho=None, qp=None, qs=None):
    """
    Write material properties at a point `i`, `j` in block `b`.

    This is a convinient function to use while looping over `i`, `j` in
    a specific block `b` for writing out material properties at each
    index `k`. At the very least `vp` should be provided. If only `vp`
    is provided, the other properties are calculated using the
    :mod:`~..material_model` module.

    Parameters
    ----------
    f : file
        Open file handle in ``'wa'`` mode for appending data to the end
        of a file in construction ot in ``'r+b'`` mode for overwriting
        existing data.

        When overwriting existing data, the user must take care to place
        the cursor in the right place in the file.

    vp : array-like
        P wave velocity at indicies of `k` in m/s.

    nc : int
        Number of components to write out. Either 3 (`rho`, `vp`, and
        `vs` if ``attenuation=0``) or 5 (also `qp` and `qs` if
        ``attenuation=1``).

    vs : array-like, optional
        S wave velocity at indicies of `k` in m/s. If not given, `vs` is
        calculated from :func:`~..material_model.get_vs`.

    rho : array-like, optional
        Density at indicies of `k` in kg/m^3. If not given, `rho` is
        calculated from :func:`~..material_model.get_rho`.

    qp : array-like, optional
        P quality factor at indicies of `k`. If not given, `qp` is
        calculated from :func:`~..material_model.get_qp`.

    qs : array-like, optional
        S quality factor at indicies of `k`. If not given, `qs` is
        calculated from :func:`~..material_model.get_qs`.
    """

    k_array = np.empty((vp.size, nc), np.float32)
    vs = vs
    rho = rho 
    qs = qs 
    qp = qp 

    k_array[:, 0] = rho 
    k_array[:, 1] = vp 
    k_array[:, 2] = vs 
    k_array[:, 3] = qp
    k_array[:, 4] = qs
    k_array.tofile(f)


