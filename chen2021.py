#!/usr/bin/env python

# Copyright (C) 2016  Gregory M. Green
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#

from __future__ import print_function, division

import numpy as np
import h5py
import os
import astropy.io.fits as fits
import astropy.coordinates as coordinates
import astropy.units as units

from .std_paths import *
from .map_base import DustMap, ensure_flat_galactic
from .unstructured_map import UnstructuredDustMap
from . import fetch_utils



class Chen2021Query(UnstructuredDustMap):
   

    def __init__(self, map_fname=None):
        """
        Args:
            map_fname (Optional[:obj:`str`]): Filename at which the map is stored.
                Defaults to ``None``, meaning that the default filename is used.
        """
        if map_fname is None:
            map_fname = os.path.join(data_dir(), 'chen2021', 'chen2021.fits')
        f=fits.open(map_fname)
        self._dists = np.array([50.0])
        self._l = f[1].data['RA']
        self._b = f[1].data['DEC']
        self._A = f[1].data['EBV']
        self._A =self._A.reshape(len(self._A),1)
        self._sigma_A = f[1].data['EBV_ERROR']
        self._sigma_A =self._sigma_A.reshape(len(self._sigma_A),1)
        self._n_pix = self._l.size
	
        # Have to filter out zero pixels
        # idx = ~np.all(self._A < 1.e-5, axis=1)
        # self._lb = self._lb[idx]
        # self._A = self._A[idx]
        # self._sigma_A = self._sigma_A[idx]

        self._n_dists = self._dists.size

        # Don't query more than this angular distance from any point
        max_pix_scale = 0.5 * units.deg

        # Tesselate the sphere
        coords = coordinates.SkyCoord(
            self._l,
            self._b,
            unit='deg',
            frame='icrs')

        super(chen2021Query, self).__init__(coords, max_pix_scale, metric_p=2)

    @ensure_flat_galactic
    def query(self, coords, return_sigma=False):
        """
        Returns E(B-V) at the given coordinates. Can also
        return uncertainties.

        Args:
            coords (:obj:`astropy.coordinates.SkyCoord`): The coordinates to query.
            return_sigma (Optional[:obj:`bool`]): If ``True``, returns the uncertainty in
                extinction as well. Defaults to ``False``.

        Returns:
            E(B-V) at the specified coordinates, in mags.
            The shape of the output depends on whether :obj:`coords` contains
            distances.

            If :obj:`coords` does not specify distance(s), then the shape of the
            output begins with :obj:`coords.shape`. If :obj:`coords` does specify
            distance(s), then the shape of the output begins with
            ``coords.shape + ([number of distance bins],)``.
        """
        n_coords_ret = coords.shape[0]

        # Determine if distance has been requested
        has_dist = hasattr(coords.distance, 'kpc')
        d = coords.distance.kpc if has_dist else None

        # Convert coordinates to pixel indices
        pix_idx = self._coords2idx(coords)

        # Determine which coordinates are out of bounds
        mask_idx = (pix_idx == self._n_pix)
        if np.any(mask_idx):
            pix_idx[mask_idx] = 0

        # Which distances to extract
        if has_dist:
            d = coords.distance.kpc
            dist_idx_ceil = np.searchsorted(self._dists, d)

            ret = np.empty((n_coords_ret,), dtype='f8')
            if return_sigma:
                sigma_ret = np.empty((n_coords_ret,), dtype='f8')

            # d < d(nearest distance slice)
            idx_near = (dist_idx_ceil == 0) & ~mask_idx
            print('d < d(nearest): {:d}'.format(np.sum(idx_near)))
            if np.any(idx_near):
                a = d[idx_near] / self._dists[0]
                ret[idx_near] = a[:] * self._A[pix_idx[idx_near], 0]
                if return_sigma:
                    sigma_ret[idx_near] = a[:] * self._sigma_A[pix_idx[idx_near], 0]

            # d > d(farthest distance slice)
            idx_far = (dist_idx_ceil == self._n_dists) & ~mask_idx
            print('d > d(farthest): {:d}'.format(np.sum(idx_far)))
            if np.any(idx_far):
                ret[idx_far] = self._A[pix_idx[idx_far], -1]
                if return_sigma:
                    sigma_ret[idx_far] = self._sigma_A[pix_idx[idx_far], -1]

            # d(nearest distance slice) < d < d(farthest distance slice)
            idx_btw = ~idx_near & ~idx_far & ~mask_idx
            print('d(nearest) < d < d(farthest): {:d}'.format(np.sum(idx_btw)))
            if np.any(idx_btw):
                d_ceil = self._dists[dist_idx_ceil[idx_btw]]
                d_floor = self._dists[dist_idx_ceil[idx_btw]-1]
                a = (d_ceil - d[idx_btw]) / (d_ceil - d_floor)
                ret[idx_btw] = (
                    (1.-a[:]) * self._A[pix_idx[idx_btw], dist_idx_ceil[idx_btw]]
                    +    a[:] * self._A[pix_idx[idx_btw], dist_idx_ceil[idx_btw]-1])
                if return_sigma:
                    w0 = (1.-a)**2
                    w1 = a**2
                    norm = 1. / (w0 + w1)
                    w0 *= norm
                    w1 *= norm
                    sigma_ret[idx_btw] = np.sqrt(
                        w0 * self._sigma_A[pix_idx[idx_btw], dist_idx_ceil[idx_btw]]**2
                        + w1 * self._sigma_A[pix_idx[idx_btw], dist_idx_ceil[idx_btw]-1]**2
                    )
        else:
            # TODO: Harmonize order of distances & samples with Bayestar.
            ret = self._A[pix_idx, :]
            if return_sigma:
                sigma_ret = self._sigma_A[pix_idx, :]

        if np.any(mask_idx):
            ret[mask_idx] = np.nan
            if return_sigma:
                sigma_ret[mask_idx] = np.nan

        if return_sigma:
            return ret, sigma_ret

        return ret

    @property
    def distances(self):
        """
        Returns the distance bins that the map uses. The return type is
        :obj:`astropy.units.Quantity`, which stores unit-full quantities.
        """
        return self._dists * units.kpc




def fetch(clobber=False):
    """
    Downloads the 3D dust map of schultheis et al. (2014).
    Args:
        clobber (Optional[bool]): If ``True``, any existing file will be
            overwritten, even if it appears to match. If ``False`` (the
            default), ``fetch()`` will attempt to determine if the dataset
            already exists. This determination is not 100\% robust against data
            corruption.
    """
    dest_dir = fname_pattern = os.path.join(data_dir(), 'chen2021')
    table_fname = os.path.join(dest_dir, 'chen2021.fits')
    
    
    # Check if the FITS table already exists
    table_md5sum = '36e0e3bdca6d1e7f85c24aa3a58031ea'

    if (not clobber) and fetch_utils.check_md5sum(table_fname, table_md5sum):
        print('File appears to exist already. Call `fetch(clobber=True)` '
              'to force overwriting of existing file.')
        return

    # Download from the server
    url = 'http://paperdata.china-vo.org/guo/dust/chen2021.zip'
    archive_fname = os.path.join(dest_dir, 'chen2021.zip')
    archive_md5sum = '30e05c166aee7eca4ed0cbaedb3b57fb'

    fetch_utils.download_and_verify(url, archive_md5sum, archive_fname)

    # Extract the FITS table
    print('Exracting FITS table from Zip archive ...')
    import zipfile

    readme_fname = os.path.join(dest_dir, 'RedMe')

    with zipfile.ZipFile(archive_fname, 'r') as f:
        f.extractall(path=dest_dir)
        
      

    # Delete the Zip archive
    print('Removing Zip archive ...')
    os.remove(archive_fname)
