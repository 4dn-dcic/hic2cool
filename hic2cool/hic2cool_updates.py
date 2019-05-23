"""
This module contains updates used with the `hic2cool update` command.
See usage in hic2cool.hic2cool_utils.hic2cool_update
"""
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals
)
import h5py
from .hic2cool_config import *


def prepare_hic2cool_updates(version_nums):
    """
    Find what must be done when actually running `hic2cool update`
    Determines what updates are necessary based off of version numbers
    Version numbers is a list of ints in form: [major, minor, release]
    """
    updates = []
    # normalization vectors were inverted before version 0.5.0
    if version_nums[0] == 0 and version_nums[1] < 5:
        updates.append(
            {
                'title': 'Invert weights',
                'effect': 'Invert cooler weights so that they match original hic normalization values',
                'detail': 'cooler uses multiplicative weights and hic uses divisive weights. Before version 0.5.0, hic2cool inverted normalization vectors for consistency with cooler behavior, but now that is no longer done for consistency with 4DN analysis pipelines.',
                'function': update_invert_weights
            }
        )
    # import cooler attributes added in version 0.6.0
    if version_nums[0] == 0 and version_nums[1] < 6:
        updates.append(
            {
                'title': 'Add cooler schema version',
                'effect': 'Add a couple important cooler schema attributes',
                'detail': 'Adds format-version and storage-mode attributes to hdf5 for compatibility with cooler schema v3.',
                'function': update_cooler_schema_v3
            }
        )
    # import mcool attributes added in version 0.7.1
    if version_nums[0] == 0 and ((version_nums[1] == 7 and version_nums[2] < 1) or version_nums[1] < 7):
        updates.append(
            {
                'title': 'Add mcool schema attributes',
                'effect': 'Adds missing schema attributes if this is a multi-resolution cooler',
                'detail': 'Adds format and format-version attributes to the "/" hdf5 collection for mcool schema v2.',
                'function': update_mcool_schema_v2
            }
        )
    return updates


def norm_convert(val):
    """
    hic2cool now just uses hic normalization vectors as-is,
    without attempting to invert them to match cooler convention. This function
    is now only used with `hic2cool update` to revert cooler weights to their
    original hic values.

    Simply invert norm vectors, since hic norms are divisive and cooler
    weights are multiplicative.
    """
    if val != 0.0:
        return 1 / val
    else:
        return np.nan


def update_invert_weights(writefile):
    """
    Invert all the weights from each resolution (if a mult-res file) or the
    top level (if a single-res file)
    """
    # helper fxn
    def update_invert_weight_for_resolution(h5_data, res=None):
        """
        Access the bins table, find the weights, and invert
        """
        found_weights = [val for val in h5_data if val not in ['chrom', 'start', 'end']]
        for weight in found_weights:
            h5_weight = h5_data[weight][:]
            h5_data[weight][:] = list(map(norm_convert, h5_weight))
        if res:
            print('... For resolution %s, inverted following weights: %s' % (res, found_weights))
        else:
            print('... Inverted following weights: %s' % found_weights)

    with h5py.File(writefile) as h5_file:
        if 'resolutions' in h5_file:
            for res in h5_file['resolutions']:
                update_invert_weight_for_resolution(h5_file['resolutions'][res]['bins'], res=res)
        else:
            update_invert_weight_for_resolution(h5_file['bins'])


def update_cooler_schema_v3(writefile):
    """
    Add format-version and storage-mode attributes to given cooler
    """
    # helper fxn
    def add_v3_attrs(h5_data, res=None):
        info = {
            'format-version': COOLER_FORMAT_VERSION,
            'storage-mode': 'symmetric-upper'
        }
        h5_data.attrs.update(info)
        if res:
            print('... For resolution %s, added format-version and storage-mode attributes' % res)
        else:
            print('... Added format-version and storage-mode attributes')

    with h5py.File(writefile) as h5_file:
        if 'resolutions' in h5_file:
            for res in h5_file['resolutions']:
                add_v3_attrs(h5_file['resolutions'][res], res=res)
        else:
            add_v3_attrs(h5_file)


def update_mcool_schema_v2(writefile):
    """
    Add format and format-version attributes to the base level of an mcool
    """
    with h5py.File(writefile) as h5_file:
        # only run if it's an mcool and 'resolutions' exist
        if 'resolutions' in h5_file:
            mcool_info = {
                'format': MCOOL_FORMAT,
                'format-version': MCOOL_FORMAT_VERSION
            }
            h5_file.attrs.update(mcool_info)
            print('... Added format and format-version attributes for the mcool')
        else:
            print('... Not a multi-res file, so will not add mcool schema attributes')
