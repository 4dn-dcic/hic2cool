# hic2cool #

CLI for a converter between .hic files (from juicer) and .cool files (for cooler).  Both hic and cool files describe Hi-C contact matrices. Intended to be lightweight, this can be used as a simple imported package or a stand-alone Python file. Written by Carl Vitzthum and directed by Soo Lee of the HMS DBMI Park lab.
Originally published 1/26/17.

The .hic parsing code is based off the straw project by Neva C. Durand and Yue Wu (https://github.com/theaidenlab/straw). The h5py structure used for .cool file writing is based off code from the cooler repository: https://github.com/mirnylab/cooler.

## Using the Python package
```
$ pip install hic2cool
```

Once the package is installed, the main method is hic2cool_convert. It takes the same parameters as hic2cool.py, described at the bottom of this README. Example usage in a Python script is shown below or in /test/test.py.
```
from hic2cool import hic2cool_convert
hic2cool_convert(<infile>, <outfile>, <optional resolution>, <optional normalization type> <optional boolean to exclude MT>)
```

Or you can use the script directly. The following Python packages are required and can both be installed using pip:

## Command line converter:
```
$ python hic2cool/hic2cool_utils.py <infile> <outfile> -r <resolution> -n <normalization type> -e
```
Arguments:
**infile** is a .hic input file.
**outfile** is a .cool output file.
**resolution** is a integer bp resolution supported by the file. If 0 is given, will use all resolutions to build a multi-resolution file. Default is 0.
**normalization type** is one of: 'KR', 'NONE', 'VC', or 'VC_SQRT'. Defaults to 'KR'.
**-e**, or --exclude_MT, ignores the mitochondrial contacts if provided.

## File structure
All the information for a complete cooler file is stored in each resolution. The hdf5 hierarchy is organized as such:
File --> 'resolutions' --> '###' (where ### is the resolution in bp).
Specific resolutions are referenced using a dictionary syntax (i.e. h5File['resolutions']['###'] or h5File['resolutions/###']).
For example, see the code below that generates a one-resolution file of 10000 bp.
```
from hic2cool import hic2cool_convert
import cooler
hic2cool_convert('my_hic.hic', 'my_cool.cool', 10000, 'KR')
h5file = h5py.File('my_cool.cool', 'r')
my_cooler = cooler.Cooler(h5file['resolutions/10000'])
### do things with the my_cooler object
```

## Changelog:
### 0.2.3
Added multi-resolution format to output cool files. Setup argparse. Improved speed. Added tests for new resolutions format.
