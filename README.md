# hic2cool #

A converter between .hic files (from juicer) and .cool files (for cooler).

Both hic and cool files describe Hi-C contact matrices. Intended to be lightweight, this can be used as a simple imported package or a stand-alone Python file. Written by Carl Vitzthum and directed by Soo Lee of the HMS DBMI Park lab.
Originally published 1/26/17.

The .hic parsing code is based off the straw project by Neva C. Durand and Yue Wu (https://github.com/theaidenlab/straw). The h5py structure used for .cool file writing is based off code from the cooler repository (https://github.com/mirnylab/cooler).

## Using the Python package
```
$ pip install hic2cool
```

Once the package is installed, the main method is hic2cool_convert. It takes the same parameters as hic2cool.py, described in the next section. Example usage in a Python script is shown below or in /test/test.py.
```
from hic2cool import hic2cool_convert
hic2cool_convert(<infile>, <outfile>, <optional resolution>, <optional normalization type> <optional boolean to exclude MT>)
```
Please note that you need to install cooler (pip install cooler) to run the test package. It is included in requirements.txt.


## Command line use:

If you install hic2cool itself using pip, you can simply type:
```
$ hic2cool <infile> <outfile> -r <resolution> -n <normalization type> -e
```

Otherwise, first ensure that dependencies are installed (done automatically 'pip installing' hic2cool package).
```
pip install -r requirements.txt
```

Then, from the root directory, run:
```
$ python -m hic2cool <infile> <outfile> -r <resolution> -n <normalization type> -e
```

Arguments:

**infile** is a .hic input file.

**outfile** is a .cool output file.

**resolution** is a integer bp resolution supported by the hic file. *Please note* that only resolutions contained within the original hic file can be used. If 0 is given, will use all resolutions to build a multi-resolution file. Default is 0.

**normalization type** is one of: 'KR', 'NONE', 'VC', or 'VC_SQRT'. Defaults to 'KR'.

**-e**, or --exclude_MT, ignores the mitochondrial contacts if provided.

Running hic2cool from the command line will cause some helpful information about the hic file to be printed to stdout.



## Output file structure
If you elect to use all resolutions, a multi-res .cool file will be produced. This changes the hdf5 structure of the file from a typical .cool file. Namely, all of the information needed for a complete cooler file is stored in separate hdf5 groups named by resolution. The hdf5 hierarchy is organized as such:

File --> 'resolutions' --> '###' (where ### is the resolution in bp).
For example, see the code below that generates a multi-res file and then accesses the specific resolution of 10000 bp.

```
from hic2cool import hic2cool_convert
import cooler
### using 0 triggers a mult-res output
hic2cool_convert('my_hic.hic', 'my_cool.cool', 0, 'KR')
### will give you the cooler object with resolution = 10000 bp
my_cooler = cooler.Cooler('out.cool::resolutions/10000')
```

When using only one resolution, the .cool file produced stores all the necessary information at the top level. Thus, organization in the multi-res format is not needed. The code below produces a file with one resolution, 10000 bp, and opens it with a cooler object.

```
from hic2cool import hic2cool_convert
import cooler
hic2cool_convert('my_hic.hic', 'my_cool.cool', 10000, 'KR')
h5file = h5py.File('my_cool.cool', 'r')
### will give you the cooler object with resolution = 10000 bp
my_cooler = cooler.Cooler(h5file)
```

## Changelog:

### 0.3.6
Simple release to fix pip execution
### 0.3.5
README updates, switched cooler syntax in test, and added helpful printing of hic file header info when using the command line tool.
### 0.3.4
Fixed issue where chromosome name was not getting properly set for 'All' vs 'all'.
### 0.3.3
Removed rounding fix. For now, allow py2 and py3 weights to have different number of significant figures (they're very close).
### 0.3.2
Changed output file structure for single resolution file. Resolved an issue where rounding for weights was different between python 2 and 3.
### 0.3.1
Added .travis.yml for automated testing. Changed command line running scheme. Python3 fix in hic2cool_utils.
### 0.3.0
Added multi-resolution format to output cool files. Setup argparse. Improved speed. Added tests for new resolutions format.
