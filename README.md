# hic2cool #

A converter between .hic files (from juicer) and .cool files (for cooler).

Both hic and cool files describe Hi-C contact matrices. Intended to be lightweight, this can be used as a simple imported package or a stand-alone Python file.

The .hic parsing code is based off the straw project by Neva C. Durand and Yue Wu (https://github.com/theaidenlab/straw). The h5py structure used for .cool file writing is based off code from the cooler repository (https://github.com/mirnylab/cooler).

## Using the Python package
```
$ pip install hic2cool
```

Once the package is installed, the main method is hic2cool_convert. It takes the same parameters as hic2cool.py, described in the next section. Example usage in a Python script is shown below or in test.py.
```
from hic2cool import hic2cool_convert
hic2cool_convert(<infile>, <outfile>, <optional resolution>, <optional show warnings>)
```
Please note that you need to install cooler (pip install cooler) to run the test package. It is included in requirements.txt.


## Command line use:

If you install hic2cool itself using pip, you can simply type:
```
$ hic2cool <infile> <outfile> -r <resolution>
```

Otherwise, first ensure that dependencies are installed (done automatically 'pip installing' hic2cool package).
```
pip install -r requirements.txt
```

Then, from the root directory, run:
```
$ python -m hic2cool <infile> <outfile> -r <resolution>
```

Arguments:

**infile** is a .hic input file.

**outfile** is a .cool output file.

**-r**, or --resolution, is a integer bp resolution supported by the hic file. *Please note* that only resolutions contained within the original hic file can be used. If 0 is given, will use all resolutions to build a multi-resolution file. Default is 0.

**-w**, or --warnings, causes warnings to be explicitly printed to the console. These are hidden by default, though there are a few cases in which hic2cool will exit with an error based on the input hic file.

**-v**, or --version, print out hic2cool package version and exit.

Running hic2cool from the command line will cause some helpful information about the hic file to be printed to stdout.



## Output file structure
If you elect to use all resolutions, a multi-res .multi.cool file will be produced. This changes the hdf5 structure of the file from a typical .cool file. Namely, all of the information needed for a complete cooler file is stored in separate hdf5 groups named by the individual resolutions. The hdf5 hierarchy is organized as such:

File --> 'resolutions' --> '###' (where ### is the resolution in bp).
For example, see the code below that generates a multi-res file and then accesses the specific resolution of 10000 bp.

```
from hic2cool import hic2cool_convert
import cooler
### using 0 triggers a multi-res output
hic2cool_convert('my_hic.hic', 'my_cool.cool', 0)
### will give you the cooler object with resolution = 10000 bp
my_cooler = cooler.Cooler('my_cool.cool::resolutions/10000')
```

When using only one resolution, the .cool file produced stores all the necessary information at the top level. Thus, organization in the multi-res format is not needed. The code below produces a file with one resolution, 10000 bp, and opens it with a cooler object.

```
from hic2cool import hic2cool_convert
import cooler
### giving a specific resolution below (e.g. 10000) triggers a single-res output
hic2cool_convert('my_hic.hic', 'my_cool.cool', 10000)
h5file = h5py.File('my_cool.cool', 'r')
### will give you the cooler object with resolution = 10000 bp
my_cooler = cooler.Cooler(h5file)
```

## Changelog:

### 0.4.1
* Fixed error in reading counts from hic files of version 6
* Chromosome names are now directly taken from hic file (with exception of 'all')
### 0.4.0
Large patch, should fix most memory issues and improve runtimes:
* Changed run parameters. Removed -n and -e; added -v (--version) and -w (--warnings)
* Improved memory usage
* Improved runtime (many thanks to Nezar Abdennur)
* hic2cool now does a 'direct' conversion of files and does not fail on missing chr-chr contacts or missing normalization vectors. Finding these issues will cause warnings to be printed (controlled by -w flag)
* No longer uses the 'weights' column, which is reserved for cooler
* No longer takes a normalization type argument. All normalization vectors from the hic file are automatically added to the bins table in the output .cool
* Many other minor bug fixes/code improvement
### 0.3.7
Fixed issue with bin1_offset not containing final entry (should be length nbins + 1).
### 0.3.6
Simple release to fix pip execution.
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

## Contributors
Written by Carl Vitzthum (1), Nezar Abdennur (2), Soo Lee (1), and Peter Kerpedjiev (3).

(1) Park lab, Harvard Medical School DBMI

(2) Mirny lab, MIT

(3) Gehlenborg lab, Harvard Medical School DBMI

Originally published 1/26/17.
