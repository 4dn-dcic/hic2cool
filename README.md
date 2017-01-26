# hic2cool #

CLI for a converter between .hic files (from juicer) and .cool files (for
cooler).  Both hic and cool files describe Hi-C contact matrices.
Intended to be lightweight, this is a stand-alone Python file written by
Carl Vitzthum and directed by Soo Lee of the HMS DBMI Park lab.
Originally published 1/26/17.

This code was originally based off of the straw project by Neva C. Durand and
Yue Wu (https://github.com/theaidenlab/straw).
The cooler file writing was based off of much of the CLI code contained in
this repo: https://github.com/mirnylab/cooler.

The following Python packages are required and can both be installed using pip:
**h5py** pip install h5py
**numpy** pip install numpy

Usage of the converter is:
**python hic2cool.py <.hic infile> <.cool outfile> <bin size in bp>**

If an invalid bin size is given (i.e. one that does not correspond to a cooler resolution,
the program will terminate and prompt you to enter a valid one)
