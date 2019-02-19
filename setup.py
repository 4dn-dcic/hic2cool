#!/usr/bin/env python
# Written by Carl Vitzthum with help from Soo Lee
# This code was originally based off of the straw project by Neva C. Durand and
# Yue Wu (https://github.com/theaidenlab/straw). The cooler file writing was
# based off of much of the CLI code contained in this
# repo: https://github.com/mirnylab/cooler.

import io
from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

with io.open(path.join(this_directory, 'requirements.txt')) as f:
    requires = f.read().splitlines()
requires = [req.strip() for req in requires]

this_version = io.open(path.join(this_directory, "hic2cool/_version.py")).readlines()[-1].split()[-1].strip("\"'")

setup(
    name = "hic2cool",
    version = this_version,
    packages = ['hic2cool'],
    description = """Converter between hic files (from juicer) and single-resolution or multi-resolution cool files (for cooler).  Both hic and cool files describe Hi-C contact matrices. Intended to be lightweight, this can be used as an imported package or a stand-alone Python tool for command line conversion.""",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url = "https://github.com/4dn-dcic/hic2cool",
    download_url = "https://github.com/4dn-dcic/hic2cool/tarball/" + this_version,
    author = "Carl Vitzthum",
    author_email = "carl.vitzthum@gmail.com",
    license = "MIT",
    keywords = ["bioinformatics", "genomics", "hi-c", "juicer", "cooler", "contact-matrix", "file-format"],
    install_requires = requires,
    setup_requires = requires,
    tests_require = requires,
    test_suite = "test",
    entry_points = {
        'console_scripts': [
             'hic2cool = hic2cool.__main__:main',
        ]
    },
    classifiers = [
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
