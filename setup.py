#!/usr/bin/env python
# Written by Carl Vitzthum with help from Soo Lee
# This code was originally based off of the straw project by Neva C. Durand and
# Yue Wu (https://github.com/theaidenlab/straw). The cooler file writing was
# based off of much of the CLI code contained in this
# repo: https://github.com/mirnylab/cooler.

from setuptools import setup

with open('requirements.txt') as f:
    requires = f.read().splitlines()
requires = [req.strip() for req in requires]

this_version = open("hic2cool/_version.py").readlines()[-1].split()[-1].strip("\"'")

setup(
    name = "hic2cool",
    version = this_version,
    packages = ['hic2cool'],
    description = """Converter between .hic files (from juicer) and single-resolution or multi-resolution .cool files (for cooler).  Both hic and cool files describe Hi-C contact matrices. Intended to be lightweight, this can be used as a simple imported package or a stand-alone Python file for command line conversion.""",
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
