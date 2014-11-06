#!/usr/bin/env python

import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="tailseq",
    version="0.0.1",
    packages=["tailseq"],
    install_requires=["pysam",
                      "ipython-cluster-helper"],
    scripts=['scripts/run_analysis.py']
)
