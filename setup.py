#! /usr/bin/env python

from setuptools import setup, find_packages
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name="greco_technologies",
    version="0.0.1",
    author="Inia Steinbach",
    author_email="inia.steinbach@rl-institut.de",
    description="This repository contains modules for modeling the greco "
    "PV technologies.",
    namespace_package=["greco_technologies"],
    long_description=read("README.md"),
    packages=find_packages(),
    package_dir={"greco_technologies": "greco_technologies"},
    package_data={
        "greco_technologies": [
            "perosi/*.exe",
            "perosi/Gases/*.dat",
            "perosi/Solar/*.dat",
            "perosi/Albedo/*.dat",
            "perosi/data/CHEN_2020_EQE_curve_pero_corrected.csv",
            "perosi/data/CHEN_2020_EQE_curve_si_corrected.csv",
        ]
    },
    extras_require={"dev": ["sphinx", "sphinx_rtd_theme", "requests"]},
    install_requires=[
        "pvlib >= 0.6.3",
        "pandas >= 0.24.1",
        "numpy",
        "matplotlib",
        "xarray",
        "feedinlib==0.1.0rc2",
    ],
)
