#! /usr/bin/env python

from setuptools import setup, find_packages
import os


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


setup(
    name='greco_technologies',
    version='0.0.1',
    author='Inia Steinbach',
    author_email='inia.steinbach@rl-institut.de',
    description='This repository contains modules for modeling the greco '
                'PV technologies.',
    namespace_package=['greco_technologies'],
    long_description=read('README.md'),
    packages=find_packages(),
    package_dir={'greco_technologies': 'greco_technologies'},
    extras_require={
          'dev': ['sphinx', 'sphinx_rtd_theme', 'requests']},
    install_requires=[
        'pvlib >= 0.6.3',
        'pandas >= 0.24.1',
        'numpy',
        'matplotlib',
        'xarray==0.15.0'])
