#!/usr/bin/env python

from setuptools import setup, find_packages


setup(
    name='isorefiner',
    version='0.0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'isorefiner = isorefiner.main:main',
        ],
    },
)
