# try:
#     from setuptools import setup, find_packages
# except ImportError:
#     from distutils.core import setup, find_packages
#
# setup(
#     name='openrdp',
#     #packages=find_packages("scripts"),
#     package_dir={'': 'scripts'}
# )


# !/usr/bin/env python3
from distutils.core import setup
from setuptools.command.install import install

setup(
    name='OpenRDP',
    version="0.0.1",
    description='Open Source implementation of RDP5',
    packages=['openrdp'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        'scipy>=1.5.0,<1.6.0',
        'numpy>=1.17.4,<1.20.0',
    ],

    python_requires='>=3.6'
)
