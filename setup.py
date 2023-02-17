# !/usr/bin/env python3
from setuptools import setup

setup(
    name='OpenRDP',
    version="0.1.0",
    description='Open Source implementation of RDP5',
    packages=['openrdp'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        'scipy>=1.5.0',
        'numpy>=1.17.4',
        'h5py>=3.8.0'
    ],
    python_requires='>=3.6',
    scripts=['bin/openrdp'],
    options={'build_scripts': {'executable': '/usr/bin/env python3'}},
    package_data={
        'openrdp': [
            'default.ini',
            'bin/3Seq/*',
            'bin/GENECONV/*'
        ]
    },
    zip_safe=False
)
