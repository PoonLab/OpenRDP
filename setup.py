# !/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name='OpenRDP',
    version="0.0.1",
    description='Open Source implementation of RDP5',
    packages=find_packages(include=['openrdp', 'openrdp.*']),
    # packages=['openrdp'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        'scipy>=1.5.0,<1.6.0',
        'numpy>=1.17.4,<1.20.0',
    ],
    python_requires='>=3.6',
    package_data={'openrdp': ['tests/*.fasta', 'tests/*.fa', 'tests/*.ini', 'bin/3Seq/*', 'bin/GENECONV/*']},
    zip_safe=False
)
