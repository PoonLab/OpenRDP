try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

setup(
    name='openrdp',
    #packages=find_packages("scripts"),
    package_dir={'': 'scripts'}
)
