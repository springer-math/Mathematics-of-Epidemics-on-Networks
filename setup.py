#!/usr/bin/env python

r'''
Setup script for EoN (Epidemics on Networks)

to install from this script, run

    python setup.py install

Alternately, you can install with pip.
    
    pip install EoN
    
If this is a "release candidate" (has an "rc" in the version name below), then
pip will download the previous version - see the download_url below.
'''

from setuptools import setup

setup(name='EoN',
      packages = ['EoN'], 
      version='1.2',  #http://semver.org/
      description = 'Epidemics on Networks',
      author = 'Joel C. Miller, Istvan Z. Kiss, and Peter Simon',
      author_email = 'joel.c.miller.research@gmail.com',
      url = 'https://springer-math.github.io/Mathematics-of-Epidemics-on-Networks/',
      download_url = 'https://github.com/springer-math/Mathematics-of-Epidemics-on-Networks/archive/1.2.tar.gz',
      keywords = ['Epidemics on Networks', 'Epidemic Sonnet Works'],
      install_requires = [
          'networkx>=2',
          'numpy',
          'scipy',
          'matplotlib'
          ],
      )
