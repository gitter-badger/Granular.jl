#!/usr/bin/env python
'''
SeaIce

:Package name: seaice
:Release date: 2017-04-13
:Authors: Anders Damsgaard:
:URL: https://github.com/anders-dc/seaice
:License: GPLv3
'''
from __future__ import absolute_import

from .icefloe import IceFloeCylindrical
from .grid import SquareGrid
# from packing import *

__version__ = '0.1b'
__author__ = 'Anders Damsgaard <anders.damsgaard@noaa.gov>'

__all__ = ['IceFloeCylindrical', 'SquareGrid']
