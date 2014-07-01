"""
Tests for hdfeos zoo example codes.
"""
import os
import re
import shutil
import struct
import sys
import tempfile
import unittest
import uuid
from xml.etree import cElementTree as ET

import warnings

import numpy as np

import zoo

def fullpath(fname):
    """
    Short cut for creating the full path.
    """
    return os.path.join(os.environ['HDFEOS_ZOO_DIR'], fname)

class TestGesdisc(unittest.TestCase):
    """
    Run GESDISC codes.
    """
    def test_omi_l2_o2_cloudfraction(self):
        """
        """
        hdffile = 'OMI-Aura_L2-OMNO2_2008m0720t2016-o21357_v003-2008m0721t101450.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.run(hdffile)

        # Now run using h5py instead.
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.USE_NETCDF4 = False
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.run(hdffile)

class TestNSIDC(unittest.TestCase):
    """
    Run NSIDC codes.
    """
    def test_modis_snowcover(self):
        """
        """
        filename = fullpath('MOD10_L2.A2000065.0040.005.2008235221207.hdf')
        zoo.nsidc.modis.MOD10_L2_SnowCover_P.run(filename)

    def test_modis_ice_surface_temperature(self):
        """
        """
        hdffile = fullpath('MOD29.A2013196.1250.005.2013196195940.hdf')
        zoo.nsidc.modis.MOD29_A2013196_1250_005_2013196195940_hdf.run(hdffile)

