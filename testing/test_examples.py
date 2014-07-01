"""
Tests for hdfeos zoo example codes.
"""
import os
import unittest

import matplotlib.pyplot as plt

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
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.clf()

    def test_airs_l2_radiances_channel567(self):
        """
        """
        hdffile = 'AIRS.2002.12.31.001.L2.CC_H.v4.0.21.0.G06100185050.hdf'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.airs.AIRS_L2_radiances_channel567.run(hdffile)

    def test_airs_l3_relhumid_a_lvls11(self):
        """
        """
        hdffile = 'AIRS.2002.08.01.L3.RetStd_H031.v4.0.21.0.G06104133732.hdf'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.airs.AIRS_L3_RelHumid_A_Lvls11.run(hdffile)

    def test_airs_l3_temperature_mw_a_lvls11(self):
        """
        """
        hdffile = 'AIRS.2002.08.01.L3.RetStd_H031.v4.0.21.0.G06104133732.hdf'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.airs.AIRS_L3_Temperature_MW_A_Lvls11.run(hdffile)

    def test_omi_l2_o2_cloudfraction_netcdf4(self):
        """
        Run using netCDF4
        """
        hdffile = 'OMI-Aura_L2-OMNO2_2008m0720t2016-o21357_v003-2008m0721t101450.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.run(hdffile)

    def test_omi_l2_o2_cloudfraction_h5py(self):
        """
        Run using h5py
        """
        hdffile = 'OMI-Aura_L2-OMNO2_2008m0720t2016-o21357_v003-2008m0721t101450.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.USE_NETCDF4 = False
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.run(hdffile)

class TestNSIDC(unittest.TestCase):
    """
    Run NSIDC codes.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.clf()

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

