"""
Tests for LPDAAC MODIS codes.
"""
import inspect
import os
import unittest

import matplotlib.pyplot as plt

import zoo

def fullpath(fname):
    """
    Short cut for creating the full path.
    """
    return os.path.join(os.environ['HDFEOS_ZOO_DIR'], fname)


class TestSwaths(unittest.TestCase):
    """
    Run LPDAAC/MOD swath examples.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.close()

    def test_MOD11_L2_LST(self):
        """
        """
        hdffile = 'MOD11_L2.A2007278.0350.005.2007280073443.hdf'
        hdffile = fullpath(hdffile)
        zoo.lpdaac.mod.MOD11_L2_LST.run(hdffile)

class TestGrids(unittest.TestCase):
    """
    Run LPDAAC/MOD grid examples.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.close()

    def test_MOD13C2_CMG_0_05_Deg_Monthly_NDVI(self):
        """
        """
        hdffile = 'MOD13C2.A2007001.005.2007108072029.hdf'
        hdffile = fullpath(hdffile)
        zoo.lpdaac.mod.MOD13C2_CMG_0_05_Deg_Monthly_NDVI.run(hdffile)

    def test_MOD11C2_LST_Night_CMG(self):
        """
        """
        hdffile = 'MOD11C2.A2007073.005.2007098050130.hdf'
        hdffile = fullpath(hdffile)
        zoo.lpdaac.mod.MOD11C2_LST_Night_CMG.run(hdffile)

    def test_MOD09GA_Range(self):
        """
        """
        hdffile = 'MOD09GA.A2007268.h10v08.005.2007272184810.hdf'
        hdffile = fullpath(hdffile)
        zoo.lpdaac.mod.MOD09GA_Range.run(hdffile)

    def test_MOD09GHK_sur_refl_01_1(self):
        """
        """
        hdffile = 'MOD09GHK.A2007001.h31v08.004.2007003192844.hdf'
        hdffile = fullpath(hdffile)
        zoo.lpdaac.mod.MOD09GHK_sur_refl_01_1.run(hdffile)

    def test_MOD13A1_500m_16_days_EVI(self):
        """
        """
        hdffile = 'MOD13A1.A2007257.h09v05.005.2007277183254.hdf'
        hdffile = fullpath(hdffile)
        zoo.lpdaac.mod.MOD13A1_500m_16_days_EVI.run(hdffile)

    def test_MOD43B4_Nadir_Reflectance_lvl5(self):
        """
        """
        hdffile = 'MOD43B4.A2006353.h15v15.004.2007006030047.hdf'
        hdffile = fullpath(hdffile)
        zoo.lpdaac.mod.MOD43B4_Nadir_Reflectance_lvl5.run(hdffile)

    def test_MOD17A2_PsnNet_1km(self):
        """
        """
        hdffile = 'MOD17A2.A2007113.h11v09.005.2007136163924.hdf'
        hdffile = fullpath(hdffile)
        zoo.lpdaac.mod.MOD17A2_PsnNet_1km.run(hdffile)

if __name__ == "__main__":
    unittest.main()
