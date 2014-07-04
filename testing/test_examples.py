"""
Tests for hdfeos zoo example codes.
"""
import inspect
import os
import unittest

import matplotlib.pyplot as plt

import zoo

from . import fixtures

def fullpath(fname):
    """
    Short cut for creating the full path.
    """
    return os.path.join(os.environ['HDFEOS_ZOO_DIR'], fname)

class TestDocstrings(unittest.TestCase):
    """
    Verify information in the docstrings of the examples.
    """
    def test_docstring(self):
        """
        Verify that the docstring in each example has the EOS contact info.
        """
        for center_name, center_module in inspect.getmembers(zoo, inspect.ismodule):
            for inst_name, inst_module in inspect.getmembers(center_module, inspect.ismodule):
                for example_name, example_module in inspect.getmembers(inst_module, inspect.ismodule):
                    msg = "Failed to verify docstring in {0}".format(example_name)
                    docstring = example_module.__doc__.replace('\n', ' ')
                    contact_info = fixtures.contact_info.replace('\n', ' ')
                    self.assertTrue(contact_info in docstring, msg)


class TestGesdiscAirs(unittest.TestCase):
    """
    Run GESDISC/AIRS codes.
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

class TestGesdiscBuv(unittest.TestCase):
    """
    Run GESDISC/BUV codes.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.clf()

    def test_BUV_Nimbus04_L3zm_v01_00_2012m0203t144121_h5_h5py(self):
        """
        Use h5py.
        """
        zoo.gesdisc.buv.BUV_Nimbus04_L3zm_v01_00_2012m0203t144121_h5.USE_NETCDF4 = False
        hdffile = 'BUV-Nimbus04_L3zm_v01-00-2012m0203t144121.h5'
        zoo.gesdisc.buv.BUV_Nimbus04_L3zm_v01_00_2012m0203t144121_h5.run(fullpath(hdffile))

    def test_BUV_Nimbus04_L3zm_v01_00_2012m0203t144121_h5_netcdf(self):
        """
        Use netcdf4.
        """
        zoo.gesdisc.buv.BUV_Nimbus04_L3zm_v01_00_2012m0203t144121_h5.USE_NETCDF4 = True
        hdffile = 'BUV-Nimbus04_L3zm_v01-00-2012m0203t144121.h5'
        zoo.gesdisc.buv.BUV_Nimbus04_L3zm_v01_00_2012m0203t144121_h5.run(fullpath(hdffile))

    def test_SBUV2_NOAA17_L2_SBUV2N17L2_2011m1231_v01_01_2012m0905t152911_h5py(self):
        """
        Use h5py.
        """
        zoo.gesdisc.buv.SBUV2_NOAA17_L2_SBUV2N17L2_2011m1231_v01_01_2012m0905t152911_h5.USE_NETCDF4 = False
        hdffile = 'SBUV2-NOAA17_L2-SBUV2N17L2_2011m1231_v01-01-2012m0905t152911.h5'
        zoo.gesdisc.buv.SBUV2_NOAA17_L2_SBUV2N17L2_2011m1231_v01_01_2012m0905t152911_h5.run(fullpath(hdffile))

    def test_SBUV2_NOAA17_L2_SBUV2N17L2_2011m1231_v01_01_2012m0905t152911_netcdf4(self):
        """
        Use netcdf4.
        """
        zoo.gesdisc.buv.SBUV2_NOAA17_L2_SBUV2N17L2_2011m1231_v01_01_2012m0905t152911_h5.USE_NETCDF4 = True
        hdffile = 'SBUV2-NOAA17_L2-SBUV2N17L2_2011m1231_v01-01-2012m0905t152911.h5'
        zoo.gesdisc.buv.SBUV2_NOAA17_L2_SBUV2N17L2_2011m1231_v01_01_2012m0905t152911_h5.run(fullpath(hdffile))

class TestGesdiscMerra(unittest.TestCase):
    """
    Run GESDISC/MERRA codes.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.clf()

    def test_MERRA_PLE_TIME1_Height72(self):
        """
        """
        hdffile = 'MERRA300.prod.assim.inst3_3d_chm_Ne.20021201.hdf'
        zoo.gesdisc.merra.MERRA_PLE_TIME1_Height72.run(fullpath(hdffile))

    def test_MERRA_MFYC_TIME4_Height42(self):
        """
        """
        hdffile = 'MERRA300.prod.assim.tavg3_3d_chm_Nv.20021201.hdf'
        zoo.gesdisc.merra.MERRA_MFYC_TIME4_Height42.run(fullpath(hdffile))

class TestGesdiscMls(unittest.TestCase):
    """
    Run GESDISC/MLS codes.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.clf()

    def test_MLS_L2GP_v01_L2gpValue_h5py(self):
        """
        Run using h5py
        """
        hdffile = 'MLS-Aura_L2GP-BrO_v01-52-c01_2007d029.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.mls.MLS_L2GP_v01_L2gpValue.USE_NETCDF4 = False
        zoo.gesdisc.mls.MLS_L2GP_v01_L2gpValue.run(hdffile)

    def test_MLS_L2GP_v01_L2gpValue_netcdf(self):
        """
        Run using netcdf
        """
        hdffile = 'MLS-Aura_L2GP-BrO_v01-52-c01_2007d029.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.mls.MLS_L2GP_v01_L2gpValue.USE_NETCDF4 = True
        zoo.gesdisc.mls.MLS_L2GP_v01_L2gpValue.run(hdffile)

    def test_MLS_L2GP_v02_L2gpValue_h5py(self):
        """
        Run using h5py
        """
        hdffile = 'MLS-Aura_L2GP-BrO_v02-23-c01_2010d255.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.mls.MLS_L2GP_v01_L2gpValue.USE_NETCDF4 = False
        zoo.gesdisc.mls.MLS_L2GP_v01_L2gpValue.run(hdffile)

    def test_MLS_L2GP_v02_L2gpValue_netcdf(self):
        """
        Run using netcdf
        """
        hdffile = 'MLS-Aura_L2GP-BrO_v02-23-c01_2010d255.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.mls.MLS_L2GP_v01_L2gpValue.USE_NETCDF4 = True
        zoo.gesdisc.mls.MLS_L2GP_v01_L2gpValue.run(hdffile)

class TestGesdiscOmi(unittest.TestCase):
    """
    Run GESDISC/OMI codes.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.clf()

    def test_omi_omcldo2g_h5py(self):
        """
        Run using h5py
        """
        hdffile = 'OMI-Aura_L2G-OMCLDO2G_2007m0129_v002-2007m0130t174603.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_OMCLDO2G.USE_NETCDF4 = False
        zoo.gesdisc.omi.OMI_OMCLDO2G.run(hdffile)

    def test_omi_omcldo2g_netcdf4(self):
        """
        Run using netCDF4
        """
        hdffile = 'OMI-Aura_L2G-OMCLDO2G_2007m0129_v002-2007m0130t174603.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_OMCLDO2G.USE_NETCDF4 = True
        zoo.gesdisc.omi.OMI_OMCLDO2G.run(hdffile)

    def test_omi_l2_o2_cloudfraction_netcdf4(self):
        """
        Run using netCDF4
        """
        hdffile = 'OMI-Aura_L2-OMNO2_2008m0720t2016-o21357_v003-2008m0721t101450.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.USE_NETCDF4 = True
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.run(hdffile)

    def test_omi_l2_o2_cloudfraction_h5py(self):
        """
        Run using h5py
        """
        hdffile = 'OMI-Aura_L2-OMNO2_2008m0720t2016-o21357_v003-2008m0721t101450.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.USE_NETCDF4 = False
        zoo.gesdisc.omi.OMI_L2_OMNO2_CloudFraction.run(hdffile)

    def test_omi_l3_columnamounto3_h5py(self):
        """
        Run using h5py
        """
        hdffile = 'OMI-Aura_L3-OMTO3e_2005m1214_v002-2006m0929t143855.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_L3_ColumnAmountO3.USE_NETCDF4 = False
        zoo.gesdisc.omi.OMI_L3_ColumnAmountO3.run(hdffile)

    def test_omi_l3_columnamounto3_netcdf4(self):
        """
        Run using netCDF4
        """
        hdffile = 'OMI-Aura_L3-OMTO3e_2005m1214_v002-2006m0929t143855.he5'
        hdffile = fullpath(hdffile)
        zoo.gesdisc.omi.OMI_L3_ColumnAmountO3.USE_NETCDF4 = True
        zoo.gesdisc.omi.OMI_L3_ColumnAmountO3.run(hdffile)


class TestGesdiscTOMS(unittest.TestCase):
    """
    Run GESDISC/TOMS codes.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.clf()

    def test_TOMS_L3_Ozone(self):
        """
        """
        hdffile = 'TOMS-EP_L3-TOMSEPL3_2000m0101_v8.HDF'
        zoo.gesdisc.toms.TOMS_L3_Ozone.run(fullpath(hdffile))

class TestGesdiscTRMM(unittest.TestCase):
    """
    Run GESDISC/TRMM codes.
    """
    def tearDown(self):
        """
        Clear any open figure windows.
        """
        plt.clf()

    def test_TRMM_1B21_binDIDHmean(self):
        """
        """
        hdffile = fullpath('1B21.071022.56609.6.HDF')
        zoo.gesdisc.trmm.TRMM_1B21_binDIDHmean.run(hdffile)

    def test_TRMM_1B21_19971208_00170_7_HDF(self):
        """
        """
        hdffile = fullpath('1B21.19971208.00170.7.HDF')
        zoo.gesdisc.trmm.TRMM_1B21_19971208_00170_7_HDF.run(hdffile)

    def test_TRMM_1B21_CSI_binDIDHmean_zoom(self):
        """
        """
        hdffile = fullpath('1B21_CSI.990906.10217.KORA.6.HDF')
        zoo.gesdisc.trmm.TRMM_1B21_CSI_binDIDHmean_zoom.run(hdffile)

    def test_TRMM_2A12_cldWater_lvl9(self):
        """
        """
        hdffile = fullpath('2A12.100402.70512.6.HDF')
        zoo.gesdisc.trmm.TRMM_2A12_cldWater_lvl9.run(hdffile)

    def test_TRMM_2A12_20140308_92894_7_HDF(self):
        """
        """
        hdffile = '2A12.20140308.92894.7.HDF'
        zoo.gesdisc.trmm.TRMM_2A12_20140308_92894_7_HDF.run(fullpath(hdffile))

    def test_TRMM_2A25_CSI_nearSurfZ_zoom(self):
        """
        """
        hdffile = '2A25_CSI.990804.9692.KORA.6.HDF'
        zoo.gesdisc.trmm.TRMM_2A25_CSI_nearSurfZ_zoom.run(fullpath(hdffile))

    def test_TRMM_2B31_CSI_dHat(self):
        """
        """
        hdffile = '2B31_CSI.990911.10296.KORA.6.HDF'
        zoo.gesdisc.trmm.TRMM_2B31_CSI_dHat_zoom.run(fullpath(hdffile))

    def test_TRMM_3B42_precipitation_scan0(self):
        """
        """
        hdffile = '3B42.100331.21.6A.HDF'
        zoo.gesdisc.trmm.TRMM_3B42_precipitation_scan0.run(fullpath(hdffile))

    def test_TRMM_3B43_precipitation_scan0(self):
        """
        """
        hdffile = '3B43.070901.6A.HDF'
        zoo.gesdisc.trmm.TRMM_3B43_precipitation_scan0.run(fullpath(hdffile))

    def test_TRMM_3A46_ssmiData(self):
        """
        """
        hdffile = '3A46.080101.2.HDF'
        zoo.gesdisc.trmm.TRMM_3A46_ssmiData.run(fullpath(hdffile))


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

