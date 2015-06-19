"""
Tests for hdfeos zoo example codes.
"""
import inspect
import os
import unittest

import matplotlib.pyplot as plt

import zoo

class TestAll(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.data_centers = {'gesdisc': ('airs', 'buv', 'gosat', 'gsstf',
                                        'hirdls', 'merra', 'mls', 'oco2',
                                        'omi', 'swdb', 'toms', 'trmm'),
                            'laads':  ('mod', 'myd', 'viirs'),
                            'larc':  ('calipso', 'ceres', 'misr', 'mopitt',
                                      'tes'),
                            'lpdaac':  ('mcd', 'mod', 'myd', 'vip', 'weld',
                                        'ged'),
                            'nsidc': ('icesat', 'modis', 'nise', 'amsre'),
                            'podaac': ('aquarius', 'avhrr', 'quikscat',
                                       'seawinds')}

    def test_all(self):
        for center_name in self.data_centers.keys():
            center = getattr(zoo, center_name)
            if center_name != 'podaac':
                continue
            for product_name in self.data_centers[center_name]:
                #if product_name != 'quikscat':
                #    continue
                product_module = getattr(center, product_name)
                pargs = (product_module, inspect.ismodule)
                for example_name, example_module in inspect.getmembers(*pargs):

                    # Force the use of gdal where it is optional,
                    if 'TEST_GDAL' in os.environ.keys():
                        if not hasattr(example_module, 'USE_GDAL'):
                            continue
                        example_module.USE_GDAL = True

                    # Force the use of netcdf where it is optional,
                    if 'TEST_NETCDF' in os.environ.keys():
                        if not hasattr(example_module, 'USE_NETCDF4'):
                            continue
                        example_module.USE_NETCDF4 = True

                    print("Testing {}.{}.{}".format(center_name,
                                                    product_name,
                                                    example_name))
                    example_module.run()
                    plt.close()
