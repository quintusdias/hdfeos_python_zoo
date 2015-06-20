"""
Tests for hdfeos zoo example codes.
"""
import inspect
import os
import unittest

import matplotlib.pyplot as plt

import zoo

from . import fixtures

class TestDocstrings(unittest.TestCase):
    """
    Verify information in the docstrings of the examples.
    """

    def test_docstring_contact_info(self):
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

    def test_docstring_run_instructions(self):
        """
        Verify instructions to run each script.
        """
        for center_name, center_module in inspect.getmembers(zoo, inspect.ismodule):
            for inst_name, inst_module in inspect.getmembers(center_module, inspect.ismodule):
                for example_name, example_module in inspect.getmembers(inst_module, inspect.ismodule):
                    msg = "Failed to verify docstring in {0}".format(example_name)
                    docstring = example_module.__doc__.replace('\n', ' ')
                    run_info = fixtures.run_info.replace('\n', ' ')
                    run_info = run_info.format(example_name)
                    self.assertTrue(run_info in docstring, msg)


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

            if 'TEST_CENTER' in os.environ.keys():
                if center_name != os.environ['TEST_CENTER']:
                    continue

            for product_name in self.data_centers[center_name]:

                if 'TEST_PRODUCT' in os.environ.keys():
                    if product_name != os.environ['TEST_PRODUCT']:
                        continue

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
