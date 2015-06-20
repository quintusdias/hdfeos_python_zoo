Instructions
============
You can run through all examples with default settings from the UNIX command
line as follows:

    $ python -m unittest discover

Some scripts have a choice as to how a file is read.  HDF-EOS grids, for
example, might be read with either GDAL or pyhdf, and pyhdf is usually the
default choice in these cases.  Those scripts can have their GDAL code tested
as follows:

    $ TEST_GDAL=1 python -m unittest discover

Finally, only those scripts for a specific data center or product may be
tested.  For example, to test just NSIDC/MODIS codes, you could use

    $ TEST_CENTER=nsidc TEST_PRODUCT=modis python -m unittest discover

You can attempt to run all examples from within ipython if you have
certain options appropriately set in your ipython configuration file
(this works for me, YMMV).

    c.InteractiveShellApp.gui = 'qt'
    c.InteractiveShellApp.pylab = 'auto'

From the top-level directory of the repository (not the package), fire
up ipython and type

    import os, unittest
    suite = unittest.defaultTestLoader.discover(os.getcwd())
    unittest.TextTestRunner(verbosity=2).run(suite)

