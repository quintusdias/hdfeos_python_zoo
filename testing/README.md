Instructions
============
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

