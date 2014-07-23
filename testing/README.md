Instructions
============
You can attempt to run all examples from within ipython if you have
certain options appropriately set in your ipython configuration file
(this works for me, YMMV).

    c.InteractiveShellApp.gui = 'qt'
    c.InteractiveShellApp.pylab = 'auto'

From within ipython then, type

    import os, unittest
    suite = unittest.defaultTestLoader.discover(os.getcwd())
    unittest.TextTestRunner(verbosity=2).run(suite)

