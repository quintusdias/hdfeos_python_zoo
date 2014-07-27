"""
This example code illustrates how to access and visualize a GESDISC GOSAT ACOS
L2 Swath HDF5 file in Python.

If you have any questions, suggestions, or comments on this example, please use
the HDF-EOS Forum (http://hdfeos.org/forums).  If you would like to see an
example of any other NASA HDF/HDF-EOS data product that is not listed in the
HDF-EOS Comprehensive Examples page (http://hdfeos.org/zoo), feel free to
contact us at eoshelp@hdfgroup.org or post it at the HDF-EOS Forum
(http://hdfeos.org/forums).

Usage:  save this script and run

    python acos_L2s_110101_02_Production_v110110_L2s2800_r01_PolB.py

The HDF file must either be in your current working directory or in a directory
specified by the environment variable HDFEOS_ZOO_DIR.

You can get the meaningful description of data from the "README Document for
ACOS Level 2 Standard Product" [1] or the XML description file that isprovided
by GES-DISC.

References
----------
[1] ftp://aurapar1u.ecs.nasa.gov/ftp/data/s4pa/GOSAT_TANSO_Level2/ACOS_L2S.002/doc/README.ACOS_L2S_v2.8.pdf
"""

import datetime
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

USE_NETCDF4 = True

def run(FILE_NAME):
    
    if USE_NETCDF4:

        from netCDF4 import Dataset

        nc = Dataset(FILE_NAME)
        data_var = nc.groups['RetrievalResults'].variables['xco2']
        lat_var = nc.groups['SoundingGeometry'].variables['sounding_latitude_geoid']
        lon_var = nc.groups['SoundingGeometry'].variables['sounding_longitude_geoid']
        lev_var = nc.groups['SoundingGeometry'].variables['sounding_altitude']
        time_var = nc.groups['SoundingHeader'].variables['sounding_time_tai93']

        # Read the data.
        data = data_var[:]
        lat = lat_var[:]
        lon = lon_var[:]
        lev = lev_var[:]
        time = time_var[:]

        # Read the needed attributes.
        data_units = data_var.Units
        lev_units = lev_var.Units

    else:
        
        import h5py

        with h5py.File(FILE_NAME, mode='r') as f:

            dset_var = f['/RetrievalResults/xco2']
            dset_lat = f['/SoundingGeometry/sounding_latitude_geoid']
            dset_lon = f['/SoundingGeometry/sounding_longitude_geoid']
            dset_lev = f['/SoundingGeometry/sounding_altitude']
            dset_time = f['/SoundingHeader/sounding_time_tai93']

            # Read the data.
            data = dset_var[:]
            lat = dset_lat[:]
            lon = dset_lon[:]
            lev = dset_lev[:]
            time = dset_time[:]

            # Read the needed attributes.
            data_units = dset_var.attrs['Units'][0]
            lev_units = dset_lev.attrs['Units'][0]

    # First subplot is an orthographic projection using the low resolution
    # coastline database.  Plot the trajectory.
    fig = plt.figure(figsize=(15, 6))
    plt.subplot(1, 2, 1)
    m = Basemap(projection='ortho', resolution='l', lat_0=-55, lon_0 = 120)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-80., -0., 20.))
    m.drawmeridians(np.arange(-180., 181., 20.))
    x, y = m(lon, lat)
    m.plot(x, y)

    # Annotate the starting point.  Offset the annotation text by 200 km.
    m.plot(x[0], y[0], marker='o', color='red')
    plt.annotate('START',
                 xy=(x[0] + 200000, y[0]),
                 xycoords='data',
                 color='red')

    plt.title('GOSAT Trajectory')

    # Second subplot will be time vs. CO2
    ax1 = plt.subplot(1, 2, 2)
    elapsed_time = time - time[0]
    ax1.plot(elapsed_time, data, '-', color='black')
    ax1.set_xlabel('Elapsed Time (seconds)')
    ax1.set_ylabel('CO2 column averaged dry air mole fraction\n({0})'.format(data_units))

    # Save some screen space by using scientific notation for the tick labels.
    formatter = plt.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-3, 4))
    ax1.yaxis.set_major_formatter(formatter)

    # Make a twin plot of time vs altitude.
    # The color of the plot, ylabel, and ticks should be the same.
    ax2 = ax1.twinx()
    ax2.plot(elapsed_time, lev, '-', color='blue')
    ax2.set_ylabel('Altitude ({0})'.format(lev_units), color='blue')
    for tick in ax2.get_yticklabels():
        tick.set_color('blue')

    start_time = datetime.datetime(1993,1,1) + datetime.timedelta(seconds=time[0])
    end_time = datetime.datetime(1993,1,1) + datetime.timedelta(seconds=time[-1])
    titlestr = 'CO2 column averaged dry air mole fraction and Altitude\n'
    titlestr += '{0} - {1}'.format(start_time.strftime('%Y-%m-%d %H:%M:%S'),
                                   end_time.strftime('%Y-%m-%d %H:%M:%S'))

    # Don't use "plt.title" in this case, as it will overlap with the ylabel
    # text.  We need finer control as provided by the axis object methods.
    ax2.text(0.5, 1.03, titlestr,
             horizontalalignment='center',
             fontsize=12,
             transform=ax2.transAxes)

    plt.show()
    
    basename = os.path.splitext(os.path.basename(FILE_NAME))[0]
    pngfile = "{0}.{1}.png".format(basename, 'xco2')
    fig.savefig(pngfile)

if __name__ == "__main__":

    # If a certain environment variable is set, look there for the input
    # file, otherwise look in the current directory.
    hdffile = 'acos_L2s_110101_02_Production_v110110_L2s2800_r01_PolB_110124184213.h5'
    try:
        hdffile = os.path.join(os.environ['HDFEOS_ZOO_DIR'], hdffile)
    except KeyError:
        pass

    run(hdffile)
