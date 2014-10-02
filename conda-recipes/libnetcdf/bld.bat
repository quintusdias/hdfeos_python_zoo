mkdir build
cd build

REM Configure step
set CMAKE_CUSTOM=
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_DAP=OFF -DENABLE_NETCDF_4=ON -DENABLE_HDF4=ON -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% -DCMAKE_INSTALL_PREFIX:PATH=%LIBRARY_PREFIX% %CMAKE_CUSTOM% %SRC_DIR%
if errorlevel 1 exit 1

REM Build step
REM devenv netCDF.sln /Build "%RELEASE_TARGET%"
msbuild netCDF.sln
if errorlevel 1 exit 1

REM Install step
REM devenv netCDF.sln /Build "%RELEASE_TARGET%" /Project INSTALL
msbuild INSTALL.vcxproj 
if errorlevel 1 exit 1
