mkdir build
cd build

REM Configure step
set CMAKE_CUSTOM=
REM cmake -G "%CMAKE_GENERATOR%" -DCMAKE_BUILD_TYPE=Release -DHDF4_BUILD_HL_LIB=ON -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% -DCMAKE_INSTALL_PREFIX:PATH=%LIBRARY_PREFIX% %CMAKE_CUSTOM% %SRC_DIR%
cmake -DCMAKE_BUILD_TYPE=Release -DHDF4_BUILD_HL_LIB=ON -DHDF4_BUILD_FORTRAN=OFF -DHDF4_ENABLE_NETCDF=OFF -DCMAKE_PREFIX_PATH=%LIBRARY_PREFIX% -DCMAKE_INSTALL_PREFIX:PATH=%LIBRARY_PREFIX% %CMAKE_CUSTOM% %SRC_DIR%
if errorlevel 1 exit 1

REM Build step
REM devenv HDF4.sln /Build "%RELEASE_TARGET%"
msbuild HDF4.sln /p:Configuration=Release /p:Platform=Win32
if errorlevel 1 exit 1

REM Install step
REM devenv HDF4.sln /Build "%RELEASE_TARGET%" /Project INSTALL
msbuild INSTALL.vcxproj /p:Configuration=Release /p:Platform=Win32
if errorlevel 1 exit 1
