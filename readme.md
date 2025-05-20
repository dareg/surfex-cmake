# compiling on belenos

In order to compile the surfex library, fiat, netcdf, eccodes and falfilfa (branch chikhiw\_fixEccodes) needs to be installed.

falfilfa can be found here: https://github.com/walidchikhi/FALFILFA/

```
module load intel
module load intelmpi/2018.5.274

cmake .. \
    -Dfiat_ROOT=PATH_TO_FIAT_INSTALL_DIR \
    -Dfalfilfa_ROOT=PATH_TO_FALFILFA_INSTALL_DIR \
    -DNetCDF_DIR=/home/gmap/mrpm/khatib/opt/i-2018.5.274/netcdf-4.7.1/ \
    -Deccodes_ROOT=/home/gmap/mrpm/khatib/opt/i-2018.5.274/eccodes-2.27.0/ \
    -DENABLE_SINGLE_PRECISION=ON \
    -DENABLE_DOUBLE_PRECISION=ON \
    -DCMAKE_INSTALL_PREFIX=PATH_TO_DESIRED_INSTALL_DIR
```
