fiat, netcdf, eccodes are needed

# compiling on belenos
On belenos, use this path to set them with tth eintel 18 compiler:
```
module load intel
module load intelmpi/2018.5.274
# tell cmake where to find th external dependencies
cmake .. -Dfiat_ROOT=PATH_TO_FIAT_INSTALL_DIR -DNetCDF_DIR=/home/gmap/mrpm/khatib/opt/i-2018.5.274/netcdf-4.7.1/ -Deccodes_ROOT=/home/gmap/mrpm/khatib/opt/i-2018.5.274/eccodes-2.27.0/

