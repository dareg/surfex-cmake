fiat, netcdf, eccodes are needed

# compiling on belenos
On belenos, use this path to set them with tth eintel 18 compiler:
```
module load intel
module load intelmpi/2018.5.274
# Ether export the variables
export fiat_DIR=PATH_TO_FIAT_INSTALL_DIR
export NetCDF_DIR=/home/gmap/mrpm/khatib/opt/i-2018.5.274/netcdf-4.7.1/
export eccodes_DIR=/home/gmap/mrpm/khatib/opt/i-2018.5.274/eccodes-2.27.0/
# Or give them in paramter when running cmake
cmake .. -Dfiat_ROOT=PATH_TO_FIAT_INSTALL_DIR -DNetCDF_DIR=/home/gmap/mrpm/khatib/opt/i-2018.5.274/netcdf-4.7.1/ -Deccodes_ROOT=/home/gmap/mrpm/khatib/opt/i-2018.5.274/eccodes-2.27.0/

