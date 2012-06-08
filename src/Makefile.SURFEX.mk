##########################################################
#                                                        #
#           Initialisation of some variables             #
#                                                        #
##########################################################
ifdef OBJDIR_PATH
OBJDIR_ROOT=${OBJDIR_PATH}/dir_obj
else
OBJDIR_ROOT=${PWD}/dir_obj
endif
LIB_OBJS_ROOT=lib
#
ARCH_XYZ=${ARCH}${MNH_INT}-${VERSION_XYZ}
##########################################################
#                                                        #
#            Source DIRECTORY                            #
#                                                        #
##########################################################

##########################################################
#           Source MYSRC                                 #
##########################################################
ifdef VER_USER
DIR_USER += ${VER_USER}
endif
##########################################################
#           Source OFFLIN                                #
##########################################################
# PRE_BUG TEST !!!
#DIR_SURFEX += ARCH_SRC/bug_surfex
# PRE_BUG TEST !!!
#
DIR_OFFLIN += OFFLIN
CPPFLAGS_OFFLIN= -DBIN -DTXT -DLFI=lfi -DOL
#
ifdef DIR_OFFLIN
DIR_MASTER += $(DIR_OFFLIN)
CPPFLAGS   += $(CPPFLAGS_OFFLIN)
#VER_SURFEX=SURFEX-7-2-0
#ARCH_XYZ    := $(ARCH_XYZ)-$(VER_MYSRC)

#OBJS_NOCB +=  spll_mode_cover_301_573.o 

$(OBJS0): OPT = $(OPT0) 

endif

##########################################################
#           Source SURFEX                                #
##########################################################
# PRE_BUG TEST !!!
#DIR_SURFEX += ARCH_SRC/bug_surfex
# PRE_BUG TEST !!!
#
DIR_SURFEX += SURFEX
CPPFLAGS_SURFEX= -DASC -DFA=fa
INC_SURFEX = -I$(B)include
#
ifdef DIR_SURFEX
DIR_MASTER += $(DIR_SURFEX)
CPPFLAGS   += $(CPPFLAGS_SURFEX)
INC         += $(INC_SURFEX)
VER_SURFEX=SURFEX-7-2-0
#ARCH_XYZ    := $(ARCH_XYZ)-$(VER_MYSRC)

#OBJS_NOCB +=  spll_mode_cover_301_573.o 

$(OBJS0): OPT = $(OPT0) 

endif
##########################################################
#           Source TRIP                                  #
##########################################################
# PRE_BUG TEST !!!
#DIR_SURFEX += ARCH_SRC/bug_surfex
# PRE_BUG TEST !!!
#
DIR_TRIP += LIB/TRIP
#CPPFLAGS_TRIP=
#
ifdef DIR_TRIP
DIR_MASTER += $(DIR_TRIP)
CPPFLAGS   += $(CPPFLAGS_TROIP)
#VER_SURFEX=SURFEX-7-2-0
#ARCH_XYZ    := $(ARCH_XYZ)-$(VER_MYSRC)

#OBJS_NOCB +=  spll_mode_cover_301_573.o 

$(OBJS0): OPT = $(OPT0) 

endif
##########################################################
#           Source LFI                                   #
##########################################################
DIR_LFIC      += LIB/LFI_COMPRESS/src
PATH_LFIC     += LIB/LFI_COMPRESS/srcc
CPPFLAGS_LFIC = -DSWAPIO -DLINUX -DBIG_endian -Df2cFortran
INC_LFIC      = -I$(B)LIB/LFI_COMPRESS/include

ifdef DIR_LFIC
#
# Management/parametrisation of size of INTEGER ofr file > 16 GO & RECL for LFI
#
DIR_MASTER          += $(DIR_LFIC)
CPPFLAGS            += $(CPPFLAGS_LFIC)
OBJS_LISTE_MASTER   += bitbuff.o ieee_is_nan.o nearestpow2.o
INC                 += $(INC_LFIC)
VPATH               += $(PATH_LFIC)
#VER_NEWLFI=
#ARCH_XYZ    := $(ARCH_XYZ)-$(VER_NEWLFI)
endif
##########################################################
#           Librairie DR_HOOK                            #
##########################################################
DIR_HOOK     +=  LIB/drhook_CY31R2.032
INC_HOOK = -I$(B)LIB/drhook_CY31R2.032
#
ifdef DIR_HOOK
LIBS       += $(DIR_HOOK)/libdrhook.a 
LIBS       += $(DIR_HOOK)/libmpi_serial.a
INC        += $(INC_HOOK)
VPATH      += $(DIR_HOOK)
endif
##########################################################
#           Source XRD                                   #
##########################################################
DIR_XRD += LIB/XRD38/FA
DIR_XRD += LIB/XRD38/FA/mt
DIR_XRD += LIB/XRD38/LFI
DIR_XRD += LIB/XRD38/LFI/mt
DIR_XRD += LIB/XRD38/grib_mf
DIR_XRD += LIB/XRD38/module
DIR_XRD += LIB/XRD38/support
DIR_XRD += LIB/XRD38/utilities
INC_XRD = -I$(B)LIB/XRD38/include -I$(B)LIB/XRD38/FA -I$(B)LIB/XRD38/LFI
#CPPFLAGS_XRD = -fdefault-real-8
#
ifdef DIR_XRD
DIR_MASTER += $(DIR_XRD)
#CPPFLAGS   += $(CPPFLAGS_XRD)
INC        += $(INC_XRD)
endif
##########################################################
#           Source FM                                    #
##########################################################
DIR_FM += LIB/FM
#CPPFLAGS_FM =
#INC_FM=
#
ifdef DIR_FM
DIR_MASTER += $(DIR_FM)
#CPPFLAGS   += $(CPPFLAGS_FM)
INC        += $(INC_FM)
endif
##########################################################
#           Librairie GRIBEX                             #
##########################################################
#ifneq "$(ARCH)" "BG"
# Gribex bypass on BG for the moment
#DIR_GRIBEX     +=  LIB/GRIBEX
#endif
#
#ifdef DIR_GRIBEX
#LIB_GRIBEX     =  $(DIR_GRIBEX)_$(ARCH)/libgribexR64.a
#LIBS          +=    $(LIB_GRIBEX)
#R64_GRIBEX=R64
#endif
##########################################################
#           Librairie GRIBAPI                            #
##########################################################
ifneq "$(ARCH)" "BG"
# Gribapi bypass on BG for the moment
DIR_GRIBAPI?=${SRC_SURFEX}/src/LIB/grib_api-${VERSION_GRIBAPI}
GRIBAPI_PATH?=${DIR_GRIBAPI}-${ARCH}${MNH_INT}
GRIBAPI_INC?=${GRIBAPI_PATH}/include/grib_api.mod
endif
#
ifdef DIR_GRIBAPI
INC_GRIBAPI   ?= -I${GRIBAPI_PATH}/include
LIB_GRIBAPI   ?= -L${GRIBAPI_PATH}/lib -L${GRIBAPI_PATH}/lib64 -lgrib_api_f90 -lgrib_api
INC           += $(INC_GRIBAPI)
LIBS          += $(LIB_GRIBAPI)
VPATH         += $(GRIBAPI_PATH)/include
R64_GRIBAPI=R64
endif
##########################################################
#           Librairie NETCDF                             #
##########################################################
#
# NetCDF  : AUTO install of netcdf-3.6.X on PC linux to avoid problem with compiler
#  
#
ifeq "$(VER_CDF)" "CDFAUTO"
DIR_CDF?=${SRC_SURFEX}/src/LIB/netcdf-${VERSION_CDF}
CDF_PATH?=${DIR_CDF}-${ARCH}${MNH_INT}
CDF_INC?=${CDF_PATH}/include/netcdf.inc
#
INC_NETCDF     ?= -I${CDF_PATH}/include
LIB_NETCDF     ?= -L${CDF_PATH}/lib -L${CDF_PATH}/lib64 -lnetcdf_c++ -lnetcdf
INC            += $(INC_NETCDF)
LIBS           += $(LIB_NETCDF)
endif
#
# NetCDF in SGI ICE
#
ifeq "$(VER_CDF)" "CDFICE"
CDF_PATH?=/opt/software/SGI/netcdf/4.0
INC_NETCDF     ?= -I${CDF_PATH}/include
LIB_NETCDF     ?= -L${CDF_PATH}/lib -lnetcdff  -lnetcdf -i_dynamic 
INC            += $(INC_NETCDF)
LIBS           += $(LIB_NETCDF)
endif
#
# NetCDF in NEC SX
#
ifeq "$(VER_CDF)" "CDFSX"
CDF_PATH?=/SXlocal/pub/netcdf/3.6.1
INC_NETCDF     ?= -I${CDF_PATH}/include
LIB_NETCDF     ?= -L${CDF_PATH}/lib -lnetcdf_c++ -lnetcdf
INC            += $(INC_NETCDF)
LIBS           += $(LIB_NETCDF)
endif
#
ifeq "$(VER_CDF)" "CDFMFSX"
CDF_PATH?=/usr/local/SX/lib/NETCDF_size_t32
INC_NETCDF     ?= -I${CDF_PATH}/include
LIB_NETCDF     ?= -L${CDF_PATH}/lib -lnetcdf
INC            += $(INC_NETCDF)
LIBS           += $(LIB_NETCDF)
endif
#
# NetCDF in AIX S
#
ifeq "$(VER_CDF)" "CDFAIX"
CDF_PATH?=/usr/local/pub/NetCDF/3.6.2
INC_NETCDF     ?= -I${CDF_PATH}/include
LIB_NETCDF     ?= -L${CDF_PATH}/lib -lnetcdf_c++ -lnetcdf
INC            += $(INC_NETCDF)
LIBS           += $(LIB_NETCDF)
endif

#
# Linux with gfortran SUSE10.3
#
ifeq "$(VER_CDF)" "CDFGFOR"
INC_NETCDF     ?=  -I/usr/include
LIB_NETCDF     ?=  -lnetcdf -lnetcdff /usr/lib64/libgfortran.so.2
#LIB_NETCDF     ?=  -lnetcdf -lnetcdff 
INC            += $(INC_NETCDF)
LIBS           += $(LIB_NETCDF)
endif

#
# Linux with netcdf CTI 3.6.3
#
ifeq "$(VER_CDF)" "CDFCTI"
CDF_PATH?=/usr
INC_NETCDF     = -I${CDF_PATH}/include
LIB_NETCDF     = -L${CDF_PATH}/lib64 -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lsz -lz
INC            += $(INC_NETCDF)
LIBS           += $(LIB_NETCDF)
endif

#
# Linux with gfortran SUSE11.1
#
ifeq "$(VER_CDF)" "CDF3GFOR"
CDF_PATH       ?=/opt/netcdf3
INC_NETCDF     ?=  -I${CDF_PATH}/include
LIB_NETCDF     ?=  -L${CDF_PATH}/lib64  -lnetcdf_c++ -lnetcdf
INC            +=  $(INC_NETCDF)
LIBS           +=  $(LIB_NETCDF)
endif

##########################################################
#           Number of NESTED MODEL                       #
##########################################################
NSOURCE=8
##########################################################
#                                                        #
# PROG_LIST : Main program liste to compile              #
#                                                        #
##########################################################
#
ifeq "$(ARCH)" "BG"
PROG_LIST += OFFLINE 
else
PROG_LIST += PGD PREP OFFLINE OI_MAIN SODA SXPOST NCPOST MAIN_CARB_SPINUP MAIN_WOOD_SPINUP
endif
#
ifeq "$(VER_USER)" "FORC"
PROG_LIST += PRE_INPUT_EXPERIMENT
endif
##########################################################
#                                                        #
# LIB_OBJS : Librarie of all *.o                         #
#                                                        #
##########################################################
#
ARCH_XYZ        := $(ARCH_XYZ)-$(OPTLEVEL)
OBJDIR_ROOT     := $(OBJDIR_ROOT)-$(ARCH_XYZ)
LIB_OBJS_ROOT   := $(LIB_OBJS_ROOT)-$(ARCH_XYZ)
#
##########################################################
#                                                        #
# IGNORE_OBJS : some *.o to ignore                       #
#       ---> unused unsupported old routines             #
#                                                        #
##########################################################
#
IGNORE_OBJS += 
IGNORE_DEP_MASTER += 
IGNORE_DEP_MASTER += 

#
#
##########################################################
#                                                        #
#  VPATH_EXCLUDE : Some sources directory to exclude     #
#                                                        #
##########################################################
#
VPATH_EXCLUDE= %/CVS
#



