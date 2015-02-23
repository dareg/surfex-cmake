###
#  arch file for SLES 10-SP3 with gfortran v4.8.2 mpich 3.1.3 & openmp
###
##########################################################
#                                                        #
# Compiler Options                                       #
#                                                        #
##########################################################
#OBJDIR_PATH=/home/escj/azertyuiopqsdfghjklm/wxcvbn/azertyuiopqsdfghjklmwxcvbn
#
OPT_BASE  = -fopenmp -fdefault-real-8 -fdefault-double-8 -g -fno-second-underscore -fpic -fbacktrace -fconvert=swap -Wimplicit-interface -Wimplicit-procedure -Wall -Wextra -Waliasing -Wampersand -Warray-temporaries -Wcharacter-truncation -Wconversion-extra -Wsurprising -Wunderflow -Wno-compare-reals
#
OPT_PERF0 = -O0
OPT_PERF2 = -O2
OPT_CHECK = -fbounds-check -finit-real=nan -ffpe-trap=overflow,zero,invalid
OPT_I8    = -fdefault-integer-8
#
#
# Integer 4/8 option
#
#MNH_INT   ?=I4
LFI_RECL  ?=512
#
ifeq "$(MNH_INT)" "I8"
OPT_BASE         += $(OPT_I8)
LFI_INT           ?=8
MNH_MPI_RANK_KIND ?=8
else
MNH_MPI_RANK_KIND ?=4
LFI_INT           ?=4
endif
#
#
OPT       = $(OPT_BASE) $(OPT_PERF2)
OPT0      = $(OPT_BASE) $(OPT_PERF0)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF2)
#
ifeq "$(OPTLEVEL)" "DEBUG"
OPT       = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT0      = $(OPT_BASE) $(OPT_PERF0) $(OPT_CHECK)
OPT_NOCB  = $(OPT_BASE) $(OPT_PERF0)
endif
#
ifeq "$(VER_MPI)" "NOMPI"
F90= gfortran
F77= gfortran
else
F90= mpif90
F77= mpif77
endif
#
FC = $(F90)
#
F90FLAGS      = $(OPT) -ffree-form -ffree-line-length-none
F77FLAGS      = $(OPT) -ffixed-form
FX90          = $(F77)
FX90FLAGS     = $(OPT) -ffixed-form
#
CC            = gcc
CFLAGS        =
#
LDFLAGS   = -fopenmp -fbacktrace -fuse-ld=gold
#
# preprocessing flags
#
CPP = cpp -P -traditional -Wcomment
#
CPPFLAGS_SURFEX    =
CPPFLAGS_SURCOUCHE = -DMNH_MPI_DOUBLE_PRECISION -DMNH_LINUX -DMNH_MPI_BSEND -DDEV_NULL  -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND)
CPPFLAGS_RAD       =
CPPFLAGS_NEWLFI    = -DSWAPIO -DLINUX -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
CPPFLAGS_MNH       = -DMNH -DAINT=INT -DAMOD=MOD
#
# Gribex flags
#
TARGET_GRIBEX=linux
CNAME_GRIBEX=_gfortran
##########################################################
#                                                        #
# Source of MESONH PACKAGE  Distribution                 #
#                                                        #
##########################################################
#
include Makefile.SURFEX.mk
#
ifeq "$(VER_MPI)" "NOMPI"
CPPFLAGS += -DNOMPI
else
ifeq "$(VER_OASIS)" "mct"
CPPFLAGS += -DSFXOASIS -DTRIPOASIS
endif
endif
#
##########################################################
#                                                        #
# extra VPATH, Compilation flag modification             #
#         systeme module , etc ...                       #
#         external precompiled module librairie          #
#         etc ...                                        #
#                                                        #
##########################################################

#RJ ifneq "$(findstring 8,$(LFI_INT))" ""
#RJ OBJS_I8=spll_NEWLFI_ALL.o
#RJ $(OBJS_I8) : OPT = $(OPT_BASE) $(OPT_PERF2) $(OPT_I8)
#RJ endif
