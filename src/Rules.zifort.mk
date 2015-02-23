###
#  arch file for SLES 10-SP3 with ifort v12.1.6 mpich 3.1.3 & openmp
###
##########################################################
#                                                        #
# Compiler Options                                       #
#                                                        #
##########################################################
#OBJDIR_PATH=/home/escj/azertyuiopqsdfghjklm/wxcvbn/azertyuiopqsdfghjklmwxcvbn
#
#RJ: if using any interrupt signal watchers, dr_hook, mpi, valgrind or ifort itself do not forget to add:
#RJ -fpe0 -fp-model precise -assume ieee_fpe_flags ; specially ieee_fpe_flags one on mixed MPI/OMP
#
#OPT_BASE  = -openmp -openmp-threadprivate=compat -r8 -g -u -assume nosource_include -assume byterecl -fpic -traceback -fp-model precise -assume ieee_fpe_flags -convert big_endian
OPT_BASE  = -openmp -openmp-threadprivate=compat -r8 -g -u -assume nosource_include -assume byterecl -traceback -fp-model precise -assume ieee_fpe_flags -convert big_endian
#
OPT_PERF0 = -O0 -fpe0 -ftz
OPT_PERF2 = -O2 -fpe0 -ftz
OPT_CHECK = -fp-stack-check -ftrapuv -fpe3 -fp-speculation=strict -check all
OPT_I8    = -i8
#
# Integer 4/8 option
#
MNH_INT   ?=I4
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
F90= ifort
F77= ifort
else
F90= mpif90
F77= mpif77
endif
#
FC = ifort
#
F90FLAGS  = $(OPT) -free
F77FLAGS  = $(OPT) -nofree
FX90      = $(F77)
FX90FLAGS = $(OPT) -nofree
#
CC        = gcc
CFLAGS    =
#
LDFLAGS   =  -Wl,-warn-once -openmp -openmp-threadprivate=compat -traceback
#
# preprocessing flags
#
CPP = cpp -P -traditional -Wcomment
#
CPPFLAGS_SURFEX    =
CPPFLAGS_SURCOUCHE = -DMNH_MPI_DOUBLE_PRECISION -DMNH_LINUX -DMNH_MPI_BSEND -DDEV_NULL -DMNH_MPI_RANK_KIND=$(MNH_MPI_RANK_KIND)
CPPFLAGS_RAD       =
CPPFLAGS_NEWLFI    = -DSWAPIO -DLINUX -DLFI_INT=${LFI_INT} -DLFI_RECL=${LFI_RECL}
CPPFLAGS_MNH       = -DMNH
#
# Gribex flags
#
TARGET_GRIBEX=linux
CNAME_GRIBEX=_ifort
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
