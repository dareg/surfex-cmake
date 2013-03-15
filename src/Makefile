##########################################################
#                                                        #
#  MACRO : dependance_flist                              #
#          generate the dependance list for ALL fortran  #
#          "*.f*" files                                  #
#          searching in PATH=$(1)                        #
##########################################################
#
dependance_flist = $(notdir $(shell find $(1) -follow -type f \
               -name '*.f' -o -name '*.f90' -o -name '*.fx90' -o -name '*.F' -o -name '*.F90'\
               | grep -v ',v' | sed -e 's/\.f.*/.D/g' -e 's/\.F.*/.D/g' ))
obj_flist        = $(notdir $(shell find $(1) -follow -type f \
               -name 'spll_*.f*' | sed -e 's/\..*/.o/g' )) 
#
##########################################################
#                                                        #
# COMPILER & ARCHITECTURE CONFIGURATION                  #
#                                                        #
##########################################################
#
include Rules.$(ARCH)$(F).mk
#
##########################################################
##########################################################
##########################################################
#                                                        #
#    DONT MODIFY LINES AFTER THIS  !!!!!!!!!!!!!         #         
#                                                        #
##########################################################
##########################################################
##########################################################

##########################################################
#                                                        #
#    AUTOMATIQUE  DEPENDANCE GENERATION !!!!!!!!!!!!!    #
#                                                        #
#   "*.D" files contain dependance deduceed              #
#   form all "*.f*" FORTRAN sources files                #
#                                                        #   
# "use" key word                                         #
#     --> dependance on "*.mod" module files             #
# "module" key word key word                             #
#     --> creation of"*.mod"  module files               #
#                                                        #
#   The "*.D" list of files is contructed from every     #
#   directories of fortran sources :                     #
#       ---> user, master , bugfix , etc...              #
#   with or without RCS version control                  #
#                                                        #
#   The "*.D" list is reconstructed automaticaly         #
#   by the commande "spll" , each time a FORTRAN source  #
#   are modified or added .                              #
#                                                        #
#   All modification are allowed !!!!!                   #
#   adding new subroutines                               #
#       or                                               #
#   adding new modules                                   #       
#                                                        #
#   REM : if during modification, you deleting some      #
#         FORTRAN subroutines you must also deleted the  #
#         "spll_*.o" files created by this one by hand   #
#         to avoid potential link problems               #
#                                                        #
##########################################################

##########################################################
#                                                        #
# DEP_MASTER : Automatique -->                           #
#    Dependance file of MASTER source                    #
#                                                        #
##########################################################
#
OBJDIR_MASTER=$(OBJDIR_ROOT)/MASTER
DEP_MASTER= $(call dependance_flist,$(DIR_MASTER))
#
#  init VPATH depending on compilation phase
#
ifdef DO_DEP_MASTER
VPATH += $(shell find $(DIR_MASTER) -follow -type d )
VPATH += $(OBJDIR_MASTER)
endif
#
ifdef DO_COMP_MASTER
VPATH += $(OBJDIR_MASTER) $(OBJDIR_MASTER)/MOD
endif
#
##########################################################
#                                                        #
# DEP_USER : Automatique -->                             #
#    Dependance file of MY source                        #
#                                                        #
##########################################################
#
ifdef VER_USER
OBJDIR_USER=$(OBJDIR_ROOT)/$(VER_USER)
#
#  init VPATH depending on compilation phase
#
ifdef DO_DEP_USER
DEP_USER= $(call dependance_flist,$(DIR_USER))
VPATH += $(shell find $(DIR_USER) -follow -type d )
VPATH += $(OBJDIR_USER)
endif
#
ifdef DO_COMP_USER
VPATH += $(OBJDIR_USER) $(OBJDIR_USER)/MOD $(OBJDIR_MASTER) $(OBJDIR_MASTER)/MOD
INC += -I $(OBJDIR_MASTER)/MOD
endif
#
endif
##########################################################
#                                                        #
# VPATH ADJUSTEMENT :                                    #
#                                                        #
##########################################################
#
#
#VPATH := $(filter-out $(VPATH_EXCLUDE),$(VPATH))
#
##########################################################
#                                                        #
# OBJS : all *.o to compile                              #
#                                                        #
##########################################################
#
# MASTER 
#
OBJS_LISTE_MASTER :=  $(OBJS_LISTE_MASTER) $(notdir $(shell find $(OBJDIR_MASTER) -follow -type f -name "spll_*.f*" | sed -e 's/\..*/.o/g' )) 
#OBJS_LISTE_MASTER :=  $(OBJS_LISTE_MASTER) $(obj_flist $(OBJDIR_MASTER)) 
#
OBJS_LISTE_MASTER := $(sort $(filter-out $(IGNORE_OBJS),$(OBJS_LISTE_MASTER))) 
#
LIB_MASTER=$(OBJDIR_MASTER)/$(LIB_OBJS_ROOT)-MASTER.a
#
# USER
#
#ifdef DO_COMP_USER
OBJS_LISTE_USER :=  $(OBJS_LISTE_USER) $(notdir $(shell find $(OBJDIR_USER) -follow -type f -name "spll_*.f*" | sed -e 's/\..*/.o/g' )) 
#OBJS_LISTE_USER :=  $(OBJS_LISTE_USER) $(obj_flist $(OBJDIR_USER)) 
#
OBJS_LISTE_USER := $(sort $(filter-out $(IGNORE_OBJS), $(OBJS_LISTE_USER) $(OBJS_LISTE_MASTER) )) 
#
LIB_USER=$(OBJDIR_USER)/$(LIB_OBJS_ROOT)-$(VER_USER).a
#endif
##########################################################
#                                                        #
# DEP_ALL : All dependance *.D to generate               #
#                                                        #
##########################################################
#
DEP_ALL_MASTER = $(sort $(filter-out $(IGNORE_DEP_MASTER),$(DEP_MASTER)) )
#
DEP_ALL_USER   = $(sort $(filter-out $(IGNORE_DEP_USER)  ,$(DEP_USER)) )
#

.PHONY : clean all objdir objdirmaster depall objs lib prog print

.DELETE_ON_ERROR:  

.SUFFIXES:
.SUFFIXES:  .o .f90 .f .fx90 .mod .D .a

.INTERMEDIATE:   $(LIB_MASTER)


##########################################################
#                                                        #
# MASTER RULES                                           #
#                                                        #
##########################################################
ifeq "$(VER_CDF)" "CDFAUTO"
all : drhook ecoclimap gribapi cdf progmaster
else
all : drhook ecoclimap gribapi progmaster
endif

objdirmaster :  $(OBJDIR_MASTER)/.dummy

$(OBJDIR_MASTER)/.dummy :
	-mkdir -p $(OBJDIR_MASTER)/.dummy
	-mkdir -p $(OBJDIR_MASTER)/MOD

depmaster : objdirmaster
	@$(MAKE) DO_DEP_MASTER=YES OBJDIR=$(OBJDIR_MASTER) filedepallmaster

filedepallmaster :  $(DEP_ALL_MASTER)
	find  $(OBJDIR_MASTER) -follow -name "*.D" >  $(OBJDIR_MASTER)/filemaster
	rm -f  $(OBJDIR_MASTER)/filedepallmaster
	cat $(OBJDIR_MASTER)/filemaster  | xargs cat  >>  $(OBJDIR_MASTER)/filedepallmaster

libmaster :  depmaster
	@$(MAKE)  -I $(OBJDIR_MASTER)  DO_COMP_MASTER=YES OBJDIR=$(OBJDIR_MASTER) objmaster


bibmaster :  libmaster
	@$(MAKE)  -I $(OBJDIR_MASTER)  DO_COMP_MASTER=YES OBJDIR=$(OBJDIR_MASTER) $(LIB_MASTER)
	mv $(LIB_MASTER)  $(subst lib-,bib-,$(LIB_MASTER))


objmaster :  $(OBJS_LISTE_MASTER)

$(LIB_MASTER) :  $(OBJS_LISTE_MASTER)
	cd $(OBJDIR_MASTER) ; find . -follow -name '*.o'| xargs $(AR) r $(LIB_MASTER) 

progmaster : libmaster
	@$(MAKE) -I $(OBJDIR_MASTER)  DO_COMP_MASTER=YES  \
	OBJDIR=$(OBJDIR_MASTER) \
	LIB_MNH=$(LIB_MASTER) \
	PROG_DIR=$(OBJDIR_MASTER) \
	$(PROG_LIST) cleanlibmaster

installmaster : libmaster
	@$(MAKE) -I $(OBJDIR_MASTER)  DO_COMP_MASTER=YES  OBJDIR=$(OBJDIR_MASTER) install

cleanlibmaster :
	rm -f  $(LIB_MASTER) 

cleanmaster : gribex_clean 
	test -d $(OBJDIR_MASTER) && rm -fr  $(OBJDIR_MASTER)

cleanobjmaster :  
	test -d $(OBJDIR_MASTER) && find $(OBJDIR_MASTER) -follow -name "*.o" -o -name "*.mod" | xargs rm 

cleanprogmaster :
	test -d $(OBJDIR_MASTER) && cd $(OBJDIR_MASTER) && rm -f $(PROG_LIST)

##########################################################
#                                                        #
#  USER RULES                                            #
#                                                        #
##########################################################

user : proguser

objdiruser :  $(OBJDIR_USER)/.dummy

$(OBJDIR_USER)/.dummy :
	-mkdir -p $(OBJDIR_USER)/.dummy
	-mkdir -p $(OBJDIR_USER)/MOD

depuser : objdiruser
	@$(MAKE) DO_DEP_USER=YES OBJDIR=$(OBJDIR_USER) filedepalluser

filedepalluser : $(DEP_ALL_USER)
	find  $(OBJDIR_USER) -follow -name "*.D" >  $(OBJDIR_USER)/fileuser 
	cat $(OBJDIR_USER)/fileuser | sed -e "s!$(OBJDIR_USER)!$(OBJDIR_MASTER)!g"  >  $(OBJDIR_USER)/filenotmaster
	grep -F -v -f   $(OBJDIR_USER)/filenotmaster  $(OBJDIR_MASTER)/filemaster >  $(OBJDIR_USER)/filemaster
	rm -f  $(OBJDIR_USER)/filedepalluser
	cat $(OBJDIR_USER)/fileuser $(OBJDIR_USER)/filemaster | xargs cat >>  $(OBJDIR_USER)/filedepalluser


libuser :   depuser 
	@$(MAKE) -I $(OBJDIR_USER) DO_COMP_USER=YES OBJDIR=$(OBJDIR_USER) objuser

objuser :  $(OBJS_LISTE_USER)

$(LIB_USER) :  $(OBJS_LISTE_USER)
	rm -f $(LIB_USER)
	cd $(OBJDIR_MASTER) ; find . -name "*.o" | xargs $(AR) r $(LIB_USER) 
	cd $(OBJDIR_USER)   ; find . -name "*.o" | xargs $(AR) r $(LIB_USER) 

proguser : libuser
	@$(MAKE) -I $(OBJDIR_USER)  DO_COMP_USER=YES \
	OBJDIR=$(OBJDIR_USER) \
	LIB_MNH=$(LIB_USER) \
	PROG_DIR="$(OBJDIR_USER) $(OBJDIR_MASTER)" \
	$(PROG_LIST) # cleanlibuser

cleanlibuser :
	rm -f  $(LIB_USER) 

installuser : 
	@$(MAKE)  -I $(OBJDIR_USER) DO_COMP_USER=YES  OBJDIR=$(OBJDIR_USER) install	

cleanproguser :
	test -d $(OBJDIR_USER) && cd $(OBJDIR_USER) && rm -f $(PROG_LIST)

##########################################################
#                                                        #
# EXTRA LIB : GRIBEX                                     #
#                                                        #
##########################################################
ifneq "$(findstring 64,$(shell uname -m))" ""
A64=A64
endif
#$(LIB_GRIBEX) : libmaster
$(LIB_GRIBEX) : 
ifneq "$(ARCH)" "SX8"
	- [ ! -d $(DIR_GRIBEX)_$(ARCH) ] && cp -Rp  $(DIR_GRIBEX) $(DIR_GRIBEX)_$(ARCH) 
	cd $(DIR_GRIBEX)_$(ARCH) && $(MAKE) -j 1 ARCH=$(TARGET_GRIBEX) R64=$(R64_GRIBEX) CNAME=$(CNAME_GRIBEX) A64=$(A64)
else
ifneq "$(findstring brodie,$(shell uname -n))" ""
#brodie
	- [ ! -d $(DIR_GRIBEX)_$(ARCH) ] && mkdir -p ${workdir}/$(DIR_GRIBEX)_$(ARCH) \
	&& ln -s ${workdir}/$(DIR_GRIBEX)_$(ARCH)  $(DIR_GRIBEX)_$(ARCH)
	cp /home/rech/mnh/rmnh007/aeronec/mesonh/binaries/libemosR64.a $(DIR_GRIBEX)_$(ARCH)/libgribexR64.a
else
#tori & yuki
	- [ ! -d $(DIR_GRIBEX)_$(ARCH) ] && mkdir -p $(DIR_GRIBEX)_$(ARCH)
	cp /usr/local/SX/lib/libgribex.a $(DIR_GRIBEX)_$(ARCH)/libgribexR64.a
endif
endif

gribex_clean :
	- [  -d $(DIR_GRIBEX)_$(ARCH) ] && rm -fr $(DIR_GRIBEX)_$(ARCH) 
##########################################################
#                                                        #
# EXTRA LIB : GRIBAPI                                    #
#                                                        #
##########################################################
gribapi : $(GRIBAPI_INC)
$(GRIBAPI_INC) : 
	cd ${DIR_GRIBAPI} && ./configure --disable-jpeg --prefix=${GRIBAPI_PATH} FC="$(FC)" && $(MAKE) -j 1 clean && $(MAKE) -j 1 && $(MAKE) -j 1 install	

gribapi_clean :
	- [  -d ${GRIBAPI_PATH} ] && rm -fr ${GRIBAPI_PATH} 
##########################################################
#                                                        #
# EXTRA LIB : NETCDF                                     #
#                                                        #
##########################################################
ifeq "$(VER_CDF)" "CDFAUTO"
#
cdf : $(CDF_INC)
$(CDF_INC) : 
	cd ${DIR_CDF} && ./configure --prefix=${CDF_PATH} FC="$(FC)" FCFLAGS="$(OPT_BASE)" CPPFLAGS="-DpgiFortran" && make && make install

cleanmaster : cleancdf
cleancdf :
	- [ -d ${CDF_PATH} ] && rm -fr  ${CDF_PATH}
endif
##########################################################
#                                                        #
# EXTRA LIB : DR_HOOK                                    #
#                                                        #
##########################################################
ifeq "$(F90)" "gfortran"
ARCHOOK="linuxgfortran"
else 
ifeq "$(F90)" "mpif90"
ARCHOOK="linuxmpif90"
else
ifeq "$(F90)" "LXpgi"
ARCHOOK="linux"
else
ARCHOOK=$(ARCH)
ifeq "$(VER_MPI)" "NOMPI"
ARCHOOK=$(ARCH)_nompi
endif
endif
endif
endif

drhook : $(DIR_HOOK)/parkind1.mod
$(DIR_HOOK)/parkind1.mod : 
	cd ${DIR_HOOK} && ./configure ARCH=$(ARCHOOK) -a $(ARCHOOK) && make ARCH=$(ARCHOOK) drhook 2>&1 | tee log.make 

cleandrhook : 
	cd ${DIR_HOOK} && rm -f *.mod *.a
##########################################################
#                                                        #
# EXTRA LIB : ECOCLIMAP                                  #
#                                                        #
##########################################################
ecoclimap:${SRC_SURFEX}/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin


ifeq "$(F90)" "ifort"
${SRC_SURFEX}/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin : 
	cd ${SRC_SURFEX}/MY_RUN/ECOCLIMAP && $(F90) $(F90FLAGS) convert_ecoclimap_param.f90 && ./a.out
else
ifeq "$(F90)" "mpiifort"
${SRC_SURFEX}/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin : 
	cd ${SRC_SURFEX}/MY_RUN/ECOCLIMAP && $(F90) $(F90FLAGS) convert_ecoclimap_param.f90 && ./a.out
else
ifeq "$(F90)" "pgf90"
${SRC_SURFEX}/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin : 
	cd ${SRC_SURFEX}/MY_RUN/ECOCLIMAP && $(F90) -byteswapio convert_ecoclimap_param.f90  && ./a.out
else
ifeq "$(ARCH)" "SX8"
${SRC_SURFEX}/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin :
	cd ${SRC_SURFEX}/MY_RUN/ECOCLIMAP && $(F90) convert_ecoclimap_param.f90 && qsub job_ecoclimap_sx8
else
ifeq "$(ARCH)" "AIX64"
${SRC_SURFEX}/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin :
	cd ${SRC_SURFEX}/MY_RUN/ECOCLIMAP && $(F90) convert_ecoclimap_param.f90 && llsubmit job_ecoclimap_aix64
else
${SRC_SURFEX}/MY_RUN/ECOCLIMAP/ecoclimapI_covers_param.bin : 
	cd ${SRC_SURFEX}/MY_RUN/ECOCLIMAP && $(F90) -fconvert=swap convert_ecoclimap_param.f90 && ./a.out
endif
endif
endif
endif
endif

cleanecoclimap : 
	cd ${SRC_SURFEX}/MY_RUN/ECOCLIMAP && rm -f *.bin

##########################################################
#                                                        #
# PROGRAM RULES                                          #
#                                                        #
##########################################################
ifdef PROG_LIST

prog :   lib $(LIB_GRIBEX)
	@$(MAKE) -I $(OBJDIR) DEP=YES  $(PROG_LIST)  

install : $(patsubst %,%-${ARCH_XYZ},$(PROG_LIST))

$(PROG_LIST) : OBJ_PROG=$(shell find $(PROG_DIR) -follow -type f -name "spll_*.f*" \
              | xargs grep -l -E -i "^[[:space:]]*program *$@" |  sed -e 's/\.f.*/.o/g' | head -1 )  

$(PROG_LIST) :  $(LIB_MNH)  $(LIB_GRIBEX)
	$(F90) $(LDFLAGS) -o $(OBJDIR)/$@ $(OBJ_PROG) $(LIB_MNH) $(LIBS)

%-${ARCH_XYZ} :%
ifndef DO_COMP_USER
	ln -sf $(OBJDIR)/$(*F) $(SRC_SURFEX)/exe/$(*F)-${ARCH_XYZ}
else
#	ln -sf $(OBJDIR)/$(*F) $(SRC_SURFEX)/exe/$(*F)-${ARCH_XYZ}-${VER_USER}
	ln -sf $(OBJDIR)/$(*F) $(SRC_SURFEX)/exe/$(*F)${XYZ}
endif
	rm -f $(SRC_SURFEX)/exe/$(*F)_${ARCH_XYZ}~

##########################################################
#                                                        #
# EXAMPLES RULES                                         #
#                                                        #
##########################################################

PWD_KTEST=$(SRC_SURFEX)/MY_RUN/KTEST
EXAMPLES=001_2Drelief 002_3Drelief 003_KW78 004_Reunion 007_16janvier
EXAMPLES_NCL=004_Reunion 007_16janvier


examples : $(EXAMPLES)

ncl1_examples :  $(patsubst %,%_ncl1,$(EXAMPLES_NCL))

ncl2_examples :  $(patsubst %,%_ncl2,$(EXAMPLES_NCL))

examplesclean :  $(patsubst %,%_clean,$(EXAMPLES))

$(EXAMPLES) :
	cd $(PWD_KTEST)/$@ ; run_$@

%_ncl1 : 
	cd $(PWD_KTEST)/$(*F) ; run_$(*F) all_ncl_extractdia

%_ncl2 : 
	cd $(PWD_KTEST)/$(*F) ; run_$(*F) all_ncl_lfi2cdf

%_clean : 
	cd $(PWD_KTEST)/$(*F) ; $(MAKE) clean

EXAMPLES_CHIMIE=009_ICARTT 011_KW78CHEM

chimie : $(EXAMPLES_CHIMIE)

$(EXAMPLES_CHIMIE) :
	cd $(PWD_KTEST)/$@ ; run_$@
#
# rules to remove spll* source on MASTER directory
#
purgelist : $(patsubst %,purge_%,$(notdir $(PL)))
purge_%.f90 :
	cd $(OBJDIR_MASTER) ; grep -l "depfile:$(*F).D" *.f* | xargs rm -f ; rm -f $(*F).D ; rm -f filedepallmaster filemaster 
purge_%.fx90 :
	cd $(OBJDIR_MASTER) ; grep -l "depfile:$(*F).D" *.f* | xargs rm -f ; rm -f $(*F).D ; rm -f filedepallmaster filemaster 
purge_%.f :
	cd $(OBJDIR_MASTER) ; grep -l "depfile:$(*F).D" *.f* | xargs rm -f ; rm -f $(*F).D ; rm -f filedepallmaster filemaster 

#
# rules to remove spll* source on USER directory
#
purgelistu : $(patsubst %,purgeu_%,$(notdir $(PL)))
purgeu_%.f90 :
	cd $(OBJDIR_USER) ; grep -l "depfile:$(*F).D" *.f* | xargs rm -f ; rm -f $(*F).D ; rm -f filedepalluse fileuser
purgeu_%.fx90 :
	cd $(OBJDIR_USER) ; grep -l "depfile:$(*F).D" *.f* | xargs rm -f ; rm -f $(*F).D ; rm -f filedepalluse fileuser
purgeu_%.f :
	cd $(OBJDIR_USER) ; grep -l "depfile:$(*F).D" *.f* | xargs rm -f ; rm -f $(*F).D ; rm -f filedepalluse fileuser

stat_f90 :
	cd $(OBJDIR_MASTER) ; find . -name "*.f90" | wc -l
stat_o :
	cd $(OBJDIR_MASTER) ; find . -name "*.o" | wc -l

endif





cleanuser :  
	test -d $(OBJDIR_USER) && rm -fr  $(OBJDIR_USER)

clean :  examplesclean gribex_clean
	rm -fr  $(OBJDIR_ROOT)
cleanprog :
	cd $(OBJDIR_ROOT) ; rm -f $(PROG_LIST)

print :
	@echo $(VAR)=$($(VAR))

ARFLAGS=r
#(%.o): %.o
#	$(AR) $(ARFLAGS)   $(OBJDIR)/$(@F) $(*F).o
#	$(RM) $(*F).o

#
#  GENERATION OF DEPENDANCE FILES : *.D
# 

%.D:%.F90
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(OBJDIR)/$(*F).f90 
	spll $(NSOURCE) $(OBJDIR)/$(*F).f90 
	rm -f $(OBJDIR)/$(*F).f90 

%.D:%.f90
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(OBJDIR)/$(*F).f90 
	spll $(NSOURCE) $(OBJDIR)/$(*F).f90 
	rm -f $(OBJDIR)/$(*F).f90 

%.D:%.fx90
	$(CPP)  $(INC) $(CPPFLAGS) $< > $(OBJDIR)/$(*F).fx90
	spll $(NSOURCE) $(OBJDIR)/$(*F).fx90 
	rm -f $(OBJDIR)/$(*F).fx90
 
%.D:%.F
	$(CPP) $(INC) $(CPPFLAGS) $< >  $(OBJDIR)/$(*F).f
	spll $(NSOURCE) $(OBJDIR)/$(*F).f
	rm -f $(OBJDIR)/$(*F).f

%.D:%.f
	$(CPP)  $(INC) $(CPPFLAGS) $< > $(OBJDIR)/$(*F).f
	spll $(NSOURCE) $(OBJDIR)/$(*F).f 
	rm -f $(OBJDIR)/$(*F).f
 
#
# GENERATION OF OBJET FILES : *.o
#

%.o:%.f90
	$(F90) -I$(OBJDIR) $(INC) -c $(F90FLAGS) $< 
	-mv  $(*F).o $(OBJDIR)/. || echo OK $(*F).o

%.o:%.f
	$(F77) -I$(OBJDIR) $(INC) -c $(F77FLAGS)  $<  
	-mv  $(*F).o $(OBJDIR)/. || echo OK $(*F).o

%.o:%.fx90
	cp $< $(OBJDIR)/$(*F).f 
	$(FX90) -I$(OBJDIR) $(INC) -c $(FX90FLAGS) $(OBJDIR)/$(*F).f   
	-mv  $(*F).o $(OBJDIR)/. || echo OK $(*F).o

%.o : %.c
	$(CC) $(INC) $(CFLAGS) $(CPPFLAGS) -c $< -o $(OBJDIR)/$(*F).o
	-mv  $(*F).o $(OBJDIR)/. || echo OK $(*F).o

ifeq "$(DO_COMP_MASTER)" "YES"
include  filedepallmaster
endif 

ifeq "$(DO_COMP_USER)" "YES"
include  filedepalluser
endif 


