# Makefile for compiling COSMO
# fuo, 07.03.2008, oliver.fuhrer@meteoswiss.ch
#
# Description: This makefile automagically searches for source files in the src
#              directory and creates dependencies. Afterwards files are compiled
#              and linked together.
#              Suitable for both, lm and int2lm
# History:
#   fuo  07.03.2008   first release


### Macros ###########################

# define shell
SHELL      := /bin/sh

# directory definitions
SRCDIR     := src
OBJDIR     := obj
BINDIR     := bin
ifeq (0,${MAKELEVEL})
  ROOT      := $(shell pwd)
  VPATH     := .:$(ROOT)/$(SRCDIR)
endif

# define targets
TARGETS    := \
  extpar_alb_to_buffer.new \
  extpar_aot_to_buffer.new \
  extpar_consistency_check.new \
  extpar_cru_to_buffer.new \
  extpar_flake_to_buffer.new \
  extpar_landuse_to_buffer.new \
  extpar_ndvi_to_buffer.new \
  extpar_soil_to_buffer.new \
  extpar_topo_to_buffer.new \
  extpar_ahf_to_buffer.new \
  extpar_isa_to_buffer.new \
  extpar_sgsl_to_buffer.new 


# generate list of source files
-include $(ROOT)/Srcfiles
ifndef SRC
  $(error Could not load Srcfiles file)
endif

# generate list of object files
-include $(ROOT)/Objfiles
ifndef OBJ
  $(error Could not load Objfiles file)
endif

# dynamically generated dependency file
DEPF       := .depend
IGN        := --ignore netcdf --ignore omp_lib --ignore iso_fortran_env

# if machine is not defined, try to determine it
ifndef MACH
  MACH       := $(shell hostname | sed 's/[-0-9]//g')
endif

# select machine dependent stuff
-include $(ROOT)/Options.$(MACH)
ifndef F90
  $(error Host $(MACH) not known. Must create platform specific Options.$(MACH) file)
endif

# handle options
ifdef VERBOSE
  FFLAGS  += $(VFLAGS)
endif
ifdef DEBUG
  FFLAGS  += $(FDBG)
  LIB     += $(DBGL)
  INC     += $(DBGI)
endif
ifdef OPT
  FFLAGS  += $(FOPT)
  LIB     += $(OPTL)
  INC     += $(OPTI)
endif

export ROOT VPATH MACH VERBOSE DEBUG OPT

### Phony targets ###########################

.PHONY : default info depend opt debug clean

default : opt

depend :
	@echo "generating dependencies"
	@echo $(SRC)
	@$(ROOT)/bin/sfmakedepend --case down --longpath $(INC) $(IGN) --file $(ROOT)/$(DEPF) $(SRC)

info :
	@echo "generating compile information"
	@-rm -f $(ROOT)/.fconfig
	@echo "Target           : $(TARGET)" > $(ROOT)/.fconfig
	@echo "Compiler command : $(F90)" > $(ROOT)/.fconfig
	@echo "Compiler version : "`$(F90) -V 2>/dev/null | grep pgf` >> $(ROOT)/.fconfig
	@echo "Compiler includes: $(INC)" >> $(ROOT)/.fconfig
	@echo "Compiler flags   : $(PFLAGS) $(FFLAGS1)" >> $(ROOT)/.fconfig
	@echo "Linker command   : $(LD)" >> $(ROOT)/.fconfig
	@echo "Linker version   : "`$(LD) -V 2>/dev/null | grep pgf` >> $(ROOT)/.fconfig
	@echo "Linker flags     : $(LFLAGS) $(FFLAGS1)" >> $(ROOT)/.fconfig
	@echo "Linker libraries : $(LIB)" >> $(ROOT)/.fconfig
	@$(ROOT)/$(BINDIR)/gen_info.sh $(ROOT)/.fconfig $(ROOT)/$(SRCDIR)
	@-rm -f $(ROOT)/.fconfig

opt :
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile OPT=1 info depend
	@for target in $(TARGETS) ; do \
     echo "generating target $$target" ; \
	   $(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile OPT=1 $$target ; \
   done

debug :
	@$(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile DEBUG=1 info depend
	@for target in $(TARGETS) ; do \
     echo "generating target $$target" ; \
	   $(MAKE) -C $(OBJDIR) -f $(ROOT)/Makefile DEBUG=1 $$target ; \
   done

clean :
	-rm -f $(DEPF) $(DEPF).old $(OBJDIR)/* $(BINDIR)/extpar_*.new

### Suffix Rules ###########################

.SUFFIXES: .new .o .mod .f90

%.new : $(notdir %.o) $(OBJ)
	@echo "linking $@"
	@$(LD) $(LFLAGS) $(FFLAGS) $(INC) $(patsubst %.new,%.o,$(notdir $@)) $(OBJ) $(LIB) -o $(ROOT)/$(BINDIR)/$@

%.o : %.f90
	@echo "compiling $(patsubst %.o,%.f90,$(notdir $@))"
	@$(F90) -c $(PFLAGS) $(FFLAGS) $(INC) -o $@ $(ROOT)/$(SRCDIR)/$(patsubst %.o,%.f90,$(notdir $@))

# include dynamically generated dependency file
-include $(ROOT)/$(DEPF)

# goodbye earthling!
