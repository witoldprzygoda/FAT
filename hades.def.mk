###############################################################################
#
#  Makefile containing the global definitions for all Hydra modules.
#
#  See the global file Makefile for additional information.
#
###############################################################################


ifndef ROOTSYS
   $(error ROOTSYS not set!)
endif


GCCVERSION := $(shell gcc -dumpversion | cut -f1 -d.)

GCC8LINKADD = 

ifeq ($(shell test $(GCCVERSION) -gt 5; echo $$?),0)
    GCC8LINKADD = -no-pie
endif


### Hydra specific flags ######################################################

DEBUG_LEVEL ?= 0

ifeq ($(strip $(USES_DEBUG)),yes)
      DEBUGOPT = -g
endif

ifeq ($(strip $(USES_DEBUG)),)
      DEBUGOPT = -g
endif


ifeq ($(strip $(USES_OPT)),O3)
      OPTLEVEL = -O3
endif

ifeq ($(strip $(USES_OPT)),)
      OPTLEVEL = -O2
endif

HYDRA_LIBS_DEP  = -lAlignment -lQA -lSts -lFRpc -lForward -lParticle -lKalman -lPionTracker -lMdcTrackG -lMdcTrackD
HYDRA_LIBS_DEP += -lWall -lShowerUtil -lShower -lTof -lMdc -lRich -lEmc -lRpc -liTof -lStart -lTools
ifeq ($(strip $(USES_RFIO)),yes)
    HYDRA_LIBS += -lRFIOtsm
endif
HYDRA_LIBS_DEP += -lHydra
ifdef PLUTODIR
HYDRA_LIBS_DEP += -lPluto
endif
ifdef DABCSYS
HYDRA_LIBS_DEP += -lDabcHydra
endif


HYDRA_LIBS = $(HYDRA_LIBS_DEP)

BASE_PCM=Correlation DataSource DataStruct EventHandling EvtServer GeantUtil Geometry RuntimeDb Util Haddef
#HYDRA_SO:= $(patsubst -l%,lib%.so,$(HYDRA_LIBS))
HYDRA_SO:= $(patsubst -l%,-rml lib%.so,$(HYDRA_LIBS_DEP))

ifeq ($(strip $(USES_ORACLE)),yes)
    HYDRA_LIBS += -lOra
endif

### variables related to third party software #################################


X11_INC_DIRS := /usr/X11/include
X11_LIB_DIRS := /usr/lib/x86_64-linux-gnu
X11_LIBS     := -lXmu -lXt -lXpm -lXext -lX11

KUIPC             := $(CERN_ROOT)/bin/kuipc
CERNLIB_INC_DIRS  := $(CERN_ROOT)/include/cfortran $(CERN_ROOT)/include
CERNLIB_LIB_DIRS  := $(CERN_ROOT)/lib
CERNLIB_CPP_FLAGS := -Wno-unused -Df2cFortran
CERNLIB_FPP_FLAGS := -DCERNLIB_TYPE -std=legacy
CERNLIB_LIBS      := -lgeant321 -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib
CERNLIB_LIBS      += -lkernlib -llapack3 -lblas -lm -lnsl -ldl -lcrypt

ifeq ($(strip $(USES_GFORTRAN)),yes) 
   CERNLIB_LIBS      += -lgfortran
endif
ifeq ($(strip $(USES_GFORTRAN)),no) 
   CERNLIB_LIBS      += -lg2c
endif


ROOTCINT       := $(ROOTSYS)/bin/rootcint
ROOTCONFIG     := $(ROOTSYS)/bin/root-config
ROOT_INC_DIRS  := $(shell $(ROOTCONFIG) --incdir)
ROOT_INC_DIRS  += $(shell $(ROOTCONFIG) --incdir)/ROOT
ROOT_LIB_FLAGS := $(shell $(ROOTCONFIG) --glibs --evelibs)
ROOT_LIB_DIRS  := $(subst -L,,$(filter -L%,$(ROOT_LIB_FLAGS)))
ROOT_LIBS      := $(filter-out -L%,$(ROOT_LIB_FLAGS))
ROOT_LIBS      += -lMinuit -lProof -lXMLParser

ifeq ($(strip $(shell gcc -dumpversion)),2.95.4) 
   ROOT_LIB_DIRS += /usr/local/globus2.0/lib
endif

ifeq ($(strip $(USES_ORACLE)),yes)
ORACLE_INC_DIRS := $(ORACLE_HOME)/precomp/public	\
                   $(ORACLE_HOME)/rdbms/public		\
                   $(ORACLE_HOME)/plsql/include         \
                   $(ORACLE_HOME)/sdk/include
ORACLE_LIB_DIRS := $(ORACLE_HOME)/lib
ORACLE_LIBS     := -lclntsh $(shell cat $(ORACLE_HOME)/lib/sysliblist) -ldl -lm
ORACLE_CPP_FLAGS := -DORACLE_SUPPORT
endif

ifeq ($(strip $(USES_RFIO)),yes)
RFIO_LIBS       := -lrawapiclin
RFIO_LIB_DIRS   := $(ADSM_BASE_NEW)/lib
RFIO_INC_DIRS   := $(ADSM_BASE_NEW)/include
RFIO_CPP_FLAGS  := -DRFIO_SUPPORT
endif

ifeq ($(strip $(USES_DABC)),yes)
DABC_LIBS       := -lDabcBase -lDabcMbs -lDabcHadaq
DABC_LIB_DIRS   := $(DABCSYS)/lib
DABC_INC_DIRS   := $(DABCSYS)/include
DABC_CPP_FLAGS  := -DDABC_SUPPORT
endif

### variables related to make and the makefiles ###############################


SHELL   := bash
CP      := cp
RM      := rm
SED     := sed
RMDIR   := rmdir
SLEEP   := sleep
ECHO    := echo
INSTALL := install -p --backup=none

# these statement substitute characters with special meaning for Make
empty :=
space := $(empty) $(empty)
comma := ,

# this is something like a make internal "which" command
pathsearch = $(firstword $(wildcard $(addsuffix /$(1),$(subst :, ,$(PATH)))))

# this removes non existing files or directories from a list
cleanlist = $(strip $(foreach file,$(1),$(wildcard $(file))))

### the following functions work on list of files with and without paths

# <new list> = $(call remove-files,<remove these files>,<from this list>)
remove-files = $(filter-out $(join $(dir $(1)),$(notdir $(1))),$(2))

# <new list> = $(call append-files,<add these files>,<to this list>)
append-files = $(2) $(join $(dir $(1)),$(notdir $(1)))

# <new list> = $(call prefix-with-path,<list of files>)
prefix-with-path = $(join $(dir $(1)),$(notdir $(1)))


### check, if necessary but externally defined variables were set #############


ifdef MODULES
   ifeq ($(MAKELEVEL),0)
      ifdef LIB_NAME
         $(error Illegal definition LIB_NAME in this global makefile!)
      endif
      ifdef APP_NAME
         $(error Illegal definition APP_NAME in this global makefile!)
      endif
      ifdef SOURCE_FILES
         $(error Illegal definition SOURCE_FILES in this global makefile!)
      endif
   endif
endif
ifdef LIB_NAME
   ifeq ($(origin MODULES),file)
      $(error Illegal definition MODULES in this module makefile!)
   endif
   ifdef APP_NAME
      $(error Illegal definition APP_NAME in this module makefile!)
   endif
   ifndef SOURCE_FILES
      $(error List of SOURCE_FILES not defined in this module makefile!)
   endif
endif
ifdef APP_NAME
   ifeq ($(origin MODULES),file)
      $(error Illegal definition MODULES in this application makefile!)
   endif
   ifdef LIB_NAME
      $(error Illegal definition LIB_NAME in this application makefile!)
   endif
   ifndef SOURCE_FILES
      $(error List of SOURCE_FILES not defined in this application makefile!)
   endif
endif

ifeq ($(strip $(MODULES)$(LIB_NAME)$(APP_NAME)),)
   $(error One of these targets MUST be defined: MODULES, LIB_NAME, APP_NAME)
endif


### directory and file definitions ############################################


# input/output file extensions
HEADER_EXT := h
INC_EXT    := inc
CXX_EXT    := cc
C_EXT      := c
F77_EXT    := F
ORACLE_EXT := pc
CDF_EXT    := cdf
LIST_EXT   := lis
DEP_EXT    := d
OBJ_EXT    := o
SO_EXT     := so
PCM_EXT    := pcm
MAP_EXT    := rootmap

# build/install directory definitions
ifdef MODULES
   ifeq ($(MAKELEVEL),0)
      export MODULES
      ifndef BUILD_DIR
         BUILD_DIR := ./build
      else
         export BUILD_DIR
      endif
      ifndef INSTALL_DIR
         ifndef MYHADDIR
            INSTALL_DIR := ./install
         else
            INSTALL_DIR := $(MYHADDIR)
         endif
      else
         export INSTALL_DIR
      endif
   endif
endif
ifdef LIB_NAME
   BUILD_DIR ?= ../build
   ifndef INSTALL_DIR
      ifdef MYHADDIR
         INSTALL_DIR := $(MYHADDIR)
      else
         INSTALL_DIR := ../install
      endif
   endif
   export BUILD_DIR INSTALL_DIR
endif
ifdef APP_NAME
   BUILD_DIR   ?= ./build
   INSTALL_DIR ?= .
   export BUILD_DIR INSTALL_DIR
endif

# output directories
OBJ_DIR := $(BUILD_DIR)/obj
LIB_DIR := $(BUILD_DIR)/lib
DEP_DIR := $(BUILD_DIR)/depend
DOC_DIR	:= $(BUILD_DIR)/doc
PC_DIR  := $(BUILD_DIR)/pc

# search path for prerequisites
SOURCE_DIRS := $(subst /,,$(sort $(dir . $(SOURCE_FILES))))
vpath %.$(HEADER_EXT) $(SOURCE_DIRS)
vpath %.$(INC_EXT)    $(SOURCE_DIRS)
vpath %.$(CDF_EXT)    $(SOURCE_DIRS)
vpath %.$(CXX_EXT)    $(SOURCE_DIRS)
vpath %.$(C_EXT)      $(SOURCE_DIRS)
vpath %.$(F77_EXT)    $(SOURCE_DIRS)
vpath %.$(ORACLE_EXT) $(SOURCE_DIRS)

# lists of header files (excluding LinkDef files) within SOURCE_DIRS - with
# and without path
HEADER_FILES := $(foreach dir,$(SOURCE_DIRS), \
   $(wildcard $(dir)/*.$(HEADER_EXT)))    
HEADER_FILES := $(filter-out %LinkDef.$(HEADER_EXT),$(HEADER_FILES))
HEADER_FILE_NAMES = $(notdir $(HEADER_FILES))

# list of CINT dictionary definition files within SOURCE_DIRS - with
# and without path
DICT_DEF_FILES := $(foreach dir,$(SOURCE_DIRS), \
   $(wildcard $(dir)/*LinkDef.$(HEADER_EXT)))
DICT_DEF_FILE_NAMES = $(notdir $(DICT_DEF_FILES))

# lists of KUIP files - with and without path
CDF_FILES      := $(foreach dir,$(SOURCE_DIRS),$(wildcard $(dir)/*.$(CDF_EXT)))
CDF_FILE_NAMES := $(notdir $(CDF_FILES))

# list of input to CINT dictionary - with and without path; a list of
# header files deduced from SOURCE_FILES + headers named *def.h
DICT_INPUT_FILES     := $(SOURCE_FILES:.$(CXX_EXT)=.$(HEADER_EXT))
DICT_INPUT_FILES     := $(DICT_INPUT_FILES:.$(ORACLE_EXT)=.$(HEADER_EXT))
DICT_INPUT_FILES     += $(foreach dir,$(SOURCE_DIRS), \
                           $(wildcard $(dir)/*def.$(HEADER_EXT)))
DICT_INPUT_FILES     := $(call cleanlist,$(DICT_INPUT_FILES))
DICT_INPUT_FILES     := $(call prefix-with-path,$(DICT_INPUT_FILES))
DICT_INPUT_FILE_NAMES = $(notdir $(DICT_INPUT_FILES))

# list of source files - with and without path, deduced from SOURCE_FILES
C_FILES           := $(filter %.$(C_EXT),$(SOURCE_FILES))
C_FILES           := $(call prefix-with-path,$(C_FILES))
C_FILE_NAMES       = $(notdir $(C_FILES))
CXX_FILES         := $(filter %.$(CXX_EXT),$(SOURCE_FILES))
CXX_FILES         := $(call prefix-with-path,$(CXX_FILES))
CXX_FILE_NAMES     = $(notdir $(CXX_FILES))
F77_FILES         := $(filter %.$(F77_EXT),$(SOURCE_FILES))
F77_FILES         := $(call prefix-with-path,$(F77_FILES))
F77_FILE_NAMES     = $(notdir $(F77_FILES))
ORACLE_FILES      := $(filter %.$(ORACLE_EXT),$(SOURCE_FILES))
ORACLE_FILES      := $(call prefix-with-path,$(ORACLE_FILES))
ORACLE_FILE_NAMES  = $(notdir $(ORACLE_FILES))

# intermediate files with paths: dependency files
DEPENDENCIES = \
   $(C_FILE_NAMES:%.$(C_EXT)=$(DEP_DIR)/%.$(C_EXT).$(DEP_EXT))
DEPENDENCIES = \
   $(CXX_FILE_NAMES:%.$(CXX_EXT)=$(DEP_DIR)/%.$(CXX_EXT).$(DEP_EXT))
DEPENDENCIES += \
   $(F77_FILE_NAMES:%.$(F77_EXT)=$(DEP_DIR)/%.$(F77_EXT).$(DEP_EXT))
DEPENDENCIES += \
   $(ORACLE_FILE_NAMES:%.$(ORACLE_EXT)=$(DEP_DIR)/%.$(ORACLE_EXT).$(DEP_EXT))
DEPENDENCIES += \
   $(DICT_DEF_FILE_NAMES:%LinkDef.$(HEADER_EXT)=$(DEP_DIR)/%Dict.$(CXX_EXT).$(DEP_EXT))

# intermediate files with paths: output files of precompilers
PRECOMPILED  = \
   $(CDF_FILE_NAMES:%.$(CDF_EXT)=$(PC_DIR)/%.$(CDF_EXT).$(F77_EXT))
PRECOMPILED += \
   $(ORACLE_FILE_NAMES:%.$(ORACLE_EXT)=$(PC_DIR)/%.$(LIST_EXT))
PRECOMPILED += \
   $(ORACLE_FILE_NAMES:%.$(ORACLE_EXT)=$(PC_DIR)/%.$(ORACLE_EXT).$(CXX_EXT))
PRECOMPILED += \
   $(ORACLE_FILE_NAMES:%.$(ORACLE_EXT)=$(PC_DIR)/%.$(LIST_EXT))
PRECOMPILED += \
   $(DICT_DEF_FILE_NAMES:%LinkDef.$(HEADER_EXT)=$(PC_DIR)/%Dict.$(CXX_EXT))
PRECOMPILED += \
   $(DICT_DEF_FILE_NAMES:%LinkDef.$(HEADER_EXT)=$(PC_DIR)/%Dict.$(HEADER_EXT))

# intermediate files with paths: object output files of compilers
OBJECTS	 = \
   $(C_FILE_NAMES:%.$(C_EXT)=$(OBJ_DIR)/%.$(C_EXT).$(OBJ_EXT))
OBJECTS	+= \
   $(CXX_FILE_NAMES:%.$(CXX_EXT)=$(OBJ_DIR)/%.$(CXX_EXT).$(OBJ_EXT))
OBJECTS += \
   $(F77_FILE_NAMES:%.$(F77_EXT)=$(OBJ_DIR)/%.$(F77_EXT).$(OBJ_EXT))
OBJECTS += \
   $(ORACLE_FILE_NAMES:%.$(ORACLE_EXT)=$(OBJ_DIR)/%.$(ORACLE_EXT).$(OBJ_EXT))
OBJECTS += \
   $(DICT_DEF_FILE_NAMES:%LinkDef.$(HEADER_EXT)=$(OBJ_DIR)/%Dict.$(OBJ_EXT))
OBJECTS += \
   $(CDF_FILE_NAMES:%.$(CDF_EXT)=$(OBJ_DIR)/%.$(CDF_EXT).$(OBJ_EXT))

# shared library path and name
SO_FILE_NAME := lib$(LIB_NAME).$(SO_EXT)
SO_FILE      := $(LIB_DIR)/$(SO_FILE_NAME)
PCM_FILE_NAME := $(LIB_NAME)Dict_rdict.$(PCM_EXT)
MAP_FILE_NAME := lib$(LIB_NAME).$(MAP_EXT)

# miscellaneous files
CREATE_DOC_MACRO := macros/makeDocs.C


### (pre-)compilers and linker settings #######################################


# default list and order of library directories (at build- and runtime)
LIB_DIRS        := $(ROOT_LIB_DIRS)
SO_LOADER_DIRS  := $(ROOT_LIB_DIRS)
APP_LOADER_DIRS := $(ROOT_LIB_DIRS)
ifdef HADDIR
   LIB_DIRS        := $(HADDIR)/lib $(LIB_DIRS)
   APP_LOADER_DIRS := $(HADDIR)/lib $(APP_LOADER_DIRS)
endif
ifdef MYHADDIR
   LIB_DIRS        := $(MYHADDIR)/lib $(LIB_DIRS)
   APP_LOADER_DIRS := $(MYHADDIR)/lib $(APP_LOADER_DIRS)
endif

# default list of linked libraries
LIBS := $(ROOT_LIBS)

# default list and order of include directories - currently, sub-modules 
# must be hard coded, here :-(
INC_DIRS := $(SOURCE_DIRS)
ifndef APP_NAME
   INC_DIRS += ../base/datasource ../base/datastruct ../base/eventhandling
   INC_DIRS += ../base/geantutil ../base/geometry ../base/runtimedb 
   INC_DIRS += ../base/util
   INC_DIRS += $(foreach module,$(MODULES),$(addprefix ../,$(module)))
endif
ifdef MYHADDIR
   INC_DIRS += $(MYHADDIR)/include
endif
ifdef HADDIR
   INC_DIRS += $(HADDIR)/include
endif
ifdef PLUTODIR
   INC_DIRS += $(PLUTODIR)/include
endif
INC_DIRS += $(ROOT_INC_DIRS)


# add-on's in case remote file IO service at GSI is used
ifeq ($(strip $(USES_RFIO)),no)
   export USES_RFIO
endif
ifeq ($(strip $(USES_RFIO)),yes)
   ifndef ADSM_BASE_NEW
      $(error ADSM_BASE_NEW not set!)
   endif
   INC_DIRS        += $(RFIO_INC_DIRS)
   LIB_DIRS        += $(RFIO_LIB_DIRS)
   SO_LOADER_DIRS  += $(RFIO_LIB_DIRS)
   APP_LOADER_DIRS += $(RFIO_LIB_DIRS)
   LIBS            += $(RFIO_LIBS)
   CPP_FLAGS       += $(RFIO_CPP_FLAGS)
endif

# add-on's in case DABC is used
ifeq ($(strip $(USES_DABC)),no)
   export USES_DABC
endif
ifeq ($(strip $(USES_DABC)),yes)
   INC_DIRS        += $(DABC_INC_DIRS)
   LIB_DIRS        += $(DABC_LIB_DIRS)
   SO_LOADER_DIRS  += $(DABC_LIB_DIRS)
   APP_LOADER_DIRS += $(DABC_LIB_DIRS)
   LIBS            += $(DABC_LIBS)
   CPP_FLAGS       += $(DABC_CPP_FLAGS)
endif



ifeq ($(strip $(USES_ORACLE)),yes)
  CPP_FLAGS += $(ORACLE_CPP_FLAGS)
endif
# add-on's in case Oracle is used 

ifneq ($(strip $(ORACLE_FILES)),)
   USES_ORACLE := yes
endif

ifeq ($(strip $(USES_ORACLE)),no)
   export USES_ORACLE
endif
ifeq ($(strip $(USES_ORACLE)),yes)
   ifndef ORACLE_HOME
      $(error ORACLE_HOME not set!)
   endif
   ifndef ORA_USER
      $(error ORA_USER not set!)
   endif
   INC_DIRS        += $(ORACLE_INC_DIRS)
   LIB_DIRS        += $(ORACLE_LIB_DIRS)
   SO_LOADER_DIRS  += $(ORACLE_LIB_DIRS)
   APP_LOADER_DIRS += $(ORACLE_LIB_DIRS)
   LIBS            += $(ORACLE_LIBS)
endif

# add-on's in case the CERN library is used
ifeq ($(strip $(USES_CERNLIB)),no)
   export USES_CERNLIB
endif
ifeq ($(strip $(USES_CERNLIB)),yes)
   ifndef CERN_ROOT
      $(error CERN_ROOT not set!)
   endif
   INC_DIRS        := $(CERNLIB_INC_DIRS) $(INC_DIRS)
   FPP_FLAGS       += $(CERNLIB_FPP_FLAGS)
   CPP_FLAGS       += $(CERNLIB_CPP_FLAGS)
   LIB_DIRS        += $(CERNLIB_LIB_DIRS)
   SO_LOADER_DIRS  += $(CERNLIB_LIB_DIRS)
   APP_LOADER_DIRS += $(CERNLIB_LIB_DIRS)
   LIBS            += $(CERNLIB_LIBS)
endif

# add-on's in case the X11 system is used
ifeq ($(strip $(USES_X11)),no)
   export USES_X11
endif
ifeq ($(strip $(USES_X11)),yes)
   INC_DIRS        += $(X11_INC_DIRS)
   LIB_DIRS        += $(X11_LIB_DIRS)
   SO_LOADER_DIRS  += $(X11_LIB_DIRS)
   APP_LOADER_DIRS += $(X11_LIB_DIRS)
   LIBS            += $(X11_LIBS)
endif

# C, C++ and Fortran preprocessor flags
CPP_FLAGS += -DR__GLIBC -DDEBUG_LEVEL=$(DEBUG_LEVEL)

ROOTMAJVERSION=$(shell $(ROOTCONFIG) --version | cut -d "." -f 1)

ifneq ($(ROOTMAJVERSION),5)
CPP_FLAGS += -std=c++14
endif

CPP_FLAGS += $(addprefix -I,$(call cleanlist,$(INC_DIRS)))
FPP_FLAGS += $(addprefix -I,$(call cleanlist,$(INC_DIRS)))

# C compiler flags
CC         := gcc
SO_C_FLAGS  = $(OPTLEVEL) -Wall -Wextra -Wno-unused-parameter -pipe -pthread -fPIC $(DEBUGOPT) $(CPP_FLAGS)
APP_C_FLAGS = $(OPTLEVEL) -Wall -Wextra -Wno-unused-parameter -pipe -pthread $(DEBUGOPT) $(CPP_FLAGS)

# C++ compiler flags
CXX          := g++

SO_CXX_FLAGS  = $(OPTLEVEL) -Wall -Wextra -Wno-unused-parameter -pipe -pthread -fPIC $(DEBUGOPT) $(CPP_FLAGS)
APP_CXX_FLAGS = $(OPTLEVEL) -Wall -Wextra -Wno-unused-parameter -pipe -pthread $(DEBUGOPT) $(CPP_FLAGS)

# flags of Root's dictionary file precompiler
ROOTCINT_FLAGS = -c -p $(CPP_FLAGS)
#ROOTCINT_FLAGS = $(filter-out -std=c++14,$(CPP_FLAGS))

ifneq ($(ROOTMAJVERSION),5)
   ROOTCINT_FLAGS = $(filter-out -std=c++14,$(CPP_FLAGS))
endif

# Fortran compiler flags
F77          := g77

ifeq ($(strip $(USES_GFORTRAN)),yes) 
   F77       = gfortran
endif

SO_F77_FLAGS  = $(DEBUGOPT) -std=legacy -fno-optimize-sibling-calls $(GCC8LINKADD) $(FPP_FLAGS)
APP_F77_FLAGS = $(DEBUGOPT) -std=legacy -fno-optimize-sibling-calls $(GCC8LINKADD) $(FPP_FLAGS)

# flags of the Oracle pre-compiler
PC      := $(ORACLE_HOME)/bin/proc
PC_FLAGS = $(addprefix include=,$(call cleanlist,$(INC_DIRS))) \
           oraca=yes code=cpp parse=partial sqlcheck=semantics \
           ireclen=130 oreclen=130 userid=$(ORA_USER)

# linker and loader flags to build Hydra shared libraries and applications
LD           := g++
SO_LD_FLAGS   = $(DEBUGOPT) -shared -Wl,-soname,lib$(LIB_NAME).$(SO_EXT)
SO_LD_FLAGS  += $(addprefix -Wl$(comma)-rpath$(comma),$(SO_LOADER_DIRS))
SO_LD_FLAGS  += $(addprefix -L,$(call cleanlist,$(LIB_DIRS)))
APP_LD_FLAGS  = $(DEBUGOPT) $(GCC8LINKADD)
APP_LD_FLAGS += $(addprefix -Wl$(comma)-rpath$(comma),$(APP_LOADER_DIRS))
APP_LD_FLAGS += $(addprefix -L,$(call cleanlist,$(LIB_DIRS)))
