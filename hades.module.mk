###############################################################################
#
# Makefile hades.module.mk - to be included from all module Makefiles.
#
# This makefile is part of the Hades build system. Basically this makefile
# contains a related target to all targets in hades.rules.mk, but on module
# level.
#
###############################################################################


.PHONY: all check depend build install deinstall clean distclean


### alias targets


all: build


### normal targets called by the user


check:
	@$(ECHO) -e \
	   "\n\033[1mChecking settings for library $(SO_FILE_NAME) ...\033[0m"

        # check relevant path settings
	@for dir in $(INC_DIRS); do		 			\
	if [ ! -d $$dir ]; then 					\
		$(ECHO) -n "WARNING: ";					\
		$(ECHO) "Include directory $$dir does not exist!";	\
	fi;								\
	done;								\
	$(ECHO) "Order of include paths:";				\
	for dir in $(INC_DIRS); do					\
	if [ -d $$dir ]; then 						\
		$(ECHO) "   $$dir";					\
	fi;								\
	done;								\
	for dir in $(LIB_DIRS); do					\
	if [ ! -d $$dir ]; then 					\
		$(ECHO) -n "WARNING: ";					\
		$(ECHO) "Library directory $$dir does not exist!";	\
	fi;								\
	done;								\
	$(ECHO) "Order of library paths:";				\
	for dir in $(LIB_DIRS); do					\
	if [ -d $$dir ]; then 						\
		$(ECHO) "   $$dir";					\
	fi;								\
	done;

        # is there a header file for each source file
	@for file in $(CXX_FILES) $(ORACLE_FILES); do			\
	if [ ! -f $${file//.$(CXX_EXT)/.$(HEADER_EXT)} ]; then		\
		$(ECHO) -n "HINT: ";					\
		$(ECHO) "$$file has no corresponding header file.";	\
	fi;								\
	done

        # is there a header file for each Oracle source file
	@for file in $(ORACLE_FILES); do				\
	if [ ! -f $${file//.$(ORACLE_EXT)/.$(HEADER_EXT)} ]; then	\
		$(ECHO) -n "HINT: ";					\
		$(ECHO) "$$file has no corresponding header file.";	\
	fi;								\
	done

        # is there a source file for each header file
	@for file in $(HEADER_FILES); do				\
	case $$file in							\
	*LinkDef.$(HEADER_EXT)) continue ;;				\
	esac;								\
	cxx_file=$${file//.$(HEADER_EXT)/.$(CXX_EXT)};			\
	oracle_file=$${file//.$(HEADER_EXT)/.$(ORACLE_EXT)};		\
	if [ ! -f $$cxx_file ]  &&  [ ! -f $$oracle_file ]; then	\
		$(ECHO) -n "HINT: ";					\
		$(ECHO) "$$file has no corresponding source file.";	\
	fi;								\
	done

        # check for orphan C source files
	@for file in $(wildcard *.$(C_EXT) */*.$(C_EXT)); do		\
		case "$(C_FILES)" in					\
		*$$file*) continue ;;					\
		*) $(ECHO) -n "WARNING: ";				\
		   $(ECHO) "$$file not included in build process." ;;	\
		esac;							\
	done

        # check for orphan C++ source files
	@for file in $(wildcard *.$(CXX_EXT) */*.$(CXX_EXT)); do	\
		case "$(CXX_FILES)" in					\
		*$$file*) continue ;;					\
		*) $(ECHO) -n "WARNING: ";				\
		   $(ECHO) "$$file not included in build process." ;;	\
		esac;							\
	done

        # check for orphan Oracle source files
	@for file in $(wildard *.$(ORACLE_EXT) */*.$(ORACLE_EXT)); do	\
		case "$(ORACLE_FILES)" in				\
		*$$file*) continue ;;					\
		*) $(ECHO) -n "WARNING: ";				\
		   $(ECHO) "$$file not included in build process." ;;	\
		esac;							\
	done

        # check for not existing source files
	@for file in $(SOURCE_FILES); do				\
	if [ ! -f $$file ]; then					\
		$(ECHO) -n "ERROR: ";					\
		$(ECHO) "Source file $$file does not exist!";		\
		exit 1;							\
	fi;								\
	done

# just create the output directories
pre-build:
	@set -e;					\
	[ -e $(BUILD_DIR) ]  ||  mkdir $(BUILD_DIR);	\
	[ -e $(DEP_DIR) ]    ||  mkdir $(DEP_DIR);	\
	[ -e $(PC_DIR) ]     ||  mkdir $(PC_DIR);	\
	[ -e $(OBJ_DIR) ]    ||  mkdir $(OBJ_DIR);	\
	[ -e $(LIB_DIR) ]    ||  mkdir $(LIB_DIR)


# just create all dependency files
depend: pre-build $(DEPENDENCIES)


# build a shared libary
build: pre-build $(SO_FILE)


# install all neccessary header and library files
install:
	@[ -e $(INSTALL_DIR) ]         ||  mkdir $(INSTALL_DIR);	\
	[ -e $(INSTALL_DIR)/lib ]      ||  mkdir $(INSTALL_DIR)/lib;	\
	[ -e $(INSTALL_DIR)/include ]  ||  mkdir $(INSTALL_DIR)/include;

	@set -e;							\
	if [ ! -f $(SO_FILE) ]; then					\
		$(ECHO) -ne "ERROR: ";					\
		$(ECHO) "Library $(LIB_NAME) was not build!";		\
		exit 1;							\
	fi;								\
	$(ECHO) "$(INSTALL) -m 755 $(SO_FILE) $(INSTALL_DIR)/lib";	\
	$(INSTALL) -m 755 $(SO_FILE) $(INSTALL_DIR)/lib

	@set -e;										\
	if [ -f $(PC_DIR)/$(PCM_FILE_NAME) ]; then						\
		$(ECHO) "$(INSTALL) -m 755 $(PC_DIR)/$(PCM_FILE_NAME) $(INSTALL_DIR)/lib";	\
		$(INSTALL) -m 755 $(PC_DIR)/$(PCM_FILE_NAME) $(INSTALL_DIR)/lib;		\
	fi

	@set -e;										\
	if [ -f $(PC_DIR)/$(MAP_FILE_NAME) ]; then						\
		$(ECHO) "$(INSTALL) -m 755 $(PC_DIR)/$(MAP_FILE_NAME) $(INSTALL_DIR)/lib";	\
		$(INSTALL) -m 755 $(PC_DIR)/$(MAP_FILE_NAME) $(INSTALL_DIR)/lib;		\
	fi

	@set -e;								\
	if [ $(findstring libHydra.so,$(SO_FILE)) ]; then     			\
	for d in $(BASE_PCM) ; do						\
	dir="$$d"Dict_rdict.pcm ;						\
	if [ -f $(PC_DIR)/$$dir ]; then						\
		$(ECHO) "$(INSTALL) -m 755 $(PC_DIR)/$$dir $(INSTALL_DIR)/lib";	\
		$(INSTALL) -m 755 $(PC_DIR)/$$dir $(INSTALL_DIR)/lib;		\
	fi;							 		\
	dir=lib"$$d".rootmap ;						\
	if [ -f $(PC_DIR)/$$dir ]; then						\
		$(ECHO) "$(INSTALL) -m 755 $(PC_DIR)/$$dir $(INSTALL_DIR)/lib";	\
		$(INSTALL) -m 755 $(PC_DIR)/$$dir $(INSTALL_DIR)/lib;		\
	fi;							 		\
	done;									\
	fi

	@set -e;							\
	if [ "x$(HEADER_FILES)" = "x" ]; then				\
		$(ECHO) -ne "ERROR: ";					\
		$(ECHO) "No header files found!";			\
		exit 1;							\
	fi;								\
	for file in $(HEADER_FILES); do					\
	case $$file in							\
		*LinkDef.$(HEADER_EXT)) continue ;;			\
		*Cint.$(HEADER_EXT)) continue ;;			\
	esac;								\
	$(ECHO) "$(INSTALL) -m 644 $$file $(INSTALL_DIR)/include";	\
	$(INSTALL) -m 644 $$file $(INSTALL_DIR)/include;		\
	done


# deinstall all header and library files of this module
deinstall:
	@set -e;							\
	if [ "x$(HEADER_FILE_NAMES)" = "x" ]; then			\
		$(ECHO) "ERROR: No header files found!";		\
		exit 1;							\
	fi;								\
	for file in $(HEADER_FILE_NAMES); do				\
		$(ECHO) "$(RM) -f $(INSTALL_DIR)/include/$$file";	\
		$(RM) -f $(INSTALL_DIR)/include/$$file;			\
	done
	$(RM) -f $(INSTALL_DIR)/lib/$(SO_FILE_NAME);			\
	$(RM) -f $(INSTALL_DIR)/lib/$(PCM_FILE_NAME)                    \
	$(RM) -f $(INSTALL_DIR)/lib/$(MAP_FILE_NAME)                    

	@set -e;								\
	if [ $(findstring libHydra.so,$(SO_FILE)) ]; then     			\
	for d in $(BASE_PCM) ; do						\
	dir="$$d"Dict_rdict.pcm ;						\
	if [ -f $(INSTALL_DIR)/lib/$$dir ]; then				\
		$(RM) -f $(INSTALL_DIR)/lib/$$dir;				\
	fi;									\
	dir=lib"$$d".rootmap ;						\
	if [ -f $(INSTALL_DIR)/lib/$$dir ]; then				\
		$(RM) -f $(INSTALL_DIR)/lib/$$dir;				\
	fi;									\
	done;									\
	fi

# delete all but the source and dependency files
clean:
	$(RM) -f $(SO_FILE) $(PC_DIR)/$(PCM_FILE_NAME) $(PC_DIR)/$(MAP_FILE_NAME) $(OBJECTS) $(PRECOMPILED) 


# delete all but the source files
distclean: clean
	$(RM) -f $(DEPENDENCIES) *~


### Pattern and implicite rules ###############################################


# delete all suffix rules - only pattern rules are used here
.SUFFIXES:

# do not delete this files automatically
.SECONDARY: $(DEPENDENCIES) $(PRECOMPILED) $(OBJECTS)

# include the dependencies only if there are files at all
ifneq ($(call cleanlist,$(DEPENDENCIES)),)
   include $(call cleanlist,$(DEPENDENCIES))
endif


# precompile Oracle input files (.pc), resulting in .cc files
$(PC_DIR)/%.$(ORACLE_EXT).$(CXX_EXT): %.$(ORACLE_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Precompiling $<"
	$(PC) $(PC_FLAGS) lname=$(subst .$(ORACLE_EXT).$(CXX_EXT),.$(LIST_EXT),$@) iname=$< oname=$@


# precompile Geant menu definition file
$(PC_DIR)/%.$(CDF_EXT).$(F77_EXT): %.$(CDF_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Precompiling $<"
	$(KUIPC) $< $@

	
# create (precompile) Cint dictionary source file based on LinkDef.h file
$(PC_DIR)/%Dict.$(CXX_EXT): %LinkDef.$(HEADER_EXT) $(DICT_INPUT_FILE_NAMES)
	@$(ECHO) "==> $(LIB_NAME) : Creating Dictionary $@"
        ifneq ($(ROOTMAJVERSION),5)
        ifeq (,$(findstring $(SO_FILE_NAME),$(HYDRA_SO)))
	$(ROOTCINT) -f $@ -rml $(SO_FILE_NAME) $(HYDRA_SO) -rmf $(PC_DIR)/lib$(subst Dict.cc,,$(notdir $@)).$(MAP_EXT) -noIncludePaths $(ROOTCINT_FLAGS) $(notdir $(filter $(dir $<)%.$(HEADER_EXT),$(DICT_INPUT_FILES))) $< 
        else
	$(ROOTCINT) -f $@ $(HYDRA_SO) -rmf $(PC_DIR)/lib$(subst Dict.cc,,$(notdir $@)).$(MAP_EXT) -noIncludePaths $(ROOTCINT_FLAGS) $(notdir $(filter $(dir $<)%.$(HEADER_EXT),$(DICT_INPUT_FILES))) $< 
        endif
        else
	$(ROOTCINT) -f $@ $(ROOTCINT_FLAGS) $(filter $(dir $<)%.$(HEADER_EXT),$(DICT_INPUT_FILES)) $<
        endif

ifeq ($(strip $(shell $(CXX) -dumpversion)),2.95.4) 

# GCC 2.95.4: create dependencies of precompiled Oracle files
$(DEP_DIR)/%.$(ORACLE_EXT).$(DEP_EXT): $(PC_DIR)/%.$(ORACLE_EXT).$(CXX_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -M $(CPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(ORACLE_EXT).$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

# create dependencies of Cint dictionary files
$(DEP_DIR)/%Dict.$(CXX_EXT).$(DEP_EXT): $(PC_DIR)/%Dict.$(CXX_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -M $(CPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(CXX_EXT).$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

# GCC 2.95.4: create dependencies of C++ source files
$(DEP_DIR)/%.$(CXX_EXT).$(DEP_EXT): %.$(CXX_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -M $(CPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(CXX_EXT).$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

# GCC 2.95.4: create dependencies of C source files
$(DEP_DIR)/%.$(C_EXT).$(DEP_EXT): %.$(C_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(CC) -M $(CPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(C_EXT).$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

# GCC 2.95.4: create dependencies of F77 source files
$(DEP_DIR)/%.$(F77_EXT).$(DEP_EXT): %.$(F77_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(F77) -M $(FPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(F77_EXT).$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

else

# create dependencies of precompiled Oracle files
$(DEP_DIR)/%.$(ORACLE_EXT).$(DEP_EXT): $(PC_DIR)/%.$(ORACLE_EXT).$(CXX_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -MM $(CPP_FLAGS) -MT $(OBJ_DIR)/$*.$(ORACLE_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

# create dependencies of Cint dictionary files
$(DEP_DIR)/%Dict.$(CXX_EXT).$(DEP_EXT): $(PC_DIR)/%Dict.$(CXX_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -MM $(CPP_FLAGS) -MT $(OBJ_DIR)/$*.$(CXX_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

# create dependencies of C++ source files
$(DEP_DIR)/%.$(CXX_EXT).$(DEP_EXT): %.$(CXX_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -MM $(CPP_FLAGS) -MT $(OBJ_DIR)/$*.$(CXX_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

# create dependencies of C source files
$(DEP_DIR)/%.$(C_EXT).$(DEP_EXT): %.$(C_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(CC) -MM $(CPP_FLAGS) -MT $(OBJ_DIR)/$*.$(C_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

# create dependencies of F77 source files
$(DEP_DIR)/%.$(F77_EXT).$(DEP_EXT): %.$(F77_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Determining Dependencies of $<"
	set -e; $(F77) -MM $(FPP_FLAGS) -MT $(OBJ_DIR)/$*.$(F77_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

endif


# compile C++ source files created by the Oracle precompiler
$(OBJ_DIR)/%.$(ORACLE_EXT).$(OBJ_EXT): \
		$(PC_DIR)/%.$(ORACLE_EXT).$(CXX_EXT) \
		$(DEP_DIR)/%.$(ORACLE_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Compiling $<"
	$(CXX) $(SO_CXX_FLAGS) -o $@ -c $<


# compile Fortran 77 files created by KUIPC
$(OBJ_DIR)/%.$(CDF_EXT).$(OBJ_EXT): \
		$(PC_DIR)/%.$(CDF_EXT).$(F77_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Compiling $<" 
	$(F77) $(SO_F77_FLAGS) -o $@ -c $<


# compile C++ source files
$(OBJ_DIR)/%.$(CXX_EXT).$(OBJ_EXT): \
		%.$(CXX_EXT) $(DEP_DIR)/%.$(CXX_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Compiling $<" 
	$(CXX) $(SO_CXX_FLAGS) -o $@ -c $<


# compile C source files
$(OBJ_DIR)/%.$(C_EXT).$(OBJ_EXT): \
		%.$(C_EXT) $(DEP_DIR)/%.$(C_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Compiling $<" 
	$(CC) $(SO_C_FLAGS) -o $@ -c $<


# compile Fortran 77 files
$(OBJ_DIR)/%.$(F77_EXT).$(OBJ_EXT): \
		%.$(F77_EXT) $(DEP_DIR)/%.$(F77_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Compiling $<" 
	$(F77) $(SO_F77_FLAGS) -o $@ -c $<


# compile Cint dictionary
$(OBJ_DIR)/%Dict.$(OBJ_EXT): $(PC_DIR)/%Dict.$(CXX_EXT) $(DEP_DIR)/%Dict.$(CXX_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(LIB_NAME) : Compiling $<"
	$(CXX) $(SO_CXX_FLAGS) -I$(PC_DIR) -o $@ -c $<


# link shared library
$(LIB_DIR)/%.$(SO_EXT): $(OBJECTS)
	@$(ECHO) "==> $(LIB_NAME) : Linking $@"
	$(LD) $(SO_LD_FLAGS) -o $@ $+ $(LIBS)
	@$(ECHO) -e "\n==> Shared library \033[1m$(LIB_NAME)\033[0m created.\n"
