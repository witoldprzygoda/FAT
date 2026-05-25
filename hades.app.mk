###############################################################################
#
# Makefile hades.app.mk - to be included from the global Makefiles.
#
# This makefile is part of the Hydra build system. It contains all default
# rules to build all Hydra modules. To use this makefile, include it in your
# global Makefile which lists the modules which should be build together.
#
# See the file Makefile for details. 
#
###############################################################################


.PHONY: all build check pre-build depend clean distclean install deinstall


### alias targets #############################################################


all:   $(APP_NAME)
build: $(APP_NAME)


### normal targets called by the user #########################################


check:
        # check relevant path and library settings
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
	done;								\
	$(ECHO) "Linked libraries:";					\
	for lib in $(LIBS); do						\
		$(ECHO) "   $${lib#-l}";				\
	done

        # check for not existing source files
	@for file in $(SOURCE_FILES); do				\
	if [ ! -f $$file ]; then					\
		$(ECHO) -ne "ERROR: ";					\
		$(ECHO) "Source file $$file for does not exist!";	\
		exit 1;							\
	fi;								\
	done


pre-build:
	@set -e;					\
	[ -e $(BUILD_DIR) ]  ||  mkdir $(BUILD_DIR);	\
	[ -e $(DEP_DIR) ]    ||  mkdir $(DEP_DIR);	\
	[ -e $(PC_DIR) ]     ||  mkdir $(PC_DIR);	\
        [ -e $(OBJ_DIR) ]    ||  mkdir $(OBJ_DIR)


depend: pre-build $(DEPENDENCIES)
	@$(ECHO) -e "\n==> Dependency files created.\n"


$(APP_NAME): pre-build $(BUILD_DIR)/$(APP_NAME)
	@$(ECHO) -e "\n==> Application \033[1m$(APP_NAME)\033[0m created.\n"


install:
	$(INSTALL) -D -m 755 $(BUILD_DIR)/$(APP_NAME) $(INSTALL_DIR)/$(APP_NAME)
	@set -e;												\
	if [ -f $(PC_DIR)/KineDict_rdict.pcm ]; then								\
	        $(ECHO) "$(INSTALL) -m 755 $(PC_DIR)/KineDict_rdict.pcm $(INSTALL_DIR)/KineDict_rdict.pcm";	\
	        $(INSTALL) -m 755 $(PC_DIR)/KineDict_rdict.pcm $(INSTALL_DIR)/KineDict_rdict.pcm;		\
	fi
	@$(ECHO) -e "\n==> Application \033[1m$(APP_NAME)\033[0m installed.\n"


uninstall:
	$(RM) $(INSTALL_DIR)/$(APP_NAME)
	@set -e;												\
	if [ -f $(INSTALL_DIR)/KineDict_rdict.pcm ]; then								\
	        $(ECHO) "rm $(INSTALL_DIR)/KineDict_rdict.pcm";	\
	        $(RM) $(INSTALL_DIR)/KineDict_rdict.pcm;		\
	fi
	@$(ECHO) -e "\n==> Application \033[1m$(APP_NAME)\033[0m removed.\n"


# delete all but the source and dependency files
clean:
	$(RM) -f $(BUILD_DIR)/$(APP_NAME) $(OBJECTS) $(PRECOMPILED)
	@$(ECHO) -e "\n==> All files but dependency files deleted.\n"


# delete all but the source files
distclean: clean
	$(RM) -f $(DEPENDENCIES) *~
	- $(RMDIR) $(OBJ_DIR) $(DEP_DIR) $(PC_DIR) $(BUILD_DIR)
	@$(ECHO) -e "\n==> All files but sources deleted.\n"


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
	@$(ECHO) "==> $(APP_NAME) : Precompiling $<"
	$(PC) $(PC_FLAGS) lname=$(subst .$(ORACLE_EXT).$(CXX_EXT),.$(LIST_EXT),$@) iname=$< oname=$@


# precompile Geant menu definition file
$(PC_DIR)/%.$(CDF_EXT).$(F77_EXT): %.$(CDF_EXT)
	@$(ECHO) "==> $(APP_NAME) : Precompiling $<"
	$(KUIPC) $< $@


ifeq ($(strip $(shell $(CXX) -dumpversion)),2.95.4) 

# GCC 2.95.4: create dependencies of precompiled Oracle files
$(DEP_DIR)/%.$(ORACLE_EXT).$(DEP_EXT): $(PC_DIR)/%.$(ORACLE_EXT).$(CXX_EXT)
	@$(ECHO) "==> $(APP_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -M $(CPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

# GCC 2.95.4: create dependencies of C++ source files
$(DEP_DIR)/%.$(CXX_EXT).$(DEP_EXT): %.$(CXX_EXT)
	@$(ECHO) "==> $(APP_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -M $(CPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(CXX_EXT).$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

# GCC 2.95.4: create dependencies of C source files
$(DEP_DIR)/%.$(C_EXT).$(DEP_EXT): %.$(C_EXT)
	@$(ECHO) "==> $(APP_NAME) : Determining Dependencies of $<"
	set -e; $(CC) -M $(CPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(C_EXT).$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

# GCC 2.95.4: create dependencies of F77 source files
$(DEP_DIR)/%.$(F77_EXT).$(DEP_EXT): %.$(F77_EXT)
	@$(ECHO) "==> $(APP_NAME) : Determining Dependencies of $<"
	set -e; $(F77) -M $(FPP_FLAGS) $< | $(SED) 's!\($*\)\.$(OBJ_EXT)[ :]*!\1.$(F77_EXT).$(OBJ_EXT) $@ : !g' > $@; [ -s $@ ] || $(RM) -f $@

else

# create dependencies of precompiled Oracle files
$(DEP_DIR)/%.$(ORACLE_EXT).$(DEP_EXT): $(PC_DIR)/%.$(ORACLE_EXT).$(CXX_EXT)
	@$(ECHO) "==> $(APP_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -MM $(CPP_FLAGS) -MT $(OBJ_DIR)/$*.$(ORACLE_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

# create dependencies of C++ source files
$(DEP_DIR)/%.$(CXX_EXT).$(DEP_EXT): %.$(CXX_EXT)
	@$(ECHO) "==> $(APP_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -MM $(CPP_FLAGS) -MT $(OBJ_DIR)/$*.$(CXX_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

# create dependencies of C source files
$(DEP_DIR)/%.$(C_EXT).$(DEP_EXT): %.$(C_EXT)
	@$(ECHO) "==> $(APP_NAME) : Determining Dependencies of $<"
	set -e; $(CXX) -MM $(CPP_FLAGS) -MT $(OBJ_DIR)/$*.$(C_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

# create dependencies of F77 source files
$(DEP_DIR)/%.$(F77_EXT).$(DEP_EXT): %.$(F77_EXT)
	@$(ECHO) "==> $(APP_NAME) : Determining Dependencies of $<"
	set -e; $(F77) -MM $(FPP_FLAGS) -MT $(OBJ_DIR)/$*.$(F77_EXT).$(OBJ_EXT) -MF $@ $<; [ -s $@ ] || $(RM) -f $@

endif


# compile C++ source files created by the Oracle precompiler
$(OBJ_DIR)/%.$(ORACLE_EXT).$(OBJ_EXT): \
		$(PC_DIR)/%.$(ORACLE_EXT).$(CXX_EXT) \
		$(DEP_DIR)/%.$(ORACLE_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(APP_NAME) : Compiling $<"
	$(CXX) $(APP_CXX_FLAGS) -o $@ -c $<


# compile Fortran 77 files created by KUIPC
$(OBJ_DIR)/%.$(CDF_EXT).$(OBJ_EXT): \
		$(PC_DIR)/%.$(CDF_EXT).$(F77_EXT)
	@$(ECHO) "==> $(APP_NAME) : Compiling $<" 
	$(F77) $(APP_F77_FLAGS) -o $@ -c $<


# compile C++ source files
$(OBJ_DIR)/%.$(CXX_EXT).$(OBJ_EXT): \
		%.$(CXX_EXT) $(DEP_DIR)/%.$(CXX_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(APP_NAME) : Compiling $<"
	$(CXX) $(APP_CXX_FLAGS) -o $@ -c $<


# compile Fortran 77 files
$(OBJ_DIR)/%.$(F77_EXT).$(OBJ_EXT): \
		%.$(F77_EXT) $(DEP_DIR)/%.$(F77_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(APP_NAME) : Compiling $<" 
	$(F77) $(APP_F77_FLAGS) -o $@ -c $<


# compile C source files
$(OBJ_DIR)/%.$(C_EXT).$(OBJ_EXT): \
		%.$(C_EXT) $(DEP_DIR)/%.$(C_EXT).$(DEP_EXT)
	@$(ECHO) "==> $(APP_NAME) : Compiling $<"
	$(CC) $(APP_C_FLAGS) -o $@ -c $<


# create Cint dictionary source file based on LinkDef.h file
$(PC_DIR)/%Dict.$(CXX_EXT): %LinkDef.$(HEADER_EXT) $(DICT_INPUT_FILE_NAMES)
	@$(ECHO) "==> $(APP_NAME) : Creating Dictionary $@"
	$(ROOTCINT) -f $@ $(ROOTCINT_FLAGS) $(filter $(dir $<)%.$(HEADER_EXT),$(DICT_INPUT_FILES)) $<


# compile Cint dictionary
$(OBJ_DIR)/%Dict.$(OBJ_EXT): $(PC_DIR)/%Dict.$(CXX_EXT)
	@$(ECHO) "==> $(APP_NAME) : Compiling $<"
	$(CXX) $(APP_CXX_FLAGS) -I$(PC_DIR) -o $@ -c $<


$(BUILD_DIR)/$(APP_NAME): $(OBJECTS)
	@$(ECHO) "==> $(APP_NAME) : Linking $@"
	$(LD) $(APP_LD_FLAGS) -o $@ -Wl,--start-group $(OBJECTS) -Wl,--end-group -Wl,--start-group $(HYDRA_LIBS) -Wl,--end-group $(LIBS)
# if one uses this line, possible circular dependencies are printed on stderr
#	$(LD) $(APP_LD_FLAGS) -o $@ `lorder $(OBJECTS) | tsort` `lorder $(HYDRA_LIBS) | tsort` $(LIBS)
