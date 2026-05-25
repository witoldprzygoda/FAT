###############################################################################
#
# Makefile hades.rules.mk - to be included from the global Makefiles.
#
# This makefile is part of the Hades build system. It contains all default
# rules to build all Hydra modules. To use this makefile, include it in your
# global Makefile which lists the modules which should be build together.
#
###############################################################################

# loop (in non-parallel mode) over all given submodules and change to the
# related directories, too
define LoopOverModules
	@set -e;				\
	for module in $(MODULES); do		\
		$(MAKE) -C $$module $@; 	\
	done
endef


.PHONY:	all check pre-build depend build install \
	doc package deinstall clean distclean

.NOTPARALLEL: check pre-build install deinstall doc clean distclean 


### alias targets


all: build


### normal targets called by the user


check:
   ifneq ($(strip $(call pathsearch,perl)),)
	@$(ECHO) -e "\033[1mChecking for duplicate source code files:\033[0m"
	@perl -e '							\
	use File::Basename;						\
	use File::Find;							\
	find( { wanted => \&processFile, no_chdir => 1 }, "." );	\
	for $$File (sort keys %Files)					\
	{								\
	   if ($$#{$$Files{$$File}} > 0)				\
	   {								\
	      print "$$File found in\n";				\
	      for $$i (0..$$#{$$Files{$$File}})				\
	      {								\
		 print "  $$Files{$$File}[$$i]\n";			\
	      }								\
	   }								\
	}								\
	sub processFile()						\
	{								\
	   return if (-d $$File::Find::name);				\
	   $$filename = basename( $$File::Find::name );			\
	   return unless ($$filename =~					\
	       /(\.$(HEADER_EXT)|\.$(CXX_EXT)|\.$(ORACLE_EXT))$$/);	\
	   $$directory = dirname( $$File::Find::name );			\
	   return if ($$directory eq "$(INSTALL_DIR)/include");		\
	   return if ($$directory eq "$(BUILD_DIR)/include");		\
	   push( @{$$Files{"$$filename"}}, $$directory );		\
	} '
   else
	@$(ECHO) -n "HINT: Skipping checks on duplicate source code ";	\
	$(ECHO) "files skipped (Perl not found)."
   endif
        # loop over sub-modules without printing the directory change
	@set -e;						\
	for module in $(MODULES); do				\
		$(MAKE) --no-print-directory -C $$module $@;	\
	done

	@$(ECHO) -e "\n==> Consistency checks done.\n"


# create directories before several jobs might are forked,
# otherwise directory creation can lead to a race condition
pre-build:
	@set -e;					\
	[ -e $(BUILD_DIR) ]  ||  mkdir $(BUILD_DIR);	\
	[ -e $(DEP_DIR) ]    ||  mkdir $(DEP_DIR);	\
	[ -e $(PC_DIR) ]     ||  mkdir $(PC_DIR);	\
	[ -e $(OBJ_DIR) ]    ||  mkdir $(OBJ_DIR);	\
	[ -e $(LIB_DIR) ]    ||  mkdir $(LIB_DIR)


# this target is used to loop over all modules to build the libraries and
# to include dependencies just here (using this method is necessary, since
# it is a must to do a direct fork of a sub-make to get "-j" working - so
# one cannot use a for-loop running in a shell)
%_depend %_build: pre-build
	@$(MAKE) -C $(subst _,$(space),$@)


# just create all dependency files
depend: $(addsuffix _depend,$(MODULES))
	@$(ECHO) -e "\n==> Dependency files created.\n"


# build shared libraries
build: $(addsuffix _build,$(MODULES))
	@$(ECHO) -e "\n==> All Hydra modules created.\n"


# create the class documentation using Root's inline code documentation
doc: $(CREATE_DOC_MACRO)
	@$(ROOT) -q -l $(CREATE_DOC_MACRO)
	@$(ECHO) -e "\n==> Hydra class documentation created.\n"


# install all neccessary header and library files
install:
	@[ -e $(INSTALL_DIR) ]         ||  mkdir $(INSTALL_DIR);	\
	[ -e $(INSTALL_DIR)/lib ]      ||  mkdir $(INSTALL_DIR)/lib;	\
	[ -e $(INSTALL_DIR)/include ]  ||  mkdir $(INSTALL_DIR)/include;\
	[ -e $(INSTALL_DIR)/macros ]   ||  mkdir $(INSTALL_DIR)/macros;
	$(LoopOverModules)
	@for file in hades.*.mk; do					\
	if [ -f $$file ]; then						\
		$(ECHO) "$(INSTALL) -m 644 $$file $(INSTALL_DIR)";      \
		$(INSTALL) -m 644 $$file $(INSTALL_DIR);                \
	fi;								\
	done
	@for file in macros/rootlog*.C; do                              \
	if [ -f $$file ]; then                                          \
		$(ECHO) "$(INSTALL) -m 644 $$file $(INSTALL_DIR)/macros";\
		$(INSTALL) -m 644 $$file $(INSTALL_DIR)/macros;         \
	fi;								\
	done
	@content=`$(ECHO) $(DOC_DIR)/*`;				\
	if [ "$$content" = "$(DOC_DIR)/*" ]; then			\
		$(ECHO) "HINT: No documents found (nothing installed).";\
	else								\
		$(CP) -r $(DOC_DIR)/* $(INSTALL_DIR)/doc;		\
	fi
	@$(ECHO) -e "\n==> Installation to $(INSTALL_DIR) complete.\n"


# deinstall all header and library files
deinstall:
	$(LoopOverModules)
	- $(RMDIR) $(INSTALL_DIR)/lib
	- $(RMDIR) $(INSTALL_DIR)/include
	- $(RM) -r $(INSTALL_DIR)/doc
	@$(ECHO) -e "\n==> Deinstallation from $(INSTALL_DIR) complete.\n"


package:
	@$(ECHO) -e "\n==> Sorry, not yet implemented!\n"


# delete all but the source and dependency files
clean:
	$(LoopOverModules)
	@$(ECHO) -e "\n==> All files but dependency files deleted.\n"


# delete all but the source files
distclean:
	$(LoopOverModules)
	$(RM) -f *~
	- $(RMDIR) $(OBJ_DIR) $(LIB_DIR) $(DEP_DIR) $(PC_DIR) $(BUILD_DIR)
	@$(ECHO) -e "\n==> All files but sources deleted.\n"
