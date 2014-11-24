#
# make -j 16 sopt popt ssmp psmp
#
# will now perform a parallel build of 4 cp2k executables
#
.SUFFIXES: .F .dep .o .f90 .cu .c .instantiation .instantiated .template
SHELL = /bin/sh
#
# the home dir is taken from the current directory: you must cd to the makefiles
# directory before executing make.
# adding a symbolic link from the src directory to the Makefile in makefiles/
# also works fine
#
CP2KHOME     := $(abspath $(PWD)/..)
ARCH         := $(shell $(CP2KHOME)/tools/build_utils/get_arch_code)
export VERSION=sopt

MAKEFILE     := $(CP2KHOME)/makefiles/Makefile
DOXYGENDIR   := $(CP2KHOME)/doc/doxygen
MAINEXEDIR   := $(CP2KHOME)/exe
MAINLIBDIR   := $(CP2KHOME)/lib
MAINOBJDIR   := $(CP2KHOME)/obj
PRETTYOBJDIR := $(CP2KHOME)/obj/prettified
DOXIFYOBJDIR := $(CP2KHOME)/obj/doxified
INSTANTIATEDDIR := $(CP2KHOME)/obj/instantiated
TOOLSRC      := $(CP2KHOME)/tools
SRCDIR       := $(CP2KHOME)/src
EXEDIR       := $(MAINEXEDIR)/$(ARCH)
REVISION     := $(shell $(CP2KHOME)/tools/build_utils/get_revision_number $(SRCDIR))


# Discover programs =========================================================
ifeq ($(ALL_EXE_FILES),)
export ALL_EXE_FILES := $(sort $(shell $(TOOLSRC)/build_utils/discover_programs.py $(SRCDIR)))
endif
EXE_NAMES := $(basename $(notdir $(ALL_EXE_FILES)))

# Once we are down to a single version ======================================
# this only happens on stage 3 and 4
ifneq ($(ONEVERSION),)
MODDEPS = "lower"
include $(CP2KHOME)/arch/$(ARCH).$(ONEVERSION)
LIBDIR  := $(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)
OBJDIR  := $(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)
ifeq ($(NVCC),)
EXE_NAMES := $(basename $(notdir $(filter-out %.cu, $(ALL_EXE_FILES))))
endif
endif

# Declare PHONY targets =====================================================
.PHONY : $(VERSION) $(EXE_NAMES) \
         dirs makedep default_target all \
         doxifydir doxify doxifyclean \
         prettydir pretty prettyclean doxygen/clean doxygen \
         install clean realclean distclean mrproper help templates instantiationdir templatesclean


# Discover files and directories ============================================
ALL_SRC_DIRS := $(shell find $(SRCDIR) -type d ! -name preprettify  ! -path "*/.svn*" | awk '{printf("%s:",$$1)}')
LIBCUSMM_DIR := $(shell find $(SRCDIR) -type d -name "libcusmm")
ALL_PREPRETTY_DIRS = $(shell find $(SRCDIR) -type d -name preprettify)

ALL_PKG_FILES  = $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "PACKAGE")
OBJ_SRC_FILES  = $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.F")
OBJ_SRC_FILES += $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.c")
ALL_INST_FILES = $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.instantiation")

ifneq ($(NVCC),)
OBJ_SRC_FILES += $(shell find $(SRCDIR) ! -path "*/preprettify/*" ! -name "libcusmm.cu" -name "*.cu")
OBJ_SRC_FILES += $(LIBCUSMM_DIR)/libcusmm.cu
endif

# Include also source files which won't compile into an object file
ALL_SRC_FILES  = $(OBJ_SRC_FILES)
ALL_SRC_FILES += $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.f90")
ALL_SRC_FILES += $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.h")

ALL_OBJECTS        = $(addsuffix .o, $(basename $(notdir $(OBJ_SRC_FILES))))
ALL_EXE_OBJECTS    = $(addsuffix .o, $(EXE_NAMES))
ALL_NONEXE_OBJECTS = $(filter-out $(ALL_EXE_OBJECTS), $(ALL_OBJECTS))


# Common Targets ============================================================
default_target: all

# stage 1: Call make recursively with each element in $(VERSION) as a target,
#          The actual target is stored in ORIG_TARGET.
#          Make can then parallelize over multiple versions targets.
ifeq ($(ONEVERSION),)
ORIG_TARGET = default_target

fes :
	@$(MAKE) --no-print-directory -f $(MAKEFILE) $(VERSION) ORIG_TARGET=graph

$(EXE_NAMES) all toolversions:
	@$(MAKE) --no-print-directory -f $(MAKEFILE) $(VERSION) ORIG_TARGET=$@

# stage 2: Store the version target in $(ONEVERSION),
#          Call make recursively with $(ORIG_TARGET) as target.
$(VERSION) :
	@$(MAKE) --no-print-directory -f $(MAKEFILE) $(ORIG_TARGET) ORIG_TARGET="" VERSION="" ONEVERSION=$@

else

# stage 3: Include arch-file, create dirs, and run makedep.py for given $(ONEVERSION).
#          Afterwards, call make recursively again with -C $(OBJDIR) and INCLUDE_DEPS=true
ifeq ($(INCLUDE_DEPS),)
$(EXE_NAMES): dirs makedep
	@$(MAKE) --no-print-directory -C $(OBJDIR) -f $(MAKEFILE) $(EXEDIR)/$@.$(ONEVERSION) INCLUDE_DEPS=true

all: dirs makedep
	@$(MAKE) --no-print-directory -C $(OBJDIR) -f $(MAKEFILE) all INCLUDE_DEPS=true

dirs:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(EXEDIR)

toolversions:
ifneq ($(CC),)
	@echo "=========== CC ($(ONEVERSION)) ==========="
	$(CC) --version
endif
ifneq ($(FC),)
	@echo "=========== FC ($(ONEVERSION)) ==========="
	$(FC) --version
endif
ifneq ($(AR),)
	@echo "=========== AR ($(ONEVERSION)) ==========="
	$(firstword $(AR)) V
	@echo ""
endif
else

# stage 4: Include $(OBJDIR)/all.dep, expand target all, and perform actual build.
all: $(foreach e, $(EXE_NAMES), $(EXEDIR)/$(e).$(ONEVERSION))

endif
endif

OTHER_HELP += "toolversions : Print versions of build tools"

#   extract help text from doxygen "\brief"-tag
help:
	@echo "=================== Binaries ===================="
	@echo "all                         Builds all executables (default target)"
	@for i in $(ALL_EXE_FILES); do \
	basename  $$i | sed 's/^\(.*\)\..*/\1/' | awk '{printf "%-28s", $$1}'; \
	grep "brief" $$i | head -n 1 | sed 's/^.*\\brief\s*//'; \
	done
	@echo ""
	@echo "===================== Tools ====================="
	@printf "%s\n" $(TOOL_HELP) | awk -F ':' '{printf "%-28s%s\n", $$1, $$2}'
	@echo ""
	@echo "================= Other Targets ================="
	@printf "%s\n" $(OTHER_HELP) | awk -F ':' '{printf "%-28s%s\n", $$1, $$2}'
	@echo "help                         Print this help text"

#
# so far CP2K does not install, but give a hint to the user
#
install:
	@echo ""
	@echo "The CP2K executable is $(foreach v, $(VERSION), $(EXEDIR)/cp2k.$(v))"
	@echo ""
OTHER_HELP += "install : Print installation help"

#
# delete the intermediate files, but not the libraries and executables. 
# Most useful to save space on the disk
#
# cleaning stuff ============================================================
clean:
	rm -rf $(foreach v, $(VERSION), $(MAINOBJDIR)/$(ARCH)/$(v)) $(LIBCUSMM_DIR)/libcusmm.cu
OTHER_HELP += "clean : Remove intermediate files, but not the libraries and executables, for given ARCH and VERSION"

#
# delete the intermediate files, the programs and libraries and anything that might be in the objdir directory
# Use this if you want to rebuild the executable (for a given compiler and or VERSION)
#
realclean: clean
	rm -rf $(foreach v, $(VERSION), $(MAINLIBDIR)/$(ARCH)/$(v) $(EXEDIR)/*.$(v))
OTHER_HELP += "realclean : Remove all file for given ARCH and VERSION"

#
# Remove all files from previous builds
#
distclean: prettyclean doxifyclean
	rm -rf $(DOXYGENDIR) $(MAINEXEDIR) $(MAINOBJDIR) $(MAINLIBDIR) $(LIBCUSMM_DIR)/libcusmm.cu
OTHER_HELP += "distclean : Remove all files from previous builds"

#
# Remove ALL files not in SVN
# NOTE: This removes really everything not committed to SVN
#
mrproper:
	cd $(CP2KHOME); ./tools/clean_cwd.sh
OTHER_HELP += "mrproper : Remove ALL files not in SVN"

#
# instantiate the templates
#
vpath %.instantiation   $(ALL_SRC_DIRS)

# since the same named .instantiation file can be present in many directories, we need to copy the full tree
# to store the instantiated file in a unique location
# TODO: proper dependencies based on .template
instantiateddir:
	-rm -fR $(INSTANTIATEDDIR)
	@mkdir -p $(dir $(addprefix $(INSTANTIATEDDIR), $(ALL_INST_FILES:.instantiation=.instantiated)))

templates: instantiateddir $(addprefix $(INSTANTIATEDDIR), $(ALL_INST_FILES:.instantiation=.instantiated))
TOOL_HELP += "templates : instantiate all templates"

$(addprefix $(INSTANTIATEDDIR), %.instantiated): %.instantiation instantiateddir
	@cd $(dir $<); $(TOOLSRC)/instantiateTemplates.py $(notdir $<)
	@touch $@

# Prettyfier stuff ==========================================================
vpath %.pretty $(PRETTYOBJDIR)

prettydir:
	@mkdir -p $(PRETTYOBJDIR)

pretty: doxify prettydir $(addprefix $(PRETTYOBJDIR)/, $(ALL_OBJECTS:.o=.pretty))
TOOL_HELP += "pretty : Reformat all source files in a pretty way."

prettyclean:
	-rm -rf $(PRETTYOBJDIR) $(ALL_PREPRETTY_DIRS)
TOOL_HELP += "prettyclean : Remove prettify marker files and preprettify directories"

$(PRETTYOBJDIR)/%.pretty: %.F
	cd $(dir $<); $(TOOLSRC)/prettify.py $(notdir $<)
	@touch $@

$(PRETTYOBJDIR)/%.pretty: %.c
#   TODO: call indent here?
	@touch $@

# Doxyifier stuff ===========================================================
vpath %.doxified $(DOXIFYOBJDIR)

doxifydir:
	@mkdir -p $(DOXIFYOBJDIR)

doxify: doxifydir $(addprefix $(DOXIFYOBJDIR)/, $(ALL_OBJECTS:.o=.doxified))
TOOL_HELP += "doxify : Autogenerate doxygen headers for subroutines"

doxifyclean:
	-rm -rf $(DOXIFYOBJDIR)
TOOL_HELP += "doxifyclean : Remove doxify marker files"

$(DOXIFYOBJDIR)/%.doxified: %.F
	$(TOOLSRC)/doxify/doxify.sh $<
	@touch $@

$(DOXIFYOBJDIR)/%.doxified: %.c
	@touch $@

# doxygen stuff =============================================================
doxygen/clean:
	-rm -rf $(DOXYGENDIR)
TOOL_HELP += "doxygen/clean : Remove the generated doxygen documentation"

# Automatic source code documentation using Doxygen
# Prerequisites:
# - stable doxygen release 1.5.4 (Oct. 27, 2007)
# - graphviz (2.16.1)
# - webdot (2.16)
#
doxygen: doxygen/clean
	@mkdir -p $(DOXYGENDIR)
	@mkdir -p $(DOXYGENDIR)/html
	@echo "<html><body>Sorry, the Doxygen documentation is currently being updated. Please try again in a few minutes.</body></html>" > $(DOXYGENDIR)/html/index.html
	cp $(ALL_SRC_FILES) $(DOXYGENDIR)
	@for i in $(DOXYGENDIR)/*.F ; do mv $${i}  $${i%%.*}.f90; done ;
	@cat $(TOOLSRC)/doxify/Doxyfile.template | sed "s/#revision#/`$(TOOLSRC)/build_utils/get_revision_number $(CP2KHOME)`/"  >$(DOXYGENDIR)/Doxyfile
	cd $(DOXYGENDIR); doxygen ./Doxyfile 2>&1 | tee ./html/doxygen.out
TOOL_HELP += "doxygen : Generate the doxygen documentation"


# Libcusmm stuff ============================================================
$(LIBCUSMM_DIR)/libcusmm.cu: $(LIBCUSMM_DIR)/generate.py $(LIBCUSMM_DIR)/parameters.txt $(wildcard $(LIBCUSMM_DIR)/kernels/*.py)
	cd $(LIBCUSMM_DIR); ./generate.py $(LIBCUSMM_FLAGS)


# automatic dependency generation ===========================================
MAKEDEPMODE = "normal"
ifeq ($(HACKDEP),yes)
MAKEDEPMODE = "hackdep"
else
 ifneq ($(MC),)
 MAKEDEPMODE = "mod_compiler"
 endif
endif

ifneq ($(LD_SHARED),)
 ARCHIVE_EXT = ".$(ONEVERSION).so"
else
 ARCHIVE_EXT = ".a"
endif

# this happens on stage 3
makedep: $(ALL_SRC_FILES) $(ALL_PKG_FILES)
ifeq ($(LD_SHARED),)
	@echo "Removing stale archives for $(ONEVERSION) ... "
	@$(TOOLSRC)/build_utils/check_archives.py $(firstword $(AR)) $(SRCDIR) $(LIBDIR)
endif
	@echo "Resolving dependencies for $(ONEVERSION) ... "
	@$(TOOLSRC)/build_utils/makedep.py $(OBJDIR)/all.dep $(MODDEPS) $(MAKEDEPMODE) $(ARCHIVE_EXT) $(OBJ_SRC_FILES)

# on stage 4, load the rules generated by makedep.py
ifeq ($(INCLUDE_DEPS), true)
include $(OBJDIR)/all.dep
endif

# ================= Stuff need for compiling (stage 4) ======================
# These rules are executed in a recursive call to make -C $(OBJDIR)
# The change of $(CURDIR) allows to find targets without abs paths and vpaths.


### Slave rules ###
vpath %.F   $(ALL_SRC_DIRS)
vpath %.h   $(ALL_SRC_DIRS)
vpath %.f90 $(ALL_SRC_DIRS)
vpath %.cu  $(ALL_SRC_DIRS)
vpath %.c   $(ALL_SRC_DIRS)

#
# Add additional dependency of cp2k_info.F to SVN-entry or git-HEAD.
# Ensuring that cp2k prints the correct source code revision number in its banner.
#
SVN_ENTRY    := $(wildcard $(SRCDIR)/.svn/entries*)
ifneq ($(strip $(SVN_ENTRY)),)
cp2k_info.o: $(SVN_ENTRY)
endif

GIT_HEAD     := $(wildcard $(CP2KHOME)/../.git/HEAD*)
ifneq ($(strip $(GIT_HEAD)),)
cp2k_info.o: $(GIT_HEAD)
endif

# some practical variables for the build
ifeq ($(CPPSHELL),)
CPPSHELL := -D__COMPILE_ARCH="\"$(ARCH)\""\
            -D__COMPILE_DATE="\"$(shell date)\""\
            -D__COMPILE_HOST="\"$(shell hostname)\""\
            -D__COMPILE_REVISION="\"$(strip $(REVISION))\""
endif

ifneq ($(CPP),)
# always add the SRCDIR to the include path (-I here might not be portable) 
CPPFLAGS += $(CPPSHELL) -I$(SRCDIR)
else
FCFLAGS += $(CPPSHELL)
endif

# the rule how to generate the .o from the .F
# only if CPP is different from null we do a step over the C preprocessor (which is slower)
# in the other case the fortran compiler takes care of this directly
%.o: %.F
ifneq ($(CPP),)
	$(CPP) $(CPPFLAGS) -D__SHORT_FILE__="\"$(notdir $<)\"" $< > $*.f90
	$(FC) -c $(FCFLAGS) $*.f90 $(FCLOGPIPE)
else
	$(FC) -c $(FCFLAGS) -D__SHORT_FILE__="\"$(notdir $<)\"" $< $(FCLOGPIPE) 
endif

%.o: %.c
	$(CC) -c $(CFLAGS) $<

ifneq ($(LIBDIR),)
$(LIBDIR)/%:
ifneq ($(LD_SHARED),)
	@echo "Creating shared library $@"
	@$(LD_SHARED) -o $(@:.a=.so) $^
else
	@echo "Updating archive $@"
	@$(AR) $@ $?
endif
ifneq ($(RANLIB),)
	@$(RANLIB) $@
endif
endif

%.o: %.cu
	$(NVCC) -c $(NVFLAGS) $<


# module compiler magic =====================================================
ifeq ($(MC),)
#
# here we cheat... this tells make that .mod can be generated from .o (this holds in CP2K) by doing nothing
# it avoids recompilation if .o is more recent than .F, but .mod is older than .F
# (because it didn't change, as e.g. g95 can do)
#
# this is problematic if the module names are uppercase e.g. KINDS.mod (because this rule expands to kinds.mod)
#
%.mod: %.o
	@true
else
#
# if MC is defined, it is our 'module compiler' which generates the .mod file from the source file
# it is useful in a two-stage compile.
#
%.mod: %.F
	$(MC) -c $(FCFLAGS) -D__SHORT_FILE__="\"$(notdir $<)\"" $<
endif

#EOF
