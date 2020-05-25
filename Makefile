#
# make -j 16 sopt popt ssmp psmp
#
# will now perform a parallel build of 4 cp2k executables
#
SHELL = /bin/sh
#
# the home dir is taken from the current directory.
ifeq ($(CP2KHOME),)
CP2KHOME     := $(abspath $(shell pwd))
export CP2KHOME
endif

ifneq ($(SPACK_COMPILER_SPEC),)
 # SPACK_COMPILER_SPEC is set when running in a Spack build-env
 ARCH         := $(shell spack arch)-$(shell echo $${SPACK_COMPILER_SPEC%%@*})
else
 ARCH         := local
endif

export VERSION=sopt

MAKEFILE     := $(CP2KHOME)/Makefile
ARCHDIR      := $(CP2KHOME)/arch
DOXYGENDIR   := $(CP2KHOME)/doc/doxygen
DATA_DIR     := $(CP2KHOME)/data
MAINEXEDIR   := $(CP2KHOME)/exe
MAINLIBDIR   := $(CP2KHOME)/lib
MAINOBJDIR   := $(CP2KHOME)/obj
MAINTSTDIR   := $(CP2KHOME)/regtesting
PRETTYOBJDIR := $(CP2KHOME)/obj/prettified
DOXIFYOBJDIR := $(CP2KHOME)/obj/doxified
TOOLSRC      := $(CP2KHOME)/tools
SRCDIR       := $(CP2KHOME)/src
EXEDIR       := $(MAINEXEDIR)/$(ARCH)
REVISION     := $(shell $(CP2KHOME)/tools/build_utils/get_revision_number $(SRCDIR))

EXTSDIR      := exts
EXTSHOME     := $(CP2KHOME)/$(EXTSDIR)
EXTSPACKAGES := $(shell cd $(EXTSHOME) ; find * -maxdepth 0 -type d )

PYTHON       := /usr/bin/env python3

# Common Targets ============================================================
default_target: all

# Discover programs =========================================================
ifeq ($(ALL_EXE_FILES),)
export ALL_EXE_FILES := $(sort $(shell $(TOOLSRC)/build_utils/discover_programs.py $(SRCDIR)))
endif
EXE_NAMES := $(basename $(notdir $(ALL_EXE_FILES)))

# Once we are down to a single version ======================================
# this only happens on stage 3 and 4
ifneq ($(ONEVERSION),)
MODDEPS = "lower"
include $(ARCHDIR)/$(ARCH).$(ONEVERSION)
LIBDIR  := $(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)
LIBEXTSDIR := $(LIBDIR)/$(EXTSDIR)
OBJDIR  := $(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)
OBJEXTSDIR := $(OBJDIR)/$(EXTSDIR)
OBJEXTSINCL := $(foreach dir,$(EXTSPACKAGES),-I'$(OBJEXTSDIR)/$(dir)')
TSTDIR     := $(MAINTSTDIR)/$(ARCH)/$(ONEVERSION)
ifeq ($(NVCC),)
EXE_NAMES := $(basename $(notdir $(filter-out %.cu, $(ALL_EXE_FILES))))
endif
ifneq ($(LD_SHARED),)
 ARCHIVE_EXT := .so
else
 ARCHIVE_EXT := .a
endif
include $(EXTSHOME)/Makefile.inc
endif

# Declare PHONY targets =====================================================
.PHONY : $(VERSION) $(EXE_NAMES) \
         dirs makedep default_target all \
         toolversions extversions extclean libcp2k cp2k_shell exts python-bindings \
         doxify doxifyclean \
         pretty prettyclean doxygen/clean doxygen \
         install clean realclean distclean mrproper help \
         test testbg testclean testrealclean \
         data \
	 $(EXTSPACKAGES)

# Discover files and directories ============================================
ALL_SRC_DIRS := $(shell find $(SRCDIR) -type d ! -name preprettify | awk '{printf("%s:",$$1)}')
ALL_PREPRETTY_DIRS = $(shell find $(SRCDIR) -type d -name preprettify)

ALL_PKG_FILES  = $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "PACKAGE")
OBJ_SRC_FILES  = $(shell cd $(SRCDIR); find . ! -path "*/preprettify/*" -name "*.F")
OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . ! -path "*/preprettify/*" ! -path "*/python*" -name "*.c")
OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . ! -path "*/preprettify/*" -name "*.cpp")
OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . ! -path "*/preprettify/*" -name "*.cxx")
OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . ! -path "*/preprettify/*" -name "*.cc")
ifneq ($(NVCC),)
OBJ_SRC_FILES += $(shell cd $(SRCDIR); find . ! -path "*/preprettify/*" -name "*.cu")
endif

# Included files used by Fypp preprocessor and standard includes
INCLUDED_SRC_FILES = $(filter-out base_uses.f90, $(notdir $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.f90")))

# Include also source files which won't compile into an object file
ALL_SRC_FILES  = $(strip $(subst $(NULL) .,$(NULL) $(SRCDIR),$(NULL) $(OBJ_SRC_FILES)))
ALL_SRC_FILES += $(filter-out base_uses.f90, $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.f90"))
ALL_SRC_FILES += $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.h")
ALL_SRC_FILES += $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.hpp")
ALL_SRC_FILES += $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.hxx")
ALL_SRC_FILES += $(shell find $(SRCDIR) ! -path "*/preprettify/*" -name "*.hcc")

ALL_OBJECTS        = $(addsuffix .o, $(basename $(notdir $(OBJ_SRC_FILES))))
ALL_EXE_OBJECTS    = $(addsuffix .o, $(EXE_NAMES))
ALL_NONEXE_OBJECTS = $(filter-out $(ALL_EXE_OBJECTS), $(ALL_OBJECTS))

# stage 1: Call make recursively with each element in $(VERSION) as a target,
#          The actual target is stored in ORIG_TARGET.
#          Make can then parallelize over multiple versions targets.
ifeq ($(ONEVERSION),)
ORIG_TARGET = default_target

fes:
	@+$(MAKE) --no-print-directory -f $(MAKEFILE) $(VERSION) ORIG_TARGET=graph

$(EXE_NAMES) all toolversions extversions extclean libcp2k cp2k_shell exts $(EXTSPACKAGES) python-bindings test testbg:
	@+$(MAKE) --no-print-directory -f $(MAKEFILE) $(VERSION) ORIG_TARGET=$@

# stage 2: Store the version target in $(ONEVERSION),
#          Call make recursively with $(ORIG_TARGET) as target.
$(VERSION) :
	@+$(MAKE) --no-print-directory -f $(MAKEFILE) $(ORIG_TARGET) ORIG_TARGET="" VERSION="" ONEVERSION=$@

else

# stage 3: Include arch-file, create dirs, and run makedep.py for given $(ONEVERSION).
#          Afterwards, call make recursively again with -C $(OBJDIR) and INCLUDE_DEPS=true
ifeq ($(INCLUDE_DEPS),)
$(EXE_NAMES): makedep | dirs exts
	@+$(MAKE) --no-print-directory -C $(OBJDIR) -f $(MAKEFILE) $(EXEDIR)/$@.$(ONEVERSION) INCLUDE_DEPS=true

all: makedep | dirs exts
	@+$(MAKE) --no-print-directory -C $(OBJDIR) -f $(MAKEFILE) all INCLUDE_DEPS=true

# foreground testing, compilation happens in do_regtest
test: dirs
	cd $(TSTDIR); $(TOOLSRC)/regtesting/do_regtest -quick -arch $(ARCH) -version $(ONEVERSION) -cp2kdir ../../../  $(TESTOPTS)

# background testing, compilation happens here
testbg: dirs makedep all
	@+$(MAKE) --no-print-directory -C $(TSTDIR) -f $(MAKEFILE) testbg INCLUDE_DEPS=true

libcp2k: makedep | dirs exts
	@+$(MAKE) --no-print-directory -C $(OBJDIR) -f $(MAKEFILE) $(LIBDIR)/libcp2k$(ARCHIVE_EXT) INCLUDE_DEPS=true

python-bindings: libcp2k
	@cd $(SRCDIR)/start/python ; \
		env CC='$(CC)' LDSHARED='$(LD) -shared' CFLAGS='$(CFLAGS)' LDFLAGS='$(LDFLAGS) $(LDFLAGS_C) $(LIBS)' \
			$(PYTHON) setup.py build_ext \
				--build-temp="$(OBJDIR)/python" \
				--build-lib="$(LIBDIR)/python" \
				--library-dirs="$(LIBDIR)" \
				--libraries="$(patsubst -l%,%,$(filter -l%,$(LIBS)))"

exts: $(EXTSPACKAGES)

dirs:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(OBJEXTSDIR)
	@mkdir -p $(LIBDIR)
	@mkdir -p $(LIBEXTSDIR)
	@mkdir -p $(EXEDIR)
	@mkdir -p $(TSTDIR)

toolversions:
ifneq ($(FC),)
	@echo "=========== FC ($(ONEVERSION)) ==========="
ifeq (Cray,$(shell $(CC) -V 2>&1 | head -n1 | cut -d' ' -f1))
	$(FC) -V
else ifeq (IBM,$(shell $(CC) -qversion 2>&1 | head -n1 | cut -d' ' -f1))
	$(FC) -qversion
else
	$(FC) --version
endif
endif
ifneq ($(CC),)
	@echo "=========== CC ($(ONEVERSION)) ==========="
ifeq (Cray,$(shell $(CC) -V 2>&1 | head -n1 | cut -d' ' -f1))
	$(CC) -V
else ifeq (IBM,$(shell $(CC) -qversion 2>&1 | head -n1 | cut -d' ' -f1))
	$(CC) -qversion
else
	$(CC) --version
endif
endif
ifneq ($(NVCC),)
	@echo "========== NVCC ($(ONEVERSION)) =========="
	$(NVCC) --version
	@echo ""
endif
ifneq ($(AR),)
	@echo "=========== AR ($(ONEVERSION)) ==========="
	$(firstword $(AR)) V
	@echo ""
endif
	@echo "========== Make ($(ONEVERSION)) =========="
	$(MAKE) --version
	@echo ""
	@echo "========= Python ($(ONEVERSION)) ========="
	$(PYTHON) --version

else

# Force CP2K recompilations if exts are updated
$(ALL_OBJECTS): $(EXTSDEPS_MOD)
$(ALL_EXE_OBJECTS): $(EXTSDEPS_LIB)

# stage 4: Include $(OBJDIR)/all.dep, expand target all and libcp2k, and perform actual build.
all: $(foreach e, $(EXE_NAMES) cp2k_shell, $(EXEDIR)/$(e).$(ONEVERSION))
$(LIBDIR)/libcp2k$(ARCHIVE_EXT) : $(ALL_NONEXE_OBJECTS)

$(EXEDIR)/cp2k_shell.$(ONEVERSION): $(EXEDIR)/cp2k.$(ONEVERSION)
	cd $(EXEDIR); ln -sf cp2k.$(ONEVERSION) cp2k_shell.$(ONEVERSION)

testbg:
	@echo "testing: $(ONEVERSION) : full log in $(TSTDIR)/regtest.log "
	@$(TOOLSRC)/regtesting/do_regtest -nobuild $(ARCH) -version $(ONEVERSION) -cp2kdir ../../../  $(TESTOPTS) >& $(TSTDIR)/regtest.log
	@cat `grep 'regtesting location error_summary file:' $(TSTDIR)/regtest.log | awk '{print $$NF}'`
	@cat `grep 'regtesting location summary file:' $(TSTDIR)/regtest.log | awk '{print $$NF}'`
	@grep "Number of FAILED  tests 0" $(TSTDIR)/regtest.log >& /dev/null
	@grep "Number of WRONG   tests 0" $(TSTDIR)/regtest.log >& /dev/null

endif
endif

OTHER_HELP += "test : run the regression tests"
OTHER_HELP += "testbg : run the regression tests in background"

OTHER_HELP += "toolversions : Print versions of build tools"

OTHER_HELP += "extversions : Print versions of external modules"
OTHER_HELP += "extclean : Clean build of external modules"

#   extract help text from doxygen "\brief"-tag
help:
	@echo "=================== Binaries ===================="
	@echo "all                         Builds all executables (default target)"
	@for i in $(ALL_EXE_FILES); do \
	basename  $$i | sed 's/^\(.*\)\..*/\1/' | awk '{printf "%-28s", $$1}'; \
	grep "brief" $$i | head -n 1 | sed 's/^.*\\brief\s*//'; \
	done
	@echo "libcp2k                     Builds CP2K as a single library archive"
	@echo "cp2k_shell                  Creates symlink for backward compatibility"
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
# delete the intermediate files, but not the libraries and executables, or created directories.
# Most useful to save space on the disk or e.g. for recompiles with PGO that still needs the .gcda files in the objdir
#
# To avoid printing out all non-gcda files within a directory with 'make clean',
# a list of all non-gcda file extensions within a directory is first created with the get_extensions macro,
# and a wildcard expression is then used to remove the appropriate files
# cleaning stuff ============================================================
define get_extensions
	$(shell test -d $(1) && find $(1) -type f -name "*.*" ! -name "*.gcda" | sed 's|.*\.||' | sort -u)
endef
clean:
	@echo rm -rf $(foreach v, $(VERSION), $(MAINOBJDIR)/$(ARCH)/$(v))
	@$(foreach v, $(VERSION), $(foreach ext, $(call get_extensions, $(MAINOBJDIR)/$(ARCH)/$(v)/), $(shell rm -rf $(MAINOBJDIR)/$(ARCH)/$(v)/*.$(ext))))
	rm -rf $(foreach v, $(VERSION), $(MAINLIBDIR)/$(ARCH)/$(v))
OTHER_HELP += "clean : Remove intermediate object and mod files, but not the libraries and executables, for given ARCH and VERSION"

execlean:
	rm -rf $(foreach v, $(VERSION), $(EXEDIR)/*.$(v))
OTHER_HELP += "execlean : Remove the executables, for given ARCH and VERSION"

#
# delete the intermediate files, the programs and libraries and anything that might be in the objdir or libdir directory
# Use this if you want to fully rebuild an executable (for a given compiler and or VERSION)
#
realclean: extclean clean execlean
	rm -rf $(foreach v, $(VERSION), $(MAINOBJDIR)/$(ARCH)/$(v))
	rm -rf $(foreach v, $(VERSION), $(MAINLIBDIR)/$(ARCH)/$(v))
OTHER_HELP += "realclean : Remove all files for given ARCH and VERSION"

testclean:
	rm -rf $(foreach v, $(VERSION), $(MAINTSTDIR)/$(ARCH)/$(v)/TEST-*)
OTHER_HELP += "testclean : Remove all TEST-* files for given ARCH and VERSION"

testrealclean: testclean
	rm -rf $(foreach v, $(VERSION), $(MAINTSTDIR)/$(ARCH)/$(v)/LAST-*)
OTHER_HELP += "testrealclean : Remove all LAST-* and TEST-* files for given ARCH and VERSION"

#
# Remove all files from previous builds
#
distclean: prettyclean doxifyclean testrealclean
	rm -rf $(DOXYGENDIR) $(MAINEXEDIR) $(MAINOBJDIR) $(MAINLIBDIR) $(MAINTSTDIR)
OTHER_HELP += "distclean : Remove all files from previous builds"

# Prettyfier stuff ==========================================================
vpath %.pretty $(PRETTYOBJDIR)

pretty: $(addprefix $(PRETTYOBJDIR)/, $(ALL_OBJECTS:.o=.pretty)) $(addprefix $(PRETTYOBJDIR)/, $(INCLUDED_SRC_FILES:.f90=.pretty_included))
TOOL_HELP += "pretty : Reformat all source files in a pretty way."

prettyclean:
	-rm -rf $(PRETTYOBJDIR) $(ALL_PREPRETTY_DIRS)
TOOL_HELP += "prettyclean : Remove prettify marker files and preprettify directories"

$(PRETTYOBJDIR)/%.pretty: %.F $(DOXIFYOBJDIR)/%.doxified
	@mkdir -p $(PRETTYOBJDIR)
	cd $(dir $<); $(TOOLSRC)/prettify/prettify.py --do-backup --backup-dir=$(PRETTYOBJDIR) $(notdir $<)
	@touch $@

$(PRETTYOBJDIR)/%.pretty_included: %.f90 $(DOXIFYOBJDIR)/%.doxified_included
	@mkdir -p $(PRETTYOBJDIR)
	cd $(dir $<); $(TOOLSRC)/prettify/prettify.py --do-backup --backup-dir=$(PRETTYOBJDIR) $(notdir $<)
	@touch $@

$(PRETTYOBJDIR)/%.pretty: %.c $(DOXIFYOBJDIR)/%.doxified
#   TODO: call indent here?
	@mkdir -p $(PRETTYOBJDIR)
	@touch $@

$(PRETTYOBJDIR)/%.pretty: %.cpp $(DOXIFYOBJDIR)/%.doxified
#   TODO: call indent here?
	@mkdir -p $(PRETTYOBJDIR)
	@touch $@

# Doxyifier stuff ===========================================================
vpath %.doxified $(DOXIFYOBJDIR)

doxify: $(addprefix $(DOXIFYOBJDIR)/, $(ALL_OBJECTS:.o=.doxified)) $(addprefix $(DOXIFYOBJDIR)/, $(INCLUDED_SRC_FILES:.f90=.doxified_included))
TOOL_HELP += "doxify : Autogenerate doxygen headers for subroutines"

doxifyclean:
	-rm -rf $(DOXIFYOBJDIR)
TOOL_HELP += "doxifyclean : Remove doxify marker files"

$(DOXIFYOBJDIR)/%.doxified: %.F
	$(TOOLSRC)/doxify/doxify.sh $<
	@mkdir -p $(DOXIFYOBJDIR)
	@touch $@

$(DOXIFYOBJDIR)/%.doxified_included: %.f90
	$(TOOLSRC)/doxify/doxify.sh $<
	@mkdir -p $(DOXIFYOBJDIR)
	@touch $@

$(DOXIFYOBJDIR)/%.doxified: %.c
	@mkdir -p $(DOXIFYOBJDIR)
	@touch $@

$(DOXIFYOBJDIR)/%.doxified: %.cpp
	@mkdir -p $(DOXIFYOBJDIR)
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
	@sed -e "s/#revision#/$(REVISION)/" $(TOOLSRC)/doxify/Doxyfile.template >$(DOXYGENDIR)/Doxyfile
	cd $(DOXYGENDIR); doxygen ./Doxyfile 2>&1 | tee ./html/doxygen.out
TOOL_HELP += "doxygen : Generate the doxygen documentation"

# data stuff ================================================================
data: data/POTENTIAL
	@:

data/POTENTIAL: data/GTH_POTENTIALS data/HF_POTENTIALS data/NLCC_POTENTIALS data/ALL_POTENTIALS
	@echo "(re-)generating $@ ..."
	@cat $^ > $@

OTHER_HELP += "data : (re-)generate merged data files (e.g. data/POTENTIALS)"

# automatic dependency generation ===========================================
MAKEDEPMODE = "normal"
ifeq ($(HACKDEP),yes)
MAKEDEPMODE = "hackdep"
else
 ifneq ($(MC),)
 MAKEDEPMODE = "mod_compiler"
 endif
endif

# this happens on stage 3
makedep: $(ALL_SRC_FILES) $(ALL_PKG_FILES) dirs
ifeq ($(LD_SHARED),)
	@echo "Removing stale archives for $(ONEVERSION) ... "
	@$(TOOLSRC)/build_utils/check_archives.py $(firstword $(AR)) $(SRCDIR) $(LIBDIR)
endif
	@echo "Resolving dependencies for $(ONEVERSION) ... "
	@$(TOOLSRC)/build_utils/makedep.py $(OBJDIR)/all.dep cp2k $(MODDEPS) $(MAKEDEPMODE) $(ARCHIVE_EXT) $(SRCDIR) $(OBJ_SRC_FILES)

# on stage 4, load the rules generated by makedep.py
ifeq ($(INCLUDE_DEPS), true)
include $(OBJDIR)/all.dep
endif


# ================= Stuff need for compiling (stage 4) ======================
# These rules are executed in a recursive call to make -C $(OBJDIR)
# The change of $(CURDIR) allows to find targets without abs paths and vpaths.


### Slave rules ###
vpath %.F     $(ALL_SRC_DIRS)
vpath %.h     $(ALL_SRC_DIRS)
vpath %.f90   $(ALL_SRC_DIRS)
vpath %.cu    $(ALL_SRC_DIRS)
vpath %.c     $(ALL_SRC_DIRS)
vpath %.cpp   $(ALL_SRC_DIRS)
vpath %.cxx   $(ALL_SRC_DIRS)
vpath %.cc    $(ALL_SRC_DIRS)

#
# Add additional dependency of cp2k_info.F to git-HEAD.
# Ensuring that cp2k prints the correct source code revision number in its banner.
#
GIT_REF := ${MAINOBJDIR}/git-ref

# use a force (fake) target to always rebuild this file but have Make consider this updated
# iff it was actually rewritten (a .PHONY target is always considered new)
$(GIT_REF): FORCE
	echo $(REVISION) > "$@.tmp"
	@cmp "$@.tmp" "$@" || mv -f "$@.tmp" "$@"

FORCE: ;

cp2k_info.o: $(GIT_REF)

# Add some practical metadata about the build.
FCFLAGS += -D__COMPILE_ARCH="\"$(ARCH)\""\
           -D__COMPILE_DATE="\"$(shell date)\""\
           -D__COMPILE_HOST="\"$(shell hostname)\""\
           -D__COMPILE_REVISION="\"$(strip $(REVISION))\""\
           -D__DATA_DIR="\"$(DATA_DIR)\""

# $(FCLOGPIPE) can be used to store compiler output, e.g. warnings, for each F-file separately.
# This is used e.g. by the convention checker.

FYPPFLAGS ?= -n

%.o: %.F
	$(TOOLSRC)/build_utils/fypp $(FYPPFLAGS) $< $*.F90
	$(FC) -c $(FCFLAGS) -D__SHORT_FILE__="\"$(subst $(SRCDIR)/,,$<)\"" -I'$(dir $<)' $(OBJEXTSINCL) $*.F90 $(FCLOGPIPE)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $<

ifneq ($(LIBDIR),)
$(LIBDIR)/%:
ifneq ($(LD_SHARED),)
	@echo "Creating shared library $@"
	@$(LD_SHARED) $(LDFLAGS) -o $(@:.a=.so) $^ $(LIBS)
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
	$(MC) -c $(FCFLAGS) -D__SHORT_FILE__="\"$(subst $(SRCDIR)/,,$<)\"" $<
endif

#EOF
