#
# make -j 16 sopt popt ssmp psmp
#
# will now perform a parallel build of 4 cp2k executables
#
.SUFFIXES: .F .dep .o .f90 .cu .c
SHELL = /bin/sh
#
# the home dir is taken from the current directory: you must cd to the makefiles
# directory before executing make.
# adding a symbolic link from the src directory to the Makefile in makefiles/ 
# also works fine
#
CP2KHOME := $(PWD)/..
ARCH     := $(shell $(CP2KHOME)/tools/get_arch_code)
LIB1      = cp2k_base_lib
LIB2      = cp2k_fft_lib
LIB3      = cp2k_dbcsr_lib
LIB4      = cp2k_ma_lib
LIBCP2K   = cp2k_lib

export VERSION=sopt
export WHAT   =_progr

### Dependent variables ###

ARCHDIR               = $(CP2KHOME)/arch
DOCDIR                = $(CP2KHOME)/doc
MAINEXEDIR            = $(CP2KHOME)/exe
MAINLIBDIR            = $(CP2KHOME)/lib
MAINOBJDIR            = $(CP2KHOME)/obj
PRETTYOBJDIR          = $(CP2KHOME)/obj/prettified
TOOLDIR               = $(OBJDIR)/tools
TOOLSRC               = $(CP2KHOME)/tools
SRCDIR                = $(CP2KHOME)/src
LIBCUSMM_DIR          = $(SRCDIR)/dbcsr_lib/cuda/libcusmm
SRCDIRS               = $(SRCDIR):$(SRCDIR)/lib:$(SRCDIR)/arch:$(SRCDIR)/dbcsr_lib:$(SRCDIR)/dbcsr_lib/cuda:$(SRCDIR)/fft_lib:$(SRCDIR)/metadyn_tools:$(SRCDIR)/cuda:$(LIBCUSMM_DIR):$(SRCDIR)/machine
MAKEFILE              = $(CP2KHOME)/makefiles/Makefile
MACHINEDEFS           = $(ARCHDIR)/$(ARCH).$(VERSION)
PROG                  = $(EXEDIR)/cp2k.$(VERSION)
CP2KSHELL             = $(EXEDIR)/cp2k_shell.$(VERSION)
FES                   = $(EXEDIR)/fes.$(VERSION)
EXEDIR                = $(MAINEXEDIR)/$(ARCH)
LIBDIR                = $(MAINLIBDIR)/$(ARCH)
OBJDIR                = $(MAINOBJDIR)/$(ARCH)
OBJECTDEFS            = $(SRCDIR)/OBJECTDEFS
LIB1_OBJECTDEFS       = $(SRCDIR)/lib/OBJECTDEFS
LIB2_OBJECTDEFS       = $(SRCDIR)/dbcsr_lib/OBJECTDEFS
LIB2CUDA_OBJECTDEFS   = $(SRCDIR)/dbcsr_lib/cuda/OBJECTDEFS
LIB2CUDA_DEPENDENCIES = $(SRCDIR)/dbcsr_lib/cuda/DEPENDENCIES
LIB3_OBJECTDEFS       = $(SRCDIR)/fft_lib/OBJECTDEFS
LIB4_OBJECTDEFS       = $(SRCDIR)/arch/OBJECTDEFS
LIBNVOBJECTDEFS       = $(SRCDIR)/cuda/NVOBJECTDEFS
LIBMA_OBJECTDEFS      = $(SRCDIR)/machine/LIBMA_OBJECTDEFS
LIBCUDA_DEPENDENCIES  = $(SRCDIR)/cuda/LIBCUDA_DEPENDENCIES
LIBCUSMM_DEPENDENCIES = $(LIBCUSMM_DIR)/LIBCUSMM_DEPENDENCIES
LIBMA_DEPENDENCIES    = $(SRCDIR)/machine/LIBMA_DEPENDENCIES
LIB1_ARCHIVE          = $(LIBDIR)/$(VERSION)/lib$(LIB1).a
LIB2_ARCHIVE          = $(LIBDIR)/$(VERSION)/lib$(LIB2).a
LIB3_ARCHIVE          = $(LIBDIR)/$(VERSION)/lib$(LIB3).a
LIB4_ARCHIVE          = $(LIBDIR)/$(VERSION)/lib$(LIB4).a
LIB_CP2K_ARCHIVE      = $(LIBDIR)/$(VERSION)/lib$(LIBCP2K).a
SVN_ENTRY             = $(wildcard $(SRCDIR)/.svn/entries*)
GIT_HEAD              = $(wildcard $(CP2KHOME)/../.git/HEAD*)
FESOBJECTDEFS         = $(SRCDIR)/metadyn_tools/OBJECTDEFS
DOXYGENDIR            = $(DOCDIR)/doxygen
DBCSR_TEST_PROG       = $(EXEDIR)/dbcsr_test_driver.$(VERSION)
DBCSR_PERF_PROG       = $(EXEDIR)/dbcsr_performance_driver.$(VERSION)
#PLUMED
PLUMEDDIR             = $(SRCDIR)/src-plumed
PLUMEDINC             = $(PLUMEDDIR)/plumed.inc 
LIBCUSMM_FILES        = $(LIBCUSMM_DEPENDENCIES) $(LIBCUSMM_DIR)/libcusmm.cu $(wildcard $(LIBCUSMM_DIR)/kernels/*.pyc)

REVISION := $(shell $(CP2KHOME)/tools/get_revision_number $(SRCDIR))

# SMART_MODULES prevents unnecessary recompilation of files depending
# on unchanged modules.
# Set SMART_MODULES=yes in your environment, make argument, or ARCH file.
# Needs the md5sum program somewhere in the $PATH.
# Unneeded for gfortran >= 4.5.  No effect with compilers that output a
# different .mod files even for unchanged input files.
# Assumes file case preservation.
SMART_MODULES ?= no

### Definition of the multiple targets ###
#
# Should really be
#VERSION_TARGTES= sopt sdbg popt pdbg topt hopt tdbg hdbg
VERSION_TARGETS = sopt sdbg popt pdbg ssmp psmp tdbg hdbg
#
# uses an ugly trick, these are disguised as file names, the directory (@D) part 
# is used as a version, the file (@F) part is used as the target.
#
LIB_TARGETS = $(addsuffix /lib ,$(VERSION_TARGETS)) 
ALL_TARGETS = $(addsuffix /all ,$(VERSION_TARGETS))
CLEAN_TARGETS = $(addsuffix /clean ,$(VERSION_TARGETS)) $(addsuffix /realclean ,$(VERSION_TARGETS))

#
# lists the targets that are not files
#
.PHONY : _all libs _progr _lib build clean realclean distclean $(CLEAN_TARGETS)\
         $(ALL_TARGETS) $(LIB_TARGETS) $(VERSION_TARGETS) fes dirs deptool\
         dbcsr_test dbcsr_perf mrproper plumed libcusmm

### Master rules ###
build: dirs deptool libcusmm
	+$(MAKE) -C $(SRCDIR) -f $(MAKEFILE) all.dep
	+cp $(LIBCUDA_DEPENDENCIES) $(OBJDIR)/$(VERSION)/libc_cuda.dep
	+if [ -e $(LIBCUSMM_DEPENDENCIES) ]; then \
	  cp $(LIBCUSMM_DEPENDENCIES) $(OBJDIR)/$(VERSION)/libc_cusmm.dep; \
	fi
	+cp $(LIBMA_DEPENDENCIES) $(OBJDIR)/$(VERSION)/libc_ma.dep
	+$(MAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) $(WHAT)

fes: dirs deptool
	+$(MAKE) -C $(SRCDIR)/metadyn_tools/  -f $(MAKEFILE) fes.dep
	+$(MAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) HACKLIBDEP=yes VERSION=$(VERSION) $(FES) 

dirs:
	@mkdir -p $(EXEDIR)
	@mkdir -p $(LIBDIR)/$(VERSION)
	@mkdir -p $(TOOLDIR)
	@mkdir -p $(OBJDIR)/$(VERSION)
	@mkdir -p $(PRETTYOBJDIR)

realpretty: prettyclean pretty

pretty: dirs
	+$(MAKE) -C $(SRCDIR) -f $(MAKEFILE) _pretty
	+$(MAKE) -C $(SRCDIR)/lib       -f $(MAKEFILE) _pretty_lib1
	+$(MAKE) -C $(SRCDIR)/fft_lib   -f $(MAKEFILE) _pretty_lib2
	+$(MAKE) -C $(SRCDIR)/dbcsr_lib -f $(MAKEFILE) _pretty_lib3
	+$(MAKE) -C $(SRCDIR)/arch  -f $(MAKEFILE) _pretty_lib4

deptool: dirs
	+$(MAKE) LDFLAGS=$(MAKDEPF90_LDFLAGS) -C $(TOOLDIR) -f $(MAKEFILE) makedepf90

libs:
	+$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_lib build
all:
	+$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_all build
progr:
	+$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_progr build
cpshell:
	+$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_cpshell build

dbcsr_test:
	+$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_dbcsr_test build

dbcsr_perf:
	+$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_dbcsr_perf build

plumed: plumedobj plumed.dep
	+$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_plumed build

plumedobj: dirs
	+$(MAKE) -C $(OBJDIR)/$(VERSION)/ -f $(MAKEFILE) $(OBJ_PLUMED)

#
# on a parallel build only one task should go a build the deptool
#
$(VERSION_TARGETS): dirs deptool
	+$(MAKE) -f $(MAKEFILE) VERSION=$@ build

$(CLEAN_TARGETS):
	+$(MAKE) -f $(MAKEFILE) VERSION=$(@D) $(@F)

$(ALL_TARGETS) $(LIB_TARGETS):
	+$(MAKE) -f $(MAKEFILE) VERSION=$(@D) WHAT=_$(@F) build

#
# here we try to include all the objects that might need considering
#
include $(OBJECTDEFS)
include $(LIB1_OBJECTDEFS)
include $(LIB2_OBJECTDEFS)
include $(LIB3_OBJECTDEFS)
include $(LIB4_OBJECTDEFS)
include $(FESOBJECTDEFS)
include $(MACHINEDEFS)
include $(LIBNVOBJECTDEFS)
include $(LIBMA_OBJECTDEFS)
include $(LIB2CUDA_OBJECTDEFS)
include $(LIB2CUDA_DEPENDENCIES)

OBJECTS = $(OBJECTS_GENERIC)
FESOBJS = $(FESOBJ)

# Cuda specific instructions
ifneq ($(NVCC),)
  LIBNV_ARCHIVE= $(LIBDIR)/$(VERSION)/libcp2k_cuda.a
  LIBCUSMM_ARCHIVE = $(LIBDIR)/$(VERSION)/libcusmm.a
  LIB2CUDA_ARCHIVE = $(LIBDIR)/$(VERSION)/libdbcsr_cuda.a
endif


# libcusmm stuff
ifneq ($(NVCC),)
  libcusmm:
	cd $(LIBCUSMM_DIR); ./generate.py $(LIBCUSMM_FLAGS)
endif


# HWTOPO
ifneq ($(HWTOPO),)
  LIBMA_ARCHIVE= $(LIBDIR)/$(VERSION)/libcp2k_ma.a
endif

ifneq ($(SMART_MODULES),no)
MD5 := $(shell which md5sum)
endif

### Logging of fortran compiler output, e.g. for examining it's warnings ###
ifneq ($(FC_SAVELOG), )
   FCLOGPIPE = &> $*.log
endif

#
#PLUMED
#
ifeq ($(PLUMED),yes)
	include $(PLUMEDINC) 
	DFLAGS += -DPLUMED_CP2K
	CPPFLAGS += -DPLUMED_CP2K
endif


# Define the whole bunch of libraries needed
ALL_LIB =  $(LIB_CP2K_ARCHIVE) $(LIB3_ARCHIVE) $(LIB2_ARCHIVE) $(LIB1_ARCHIVE) $(LIB4_ARCHIVE) $(LIBMA_ARCHIVE) $(LIBNV_ARCHIVE) $(LIB2CUDA_ARCHIVE) $(LIBCUSMM_ARCHIVE)
ALL_SOURCE_FILES=$(notdir $(OBJECTS:.o=.F) $(LIB1_OBJECTS:.o=.F) $(LIB2_OBJECTS:.o=.F) $(LIB4_OBJECTS:.o=.F) $(LIB3_OBJECTS:.o=.F) cp_common_uses.h cp2k.F cp2k_shell.F)

#
# some practical variables for the build
#
ifeq ($(CPPSHELL),)
CPPSHELL := -D__COMPILE_ARCH="\"$(ARCH)\""\
            -D__COMPILE_DATE="\"$(shell date)\""\
            -D__COMPILE_HOST="\"$(shell hostname)\""\
	    -D__COMPILE_REVISION="\"$(strip $(REVISION))\""
endif

ifneq ($(CPP),)
#
# always add the SRCDIR to the include path (-I here might not be portable) 
#
CPPFLAGS += $(CPPSHELL) -I$(SRCDIR)
else
FCFLAGS += $(CPPSHELL)
endif


### Slave rules ###
vpath %.F $(SRCDIRS)
vpath %.h $(SRCDIRS)
vpath %.f90 $(SRCDIRS)
vpath %.dep $(OBJDIR)/$(VERSION)
vpath %.pretty $(PRETTYOBJDIR)
vpath %.cu $(SRCDIRS)
vpath %.c $(SRCDIRS)

#
# if the modules generated are not of the lower case form (e.g. kinds.F does not lead to kinds.mod)
# use MODDEPS = no in the arch file
# MODDEPS = upper does not work because of the rule '%.mod: %.o' is incorrect (see there), but if that is
# fixed, upper should work as well.
#
ifeq ($(MODDEPS),no)
MODSTRING=
else
ifeq ($(MODDEPS),upper)
MODSTRING=-m "%M.mod"
else
MODSTRING=-m "%m.mod"
endif
endif

#
# generation of the dependencies
#
all.dep: $(ALL_SOURCE_FILES)
ifeq ($(HACKDEP),yes)
# if HACKDEP=="yes" we hack the dependencies so that only modified source files are recompiled
# this is obviously not correct (i.e. can lead to a miscompiled CP2K), 
# but if you know what you are doing (i.e. you are sure that only information private to this module will be changed) might be useful nevertheless.
	echo $(ALL_SOURCE_FILES) | awk '{for (i=1;i<=NF;i++) {obj=$$(i) ; sub(".F",".o",obj); if ($$(i)!=obj) print obj,":",$$(i) }}' > $(OBJDIR)/$(VERSION)/all.dep
else
   ifeq ($(MC),)
	time $(TOOLDIR)/makedepf90 -I $(SRCDIRS) $(MODSTRING) -free $(ALL_SOURCE_FILES) > $(OBJDIR)/$(VERSION)/all.dep || ( rm -f $(OBJDIR)/$(VERSION)/all.dep  ; exit 1 )
   else
# if MC is defined, it is our 'module compiler' which generates the .mod file from the source file
# it is useful in a two-stage compile. Therefore we split the dependency in two parts the .o and the .mod.
# the awk below assumes our current rule of one mod per file
# The current module compiler is 'gfortran -fsyntax-only ' which should work fine for gfortran 4.8.0 and higher.
	$(TOOLDIR)/makedepf90 -I $(SRCDIRS) $(MODSTRING) -free $(ALL_SOURCE_FILES) | \
            awk '{deps=substr($$0,index($$0,":"),length($$0)-index($$0,":")+1); if ($$2==":") { print $$1,deps } else {print $$2,deps ; print $$1,deps,$$2 }}' > \
            $(OBJDIR)/$(VERSION)/all.dep || ( rm -f $(OBJDIR)/$(VERSION)/all.dep  ; exit 1 )
   endif
endif

# this is kind of a hack, fes depends on files in other directories, hence the -I and the -b.
fes.dep: $(FESOBJS:.o=.F) 
	$(TOOLDIR)/makedepf90 -I$(SRCDIR) -b. $(MODSTRING)  -free $^ > $(OBJDIR)/$(VERSION)/fes.dep || ( rm -f $(OBJDIR)/$(VERSION)/fes.dep  ; exit 1 )

plumed.dep: $(OBJECTS:.o=.F) $(LIB1_OBJECTS:.o=.F) $(LIB2_OBJECTS:.o=.F) $(LIB4_OBJECTS:.o=.F) $(LIB3_OBJECTS:.o=.F) cp_common_uses.h src-plumed/metadyn.h cp2k.F 

#ifeq ($(HACKDEP),yes)
#	echo $^ | awk '{for (i=1;i<=NF;i++) {obj=$$(i) ; sub(".F",".o",obj); if ($$(i)!=obj) print obj,":",$$(i) }}' > $(OBJDIR)/$(VERSION)/plumed.dep
#else
#	$(TOOLDIR)/makedepf90 $(MODSTRING)  -free $^ > $(OBJDIR)/$(VERSION)/plumed.dep || ( rm -f $(OBJDIR)/$(VERSION)/plumed.dep  ; exit 1 )
#endif
#
# and the rules doing the actual work
#
$(OBJ_PLUMED):
	$(CC) $(CPPFLAGS) -c $(PLUMEDDIR)/$(@:.o=.cpp)

_all: $(ALL_LIB) $(PROG) fes $(CP2KSHELL)
_progr: $(PROG)
_lib: $(ALL_LIB)
_cpshell: $(CP2KSHELL)
_dbcsr_test: $(DBCSR_TEST_PROG)
_dbcsr_perf: $(DBCSR_PERF_PROG)
_plumed: $(ALL_LIB) PROG-PLUMED

_pretty:$(OBJECTS:.o=.pretty)
_pretty_lib1:$(LIB1_OBJECTS:.o=.pretty)
_pretty_lib2:$(LIB2_OBJECTS:.o=.pretty)
_pretty_lib3:$(LIB3_OBJECTS:.o=.pretty)
_pretty_lib4:$(LIB4_OBJECTS:.o=.pretty)


$(LIB1_ARCHIVE): $(LIB1_OBJECTS)
	$(AR) $(LIB1_ARCHIVE) $?
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB1_ARCHIVE)
endif

$(LIB2_ARCHIVE): $(LIB2_OBJECTS)
	$(AR) $(LIB2_ARCHIVE) $?
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB2_ARCHIVE)
endif

$(LIB3_ARCHIVE): $(LIB3_OBJECTS)
	$(AR) $(LIB3_ARCHIVE) $?
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB3_ARCHIVE)
endif

$(LIB4_ARCHIVE): $(LIB4_OBJECTS)
	$(AR) $(LIB4_ARCHIVE) $?
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB4_ARCHIVE)
endif

#
# a hack to first compile the src/lib/ stuff before files from src/
# otherwise we have missing dependencies
#
ifneq ($(HACKLIBDEP),yes)
LIB_OBJECTS = $(LIB4_ARCHIVE) $(LIB3_ARCHIVE) $(LIB2_ARCHIVE) $(LIB1_ARCHIVE)
else
LIB_OBJECTS =
endif

$(LIB_CP2K_ARCHIVE): $(OBJECTS) sockets.o
	$(AR) $(LIB_CP2K_ARCHIVE) $?
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB_CP2K_ARCHIVE)
endif

ifneq ($(NVCC),)
  $(LIBNV_ARCHIVE): $(LIBNV_OBJECTS)
	$(AR) $(LIBNV_ARCHIVE) $?
 $(LIBCUSMM_ARCHIVE): libcusmm.o
	$(AR) $(LIBCUSMM_ARCHIVE) $?
 $(LIB2CUDA_ARCHIVE): $(DBCSR_CUDA_OBJECTS)
	$(AR) $(LIB2CUDA_ARCHIVE) $?
endif

ifneq ($(HWTOPO),)
  $(LIBMA_ARCHIVE): $(LIBMA_OBJECTS)
	$(AR) $(LIBMA_ARCHIVE) $?
endif


PROG-PLUMED: $(ALL_LIB) cp2k.o $(OBJ_PLUMED) 
	+$(LD) $(LDFLAGS) -o $(PROG) $(OBJ_PLUMED) cp2k.o $(ALL_LIB) $(LIBS) -lstdc++

$(PROG): $(ALL_LIB) cp2k.o
	+$(LD) $(LDFLAGS) -o $(PROG) cp2k.o $(ALL_LIB) $(LIBS)

$(CP2KSHELL): $(ALL_LIB) cp2k_shell.o
	+$(LD) $(LDFLAGS) -o $(CP2KSHELL) cp2k_shell.o $(ALL_LIB) $(LIBS)

$(FES): $(FESOBJS) graph.o 
	+$(LD) $(LDFLAGS) -o $(FES) graph.o $(FESOBJS)

# Allow the standalone DBCSR programs to compile with CUDA and/or HWTOPO

DBCSR_STANDALONE_DEPS += $(LIB3_ARCHIVE)

ifneq ($(HWTOPO),)
DBCSR_STANDALONE_DEPS += $(LIB4_ARCHIVE) $(LIBMA_ARCHIVE)
endif

DBCSR_STANDALONE_DEPS += $(LIB_CP2K_ARCHIVE)

ifneq ($(NVCC),)
DBCSR_STANDALONE_DEPS += $(LIB2CUDA_ARCHIVE)
DBCSR_STANDALONE_DEPS += $(LIBCUSMM_ARCHIVE)
endif


$(DBCSR_TEST_PROG): dbcsr_test_driver.o $(DBCSR_STANDALONE_DEPS)
	+$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

$(DBCSR_PERF_PROG): dbcsr_performance_driver.o $(DBCSR_STANDALONE_DEPS)
	+$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)


dbcsr_test_driver.o: $(LIB3_ARCHIVE)

dbcsr_performance_driver.o: $(LIB3_ARCHIVE)

#
# these are tools, but their integration in the build is ugly at least (e.g. see dependencies).
#
makedepf90:
	-test -d $(TOOLDIR)/makedepf90-build || ( mkdir makedepf90-build ; cp $(TOOLSRC)/makedepf90/* makedepf90-build )
	{ cd makedepf90-build ; ./configure --prefix=$(TOOLDIR) --bindir=$(TOOLDIR) ; $(MAKE) VERSION="2.8.8cp2k" ; $(MAKE) install ; }

#
# the rule how to generate the .o from the .F
# only if CPP is different from null we do a step over the C preprocessor (which is slower)
# in the other case the fortran compiler takes care of this directly
#
# With appropriate SMART_MODULES and MD5 settings an attempt is made
# to prevent unnecessary recompiling.  A generated .mod file is
# compared to a previously seen version.  If the module has not
# changed then its timestamp is reset to the past known version;
# otherwise the current version is remembered.  If no past known
# version is known then the current version is remembered.
%.o: %.F
ifneq ($(CPP),)
	$(CPP) $(CPPFLAGS) -D__SHORT_FILE__="\"$(notdir $<)\"" $< > $*.f90
	$(FC) -c $(FCFLAGS) $*.f90 $(FCLOGPIPE)
else
	$(FC) -c $(FCFLAGS) -D__SHORT_FILE__="\"$(notdir $<)\"" $< $(FCLOGPIPE)
endif
ifneq ($(SMART_MODULES),no)
ifneq ($(MD5),)
	@$(MD5) $*.mod &> /dev/null && ( \
	$(MD5) --status -c $*.mod.md5 2> /dev/null && touch -r $*.mod.md5 $*.mod \
	|| $(MD5) $*.mod > $*.mod.md5 ) || true
endif
endif

%.pretty: %.F 
	$(TOOLSRC)/prettify.py $<
	@touch $(PRETTYOBJDIR)/$*.pretty

%.o: %.cu
	$(NVCC) -c $(NVFLAGS) $<

%.o: %.c
	echo $<
	$(CC) -c $(CFLAGS) $<

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

#
# clean prettified files
#
prettyclean:
	-rm -rf $(PRETTYOBJDIR)

#
# delete the intermediate files, but not the libraries and executables. 
# Most useful to save space on the disk
#
clean:
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f *.mod ; }
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f *.o ; } 
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f *.f90 ; }
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f *.dep ; }

#
# delete the intermediate files, the programs and libraries and anything that might be in the objdir directory
# Use this if you want to rebuild the executable (for a given compiler and or VERSION)
#
realclean: clean
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f * ; }
	-rm -f $(PROG)
	-rm -f $(LIB1_ARCHIVE)
	-rm -f $(LIB2_ARCHIVE)
	-rm -f $(LIB3_ARCHIVE)
	-rm -f $(LIB4_ARCHIVE)
	-rm -f $(LIB_CP2K_ARCHIVE)
	-rm -f $(LIBCUSMM_FILES)

#
# remove plumed objectfiles
#
plumedclean:
	{ cd $(OBJDIR)/$(VERSION) ; rm $(OBJ_PLUMED) ; }
#
# Remove all files from previous builds
#
distclean:
	rm -rf $(DOXYGENDIR) $(MAINEXEDIR) $(MAINLIBDIR) $(MAINOBJDIR) $(LIBCUSMM_FILES)
#
# Remove ALL files not in SVN
# NOTE: This removes really everything not committed to SVN
#
mrproper:
	cd $(CP2KHOME); ./tools/clean_cwd.sh
#
# so far CP2K does not install, but give a hint to the user
#
install:
	@echo ""
	@echo "The CP2K executable is $(PROG)"
	@echo ""
#
doxygen/clean:
	-rm -rf $(DOXYGENDIR)
#
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
	{ cd $(SRCDIR) ; cp *.f90 *.h $(DOXYGENDIR) ; for f in $(OBJECTS:.o=) ; do cp $${f}.F $(DOXYGENDIR)/$${f}.f90; done ; }
	{ cd $(SRCDIR)/lib ; for f in $(LIB1_OBJECTS:.o=) ; do cp $${f}.F $(DOXYGENDIR)/$${f}.f90; done ; }
	{ cd $(SRCDIR)/fft_lib ; for f in $(LIB2_OBJECTS:.o=) ; do cp $${f}.F $(DOXYGENDIR)/$${f}.f90; done ; }
	{ cd $(SRCDIR)/dbcsr_lib ; for f in $(LIB3_OBJECTS:.o=) ; do cp $${f}.F $(DOXYGENDIR)/$${f}.f90; done ; }
	{ cd $(SRCDIR)/arch ; for f in $(LIB4_OBJECTS:.o=) ; do cp $${f}.F $(DOXYGENDIR)/$${f}.f90; done ; }
	{ cd $(DOCDIR) ; sed "s/#revision#/`$(TOOLSRC)/get_revision_number $(CP2KHOME)`/" Doxyfile.template >Doxyfile ; }
	{ cd $(DOXYGENDIR) ; doxygen $(DOCDIR)/Doxyfile ; }
# 
# Add additional dependency of cp2k_info.F to SVN-entry or git-HEAD.
# Ensuring that cp2k prints the correct source code revision number in its banner.
#
ifneq ($(strip $(SVN_ENTRY)),)
cp2k_info.o: $(SVN_ENTRY)
endif

ifneq ($(strip $(GIT_HEAD)),)
cp2k_info.o: $(GIT_HEAD)
endif

#
# Load the automatically generated rules of makedepf90 
#
ifneq ($(strip $(wildcard *.dep)),)
include $(wildcard *.dep)
endif
