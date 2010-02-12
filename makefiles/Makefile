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
CP2KHOME  := $(PWD)/..
ARCH      := $(shell $(CP2KHOME)/tools/get_arch_code)
LIB1      = cp2k_base_lib
LIB2      = cp2k_fft_lib
LIB3      = dbcsr
LIBCP2K   = cp2k_lib
export VERSION=sopt
export WHAT   =_progr

### Dependent variables ###

ARCHDIR      = $(CP2KHOME)/arch
DOCDIR       = $(CP2KHOME)/doc
MAINEXEDIR   = $(CP2KHOME)/exe
MAINLIBDIR   = $(CP2KHOME)/lib
MAINOBJDIR   = $(CP2KHOME)/obj
PRETTYOBJDIR = $(CP2KHOME)/obj/prettified
TOOLDIR      = $(OBJDIR)/tools
TOOLSRC      = $(CP2KHOME)/tools
SRCDIR       = $(CP2KHOME)/src
SRCDIRS      = $(SRCDIR):$(SRCDIR)/lib:$(TOOLSRC)/metadyn:$(CP2KHOME)/cuda_tools
MAKEFILE     = $(CP2KHOME)/makefiles/Makefile
MACHINEDEFS  = $(ARCHDIR)/$(ARCH).$(VERSION)
PROG         = $(EXEDIR)/cp2k.$(VERSION)
CP2KSHELL    = $(EXEDIR)/cp2k_shell.$(VERSION)
FES          = $(EXEDIR)/fes.$(VERSION)
EXEDIR       = $(MAINEXEDIR)/$(ARCH)
LIBDIR       = $(MAINLIBDIR)/$(ARCH)
OBJDIR       = $(MAINOBJDIR)/$(ARCH)
OBJECTDEFS   = $(SRCDIR)/OBJECTDEFS
LIBOBJECTDEFS= $(SRCDIR)/lib/OBJECTDEFS
LIBNVOBJECTDEFS= $(CP2KHOME)/cuda_tools/NVOBJECTDEFS
LIB1_ARCHIVE = $(LIBDIR)/$(VERSION)/lib$(LIB1).a
LIB2_ARCHIVE = $(LIBDIR)/$(VERSION)/lib$(LIB2).a
LIB3_ARCHIVE = $(LIBDIR)/$(VERSION)/lib$(LIB3).a
LIB_CP2K_ARCHIVE  = $(LIBDIR)/$(VERSION)/lib$(LIBCP2K).a
SRCENTRY     = $(wildcard $(SRCDIR)/CVS/Entr*)
FESOBJECTDEFS  = $(TOOLSRC)/metadyn/OBJECTDEFS
DOXYGENDIR   = $(DOCDIR)/doxygen

### Definition of the multiple targets ###
#
VERSION_TARGETS = sopt sdbg popt pdbg ssmp psmp
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
         $(ALL_TARGETS) $(LIB_TARGETS) $(VERSION_TARGETS) fes dirs deptool

### Master rules ###
build: dirs deptool
	$(MAKE) -C $(SRCDIR) -f $(MAKEFILE) all.dep
	$(MAKE) -C $(SRCDIR)/lib -f $(MAKEFILE) lib.dep
	$(MAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) $(WHAT)

fes: dirs deptool
	$(MAKE) -C $(TOOLSRC)/metadyn/  -f $(MAKEFILE) fes.dep
	$(MAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) VERSION=$(VERSION) $(FES) 

dirs:
	@mkdir -p $(EXEDIR)
	@mkdir -p $(LIBDIR)/$(VERSION)
	@mkdir -p $(TOOLDIR)
	@mkdir -p $(OBJDIR)/$(VERSION)
	@mkdir -p $(PRETTYOBJDIR)

realpretty: prettyclean pretty

pretty: dirs
	$(MAKE) -C $(SRCDIR) -f $(MAKEFILE) _pretty
	$(MAKE) -C $(SRCDIR)/lib -f $(MAKEFILE) _pretty_lib3

deptool: dirs
	$(MAKE) -C $(TOOLDIR) -f $(MAKEFILE) makedepf90

libs:
	$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_lib build
all:
	$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_all build
progr:
	$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_progr build
cpshell:
	$(MAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_cpshell build

#
# on a parallel build only one task should go a build the deptool
#
$(VERSION_TARGETS): dirs deptool
	$(MAKE) -f $(MAKEFILE) VERSION=$@ build

$(CLEAN_TARGETS):
	$(MAKE) -f $(MAKEFILE) VERSION=$(@D) $(@F)

$(ALL_TARGETS) $(LIB_TARGETS):
	$(MAKE) -f $(MAKEFILE) VERSION=$(@D) WHAT=_$(@F) build
#
# here we try to include all the objects that might need considering
#

include $(OBJECTDEFS)
include $(LIBOBJECTDEFS)
include $(FESOBJECTDEFS)
include $(MACHINEDEFS)
include $(LIBNVOBJECTDEFS)

OBJECTS = $(OBJECTS_GENERIC) $(OBJECTS_ARCHITECTURE)
FESOBJS = $(FESOBJ) $(OBJECTS_ARCHITECTURE)

# Cuda specific instructions
ifneq ($(NVCC),)
  LIBNV_ARCHIVE= $(LIBDIR)/$(VERSION)/libcp2k_cuda.a
endif

# Define the whole bunch of libraries needed
ALL_LIB =  $(LIB_CP2K_ARCHIVE) $(LIB3_ARCHIVE) $(LIB2_ARCHIVE) $(LIB1_ARCHIVE) $(LIBNV_ARCHIVE)

#
# some practical variables for the build.
# note that the last cvs entry is not necessarily the most recent
#
ifeq ($(CPPSHELL),)
CPPSHELL := -D__COMPILE_ARCH="\"$(ARCH)\""\
            -D__COMPILE_DATE="\"$(shell date)\""\
            -D__COMPILE_HOST="\"$(shell hostname)\""\
            -D__COMPILE_LASTCVS="\"$(shell tail -n1 $(SRCDIR)/CVS/Entries)\""
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
# if HACKDEP=="yes" we hack the dependencies so that only modified source files are recompiled
# this is obviously not correct (i.e. can lead to a miscompiled CP2K), 
# but if you know what you are doing (i.e. you are sure that only information private to this module will be changed) might be useful nevertheless.
#
all.dep: $(OBJECTS:.o=.F) cp_common_uses.h cp2k.F cp2k_shell.F
ifeq ($(HACKDEP),yes)
	echo $^ | awk '{for (i=1;i<=NF;i++) {obj=$$(i) ; sub(".F",".o",obj); if ($$(i)!=obj) print obj,":",$$(i) }}' > $(OBJDIR)/$(VERSION)/all.dep
else
	$(TOOLDIR)/makedepf90 $(MODSTRING)  -free $^ > $(OBJDIR)/$(VERSION)/all.dep || ( rm -f $(OBJDIR)/$(VERSION)/all.dep  ; exit 1 )
endif

lib.dep: $(LIB1_OBJECTS:.o=.F) $(LIB2_OBJECTS:.o=.F)  $(LIB3_OBJECTS:.o=.F)
	$(TOOLDIR)/makedepf90 $(MODSTRING)  -free $^ > $(OBJDIR)/$(VERSION)/lib.dep || ( rm -f $(OBJDIR)/$(VERSION)/lib.dep ; exit 1 )

# this is kind of a hack, fes depends on files in other directories, hence the -I and the -b.
fes.dep: $(FESOBJS:.o=.F) 
	$(TOOLDIR)/makedepf90 -I$(SRCDIR) -b. $(MODSTRING)  -free $^ > $(OBJDIR)/$(VERSION)/fes.dep || ( rm -f $(OBJDIR)/$(VERSION)/fes.dep  ; exit 1 )

#
# and the rules doing the actual work
#
_all: $(ALL_LIB) $(PROG) fes $(CP2KSHELL)
_progr: $(PROG)
_lib: $(ALL_LIB)
_cpshell: $(CP2KSHELL)

_pretty:$(OBJECTS:.o=.pretty)
_pretty_lib3:$(LIB3_OBJECTS:.o=.pretty)

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

#
# a hack to first compile the src/lib/ stuff before files from src/
# otherwise we have missing dependencies
#
$(OBJECTS): $(LIB3_ARCHIVE) $(LIB2_ARCHIVE) $(LIB1_ARCHIVE)

$(LIB_CP2K_ARCHIVE): $(OBJECTS)
	$(AR) $(LIB_CP2K_ARCHIVE) $?
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB_CP2K_ARCHIVE)
endif

ifneq ($(NVCC),)
  $(LIBNV_ARCHIVE): $(LIBNV_OBJECTS)
	$(AR) $(LIBNV_ARCHIVE) $?
endif

$(PROG): $(ALL_LIB) cp2k.o
	$(LD) $(LDFLAGS) -o $(PROG) cp2k.o $(ALL_LIB) $(LIBS)

$(CP2KSHELL): $(ALL_LIB) cp2k_shell.o
	$(LD) $(LDFLAGS) -o $(CP2KSHELL) cp2k_shell.o $(ALL_LIB) $(LIBS)

#
# these are tools, but their integration in the build is ugly at least (e.g. see dependencies).
#
makedepf90:
	-test -d $(TOOLDIR)/makedepf90-build || cp -r $(TOOLSRC)/makedepf90 makedepf90-build
	{ cd makedepf90-build ; ./configure --prefix=$(TOOLDIR) --bindir=$(TOOLDIR) ; $(MAKE) VERSION="2.8.8cp2k" ; $(MAKE) install ; }

$(FES): $(FESOBJS)
	$(LD) $(LDFLAGS) -o $(FES) $(FESOBJS)


#
# the rule how to generate the .o from the .F
# only if CPP is different from null we do a step over the C preprocessor (which is slower)
# in the other case the fortran compiler takes care of this directly
#
%.o: %.F
ifneq ($(CPP),)
	$(CPP) $(CPPFLAGS) $< > $*.f90
	$(FC) -c $(FCFLAGS) $*.f90
else
	$(FC) -c $(FCFLAGS) $<
endif

%.pretty: %.F 
	$(TOOLSRC)/prettify.py $<
	@touch $(PRETTYOBJDIR)/$*.pretty

%.o: %.cu
	$(NVCC) -c $(NVFLAGS) $<

#
# here we cheat... this tells make that .mod can be generated from .o (this holds in CP2K)
# by doing nothing
# it avoids recompilation if .o is more recent than .F, but .mod is older than .F
# (because it didn't change, as e.g. g95 can do)
#
# this is problematic if the module names are uppercase e.g. KINDS.mod (because this rule expands to kinds.mod)
#
%.mod: %.o
	@true

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
	-rm -f $(LIB_CP2K_ARCHIVE)
#
# really remove all build files not in CVS
#
distclean:
	rm -rf $(DOXYGENDIR) $(MAINEXEDIR) $(MAINLIBDIR) $(MAINOBJDIR)

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
	{ cd $(SRCDIR) ; cp *.f90 *.h $(DOXYGENDIR) ; for f in $(OBJECTS:.o=) ; do cp $${f}.F $(DOXYGENDIR)/$${f}.f90; done ; }
	{ cd $(SRCDIR)/lib ; for f in $(LIB1_OBJECTS:.o=) $(LIB2_OBJECTS:.o=) $(LIB3_OBJECTS:.o=) ; do cp $${f}.F $(DOXYGENDIR)/$${f}.f90; done ; }
	{ cd $(DOXYGENDIR) ; doxygen $(DOCDIR)/Doxyfile ; }
# 
# an additional dependency of the CVS entries in the source directory
# to recpp cp2k_info.F whenever the CVS/Entries file is updated (and exists)
#
ifneq ($(strip $(SRCENTRY)),)
cp2k_info.o:$(SRCENTRY)
endif

#
# Load the automatically generated rules of makedepf90 
#
ifneq ($(strip $(wildcard *.dep)),)
include $(wildcard *.dep)
endif
