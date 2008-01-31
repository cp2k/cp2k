.SUFFIXES: .F .dep .o .f90 
SHELL = /bin/sh
#
# the home dir is taken from the current directory: you must cd to the makefiles
# directory before executing make.
# adding a symbolic link from the src directory to the Makefile in makefiles/ 
# also works fine
#
CP2KHOME  := $(PWD)/..
ARCH      := $(shell $(CP2KHOME)/tools/get_arch_code)
SMAKE     = $(MAKE) -r
PMAKE     = $(SMAKE) -j 4
LIB1      = cp2k_base_lib
LIB2      = cp2k_fft_lib
LIBCP2K   = cp2k_lib
export VERSION=sopt
export WHAT   =_progr

### Dependent variables ###

ARCHDIR      = $(CP2KHOME)/arch
MAINEXEDIR   = $(CP2KHOME)/exe
MAINLIBDIR   = $(CP2KHOME)/lib
MAINOBJDIR   = $(CP2KHOME)/obj
PRETTYOBJDIR = $(CP2KHOME)/obj/prettified
TOOLDIR      = $(OBJDIR)/tools
TOOLSRC      = $(CP2KHOME)/tools
SRCDIR       = $(CP2KHOME)/src
SRCDIRS      = $(SRCDIR):$(SRCDIR)/lib:$(TOOLSRC)/Fist/topology_converter::$(TOOLSRC)/metadyn
MAKEFILE     = $(CP2KHOME)/makefiles/Makefile
MACHINEDEFS  = $(ARCHDIR)/$(ARCH).$(VERSION)
PROG         = $(EXEDIR)/cp2k.$(VERSION)
CP2KSHELL    = $(EXEDIR)/cp2k_shell.$(VERSION)
LEAP2FIST    = $(EXEDIR)/leap2fist.$(VERSION)
FES          = $(EXEDIR)/fes.$(VERSION)
EXEDIR       = $(MAINEXEDIR)/$(ARCH)
LIBDIR       = $(MAINLIBDIR)/$(ARCH)
OBJDIR       = $(MAINOBJDIR)/$(ARCH)
OBJECTDEFS   = $(SRCDIR)/OBJECTDEFS
LIBOBJECTDEFS= $(SRCDIR)/lib/OBJECTDEFS
LIB1_ARCHIVE = $(LIBDIR)/$(VERSION)/lib$(LIB1).a
LIB2_ARCHIVE = $(LIBDIR)/$(VERSION)/lib$(LIB2).a
LIB_CP2K_ARCHIVE  = $(LIBDIR)/$(VERSION)/lib$(LIBCP2K).a
SRCENTRY     = $(wildcard $(SRCDIR)/CVS/Entr*)
LEAPOBJECTDEFS = $(TOOLSRC)/Fist/topology_converter/OBJECTDEFS
FESOBJECTDEFS  = $(TOOLSRC)/metadyn/OBJECTDEFS

#
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
         $(ALL_TARGETS) $(LIB_TARGETS) $(VERSION_TARGETS) leap2fist fes dirs deptool

### Master rules ###
build: dirs deptool
	$(SMAKE) -C $(SRCDIR) -f $(MAKEFILE) all.dep
	$(SMAKE) -C $(SRCDIR)/lib -f $(MAKEFILE) lib.dep
	$(PMAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) $(WHAT)

leap2fist: dirs deptool
	$(SMAKE) -C $(TOOLSRC)/Fist/topology_converter/ -f $(MAKEFILE) leap2fist.dep
	$(PMAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) VERSION=$(VERSION) $(LEAP2FIST) 

fes: dirs deptool
	$(SMAKE) -C $(TOOLSRC)/metadyn/  -f $(MAKEFILE) fes.dep
	$(PMAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) VERSION=$(VERSION) $(FES) 

dirs:
	@mkdir -p $(EXEDIR)
	@mkdir -p $(LIBDIR)/$(VERSION)
	@mkdir -p $(TOOLDIR)
	@mkdir -p $(OBJDIR)/$(VERSION)
	@mkdir -p $(PRETTYOBJDIR)

pretty: dirs
	$(PMAKE) -C $(SRCDIR) -f $(MAKEFILE) _pretty

deptool: dirs
	$(SMAKE) -C $(TOOLDIR) -f $(MAKEFILE) makedepf90

libs:
	$(SMAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_lib build
all:
	$(SMAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_all build
progr:
	$(SMAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_progr build
cpshell:
	$(SMAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_cpshell build

$(VERSION_TARGETS):
	$(SMAKE) -f $(MAKEFILE) VERSION=$@ build

$(CLEAN_TARGETS):
	$(SMAKE) -f $(MAKEFILE) VERSION=$(@D) $(@F)

$(ALL_TARGETS) $(LIB_TARGETS):
	$(SMAKE) -f $(MAKEFILE) VERSION=$(@D) WHAT=_$(@F) build
#
# here we try to include all the objects that might need considering
#

include $(OBJECTDEFS)
include $(LIBOBJECTDEFS)
include $(LEAPOBJECTDEFS)
include $(FESOBJECTDEFS)
include $(MACHINEDEFS)

OBJECTS = $(OBJECTS_GENERIC) $(OBJECTS_ARCHITECTURE)
LEAPOBJS = $(LEAPOBJ) $(OBJECTS_ARCHITECTURE)
FESOBJS = $(FESOBJ) $(OBJECTS_ARCHITECTURE)

#
# some practical variables for the build.
# note that the last cvs entry is not necessarily the most recent
#
CPPSHELL := -D__COMPILE_ARCH="\"$(ARCH)\""\
            -D__COMPILE_DATE="\"$(shell date)\""\
            -D__COMPILE_HOST="\"$(shell hostname)\""\
            -D__COMPILE_LASTCVS="\"$(shell tail -n1 $(SRCDIR)/CVS/Entries)\""

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
	$(TOOLDIR)/makedepf90 -m "%m.mod" -free $^ > $(OBJDIR)/$(VERSION)/all.dep
endif

lib.dep: $(LIB1_OBJECTS:.o=.F) $(LIB2_OBJECTS:.o=.F)
	$(TOOLDIR)/makedepf90 -m "%m.mod" -free $^ > $(OBJDIR)/$(VERSION)/lib.dep

# this is kind of a hack, fes and leap2fist depend on files in other directories, hence the -I and the -b.
fes.dep: $(FESOBJS:.o=.F) 
	$(TOOLDIR)/makedepf90 -I$(SRCDIR) -b. -m "%m.mod" -free $^ > $(OBJDIR)/$(VERSION)/fes.dep

leap2fist.dep: $(LEAPOBJS:.o=.F) 
	$(TOOLDIR)/makedepf90 -I$(SRCDIR) -b. -m "%m.mod" -free $^ > $(OBJDIR)/$(VERSION)/leap2fist.dep

#
# and the rules doing the actual work
#
_all: $(LIB1_ARCHIVE) $(LIB2_ARCHIVE) $(LIB_CP2K_ARCHIVE) $(PROG) leap2fist fes $(CP2KSHELL)
_progr: $(PROG)
_lib: $(LIB1_ARCHIVE) $(LIB2_ARCHIVE) $(LIB_CP2K_ARCHIVE)
_cpshell: $(CP2KSHELL)

_pretty:$(OBJECTS:.o=.pretty)

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

$(LIB_CP2K_ARCHIVE): $(OBJECTS)
	$(AR) $(LIB_CP2K_ARCHIVE) $?
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB_CP2K_ARCHIVE)
endif

$(PROG): $(LIB1_ARCHIVE) $(LIB2_ARCHIVE) $(LIB_CP2K_ARCHIVE) cp2k.o
	$(LD) $(LDFLAGS) -o $(PROG) cp2k.o $(LIB_CP2K_ARCHIVE) $(LIB2_ARCHIVE) $(LIB1_ARCHIVE) $(LIBS)

$(CP2KSHELL): $(LIB1_ARCHIVE) $(LIB2_ARCHIVE) $(LIB_CP2K_ARCHIVE) cp2k_shell.o
	$(LD) $(LDFLAGS) -o $(CP2KSHELL) cp2k_shell.o $(LIB_CP2K_ARCHIVE) $(LIB2_ARCHIVE) $(LIB1_ARCHIVE) $(LIBS)

#
# these are tools, but their integration in the build is ugly at least (e.g. see dependencies).
#
makedepf90:
	-test -d $(TOOLDIR)/makedepf90-build || cp -r $(TOOLSRC)/makedepf90 makedepf90-build
	{ cd makedepf90-build ; ./configure --prefix=$(TOOLDIR) --bindir=$(TOOLDIR) ; $(MAKE) VERSION="2.8.8" ; $(MAKE) install ; }

$(LEAP2FIST): $(LEAPOBJS)
	$(LD) $(LDFLAGS) -o $(LEAP2FIST) $(LEAPOBJS)

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

%.pretty: %.F $(TOOLSRC)/prettify.py $(TOOLSRC)/normalizeFortranFile.py
	$(TOOLSRC)/prettify.py $<
	touch $(PRETTYOBJDIR)/$*.pretty

#
# here we cheat... this tells make that .mod can be generated from .o (this holds in CP2K)
# by doing nothing
# it avoids recompilation if .o is more recent than .F, but .mod is older than .F
# (because it didn't change, as e.g. g95 can do)
#
%.mod: %.o
	@true

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
#
# really remove all build files not in CVS
#
distclean:
	rm -rf $(MAINEXEDIR) $(MAINLIBDIR) $(MAINOBJDIR)

#
# so far CP2K does not install, but give a hint to the user
#
install:
	@echo ""
	@echo "The CP2K executable is $(PROG)"
	@echo ""

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
