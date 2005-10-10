.SUFFIXES: .F .dep .o .f90 .for
SHELL = /bin/sh
# the home dir is taken from the current directory: you must cd to the makefiles
# directory before executing make
CP2KHOME  := $(PWD)/..
ARCH      := $(shell $(CP2KHOME)/tools/get_arch_code)
SMAKE     = $(MAKE) -r
PMAKE     = $(SMAKE) -j 1
LIB1      = cp2k_base_lib
LIBCP2K   = cp2k_lib
export VERSION=sopt
export WHAT   =_progr

### Dependent variables ###

ARCHDIR      = $(CP2KHOME)/arch
MAINEXEDIR   = $(CP2KHOME)/exe
MAINLIBDIR   = $(CP2KHOME)/lib
MAINOBJDIR   = $(CP2KHOME)/obj
SRCDIR       = $(CP2KHOME)/src
SRCDIRS      = $(SRCDIR):$(SRCDIR)/lib
TOOLDIR      = $(OBJDIR)/tools
TOOLSRC      = $(CP2KHOME)/tools
FORPAR       = $(TOOLDIR)/forpar.x -chkint
MAKEFILE     = $(CP2KHOME)/makefiles/Makefile
SFMAKEDEPEND = $(CP2KHOME)/tools/sfmakedepend -m int -s -f
MACHINEDEFS  = $(ARCHDIR)/$(ARCH).$(VERSION)
PROG         = $(EXEDIR)/cp2k.$(VERSION)
EXEDIR       = $(MAINEXEDIR)/$(ARCH)
LIBDIR       = $(MAINLIBDIR)/$(ARCH)
OBJDIR       = $(MAINOBJDIR)/$(ARCH)
OBJECTDEFS   = $(SRCDIR)/OBJECTDEFS
LIB1_ARCHIVE = $(LIBDIR)/$(VERSION)/lib$(LIB1).a
LIB_CP2K_ARCHIVE  = $(LIBDIR)/$(VERSION)/lib$(LIBCP2K).a
SRCENTRY     = $(wildcard $(SRCDIR)/CVS/Entr*)

### Definition of the multiple targets ###

VERSION_TARGETS = sopt sdbg popt pdbg ssmp psmp
LIB_TARGETS = sopt/lib sdbg/lib popt/lib pdbg/lib ssmp/lib psmp/lib
ALL_TARGETS = sopt/all sdbg/all popt/all pdbg/all ssmp/all psmp/all
CLEAN_TARGETS   = sopt/clean sdbg/clean popt/clean pdbg/clean\
                  sopt/realclean sdbg/realclean popt/realclean pdbg/realclean
.PHONY : _all libs _progr _lib build clean realclean distclean $(CLEAN_TARGETS)\
         $(ALL_TARGETS) $(LIB_TARGETS) $(VERSION_TARGETS) dependencies
### Master rules ###

build:
	@mkdir -p $(EXEDIR)
	@mkdir -p $(LIBDIR)/$(VERSION)
	@mkdir -p $(TOOLDIR)
	@mkdir -p $(OBJDIR)/$(VERSION)
	@echo
	@echo "*** machine file = $(MACHINEDEFS)"
	@echo
	$(SMAKE) -C $(TOOLDIR) -f $(TOOLSRC)/Makefile SRCDIR=$(TOOLSRC) CC=$(CC)
	$(SMAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) dependencies
	$(PMAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) $(WHAT)

libs:
	$(SMAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_lib build
all:
	$(SMAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_all build
progr:
	$(SMAKE) -f $(MAKEFILE) VERSION=$(VERSION) WHAT=_progr build

$(VERSION_TARGETS):
	$(SMAKE) -f $(MAKEFILE) VERSION=$@ build

$(CLEAN_TARGETS):
	$(SMAKE) -f $(MAKEFILE) VERSION=$(@D) $(@F)

$(ALL_TARGETS) $(LIB_TARGETS):
	$(SMAKE) -f $(MAKEFILE) VERSION=$(@D) WHAT=_$(@F) build

include $(OBJECTDEFS)

include $(MACHINEDEFS)

CPPSHELL := -D__COMPILE_ARCH="\"$(ARCH)\""\
            -D__COMPILE_DATE="\"$(shell date)\""\
            -D__COMPILE_HOST="\"$(shell hostname)\""\
            -D__COMPILE_LASTCVS="\"$(shell tail -1 $(SRCDIR)/CVS/Entries)\""

CPPFLAGS += $(CPPSHELL)\
	    -I$(SRCDIR)

OBJECTS = $(OBJECTS_GENERIC) $(OBJECTS_ARCHITECTURE)

LIB1_OBJECTS = ctrig.o fftpre.o fftrot.o fftstp.o fast.o mltfftsg_tools.o\
               mltfftsg.o grid_fast.o

DEPENDENCIES = $(OBJECTS:.o=.dep) cp2k.dep

LIBRARIES = -L$(LIBDIR)/$(VERSION) -l$(LIB1) $(LIBS)

.SECONDARY: $(OBJECTS:.o=.f90) $(LIB1_OBJECTS:.o=.f90)

### Slave rules ###

dependencies: $(DEPENDENCIES)
	@echo
	@echo "** DFLAGS   = $(DFLAGS)"
	@echo "** CPPFLAGS = $(CPPFLAGS)"
	@echo "** FCFLAGS  = $(FCFLAGS)"
	@echo "** LDFLAGS  = $(LDFLAGS)"
	@echo "** LIBS     = $(LIBS)"
	@echo

_all: $(LIB1_ARCHIVE) $(LIB_CP2K_ARCHIVE) $(PROG)
_progr: $(PROG)
_lib: $(LIB1_ARCHIVE) $(LIB_CP2K_ARCHIVE)

$(LIB1_ARCHIVE): $(LIB1_OBJECTS)
	$(AR) $(LIB1_ARCHIVE) $(LIB1_OBJECTS)
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB1_ARCHIVE)
endif

$(PROG): $(LIB1_ARCHIVE) $(OBJECTS) cp2k.o
	$(LD) $(LDFLAGS) -o $(PROG) $(OBJECTS) cp2k.o $(LIBRARIES)

$(LIB_CP2K_ARCHIVE): $(OBJECTS) f77_int_low.o
	$(AR) $(LIB_CP2K_ARCHIVE) $(OBJECTS) f77_int_low.o
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB_CP2K_ARCHIVE)
endif

vpath %.F $(SRCDIRS)

%.f90: %.F
ifneq ($(POST_PROC),)
	$(CPP) $(CPPFLAGS) $< > $*.tmp2
	$(POST_PROC)  $*.tmp2 > $*.f90
	rm $*.tmp2
else
	$(CPP) $(CPPFLAGS) $< > $*.f90
endif

%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

%.dep: %.F
ifneq ($(POST_PROC),)
	$(CPP) $(CPPFLAGS) -P $< > $*.tmp
	$(POST_PROC)  $*.tmp > $*.for
	rm $*.tmp
else
	$(CPP) $(CPPFLAGS) -P $< > $*.for
endif
	$(FORPAR) $*.for
	$(PERL) $(SFMAKEDEPEND) $*.dep $*.for
	@rm -f $*.dep.old $*.for

# Special hack for parallel_include since $(FORPAR) cannot digest
# the INCLUDE statement which is active for parallel compilations

parallel_include.dep:
	touch -a parallel_include.dep parallel_include.int
parallel_include.int: parallel_include.o

clean:
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f *.mod ; }
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f *.o ; } 

realclean: clean
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f *.f90 ; }
	-test ! -d $(OBJDIR)/$(VERSION) || { cd $(OBJDIR)/$(VERSION) ; rm -f * ; }
	-rm -f $(PROG)
	-rm -f $(LIB1_ARCHIVE)

distclean:
	rm -rf $(MAINEXEDIR) $(MAINLIBDIR) $(MAINOBJDIR)
	$(SMAKE) -C $(TOOLSRC) realclean

# lib dependencies
mltfftsg.o: mltfftsg_tools.o
f77_int_low.o: f77_interface.o

# to recpp environment.F whenever the CVS/Entries file is updated (and exists)
ifneq ($(strip $(SRCENTRY)),)
environment.f90:$(SRCENTRY)
endif

### Load the automatically generated rules of sfmakedepend ###

ifneq ($(strip $(wildcard *.dep)),)
include $(wildcard *.dep)
endif
