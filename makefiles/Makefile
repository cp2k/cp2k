.SUFFIXES: .F .dep .o .f90 .for

# the home dir is taken from the current directory: you must cd to the makefiles
# directory before executing make
CP2KHOME  := $(PWD)/..
ARCH      = $(shell $(CP2KHOME)/tools/get_arch_code)
SMAKE     = $(MAKE) -r
PMAKE     = $(SMAKE) -j 4
LIB1      = cp2klib
export VERSION=sopt

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

### Definition of the multiple targets ###

VERSION_TARGETS = sopt sdbg popt pdbg
CLEAN_TARGETS   = sopt/clean sdbg/clean popt/clean pdbg/clean\
                  sopt/realclean sdbg/realclean popt/realclean pdbg/realclean
.PHONY : build clean realclean distclean $(CLEAN_TARGETS)
### Master rules ###

build:
	@mkdir -p $(EXEDIR)
	@mkdir -p $(LIBDIR)/$(VERSION)
	@mkdir -p $(TOOLDIR)
	@mkdir -p $(OBJDIR)/$(VERSION)
	@echo
	@echo "*** machine file = $(MACHINEDEFS)"
	@echo
	$(SMAKE) -C $(TOOLDIR) -f $(TOOLSRC)/Makefile SRCDIR=$(TOOLSRC)
	$(SMAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) dependencies
	$(PMAKE) -C $(OBJDIR)/$(VERSION) -f $(MAKEFILE) all

$(VERSION_TARGETS):
	$(SMAKE) -f $(MAKEFILE) VERSION=$@ build

$(CLEAN_TARGETS):
	$(SMAKE) -f $(MAKEFILE) VERSION=$(@D) $(@F)

include $(OBJECTDEFS)

include $(MACHINEDEFS)

CPPFLAGS += -D__COMPILE_ARCH="\"$(ARCH)\""\
            -D__COMPILE_DATE="\"$(shell date)\""\
            -D__COMPILE_HOST="\"$(shell hostname)\""\
	    -I$(SRCDIR)

OBJECTS = $(OBJECTS_GENERIC) $(OBJECTS_ARCHITECTURE)

LIB1_OBJECTS = ctrig.o fftpre.o fftrot.o fftstp.o fast.o mltfftsg_tools.o\
               mltfftsg.o grid_fast.o

DEPENDENCIES = $(OBJECTS:.o=.dep)

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

all: $(PROG)

$(LIB1_ARCHIVE): $(LIB1_OBJECTS)
	$(AR) $(LIB1_ARCHIVE) $(LIB1_OBJECTS)
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB1_ARCHIVE)
endif

$(PROG): $(LIB1_ARCHIVE) $(OBJECTS)
	$(LD) $(LDFLAGS) -o $(PROG) $(OBJECTS) $(LIBRARIES)

vpath %.F $(SRCDIRS)

%.f90: %.F
	$(CPP) $(CPPFLAGS) $< > $*.f90

%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

%.dep: %.F
	$(CPP) $(CPPFLAGS) -P $< > $*.for
	$(FORPAR) $*.for
	$(PERL) $(SFMAKEDEPEND) $*.dep $*.for
	@rm -f $*.dep.old $*.for

# Special hack for parallel_include since $(FORPAR) cannot digest
# the INCLUDE statement which is active for parallel compilations

parallel_include.dep:
	touch -a parallel_include.dep parallel_include.int
parallel_include.int: parallel_include.o

clean:
	cd $(OBJDIR)/$(VERSION) ; rm -f *.o *.mod

realclean: clean
	cd $(OBJDIR)/$(VERSION) ; rm -f $(PROG) *

distclean:
	rm -rf $(MAINEXEDIR) $(MAINLIBDIR) $(MAINOBJDIR)
	$(SMAKE) -C $(TOOLSRC) realclean

# lib dependencies
mltfftsg.o: mltfftsg_tools.o

### Load the automatically generated rules of sfmakedepend ###

include $(wildcard *.dep)
