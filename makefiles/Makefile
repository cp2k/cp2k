.SUFFIXES: .F .d .o

CP2KHOME  = $(HOME)/cp2k
ARCH      = $(shell $(CP2KHOME)/tools/get_arch_code)
VERSION   = sopt
SMAKE     = gmake -r
PMAKE     = $(SMAKE) -j 1
LIB1      = cp2klib

### Dependent variables ###

ARCHDIR      = $(CP2KHOME)/arch
MAINEXEDIR   = $(CP2KHOME)/exe
MAINLIBDIR   = $(CP2KHOME)/lib
MAINOBJDIR   = $(CP2KHOME)/obj
SRCDIR       = $(CP2KHOME)/src
TOOLDIR      = $(CP2KHOME)/tools
MAKEFILE     = $(CP2KHOME)/makefiles/Makefile
FORPAR       = $(CP2KHOME)/tools/forpar.x -chkint
SFMAKEDEPEND = $(CP2KHOME)/tools/sfmakedepend -m int -s -f
MACHINEDEFS  = $(ARCHDIR)/$(ARCH).$(VERSION)
PROG         = $(EXEDIR)/cp2k.$(VERSION)
EXEDIR       = $(MAINEXEDIR)/$(ARCH)
LIBDIR       = $(MAINLIBDIR)/$(ARCH)
OBJDIR       = $(MAINOBJDIR)/$(ARCH)
VPATH        = $(SRCDIR):$(SRCDIR)/lib
OBJECTDEFS   = $(SRCDIR)/OBJECTDEFS
LIB1_ARCHIVE = $(LIBDIR)/$(VERSION)/lib$(LIB1).a

### Definition of the multiple targets ###

VERSION_TARGETS = sopt sdbg popt pdbg
CLEAN_TARGETS   = sopt/clean sdbg/clean popt/clean pdbg/clean\
                  sopt/realclean sdbg/realclean popt/realclean pdbg/realclean

### Master rules ###

$(VERSION_TARGETS):
	@cd $(EXEDIR) || mkdir -p $(EXEDIR)
	@cd $(OBJDIR)/$@ || mkdir -p $(OBJDIR)/$@
	@cd $(LIBDIR)/$@ || mkdir -p $(LIBDIR)/$@
	@$(SMAKE) -C $(TOOLDIR)
	@$(SMAKE) -C $(OBJDIR)/$@ -f $(MAKEFILE) VERSION=$@ dependencies
	@$(PMAKE) -C $(OBJDIR)/$@ -f $(MAKEFILE) VERSION=$@ all

$(CLEAN_TARGETS):
	$(SMAKE) -C $(OBJDIR)/$(@D) -f $(MAKEFILE) VERSION=$(@D) $(@F)

include $(OBJECTDEFS)

include $(MACHINEDEFS)

OBJECTS = $(OBJECTS_GENERIC) $(OBJECTS_ARCHITECTURE)

LIB1_OBJECTS = ctrig.o fftpre.o fftrot.o fftstp.o mltfftsg.o fast.o

DEPENDENCIES = $(OBJECTS:.o=.d)

LIBRARIES = -L$(LIBDIR)/$(VERSION) -l$(LIB1) $(LIBS)

### Slave rules ###

dependencies: $(DEPENDENCIES)

all: $(PROG)

$(PROG): $(OBJECTS) $(LIB1_ARCHIVE)
	$(LD) $(LDFLAGS) -o $(PROG) $(OBJECTS) $(LIBRARIES)

$(LIB1_ARCHIVE): $(LIB1_OBJECTS)
	$(AR) $(LIB1_ARCHIVE) $(LIB1_OBJECTS)
ifneq ($(RANLIB),)
	$(RANLIB) $(LIB1_ARCHIVE)
endif

%.o: %.F
	$(FC) -c $(FCFLAGS) $<

%.d: %.F
	$(CPP) $(CPPFLAGS) $< > $*.f
	$(FORPAR) $*.f
	$(PERL) $(SFMAKEDEPEND) $*.d $*.f
	@rm -f $*.d.old

parallel_include.o: parallel_include.F
	$(FC_fixed) -c $(FCFLAGS) $<

parallel_include.d: parallel_include.F
	$(CPP) $(CPPFLAGS) $< > $*.f
	-$(FORPAR) -fix $*.f
	$(PERL) $(SFMAKEDEPEND) $*.d $*.f
	@rm -f $*.d.old

clean:
	rm -f *.o *.mod

realclean: clean
	rm -f $(PROG) *.f *.d *.int *~ *.lst

distclean:
	rm -rf $(MAINEXEDIR) $(MAINLIBDIR) $(MAINOBJDIR)
	$(SMAKE) -C $(TOOLDIR) realclean

### Load the automatically generated rules of sfmakedepend ###

include $(wildcard *.d)
