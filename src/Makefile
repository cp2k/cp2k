.SUFFIXES: .o .F .d

CP2KHOME= ${HOME}/CP2K
FORPAR	= $(CP2KHOME)/tools/forpar.x -chkint

include OBJECTDEFS

include MACHINEDEFS

OBJECTS	= $(OBJECTS_GENERIC) $(OBJECTS_ARCHITECTURE)

#################################

PROG	= cp2k.x

all: $(PROG)

$(PROG): $(OBJECTS)
	$(LDR) $(LFLAGS) -o $(PROG) $(OBJECTS) $(LIBS)

%.o : %.F
	$(FC) -c $(FFLAGS) $*.F
	$(CPP) -P $(CDEFS) $*.F $*.f
	-$(FORPAR) $*.f
	rm -f $*.f

parallel_include.o : parallel_include.F
	$(FCfixed) -c $(FFLAGS) $*.F
	$(CPP) -P $(CDEFS) $*.F $*.f
	-$(FORPAR) -fix $*.f
	rm -f $*.f

%.d : %.F
	$(CPP) -P $(CDEFS) $*.F $*.f
	$(PERL) $(CP2KHOME)/tools/sfmakedepend -m int -s -f $*.d $*.f
	rm -f $*.f $*.d.old

depend:
	../tools/deps

clean:
	rm -f *.o *.mod

realclean:
	rm -f $(PROG) *.o *.mod *.d *.int *~ F*.f *.lst

####################

include $(OBJECTS:.o=.d)
