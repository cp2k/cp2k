.SUFFIXES: .o .F

include ../MACHINEDEFS

OBJECTS_fftsg = mltfftsg.o mltfftsg_tools.o ctrig.o fftstp.o fftpre.o fftrot.o

#################################

LIB1	= fftsg

all: $(LIB1)

$(LIB1): $(OBJECTS_fftsg)
	$(AR) libfftsg.a $(OBJECTS_fftsg) 

%.o : %.F
	$(FC) -c $(FFLAGS) $*.F

clean:
	rm -f *.o *.mod

realclean:
	rm -f $(LIB) *.o *.mod *~ F*.f *.lst

####################
