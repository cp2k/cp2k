PROG  =	leap2fist

OBJS  =	leap2fist.o

LIBS =	

FC = f95
ifeq ($(MAKECMDGOALS),debug)
	FFLAGS  = -O0 -c -g -C=all  -mismatch_all -f77 -w=uei -w=uda -nan -gline
	LDFLAGS = -O0 -g -C=all  -mismatch_all -f77 -w=uei -w=uda -nan -gline
else
	FFLAGS  = -O3 -c
	LDFLAGS = -O3
endif

all: $(PROG)
debug: $(PROG)

$(PROG): $(OBJS)
	 $(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS)

#----------------------------------------------------------------------------
# Explicit rules
#----------------------------------------------------------------------------
.SUFFIXES:
.SUFFIXES: .o .f90

$(OBJS) : 
	@echo compiling file: $<
	$(FC) $(FFLAGS) $(@:.o=.f90)

