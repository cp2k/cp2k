ifeq ($(origin FC), default)
FC = gfortran
endif
FFLAGS ?= -O2 -std=f2008

PROGRAM = spatial_distribution.x
SOURCE = spatial_distribution.f90

.PHONY: all clean test

all: $(PROGRAM)

$(PROGRAM): $(SOURCE)
	$(FC) $(FFLAGS) -o $@ $<

test: $(PROGRAM)
	python3 test/test_spatial_distribution.py ./$(PROGRAM)

clean:
	rm -f $(PROGRAM) *.mod
	rm -rf $(PROGRAM).dSYM
