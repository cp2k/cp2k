FC = gfortran
FFLAGS ?= -O2 -std=f2008

COMMON_OBJECTS = \
        kinds.o \
        text_utils.o \
        trajectory_types.o \
        geometry.o \
        trajectory_io.o \
        command_line.o \
        cell_source.o \
        frame_controls.o \
        selections.o \
        groups.o \
        fft.o \
        time_series.o \
        spectral_tools.o

ANALYSIS_OBJECTS = \
        geometry_analysis.o \
        rdf_analysis.o \
        vacf_analysis.o \
        hbond_analysis.o \
        dsf_analysis.o

.PHONY: all clean test check

all: trajana.x trajconvert.x twopt.x

trajana.x: $(COMMON_OBJECTS) $(ANALYSIS_OBJECTS) trajana_main.o
	$(FC) $(FFLAGS) -o $@ $^

trajconvert.x: $(COMMON_OBJECTS) trajconvert.o
	$(FC) $(FFLAGS) -o $@ $^

twopt.x: $(COMMON_OBJECTS) twopt_thermodynamics.o twopt_analysis.o twopt_main.o
	$(FC) $(FFLAGS) -o $@ $^

trajectory_types.o: kinds.o
geometry.o: kinds.o trajectory_types.o
trajectory_io.o: kinds.o text_utils.o trajectory_types.o
command_line.o: kinds.o
cell_source.o: trajectory_types.o trajectory_io.o text_utils.o
selections.o: text_utils.o
groups.o: kinds.o text_utils.o
fft.o: kinds.o
time_series.o: kinds.o fft.o
spectral_tools.o: kinds.o
geometry_analysis.o: kinds.o text_utils.o trajectory_types.o trajectory_io.o cell_source.o command_line.o frame_controls.o geometry.o
rdf_analysis.o: kinds.o trajectory_types.o trajectory_io.o cell_source.o command_line.o frame_controls.o geometry.o selections.o
vacf_analysis.o: kinds.o text_utils.o trajectory_types.o trajectory_io.o command_line.o frame_controls.o selections.o time_series.o spectral_tools.o
hbond_analysis.o: kinds.o trajectory_types.o trajectory_io.o cell_source.o command_line.o frame_controls.o geometry.o groups.o
dsf_analysis.o: kinds.o text_utils.o trajectory_types.o trajectory_io.o command_line.o frame_controls.o selections.o time_series.o spectral_tools.o
trajana_main.o: text_utils.o command_line.o $(ANALYSIS_OBJECTS)
trajconvert.o: kinds.o text_utils.o trajectory_types.o trajectory_io.o cell_source.o command_line.o frame_controls.o geometry.o groups.o
twopt_thermodynamics.o: kinds.o
twopt_analysis.o: kinds.o text_utils.o trajectory_types.o trajectory_io.o cell_source.o command_line.o frame_controls.o geometry.o groups.o time_series.o twopt_thermodynamics.o
twopt_main.o: command_line.o twopt_analysis.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

test: all
	python3 tests/test_trajana.py

check:
	$(MAKE) clean
	$(MAKE) FFLAGS="-O0 -g -std=f2008 -Wall -Wextra -Werror -fcheck=all" all
	python3 tests/test_trajana.py

clean:
	rm -f *.o *.mod *.MOD trajana.x trajconvert.x twopt.x
