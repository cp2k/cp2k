


include ../config.mk

PYTHON=python

default: xyz_to_vab.o

lib_kinds.mod:
	$(FC_comp) $(FCFLAGS) -c lib_kinds.F

xyz_to_vab.o: lib_kinds.mod
	${PYTHON} template_parser.py xyz_to_vab_optimised.mod.template xyz_to_vab_optimised.F
	$(FC_comp) $(FCFLAGS)  -c xyz_to_vab_optimised.F



clean: 
	rm xyz_to_vab_optimised.mod xyz_to_vab_optimised.o lib_kinds.o lib_kinds.mod
