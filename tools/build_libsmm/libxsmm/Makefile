
# Use ROW_MAJOR matrix representation if set to 1, COL_MAJOR otherwise 
ROW_MAJOR ?= 1

# M, N, K values of the generated matrices
INDICES_M ?= $(shell seq 1 8)
INDICES_N ?= $(shell seq 1 8)
INDICES_K ?= $(shell seq 1 8)

DIR_KNC   ?= .
SCRDIR_KNC = $(DIR_KNC)/scripts
OBJDIR_KNC = $(DIR_KNC)/build
INCDIR_KNC = $(DIR_KNC)/include
SRCDIR_KNC = $(DIR_KNC)/src
LIBDIR_KNC = $(DIR_KNC)/lib

INDICES ?= $(foreach m,$(INDICES_M),$(foreach n,$(INDICES_N),$(foreach k,$(INDICES_K),$m_$n_$k)))

TARGET_COMPILE_C_KNC := icc -offload-attribute-target=mic -mkl=sequential -std=c99 -openmp
AR := xiar -qoffload-build

SRCFILES_KNC = $(patsubst %,dc_small_dnn_%.c,$(INDICES))
OBJFILES_KNC = $(patsubst %,$(OBJDIR_KNC)/dc_small_dnn_%.o,$(INDICES))

LIB_KNC  ?= $(LIBDIR_KNC)/libxsmm.a
INC_KNC   = $(INCDIR_KNC)/xsmm_knc.h
MAIN_KNC  = $(SRCDIR_KNC)/xsmm_knc.c


lib_knc: $(LIB_KNC)
ifeq ($(origin NO_MAIN), undefined)
$(LIB_KNC): $(OBJFILES_KNC) $(patsubst $(SRCDIR_KNC)/%.c,$(OBJDIR_KNC)/%.o,$(MAIN_KNC))
else
$(LIB_KNC): $(OBJFILES_KNC)
endif
	@mkdir -p $(LIBDIR_KNC)
	$(AR) -rs $@ $^

compile_knc: $(OBJFILES_KNC)
$(OBJDIR_KNC)/%.o: $(SRCDIR_KNC)/%.c
	@mkdir -p $(OBJDIR_KNC)
	$(TARGET_COMPILE_C_KNC) -I$(INCDIR_KNC) -c $< -o $@

source_knc: $(addprefix $(SRCDIR_KNC)/,$(SRCFILES_KNC))
$(SRCDIR_KNC)/%.c:
	@mkdir -p $(SRCDIR_KNC)
	@python $(SCRDIR_KNC)/xsmm_knc_gensrc.py `echo $* | awk -F_ '{ print $$4" "$$5" "$$6 }'` $(ROW_MAJOR) > $@

main_knc: $(MAIN_KNC)
$(MAIN_KNC): $(INC_KNC)
	@mkdir -p $(SRCDIR_KNC)
	@python $(SCRDIR_KNC)/xsmm_knc_genmain.py $(words $(INDICES_M)) $(words $(INDICES_N)) $(ROW_MAJOR) $(INDICES_M) $(INDICES_N) $(INDICES_K) > $@

header_knc: $(INC_KNC)
$(INC_KNC):
	@echo "#ifndef XSMM_KNC_H" > $@
	@echo "#define XSMM_KNC_H" >> $@
	@echo >> $@
	@python $(SCRDIR_KNC)/xsmm_knc_geninc.py >> $@
	@python $(SCRDIR_KNC)/xsmm_knc_geninc.py $(INDICES) >> $@
	@echo >> $@
	@echo "#endif // XSMM_KNC_H" >> $@

clean:
	rm -rf $(SRCDIR_KNC) $(OBJDIR_KNC) $(DIR_KNC)/*~ $(DIR_KNC)/*/*~

realclean: clean
	rm -rf $(LIBDIR_KNC) $(INC_KNC)
