DIR_MIC ?= ../mic

SRCDIR_MIC:=$(DIR_MIC)/src
OBJDIR_MIC:=$(DIR_MIC)/obj
LIBDIR_MIC:=$(DIR_MIC)/lib

SRCFILES_MIC=$(patsubst %,dc_small_dnn_%.c,$(INDICES))
OBJFILES_MIC=$(patsubst %,$(OBJDIR_MIC)/dc_small_dnn_%.o,$(INDICES))
LIB_MIC=libmic_$(firstword $(INDICES))__$(lastword $(INDICES)).a

lib_mic: $(LIBDIR_MIC)/$(LIB_MIC)
$(LIBDIR_MIC)/$(LIB_MIC): $(OBJFILES_MIC)
	@mkdir -p $(LIBDIR_MIC)
	xiar -rs -qoffload-build $@ $^

compile_mic: $(OBJFILES_MIC)
$(OBJDIR_MIC)/%.o: $(SRCDIR_MIC)/%.c
	@mkdir -p $(OBJDIR_MIC)
	${TARGET_COMPILE_C_MIC} -I$(DIR_MIC) -c $< -o $@

source_mic: $(addprefix $(SRCDIR_MIC)/,$(SRCFILES_MIC))
$(SRCDIR_MIC)/%.c:
	@mkdir -p $(SRCDIR_MIC)
	python $(DIR_MIC)/generate_singleigemm.py `echo $* | awk -F_ '{ print $$4" "$$5" "$$6 }'` > $@

clean:
	rm -rf $(SRCDIR_MIC) $(OBJDIR_MIC) $(LIBDIR_MIC) $(DIR_MIC)/*~
