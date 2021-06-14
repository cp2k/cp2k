ifeq ("","$(wildcard $(EXTSHOME)/dbcsr/.cp2k/Makefile)")
$(error "No DBCSR submodule available, please run 'git submodule update --init --recursive'")
endif

LIBS     := -L$(LIBEXTSDIR)/dbcsr -ldbcsr $(LIBS)
EXTSDEPS_LIB  = $(LIBEXTSDIR)/dbcsr/libdbcsr$(ARCHIVE_EXT)
EXTSDEPS_MOD = $(OBJEXTSDIR)/dbcsr/dbcsr_api.mod $(OBJEXTSDIR)/dbcsr/dbcsr_tensor_api.mod
$(EXTSDEPS_MOD) : ; # override builtin .mod rule to prevent circular dependency

ifneq (,$(findstring __DBCSR_ACC,$(OFFLOAD_FLAGS)))
ACC := $(OFFLOAD_CC)
ACCFLAGS := $(OFFLOAD_FLAGS)
else
ACC := ""
ACCFLAGS := ""
endif

extversions: dbcsrversion

extclean: dbcsrclean
	@echo "Clean EXT"

dbcsr:
	+$(MAKE) -C $(EXTSHOME)/$@ -f .cp2k/Makefile \
	   INCLUDEMAKE=$(ARCHDIR)/$(ARCH).$(ONEVERSION) \
	   LIBDIR=$(LIBEXTSDIR)/$@ \
	   OBJDIR=$(OBJEXTSDIR)/$@ \
	   ACC="$(ACC)" \
	   ACCFLAGS="$(ACCFLAGS)"

dbcsrversion:
	@$(MAKE) -C $(EXTSHOME)/dbcsr -f .cp2k/Makefile \
	   ACC="$(ACC)" \
	   ACCFLAGS="$(ACCFLAGS)" \
	   version

dbcsrclean:
	@echo "Clean DBCSR"
	@$(MAKE) -C $(EXTSHOME)/dbcsr -f .cp2k/Makefile \
	   ACC="$(ACC)" \
	   ACCFLAGS="$(ACCFLAGS)" \
	   clean
