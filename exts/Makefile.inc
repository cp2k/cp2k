ifeq ("","$(wildcard $(EXTSHOME)/dbcsr/src)")
$(error "No DBCSR submodule available, please run 'git submodule update --init --recursive'")
endif

LIBS     := -L$(LIBEXTSDIR)/dbcsr -ldbcsr $(LIBS)
EXTSDEPS_LIB  = $(LIBEXTSDIR)/dbcsr/libdbcsr$(ARCHIVE_EXT)
EXTSDEPS_MOD = $(OBJEXTSDIR)/dbcsr/dbcsr_api.mod $(OBJEXTSDIR)/dbcsr/dbcsr_tensor_api.mod
$(EXTSDEPS_MOD) : ; # override builtin .mod rule to prevent circular dependency

# Set Acceleration flags
DBCSR_USE_ACCEL := ""
DBCSR_ACC       := ""
DBCSR_ACCFLAGS  := ""
ifneq (,$(findstring __DBCSR_ACC,$(OFFLOAD_FLAGS)))
# By default, we use CUDA if OFFLOAD_TARGET is not specified
ifeq (,$(OFFLOAD_TARGET))
DBCSR_USE_ACCEL=cuda
else
DBCSR_USE_ACCEL = $(OFFLOAD_TARGET)
endif
ifeq (cuda,$(DBCSR_USE_ACCEL))
DBCSR_ACC      := $(OFFLOAD_CC)
DBCSR_ACCFLAGS := $(OFFLOAD_FLAGS)
else ifeq (hip,$(DBCSR_USE_ACCEL))
DBCSR_ACC := $(CXX)
DBCSR_ACCFLAGS := $(CXXFLAGS)
endif
endif

extversions: dbcsrversion

extclean: dbcsrclean
	@echo "Clean EXT"

dbcsr:
	+$(MAKE) -C $(EXTSHOME)/$@ -f $(EXTSHOME)/build_$@/Makefile \
	   ARCHFILE=$(ARCHDIR)/$(ARCH).$(ONEVERSION) \
	   LIBDIR=$(LIBEXTSDIR)/$@ \
	   OBJDIR=$(OBJEXTSDIR)/$@ \
	   USE_ACCEL="$(DBCSR_USE_ACCEL)" \
	   ACC="$(DBCSR_ACC)" \
	   ACCFLAGS="$(DBCSR_ACCFLAGS)"

dbcsrversion:
	@$(MAKE) -C $(EXTSHOME)/dbcsr -f $(EXTSHOME)/build_dbcsr/Makefile \
	   ARCHFILE=$(ARCHDIR)/$(ARCH).$(ONEVERSION) \
	   USE_ACCEL="$(DBCSR_USE_ACCEL)" \
	   ACC="$(DBCSR_ACC)" \
	   ACCFLAGS="$(DBCSR_ACCFLAGS)" \
	   version

dbcsrclean:
	@echo "Clean DBCSR"
	@$(MAKE) -C $(EXTSHOME)/dbcsr -f $(EXTSHOME)/build_dbcsr/Makefile \
	   ARCHFILE=$(ARCHDIR)/$(ARCH).$(ONEVERSION) \
	   LIBDIR=$(LIBEXTSDIR)/$@ \
	   OBJDIR=$(OBJEXTSDIR)/$@ \
	   USE_ACCEL="$(DBCSR_USE_ACCEL)" \
	   ACC="$(DBCSR_ACC)" \
	   ACCFLAGS="$(DBCSR_ACCFLAGS)" \
	   clean
