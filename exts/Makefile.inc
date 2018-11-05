
LIBS     := -L$(LIBEXTSDIR)/dbcsr -ldbcsr $(LIBS)
EXTSDEPS  = $(LIBEXTSDIR)/dbcsr/libdbcsr.a
dbcsr: 
	$(MAKE) -C $(EXTSHOME)/$@ \
	   INCLUDEMAKE=$(CP2KHOME)/arch/$(ARCH).$(ONEVERSION) \
	   LIBDIR=$(LIBEXTSDIR)/$@ \
	   OBJDIR=$(OBJEXTSDIR)/$@ \
	   FYPPBIN=$(TOOLSRC)/build_utils/fypp/bin/fypp
	@cp $(OBJEXTSDIR)/$@/*.mod $(OBJDIR)/
