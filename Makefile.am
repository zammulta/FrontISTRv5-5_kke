PREFIX     = @prefix@
BINDIR     = @bindir@
LIBDIR     = @libdir@
INCLUDEDIR = @includedir@

MAKE       = @make@
AR         = @ar@
CP         = @cp@
RM         = @rm@
MV         = @mv@
RANLIB     = @ranlib@
MKDIR      = @mkdir@

all:
	@cd hecmw1 && $(MAKE)
	@cd fistr1 && $(MAKE)

clean:
	@cd hecmw1 && $(MAKE) clean
	@cd fistr1 && $(MAKE) clean

distclean:
	@cd hecmw1 && $(MAKE) distclean
	@cd fistr1 && $(MAKE) distclean
	- $(RM) ./bin/*

install:
	if [ ! -d $(PREFIX) ]; then \
		$(MKDIR) $(PREFIX); \
	fi
	if [ ! -d $(BINDIR) ]; then \
		$(MKDIR) $(BINDIR); \
	fi
	@cd hecmw1 && $(MAKE) install
	@cd fistr1 && $(MAKE) install

depend:
	@cd hecmw1 && $(MAKE) depend
	@cd fistr1 && $(MAKE) depend
