noinst_PROGRAMS= \
step_po \
step_van 

#------------------------------------------------------------------------
step_po_SOURCES = step_po.cc

step_po_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------
step_van_SOURCES = step_van.cc

step_van_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

# Include path for library headers: All library headers live in 
# the include directory which we specify with -I
INCLUDES = -I@includedir@  

# Extra actions for clean and distclean
distclean-local:
	touch Validation
	rm -r -f Validation
	rm -f *~
clean-local:
	touch Validation
	rm -r -f Validation
	rm -f *~
