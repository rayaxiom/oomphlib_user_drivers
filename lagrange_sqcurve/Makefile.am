noinst_PROGRAMS= \
sqcurve_po \
sqcurve_van 

#------------------------------------------------------------------------
sqcurve_po_SOURCES = sqcurve_po.cc

sqcurve_po_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------
sqcurve_van_SOURCES = sqcurve_van.cc

sqcurve_van_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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

# Print top_builddir
print-top_builddir: 
	echo $(top_builddir)
