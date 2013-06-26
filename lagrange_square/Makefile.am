noinst_PROGRAMS= \
sq_po \
sq_van 

#------------------------------------------------------------------------
sq_po_SOURCES = sq_po.cc

sq_po_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------
sq_van_SOURCES = sq_van.cc

sq_van_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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
