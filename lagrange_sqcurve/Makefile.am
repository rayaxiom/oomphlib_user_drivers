noinst_PROGRAMS= \
sqcurve_po_dirichletall \
sqcurve_po_dirichletside \
sqcurve_po_neumannside \
milan_quarter_annulus 

#------------------------------------------------------------------------
sqcurve_po_dirichletall_SOURCES = sqcurve_po_dirichletall.cc

sqcurve_po_dirichletall_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------

sqcurve_po_dirichletside_SOURCES = sqcurve_po_dirichletside.cc

sqcurve_po_dirichletside_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------

sqcurve_po_neumannside_SOURCES = sqcurve_po_neumannside.cc

sqcurve_po_neumannside_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------
milan_quarter_annulus_SOURCES = milan_quarter_annulus.cc

milan_quarter_annulus_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)


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
