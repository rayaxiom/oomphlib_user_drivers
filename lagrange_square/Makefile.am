noinst_PROGRAMS= \
sq_lgr \
sq_tf_po \
sq_tf \
sq_van_po_semipara \
sq_po_semipara \
sq_tf_semipara \
sq_po_tf_semipara \
sq_van 

#------------------------------------------------------------------------
sq_lgr_SOURCES = sq_lgr.cc

sq_lgr_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
sq_tf_SOURCES = sq_tf.cc

sq_tf_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
sq_van_po_semipara_SOURCES = sq_van_po_semipara.cc

sq_van_po_semipara_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
sq_po_semipara_SOURCES = sq_po_semipara.cc

sq_po_semipara_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
sq_tf_semipara_SOURCES = sq_tf_semipara.cc

sq_tf_semipara_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
#------------------------------------------------------------------------
sq_po_tf_semipara_SOURCES = sq_po_tf_semipara.cc

sq_po_tf_semipara_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#------------------------------------------------------------------------
sq_tf_po_SOURCES = sq_tf_po.cc

sq_tf_po_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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
