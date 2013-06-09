noinst_PROGRAMS= \
sq_curve

#------------------------------------------------------------------------
sq_curve_SOURCES = sq_curve.cc

sq_curve_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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
