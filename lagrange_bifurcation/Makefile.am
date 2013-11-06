noinst_PROGRAMS= \
unstructured_three_d_fluid 

#------------------------------------------------------------------------
unstructured_three_d_fluid_SOURCES = unstructured_three_d_fluid.cc

unstructured_three_d_fluid_LDADD = -L@libdir@ -lnavier_stokes -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

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
