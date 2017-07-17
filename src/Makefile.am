AUTOMAKE_OPTIONS = foreign

ACLOCAL_AMFLAGS = -I m4 --install
SUFFIXES = .cpp

INCD = -I $(top_builddir)/inc -I $(top_builddir)/printcolor

noinst_LIBRARIES = libevolution.a
libevolution_a_SOURCES = CuResources.cpp subroutines.cpp causet.cpp

.cpp.$(OBJEXT):
	$(COMPILE_STR) $(INCD) -o $@ $<
