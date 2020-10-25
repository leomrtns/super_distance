AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib -I$(srcdir)/../@BIOMCMCLIB@/lib \
							@OPENMP_CPPFLAGS@  @ZLIB_CFLAGS@ @LZMA_CFLAGS@ 
AM_CFLAGS = @AM_CFLAGS@ @OPENMP_CFLAGS@ @ZLIB_CFLAGS@ @LZMA_CFLAGS@ @CHECK_CFLAGS@
LDADD = ../lib/libsuper_distance.la ../biomcmc-lib/lib/libbiomcmc_static.la \
				$(GTKDEPS_LIBS) @CHECK_LIBS@ @ZLIB_LIBS@ @LZMA_LIBS@ $(AM_LDFLAGS)

LIST_OF_TEST_PROGS= check_unit
# duplicate, since we want all to be compiled only with 'make check'
TESTS = $(LIST_OF_TEST_PROGS)           # list of test programs
check_PROGRAMS = $(LIST_OF_TEST_PROGS)  # list of programs to be compiled only with 'make check'

check_unit_SOURCES = check_unit.c
