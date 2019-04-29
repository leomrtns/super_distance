AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib -I$(srcdir)/../biomcmc-lib/lib  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../lib/libsuper_distance.la ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)

LIST_OF_TEST_PROGS= check_unit
# duplicate, since we want all to be compiled only with 'make check'
TESTS = $(LIST_OF_TEST_PROGS)           # list of test programs
check_PROGRAMS = $(LIST_OF_TEST_PROGS)  # list of programs to be compiled only with 'make check'

check_unit_SOURCES = check_unit.c
