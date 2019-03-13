AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../lib -I$(srcdir)/../biomcmc-lib/lib  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../lib/libsuper_sptree.la ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)

bin_PROGRAMS = supersp 

supersp_SOURCES = supersp_main.c
