AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ -I$(srcdir)/../biomcmc-lib/lib
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS) ../biomcmc-lib/lib/libbiomcmc.la

common_headers = super_sptree.h distance_supertree.h parsimony_supertree.h
common_src     = distance_supertree.c parsimony_supertree.c

lib_LTLIBRARIES = libsuper_sptree.la   ## change to noinst_LT... if should be linked statically (not installed globally)
## make sure all headers are below this subdirectory (be careful with spaces at the end of line!)
otherincludedir = "$(includedir)/super_sptree"
otherinclude_HEADERS = config.h $(common_headers)
libsuper_sptree_la_SOURCES = $(common_headers) $(common_src) 

