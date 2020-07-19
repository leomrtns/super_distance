AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  @OPENMP_CPPFLAGS@ -I$(srcdir)/../@BIOMCMCLIB@/lib
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) $(AM_LDFLAGS) ../biomcmc-lib/lib/libbiomcmc_static.la

common_headers = super_distance.h distance_supertree.h 
common_src     = distance_supertree.c 

noinst_LTLIBRARIES = libsuper_distance.la   ## change to noinst_LT... if should be linked statically (not installed globally)
## make sure all headers are below this subdirectory (be careful with spaces at the end of line!)
otherincludedir = "$(includedir)/super_distance"
otherinclude_HEADERS = config.h $(common_headers)
libsuper_distance_la_SOURCES = $(common_headers) $(common_src) 
