EXTRA_DIST = recipe 
ACLOCAL_AMFLAGS = -I m4
SUBDIRS = @BIOMCMCLIB@ lib src tests
# https://stackoverflow.com/questions/25065337/building-nested-autotools-package-without-installing
install installdirs: SUBDIRS = lib src tests # do not install biomcmc-lib (only works with --disable-shared atm) 
