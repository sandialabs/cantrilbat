dnl @synopsis ACX_NETCDF([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a SEAMS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_NETCDF.
dnl
dnl @version $Id: acx_netcdf.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl
AC_DEFUN([ACX_NETCDF], [dnl
# variables flag if user has provided this option
acx_netcdf_inc_provided=no
acx_netcdf_lib_provided=no
# variables flag if provided options are valid
acx_netcdf_inc_ok=no
acx_netcdf_lib_ok=no
# variable flags if enough of all provided options are valid
acx_netcdf_ok=no
# SEAMS can provide netCDF...
AC_REQUIRE([ACX_SEAMS])
NETCDF_LIBS=""
NETCDF_INCS=""
# User tells explicitly where to look for the netCDF include file.
AC_ARG_WITH([netcdf-inc],
	    [AC_HELP_STRING([--with-netcdf-inc=PATH],
	                    [find netCDF includes in PATH])],
	    [acx_netcdf_inc_provided=yes
	     case $withval in
		-*) NETCDF_INCS=${withval} ;;
		*) NETCDF_INCS=-I${withval} ;;
             esac
	     CPPFLAGS="${NETCDF_INCS} ${CPPFLAGS}"
	     AC_CHECK_HEADERS(netcdf.h,acx_netcdf_inc_ok=yes)],
	    [acx_netcdf_inc_provided=no
             if test x"$acx_seams_incroot" != x"" ; then
		AC_MSG_CHECKING([for netcdf.h in])
		AC_MSG_RESULT([$acx_seams_incroot])
		NETCDF_INCS="-I$acx_seams_incroot"
		CPPFLAGS="${NETCDF_INCS} ${CPPFLAGS}"
		AC_CHECK_HEADERS(netcdf.h,acx_netcdf_inc_ok=yes)
             else
		AC_CHECK_HEADERS(netcdf.h,acx_netcdf_inc_ok=yes)
             fi])
			
if test x"$acx_netcdf_inc_ok" = x"no" ; then
	AC_MSG_WARN([No netcdf.h include file found.])
fi

# Now, say user specifies where to look for the netCDF library...

AC_ARG_WITH([netcdf-lib],
	    [AC_HELP_STRING([--with-netcdf-lib=PATH], 
	                    [use netCDF library in PATH])],
	    [acx_netcdf_lib_provided=yes
             AC_MSG_CHECKING([user-specified netCDF library in])
	     AC_MSG_RESULT([$withval])
	     case $withval in
                -* | *.a | *.so | *.so.* | *.o) NETCDF_LIBS="$withval" ;;
	        *) NETCDF_LIBS="-L$withval" ;;
             esac
	     LIBS="$NETCDF_LIBS $LIBS"
	     AC_CHECK_LIB(netcdf, nc_open, acx_netcdf_lib_ok=yes,
			     acx_netcdf_lib_ok=no)],
	    [acx_netcdf_lib_provided=no])

# Library not specified, but SEAMS was specified or found. We'll check to see
# if that one is valid. If not, warn.

if test x"$acx_netcdf_lib_provided" = x"no" ; then
	if test x"$acx_seams_usrlibdir" != x"" ; then
		AC_MSG_CHECKING([for netCDF library in])
		AC_MSG_RESULT([$acx_seams_usrlibdir])
		NETCDF_LIBS="-L$acx_seams_usrlibdir"
		LIBS="$NETCDF_LIBS $LIBS"
		AC_CHECK_LIB(netcdf, nc_open, acx_netcdf_lib_ok=yes,
			     acx_netcdf_lib_ok=no)
	fi
fi

# Finally, check if the specified netCDF or SEAMS locations work...

AC_SUBST(NETCDF_LIBS)
AC_SUBST(NETCDF_INCS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_netcdf_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_NETCDF,1,[Define if you have NETCDF library.]),[$1])
        :
else
        acx_netcdf_ok=no
        $2
fi
])dnl ACX_NETCDF
