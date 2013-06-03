dnl @synopsis ACX_EXODUSII([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a SEAMS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_EXODUSII.
dnl
dnl @version $Id: acx_exodusii.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_EXODUSII], [

# variables flag if user has provided this option
acx_exodusii_inc_provided=no
acx_exodusii_lib_provided=no

# variables flag if provided options are valid
acx_exodusii_inc_ok=no
acx_exodusii_lib_ok=no

# variable flags if enough of all provided options are valid
acx_exodusii_ok=no

# Exodus II absolutely requires a netCDF library in order to work.

AC_REQUIRE([ACX_NETCDF])

# If the user wants to use SEAMS to provide EXODUS II, find out...
AC_REQUIRE([ACX_SEAMS])

EXODUSII_LIBS=""
EXODUSII_INCS=""

# User tells where to look for the Exodus II includes...

AC_ARG_WITH([exodusii-inc],
	    [AC_HELP_STRING([--with-exodusii-inc=PATH],
	                    [find ExodusII includes in PATH])
	    ],
	    [acx_exodusii_inc_provided=yes
	     case $withval in
		-*) EXODUSII_INCS=${withval} ;;
		*) EXODUSII_INCS=-I${withval} ;;
             esac
	     CPPFLAGS="${EXODUSII_INCS} ${CPPFLAGS}"
	     AC_CHECK_HEADERS(exodusII.h,acx_exodusii_inc_ok=yes)],
	    [acx_exodus_inc_provided=no
	     if test x"$acx_seams_incroot" != x"" ; then
		AC_MSG_CHECKING([for exodusII.h in])
		AC_MSG_RESULT([$acx_seams_incroot])
		EXODUSII_INCS="-I$acx_seams_incroot"
		CPPFLAGS="${EXODUSII_INCS} ${CPPFLAGS}"
		AC_CHECK_HEADERS(exodusII.h,acx_exodusii_inc_ok=yes)
             else
		AC_CHECK_HEADERS(exodusII.h,acx_exodusii_inc_ok=yes)
             fi
])

if test "$acx_exodusii_inc_ok" = "no" ; then
	AC_MSG_WARN([No exodusII.h include file found.])
fi		

# User says where to look for the Exodus II library...

AC_ARG_WITH([exodusii-lib],
	    [AC_HELP_STRING([--with-exodusii-lib=PATH], 
	                    [use ExodusII library in PATH])],
	    [acx_exodusii_lib_provided=yes
  	     AC_MSG_CHECKING([for EXODUS II library in])
	     AC_MSG_RESULT([$withval])
	     case $withval in
                -* | *.a | *.so | *.so.* | *.o) EXODUSII_LIBS="$withval" ;;
	        *) EXODUSII_LIBS="-L$withval" ;;
             esac
	     LIBS="$EXODUSII_LIBS $LIBS"
	     AC_CHECK_LIB(exoIIv2c, ex_open, [acx_exodusii_lib_ok=yes],
			     [acx_exodusii_lib_ok=no],[-lnetcdf])],
	    [acx_exodus_lib_provided=no
	     if test x"$acx_seams_usrlibdir" != x"" ; then
		AC_MSG_CHECKING([for EXODUS II library in])
		AC_MSG_RESULT([$acx_seams_usrlibdir])
		EXODUSII_LIBS="-L$acx_seams_usrlibdir"
		LIBS="$EXODUSII_LIBS $LIBS"
		echo libs are $LIBS
		AC_CHECK_LIB(exoIIv2c, ex_open, [acx_exodusii_lib_ok=yes],[],
		[-lnetcdf])
	     fi])

#
# Just in case things have fallen through this far without flagging
# problems of missing includes or libraries for Exodus II...
#

if test x"$acx_exodusii_inc_ok" = x"no" ; then
	AC_MSG_ERROR([Can not build without a valid EXODUS II includes.])
fi

if test x"$acx_exodusii_lib_ok" = x"no" ; then
	AC_MSG_ERROR([Can not build without a valid EXODUS II library.])
fi

AC_SUBST(EXODUSII_LIBS)
AC_SUBST(EXODUSII_INCS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_exodusii_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_EXODUSII,1,[Define if you have EXODUSII library.]),[$1])
        :
else
        acx_exodusii_ok=no
        $2
fi
])dnl ACX_EXODUSII
