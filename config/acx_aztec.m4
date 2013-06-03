dnl @synopsis ACX_AZTEC([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the AZTEC
dnl linear-algebra interface (see http://www.netlib.org/aztec/).
dnl On success, it sets the AZTEC_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with AZTEC, you should link with:
dnl
dnl 	$AZTEC_LIBS $LAPACK_LIBS $BLAS_LIBS $Y12M_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-aztec=<lib> in order to use some
dnl specific AZTEC library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the AZTEC and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a AZTEC
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_AZTEC.
dnl
dnl @version $Id: acx_aztec.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_AZTEC], [
AC_REQUIRE([ACX_BLAS])
AC_REQUIRE([ACX_Y12M])
AC_REQUIRE([ACX_LAPACK])

acx_aztec_ok=no
acx_aztec_inc_ok=no
acx_aztec_lib_ok=no

AC_ARG_WITH(aztec-lib,
	[AC_HELP_STRING([--with-aztec-lib=<lib>], [use AZTEC library <lib>])])
case $with_aztec_lib in
	yes | "") ;;
	no) acx_aztec_lib_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) AZTEC_LIBS="$with_aztec_lib" ;;
	*) AZTEC_LIBS="-l$with_aztec_lib" ;;
esac

AC_ARG_WITH(aztec-inc,
	[AC_HELP_STRING([--with-aztec-inc=<inc>],[use AZTEC include path <inc>])])
case $with_aztec_inc in
	yes | "") ;;
	no) acx_aztec_inc_ok=disable ;;
	-* ) AZTEC_INCS="$with_aztec_inc" ;;
	*) AZTEC_INCS="-I$with_aztec_inc" ;;
esac

# We cannot use AZTEC if BLAS is not found
if test "x$acx_blas_ok" != "xyes"; then
	acx_aztec_ok=noblas
	AC_MSG_WARN([The Aztec solver package requires the BLAS library.])
fi

# We cannot use AZTEC if LAPACK is not found
if test "x$acx_lapack_ok" != "xyes"; then
	acx_aztec_ok=nolapack
	AC_MSG_WARN([The Aztec solver package requires the LAPACK library.])
fi

# We cannot use AZTEC if Y12M is not found
if test "x$acx_y12m_ok" != "xyes"; then
	acx_aztec_ok=nolapack
	AC_MSG_WARN([The Aztec solver package requires the Y12M library.])
fi

# Check for the main AZTEC header file
az_hdr="az_aztec.h"

# Check the AZTEC_INCS environment variable
if test "x$AZTEC_INCS" != "x"; then
	CPPFLAGS_SAVE="$CPPFLAGS"
	CPPFLAGS="$CPPFLAGS $AZTEC_INCS"
	AC_CHECK_HEADERS([$az_hdr],
		[acx_aztec_inc_ok=yes],
		[acx_aztec_inc_ok=no])
	if test "$acx_aztec_inc_ok" = "no" ; then
		AZTEC_INCS=""
	fi
	CPPFLAGS="$CPPFLAGS_SAVE"
fi

# Next, check the standard CPP paths
if test "$acx_aztec_inc_ok" = "no" ; then
	AZTEC_INCS=""
	AC_CHECK_HEADERS([$az_hdr],
		[acx_aztec_inc_ok=yes],
		[acx_aztec_inc_ok=no])
fi


# Set linker name of AZTEC function to check for.
#az_func="AZ_set_proc_config"
az_func="AZ_solve"

EXTRA_LIBS="$LAPACK_LIBS $Y12M_LIBS $BLAS_LIBS"

# First, check the AZTEC_LIBS environment variable
if test "$acx_aztec_lib_ok" = "no" ; then
    if test "x$AZTEC_LIBS" != "x"; then
	save_LIBS="$LIBS"; LIBS="$AZTEC_LIBS $EXTRA_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $az_func in $AZTEC_LIBS])
	AC_TRY_LINK_FUNC($az_func, [acx_aztec_lib_ok=yes], [AZTEC_LIBS=""])
	AC_MSG_RESULT($acx_aztec_lib_ok)
	LIBS="$save_LIBS"
	if test "$acx_aztec_lib_ok" = "no"; then
		AZTEC_LIBS=""
	fi
    fi
fi

# AZTEC linked to by default?  (is sometimes included in BLAS lib)
if test "$acx_aztec_lib_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$EXTRA_LIBS $LIBS $FLIBS"
	AC_CHECK_FUNC($az_func, [acx_aztec_lib_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic AZTEC library?
if test "$acx_aztec_lib_ok" = "no"; then
	AC_CHECK_LIB(aztec$COMM_TAG, $az_func,
		[acx_aztec_lib_ok=yes; AZTEC_LIBS="-laztec$COMM_TAG"], [], [$EXTRA_LIBS $FLIBS])
fi

# OK, we really need both the includes and the libs
if test "$acx_aztec_inc_ok" = "yes" -a "$acx_aztec_lib_ok" = "yes" ; then
	acx_aztec_ok="yes"
else
	AC_MSG_WARN([The Aztec solver package will not be included.])
	AZTEC_LIBS=""
	AZTEC_INCS=""
fi

AC_SUBST(AZTEC_LIBS)
AC_SUBST(AZTEC_INCS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_aztec_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_AZTEC,1,[Define if you have AZTEC library.]),[$1])
        :
else
        acx_aztec_ok=no
        $2
fi
])dnl ACX_AZTEC
