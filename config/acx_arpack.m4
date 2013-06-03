dnl @synopsis ACX_ARPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the ARPACK
dnl linear-algebra interface (see http://www.netlib.org/arpack/).
dnl On success, it sets the ARPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with ARPACK, you should link with:
dnl
dnl 	$ARPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-arpack=<lib> in order to use some
dnl specific ARPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the ARPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a ARPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_ARPACK.
dnl
dnl @version $Id: acx_arpack.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_ARPACK], [
AC_REQUIRE([ACX_BLAS])
AC_REQUIRE([ACX_LAPACK])
acx_arpack_ok=no

AC_ARG_WITH(arpack,
	[AC_HELP_STRING([--with-arpack=<lib>], [use ARPACK library <lib>])])
case $with_arpack in
	yes | "") ;;
	no) acx_arpack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) ARPACK_LIBS="$with_arpack" ;;
	*) ARPACK_LIBS="-l$with_arpack" ;;
esac

# We cannot use ARPACK if BLAS is not found
if test "x$acx_blas_ok" != "xyes"; then
	acx_arpack_ok=noblas
fi

# We cannot use ARPACK if LAPACK is not found
if test "x$acx_lapack_ok" != "xyes"; then
	acx_arpack_ok=nolapack
fi

# Get fortran linker name of ARPACK function to check for.
AC_F77_FUNC(dstats)

# First, check ARPACK_LIBS environment variable
if test "x$ARPACK_LIBS" != "x"; then
	save_LIBS="$LIBS"; LIBS="$ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $dstats in $ARPACK_LIBS])
	AC_TRY_LINK_FUNC($dstats, [acx_arpack_ok=yes], [ARPACK_LIBS=""])
	AC_MSG_RESULT($acx_arpack_ok)
	LIBS="$save_LIBS"
	if test "$acx_arpack_ok" = "no"; then
		ARPACK_LIBS=""
	fi
fi

# ARPACK linked to by default?  (is sometimes included in BLAS lib)
if test "$acx_arpack_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
	AC_CHECK_FUNC($dstats, [acx_arpack_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic ARPACK library?
if test "$acx_arpack_ok" = "no"; then
	AC_CHECK_LIB(arpack, $dstats,
		[acx_arpack_ok=yes; ARPACK_LIBS="-larpack"], [],
		[$LAPACK_LIBS $BLAS_LIBS $FLIBS])
fi

AC_SUBST(ARPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_arpack_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_ARPACK,1,[Define if you have ARPACK library.]),[$1])
        :
else
        acx_arpack_ok=no
        $2
fi
])dnl ACX_ARPACK
