dnl @synopsis ACX_PARPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the PARPACK
dnl linear-algebra interface (see http://www.netlib.org/parpack/).
dnl On success, it sets the PARPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with PARPACK, you should link with:
dnl
dnl 	$PARPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-parpack=<lib> in order to use some
dnl specific PARPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the PARPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a PARPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_PARPACK.
dnl
dnl @version $Id: acx_parpack.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_PARPACK], [
AC_REQUIRE([ACX_ARPACK])
AC_REQUIRE([ACX_BLAS])
AC_REQUIRE([ACX_LAPACK])
#AC_REQUIRE([ACX_MPI])
acx_parpack_ok=no

AC_ARG_WITH(parpack,
	[AC_HELP_STRING([--with-parpack=<lib>], [use PARPACK library <lib>])])
case $with_parpack in
	yes | "") ;;
	no) acx_parpack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) PARPACK_LIBS="$with_parpack" ;;
	*) PARPACK_LIBS="-l$with_parpack" ;;
esac

# We cannot use PARPACK if ARPACK is not found
if test "x$acx_blas_ok" != "xyes"; then
	acx_parpack_ok=noarpack
fi

# We cannot use PARPACK if BLAS is not found
if test "x$acx_blas_ok" != "xyes"; then
	acx_parpack_ok=noblas
fi

# We cannot use PARPACK if LAPACK is not found
if test "x$acx_lapack_ok" != "xyes"; then
	acx_parpack_ok=nolapack
fi

# We cannot use PARPACK if MPI is not found
#if test "x$acx_mpi_ok" != "xyes"; then
#	acx_parpack_ok=nompi
#fi

# Get fortran linker name of PARPACK function to check for.
AC_F77_FUNC(pdgetv0)

# First, check PARPACK_LIBS environment variable
if test "x$PARPACK_LIBS" != "x"; then
	save_LIBS="$LIBS"; LIBS="$PARPACK_LIBS $ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
	#save_LIBS="$LIBS"; LIBS="$PARPACK_LIBS $ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $MPI_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $pdgetv0 in $PARPACK_LIBS])
	AC_TRY_LINK_FUNC($pdgetv0, [acx_parpack_ok=yes], [PARPACK_LIBS=""])
	AC_MSG_RESULT($acx_parpack_ok)
	LIBS="$save_LIBS"
	if test "$acx_parpack_ok" = "no"; then
		PARPACK_LIBS=""
	fi
fi

# PARPACK linked to by default?
if test "$acx_parpack_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
	#save_LIBS="$LIBS"; LIBS="$ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $MPI_LIBS $LIBS $FLIBS"
	AC_CHECK_FUNC($pdgetv0, [acx_parpack_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic PARPACK library?
if test "$acx_parpack_ok" = "no"; then
	AC_CHECK_LIB(parpack, $pdgetv0,
		[acx_parpack_ok=yes; PARPACK_LIBS="-lparpack"], [],
		[$ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $FLIBS])
		#[$ARPACK_LIBS $LAPACK_LIBS $BLAS_LIBS $MPI_LIBS $FLIBS])
fi

AC_SUBST(PARPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_parpack_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_PARPACK,1,[Define if you have PARPACK library.]),[$1])
        :
else
        acx_parpack_ok=no
        $2
fi
])dnl ACX_PARPACK
