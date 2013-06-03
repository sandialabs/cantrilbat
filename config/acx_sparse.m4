dnl @synopsis ACX_SPARSE([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the SPARSE
dnl linear-algebra interface (see http://www.netlib.org/sparse/).
dnl On success, it sets the SPARSE_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with SPARSE, you should link with:
dnl
dnl 	$SPARSE_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-sparse=<lib> in order to use some
dnl specific SPARSE library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the SPARSE and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a SPARSE
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_SPARSE.
dnl
dnl @version $Id: acx_sparse.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
AC_DEFUN([ACX_SPARSE], [
acx_sparse_ok=no

AC_ARG_WITH(sparse,
	[AC_HELP_STRING([--with-sparse=<lib>], [use SPARSE library <lib>])])
case $with_sparse in
	yes | "") SPARSE_LIBS="-lsparse" ;;
	no) acx_sparse_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) SPARSE_LIBS="$with_sparse" ;;
	*) SPARSE_LIBS="-l$with_sparse" ;;
esac

# Set linker name of SPARSE function to check for.
sparse_func="spSolve"

# First, check SPARSE_LIBS environment variable
if test "x$SPARSE_LIBS" != "x"; then
	save_LIBS="$LIBS"; LIBS="$SPARSE_LIBS $LIBS $FLIBS"
	AC_MSG_CHECKING([for $sparse_func in $SPARSE_LIBS])
	AC_TRY_LINK_FUNC($sparse_func, [acx_sparse_ok=yes], [SPARSE_LIBS=""])
	AC_MSG_RESULT($acx_sparse_ok)
	LIBS="$save_LIBS"
	if test "$acx_sparse_ok" = "no"; then
		SPARSE_LIBS=""
	fi
fi

# SPARSE linked to by default?
if test "$acx_sparse_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$LIBS $FLIBS"
	AC_CHECK_FUNC($sparse_func, [acx_sparse_ok=yes])
	LIBS="$save_LIBS"
fi

# Generic SPARSE library?
if test "$acx_sparse_ok" = "no"; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_LIB(sparse, $sparse_func,
		[acx_sparse_ok=yes; SPARSE_LIBS="-lsparse"], [], [$FLIBS])
	LIBS="$save_LIBS"
fi

AC_SUBST(SPARSE_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$acx_sparse_ok" = "xyes"; then
        ifelse([$1],,AC_DEFINE(HAVE_SPARSE,1,[Define if you have SPARSE library.]),[$1])
        :
else
        acx_sparse_ok=no
        $2
fi
])dnl ACX_SPARSE
