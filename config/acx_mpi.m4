dnl @synopsis ACX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a SEAMS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_MPI.
dnl
dnl @version $Id: acx_mpi.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Patrick Notz <pknotz@sandia.gov>
dnl
AC_DEFUN([ACX_MPI], [dnl

AC_REQUIRE([ACX_MPI_HOME])
MPI_INCS="$MPI_HOME/include"

MPICH_LIBS="-lmpich -lpmpich -lnsl"
MPICH_LIBS2="-lmpe -lmpich -lpmpich -lthread -lsocket -lnsl"
LAM_LIBS="-llammpi++ -lmpi -llammpio -lpmpi -llam -lnsl"

acx_mpi_ok=no
acx_mpi_inc_ok=no
acx_mpi_lib_ok=no

AC_ARG_WITH(mpi-lib,
	[AC_HELP_STRING([--with-mpi-lib=<lib>], [use MPI library <lib>])])
case $with_mpi_lib in
	yes | "") ;;
	no) acx_mpi_lib_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) MPI_LIBS="$with_mpi_lib" ;;
	*) MPI_LIBS="-l$with_mpi_lib" ;;
esac

# Check for the main MPI header file
mpi_hdr="mpi.h"

# Check the MPI_INCS environment variable
if test "$acx_mpi_inc_ok" = "no" ; then
    if test "x$MPI_INCS" != "x"; then
	save_CPPFLAGS="$CPPFLAGS"
	case $MPI_INCS in
	    -* ) MINC="$MPI_INCS";;
	    *) MINC="-I$MPI_INCS";;
	esac
	CPPFLAGS="$CPPFLAGS $MINC"
	AC_CHECK_HEADERS([$mpi_hdr],
		[acx_mpi_inc_ok=yes],
		[acx_mpi_inc_ok=no])
	CPPFLAGS="$save_CPPFLAGS"
	if test "$acx_mpi_inc_ok" = "no" ; then
		MPI_INCS=""
	else
		MPI_INCS="$MINC"
	fi
    fi
fi

# Next, check the standard CPP paths
if test "$acx_mpi_inc_ok" = "no" ; then
	MPI_INCS=""
	AC_CHECK_HEADERS([$mpi_hdr],
		[acx_mpi_inc_ok=yes],
		[acx_mpi_inc_ok=no])
fi

# Set linker name of Mpi function to check for.
mpi_func="MPI_Init"

# First, check the MPI_LIBS environment variable
if test "$acx_mpi_lib_ok" = "no" ; then
    if test "x$MPI_LIBS" != "x"; then
	save_LIBS="$LIBS"; LIBS="$MPI_LIBS $LIBS"
	AC_MSG_CHECKING([for $mpi_func in $MPI_LIBS])
	AC_TRY_LINK_FUNC($mpi_func, [acx_mpi_lib_ok=yes], [MPI_LIBS=""])
	AC_MSG_RESULT($acx_mpi_lib_ok)
	LIBS="$save_LIBS"
	if test "$acx_mpi_lib_ok" = "no"; then
		MPI_LIBS=""
	fi
    fi
fi

# Look for an MPICH implementation under MPI_HOME
if test "$acx_mpi_lib_ok" = "no"; then
	MPI_LIBS="-L$MPI_HOME/lib $MPICH_LIBS"
	save_LIBS="$LIBS"; LIBS="$MPI_LIBS $LIBS"
	AC_MSG_CHECKING([for $mpi_func in $MPI_LIBS])
	AC_TRY_LINK_FUNC($mpi_func, [acx_mpi_lib_ok=yes], [MPI_LIBS=""])
	AC_MSG_RESULT($acx_mpi_lib_ok)
	LIBS="$save_LIBS"
	if test "$acx_mpi_lib_ok" = "no"; then
		MPI_LIBS=""
	fi
fi

# Look for an MPICH implementation under MPI_HOME
if test "$acx_mpi_lib_ok" = "no"; then
	MPI_LIBS="-L$MPI_HOME/lib $MPICH_LIBS2"
	save_LIBS="$LIBS"; LIBS="$MPI_LIBS $LIBS"
	AC_MSG_CHECKING([for $mpi_func in $MPI_LIBS])
	AC_TRY_LINK_FUNC($mpi_func, [acx_mpi_lib_ok=yes], [MPI_LIBS=""])
	AC_MSG_RESULT($acx_mpi_lib_ok)
	LIBS="$save_LIBS"
	if test "$acx_mpi_lib_ok" = "no"; then
		MPI_LIBS=""
	fi
fi

# Look for a LAM implementation under MPI_HOME
if test "$acx_mpi_lib_ok" = "no"; then
	MPI_LIBS="-L$MPI_HOME/lib $LAM_LIBS"
	save_LIBS="$LIBS"; LIBS="$MPI_LIBS $LIBS"
	AC_MSG_CHECKING([for $mpi_func in $MPI_LIBS])
	AC_TRY_LINK_FUNC($mpi_func, [acx_mpi_lib_ok=yes], [MPI_LIBS=""])
	AC_MSG_RESULT($acx_mpi_lib_ok)
	LIBS="$save_LIBS"
	if test "$acx_mpi_lib_ok" = "no"; then
		MPI_LIBS=""
	fi
fi

# MPI linked to by default?  This could happen with mpicc, etc.
if test "$acx_mpi_lib_ok" = "no"; then
	AC_CHECK_FUNC($mpi_func, [acx_mpi_lib_ok=yes])
fi

# OK, we really need both the includes and the libs
if test "$acx_mpi_inc_ok" = "yes" -a "$acx_mpi_lib_ok" = "yes" ; then
	acx_mpi_ok="yes"
	AC_SUBST(MPI_LIBS)
	AC_SUBST(MPI_INCS)
	AC_DEFINE(MPI,1,[Define if you want to compile for MPI parallel operation])
	AC_DEFINE(HAVE_MPI_H,1,[Define if mpi.h is available])
else
	AC_MSG_WARN([Unable to find an MPI installation])
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$acx_mpi_ok" = "xyes"; then
        #ifelse([$1],,AC_DEFINE(HAVE_MPI,1,[Define if you have MPI library.]),[$1])
	#:
	AC_DEFINE(HAVE_MPI,1,[Define if you have MPI library.])
        $1
else
        acx_mpi_ok=no
        $2
fi
])dnl ACX_MPI
