dnl @synopsis ACX_MPI_HOME([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a SEAMS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_MPI_HOME.
dnl
dnl @version $Id: acx_mpi_home.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Patrick Notz <pknotz@sandia.gov>
dnl
AC_DEFUN([ACX_MPI_HOME], [

acx_mpi_home_ok=no

mpi_dirs="/usr/local/mpich /usr/local/mpi /usr/mpich /usr/mpi /usr /usr/local"

AC_ARG_WITH(mpi,
	[AC_HELP_STRING([--with-mpi=<path>], [use MPI installation at <path>])])
case $with_mpi in
	yes | "") MPI_HOME="";;
	no) acx_mpi_home_ok=disable ;;
	*) MPI_HOME="$with_mpi" ;;
esac

# Check for the main MPI header file
check_file="include/mpi.h"

# Check the MPI_INCS environment variable
if test "$acx_mpi_home_ok" = "no" ; then
    if test "x$MPI_HOME" != "x"; then
	AC_MSG_CHECKING([for MPI in $MPI_HOME])
	if test -r "$MPI_HOME/$check_file" ; then
	    acx_mpi_home_ok="yes"
	else
	    MPI_HOME=""
	fi
	AC_MSG_RESULT([$acx_mpi_home_ok])
    fi
fi

# Next, check our list of common locations
if test "$acx_mpi_home_ok" = "no" ; then
    for d in $mpi_dirs ; do
	AC_MSG_CHECKING([for a MPI home in $d])
	if test -r "$d/$check_file" ; then
	    acx_mpi_home_ok="yes"
	    MPI_HOME="$d"
	else
	    acx_mpi_home_ok="no"
	    MPI_HOME=""
	fi
	AC_MSG_RESULT([$acx_mpi_home_ok])
	if test "$acx_mpi_home_ok" = "yes" ; then break; fi
    done
fi

AC_SUBST(MPI_HOME)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test "x$acx_mpi_home_ok" = "xyes"; then
	AC_MSG_NOTICE([Using MPI installation under $MPI_HOME])
	$1
else
	AC_MSG_NOTICE([Unable to find an MPI installation])
	MPI_HOME=""
        acx_mpi_home_ok=no
        $2
fi
])dnl ACX_MPI_HOME
