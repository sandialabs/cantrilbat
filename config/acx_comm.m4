dnl @synopsis ACX_SET_COMM()
dnl
dnl This macro sets variables describing the preferred communication
dnl protocol.  Specificially, it sets SERIAL or PARALLEL.  For PARALLEL
dnl it also checks for a set of MPI libraries and mpi.h and determines
dnl the MPI implementation (MPICH or LAMMPI).
dnl
dnl Use --with-comm=<comm> to specify.  SERIAL is the default.
dnl 
dnl See also ACX_SET_MPICH, ACX_FIND_MPICH, ACX_SET_LAMMPI, ACX_FIND_LAMMPI.
dnl
dnl Use the environment variable COMM as an alternative to --with-comm.
dnl
dnl Command line options:
dnl   --with-comm=[LAMMPI,MPICH,MPI,PARALLEL,SERIAL]  default=SERIAL
dnl 
dnl Environment variables:
dnl   COMM=[LAMMPI,MPICH,MPI,PARALLEL,SERIAL]
dnl
dnl @version $Id: acx_comm.m4,v 1.1 2010/01/12 22:35:36 hkmoffa Exp $
dnl @author Patrick Notz <pknotz@sandia.gov>
dnl
AC_DEFUN([ACX_COMM], [

acx_comm_specified="no"
COMM="SERIAL"
COMM_TAG="_serial"

AC_ARG_WITH([comm],
	    [AC_HELP_STRING([--with-comm=PROTOCOL], 
	                    [use specified communication protocol])],
	    [acx_comm_specified="yes"
	     AC_MSG_CHECKING([for parallel or serial operation])
	     case $withval in
		serial | Serial | SERIAL )
		    COMM="SERIAL"
		    COMM_TAG="_serial"
		;;
		mpi | MPI | parallel | Parallel | PARALLEL)
			COMM=PARALLEL
			COMM_TAG="_mpi"
		;;
		*)
			AC_MSG_ERROR([serial or mpi are valid])	
		;;
	     esac
	     AC_MSG_RESULT([$COMM])],
	     [acx_comm_specified="no"])
			
AC_SUBST(COMM_TAG)

if test x"$COMM" = x"SERIAL" ; then
    AC_DEFINE(SERIAL,1,[Define if compiled for SERIAL execution])
else
    AC_DEFINE(PARALLEL,1,[Define if compiled for PARALLEL execution])
    AC_DEFINE(HAVE_MPI_H,1,[Define if mpi.h is available])
    AC_DEFINE(HAVE_MPI,1,[Define if you have MPI library.])
    AC_DEFINE(MPI,1,[Define if you want to compile for MPI parallel operation])
fi
])dnl ACX_COMM
