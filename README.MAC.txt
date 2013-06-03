===============================================================================================================================

                      Directions for Cantera_apps for the mac

                      Harry Moffat                                                                                  10/1/2012
                      John Hewson 


===============================================================================================================================

Directions for Compiling
----------------------------------------------------------------------------------- ----

 1) Set up the autoconf environment
    autocconf  (makes configure from configure.in)

 2) Setup the makefile environment. Type:

    mac_sierra-gcc-mp-4.4

    This sets up environmental variables that then runs preconfig. preconfig will then run configure.
    configure will create all of the Makefile files from Makefile.in inputs.

 3) Then, type 

    make

    You do not need to install anything in the cantera_apps environment. You will run all the 
    code in place.

 4) There is a small test suite located in cantera_apps_trunk itself. The main test suite
    is located in cantera_apps_test. To run the test suite:

    make test

    Some of the tests use ascii file comparisons only. Some of the tests do an rtol/atol comparison
    against csv or xml solutions containing double data.  All of the tests use some form of
    asci comparison. Therefore, diffs are to be expected.

    Comparisons within the mac environment are difficult, because all of the blessed files are
    created on the linux environment. And, there is numerical roundoff created on the mac
    from all calculations involving a matrix inversion or solution of a nonlinear system.
   


Compilers
----------------------------------------------------------------------------------------

    The sierra-devel compilers were chosen. This is based on gcc v. 4.4 and open mpi 1.4.2

The compilers without the mpi wrapper are located on the lan in the directories:

          /opt/local/var/macports/software/gcc44/4.4.4_1+darwin_10/opt/local/bin/g++-mp-4.4

The compiler names were labeled differently than their defaults. This I think caused problems with
trying to link in Cantera's python libraries. 


The compilers with the mpi wrapper are located in :

        /projects/openmpi/1.4.2/gcc-4.4/bin


The fortran environment

    The fortran program is gfortran-mp-4.4.

C linkage is obtained by lower casing the externals and prepending and suffixing with a single underscore.
I also used the "-g -fno-second-underscore" compiler flags to ensure that a second underscore wasn't used.


   
Directories
------------------------------------------

  CANTERA_ROOT          = Dierctory that Cantera, itself, compiles 
  CANTERA_INSTALL_DIR   = Directory which Cantera installs itself into

  CANTERA_APPS_ROOT     = Directory that cantera_apps compiles in and where the root of the source code is


  CANTERA_APPS_TEST     = Directory taht cantera_apps_test compiles in and where the root of the test suite is

  TRILINOS_INSTALL_DIR  = Installation directory for Trilinos. Yes, we are using v 9.0 in this project. and plan
                          to go to the new version soon.




Environmental Variables
------------------------------------------


CANTERA_BIN  = Binary directories for some Cantera utilities
CANTERA_DATA = Data file directory for Cantera

MP_MACHINEFILE         = file location for the machine file to use in mpi environments
                         The scripts attempt to use the current machine for the contents of this machine file.
                         The name of the machine is listed 8 times in that file.

MP_MACHINEFILE_DEFAULT = default location for the machine file 
                         The scripts attempt to use the current machine for the contents of this machine file.


Linking Environment
----------------------------------------------------------------------------------------------------------


The mac openmpi environment required the following flags as end libraries on C++ programs.

       -lmpi_cxx -lopen-rte -lopen-pal -lgfortran

This has not been rigorously determined, other than the mpi_cxx library.
Note, using mpi_cxx on non mpi programs, i.e., ones that are linked with the gcc-mp-4.4 compiler caused unsatisfied
externals to appear for some reason. Again, I don't really have an idea why.


We have found it necessary to use Cantera ctf2c library during the mac builds involving f77 code.
Therefore, we add:

	-lctf2c  -lgfortran

to the link line. I don't know of a way around that yet.


Python Interface
--------------------------------------------------------------------------------------------------------

I made an attempt to get numpy working with python 2.7. There were outstanding issues that arose during the linking
within Cantera, and I ended up punting on this. I compiled Cantera in a minimal python environment.






















