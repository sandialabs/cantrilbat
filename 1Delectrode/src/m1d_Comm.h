/**
 * @file m1d_Comm.h
 *   Communication routines that are part of the m1d Package
 */

/*
 *  $Id: m1d_Comm.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef _M1D_COMM_H
#define _M1D_COMM_H

#include "Epetra_ConfigDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#define M1D_MSG_TYPE_BASE 1056
#define M1D_MSG_TYPE_RANGE 20

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Comm.h"
#include "Epetra_VbrMatrix.h"

#include "m1d_defs.h"
#include "md_wrap_mpi.h"

#include <ostream>
#include <iomanip>
#include <cstdio>

namespace m1d
{

//! Returns the processor which has the maximum value of the first parameter
/*!
 *  This routine will figure out a consistent way to handle ties.
 *
 * @param pmax       local processor contribution
 * @param cc         Epetra Communicator
 * @param gmax       value of the maximum across all processors.
 * @return Global processor number
 */
int
procChoice_Max(double pmax, const Epetra_Comm& cc_ptr, double &gmax);

//! Returns the processor which has the minimum value of the first parameter
/*!
 *  This routine will figure out a consistent way to handle ties.
 *
 * @param pmin       local processor contribution
 * @param cc         Epetra Communicator
 * @param gmin       value of the maximum across all processors.
 * @return Global processor number
 */
int
procChoice_Min(const double pmin, const Epetra_Comm& cc_ptr, double &gmin);


//! Returns the processor which has the maximum value of the first parameter
//! and broadcasts an int from the winning processor.
/*!
 *  This routine will figure out a consistent way to handle ties.
 *
 * @param pmax       local processor contribution
 * @param cc         Epetra Communicator
 * @param gmax       value of the maximum across all processors.
 * @param pint       On entry it contains an int. On exit it returns the int
 *                   from the processor that won.
 * @return Global processor number
 */
int
procChoice_Max_Brcst1Int(double pmax, const Epetra_Comm& cc, double &gmax, int &pint);


//! Returns the processor which has the minimum value of the first parameter
//! and broadcasts an int from the winning processor.
/*!
 *  This routine will figure out a consistent way to handle ties.
 *
 * @param pmin       local processor contribution
 * @param cc         Epetra Communicator
 * @param gmin       value of the maximum across all processors.
 * @param pint       On entry it contains an int. On exit it returns the int
 *                   from the processor that won.
 * @return Global processor number
 */
int
procChoice_Min_Brcst1Int(const double pmin, const Epetra_Comm& cc, double &gmin, int &pint);


//! Print blocks
/*!
 *  Routine to allow IO between print_sync_start and print_sync_end to be printed
 * by each processor entirely before the next processor begins its IO.  The
 *printing sequence is from proc = 0 to the last processor,
 *number_of_procs = nprocs - 1.
 *
 *NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.
 *
 *
 *=======
 *
 *Return code:     void
 *============
 *
 *Parameter list:
 *===============
 *
 *proc:            Current processor number.
 *
 *do_print_line:   Boolean variable.  If true, a line of # is printed to
 *indicate the start of a print_sync I/O block.
 *
 */
void
print_sync_start(int do_print_line, int proc_config[]);

//!
/*!
 Routine to allow IO between print_sync_start and print_sync_end to be printed
 by each processor entirely before the next processor begins its IO. The
 printing sequence is from proc = 0 to the last processor,
 number_of_procs = nprocs - 1.

 NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

 Author:          John N. Shadid, SNL, 1421
 =======

 Return code:     void
 ============

 Parameter list:
 ===============

 proc:            Current processor number.

 nprocs:          Number of processors in the current machine configuration.

 do_print_line:   Boolean variable.  If true, a line of # is printed to
 indicate the start of a print_sync I/O block.

 */
void
print_sync_end(int do_print_line, Epetra_Comm &comm);

//==============================================================================
//! This is a small class that inherits from ostringstream that allows for
//! smoother handling of printing Epetra Objects from proc 0
/*!
 *  Suggested functionality:
 *      make it a singelton
 *      override endl-> tried to do this once but it didn't work
 *
 *   The way this class works is to first port all information to processor 0.
 *   Then, it will get printed from proc 0.
 *
 *   The class prints to a FILE stream, specifically. At initialization,
 *   the file stream is assigned. The default file stream is stdout. However,
 *   it can be changed to any open FILE *. The class does not  handle
 *   open or close operations on the FILE stream.
 *
 */
class stream0 : public std::ostringstream {
public:
  //! Default constructor
  stream0(FILE * ff = stdout);

  //! Override the flush() function
  /*!
   *  Flushes the current output to the file or stream
   */
  stream0 & flush();

  //! Returns a file pointer where the data is going to be written to.
  /*!
   * 
   */
  FILE *outFile() const {
    return outFILE_;
  }

  //! General print function that works just like the C function
  //! printf()
  /*!
   *
   * @param format  Free form formatting string that follows the same
   *                rules as printf()
   */
  void
  print0(const char *format, ...);

  int fprintf0only(const char *format, ...);

  void drawline0(int indentLen, int num, char linChar = '-');


  //! Counter keeping track of how many buffers are sent to processor 0
  int nBufsSent;


  //! Output file designator
  /*!
   *  defaults to stdout 
   */
  FILE *outFILE_;
};

//==============================================================================
//! Start a print0_sync block for neat printing of mpi jobs that are
//! being run on a network of workstations
/*!
 * The print0_sync block seems to get around a very sticky situation
 * that occurs when you are printing on a network of workstations.
 * By direct trial, we have found that the only way to get ordered
 * output from stderr or stdout is to print everything from processor 0.
 * We think this is due to networking issues outside of the control of mpi.
 * There is additional buffering of stderr and stdout that occurs on
 * networks of workstations and not a dedicated mp machine.
 *
 * Note between print0_sync blocks no other mpi operations should be
 * attempted.
 *
 * Between the print0_sync blocks, only certain printing-to-proc0 operations
 * should be done. Only these should be done, and operations that
 * are derivatives of these operations:
 *
 *      sprint0(const char *ibuf);
 *      ssprint0(stream0 &ss);
 *      stream0 ss;       ss.print0();
 *
 *
 * @param do_print_line    Print a line
 * @param ss               stream0 being used for printing.
 *                         Note, I don't think this is really needed.
 * @param comm             reference to the Epetra_Comm object
 */
void
print0_sync_start(int do_print_line, stream0 &ss, const Epetra_Comm &comm);

//! End of a print0_sync block for neat printing of mpi jobs that are
//! being run on a network of workstations
/*!
 *
 * @param do_print_line    Print a line
 * @param ss               stream0 being used for printing.
 *                         Note, I don't think this is really needed.
 * @param comm             reference to the Epetra_Comm object
 */
void
print0_sync_end(int do_print_line, stream0 &ss, const Epetra_Comm &comm);

//! Print from processor 0, while in between the print0_sync blocks
/*!
 *   This will cause a hang condition if this doesn't occur between
 *   print0_sync blocks
 *
 * @param ibuf     Buffer that needs to be printed. Must be null terminated
 *                 like any c string.
 * @param ff       FILE stream to print to. This defaults to stdout
 * @return         Returns the number of messages sent
 */
int
sprint0(const char *ibuf, FILE * ff = stdout);

//! Print stream0 information from processor 0, while in between the print0_sync blocks
/*!
 *   This will cause a hang condition if this doesn't occur between
 *   print0_sync blocks. The stream0 object is zeroed at the end
 *   of the operation.
 *
 * @param ss    stream0 object that will be printed. The stream0 object
 *              is zeroed at the end of the operation.
 */
void
ssprint0(stream0 &ss);

//! Print an Epetra_MultiVector using the stream0 process
/*!
 * @param oss   stream0 stream
 * @param v     Epetra_MultiVector reference
 */
void
Print0_epMultiVector(stream0 &oss, const Epetra_MultiVector &v);

void
print0_epIntVector(const Epetra_IntVector &v, const char *label, FILE *of = stdout);


void
Print0_epIntVector(stream0 &os, const Epetra_IntVector &v);

//! Prints a Multivector between print0_sync blocks
/*!
 * The output is dumped onto stdout on processor 0
 *
 * @param v  Epetra MultiVector to be printed
 * @param label optional argument will create an extra line with text
 */
void
print0_epMultiVector(const Epetra_MultiVector &v, const char *label = 0, FILE *of = stdout);

//! Print an Epetra_VbrMatrix using the stream0 process
/*!
 * @param oss   stream0 stream
 * @param v     Epetra_VbrMatrix reference
 */
void
Print0_epVbrMatrix(stream0 &os, const Epetra_VbrMatrix &mat);

//! Prints a VbrMatrix between print0_sync blocks
/*!
 * The output is dumped onto stdout on processor 0
 * @param mat  VbrMatrix to be printed
 * @param label optional argument will create an extra line with text
 */
void
print0_epVbrMatrix(const Epetra_VbrMatrix &mat, const char *label = 0, std::FILE *of = stdout);

//! Prints an Epetra_BlockMap using the stream0 process
/*!
 * @param oss   stream0 stream
 * @param v     Epetra_BlockMap reference
 */
void
Print0_epBlockMap(stream0 &os, const Epetra_BlockMap &map);

//! Prints a map between print0_sync blocks
/*!
 * The output is dumped onto stdout on processor 0
 * @param map  BlockMap to be printed
 * @param label optional argument will create an extra line with text
 */
void
print0_epBlockMap(const Epetra_BlockMap &map, const char *label = 0, FILE *of = stdout);

//! Prints an Epetra_MultiVector using the stream0 process
/*!
 * @param oss          stream0 stream
 * @param v            Epetra_MultiVector reference
 * @return             Returns a changeable reference to the input stream0
 *                     object
 */
stream0 & operator<<(stream0 &oss, const Epetra_MultiVector &v);



}

#endif
