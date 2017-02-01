/**
 *  @file md_wrap_mpi.h  Low level mpi wrappers
 *
 */

/*
 *  $Id: md_wrap_mpi.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef _MD_WRAP_MPI_H
#define _MD_WRAP_MPI_H

#include "m1d_defs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_MPI
#define M1D_MPI_Request MPI_Request
#define M1D_MPI_Comm    MPI_Comm
#else
#define M1D_MPI_Request int
#define M1D_MPI_Comm    int
#endif
//----------------------------------------------------------------------------------------------------------------------------------
namespace M1D
{
//==================================================================================================================================
int md_read_specific(void * const buf, int bytes, const int source, const int type);

//==================================================================================================================================
//! Blocking read to the current processor from a source processor
/*!
 *  Sends bytes bytes  of data from the address buf to the destination processor, dest.
 *  The processor will hang on this send until it is finished.
 *
 *  @param buf     Location of the data
 *  @param bytes   Number of bytes of data to be read
 *  @param source  Source Processor. Either this is a positive number
 *                specifying the source processor or it can be equal
 *                to -1, in which case MPI_ANY_SOURCE is used.
 *                On return, source is equal to the procID of the
 *                sending processor.
 *  @param type    type ID of the sending process.  Either this is a positive number
 *                specifying the type or it can be equal
 *                to -1, in which case MPI_ANY_TAG is used.
 *                On return, type is equal to the type ID of the message.
 *
 *  @return       The number of bytes read
 */
int md_read(void* const buf, int bytes, int* const source, int* const type);

//==================================================================================================================================
//! Blocking send from the current processor to a destination processor
/*!
 *  Sends bytes bytes  of data from the address buf to the destination processor, dest. 
 *  The processor will hang on this send until it is finished.
 *
 *  @param[in,out]           buf                 Location of the data
 *  @param[in]               bytes               Number of bytes of data to be sent
 *  @param[in]               dest                Destination processor
 *  @param[in]               type                type ID of the sending process.
 *
 *  @return                                      returns the number of bytes sent
 */
int md_write(void* const buf, const int bytes, const int dest, const int type);

//==================================================================================================================================
/*
 *
 Machine dependent wrapped message-reading communication routine for MPI.
 This call does not block


 Return code:     int
 ============

 Parameter list:
 ===============

 buf:             Beginning address of data to be sent.
 bytes:           Length of message in bytes.
 source:          Source processor number.
 type:            Message type
 */
int
md_wrap_iread(void *buf,
              int bytes,
              int *source,
              int *type,
              M1D_MPI_Request *request);

//==================================================================================================================================
int
md_wrap_iread_specific(void *buf,
                       const int bytes,
                       const int sourceC,
                       const int typeC,
                       M1D_MPI_Request *request);

//==================================================================================================================================
int
md_wrap_write(void *buf, int bytes, int dest, int type, int *flag);

//==================================================================================================================================
int
md_wrap_wait(int * const source, int * const type, M1D_MPI_Request * const request);

//==================================================================================================================================

int md_wrap_wait_specific(const int sourceC, const int typeC, M1D_MPI_Request * const request);

//==================================================================================================================================
int
md_wrap_iwrite(void * const buf,
               const int bytes,
               const int dest,
               const int type,
               M1D_MPI_Request * const request);

//==================================================================================================================================
void
md_parallel_info(int *proc, int *nprocs, int *dim, M1D_MPI_Comm comm);

//==================================================================================================================================
int
md_mpi_iread(void *buf,
             int bytes,
             int *source,
             int *type,
             M1D_MPI_Request *request,
             int *icomm);

//==================================================================================================================================
int
md_mpi_write(void *buf, int bytes, int dest, int type, int *flag, int *icomm);

//==================================================================================================================================
int
md_mpi_wait(int * const source, int * const type, M1D_MPI_Request * const request);

//==================================================================================================================================
int
md_mpi_iwrite(void *buf,
              int bytes,
              int dest,
              int type,
              M1D_MPI_Request *request,
              int *icomm);

//==================================================================================================================================
int
md_wrap_request_cancel(M1D_MPI_Request *request);

//==================================================================================================================================
void
md_throw(const char * const rout, const int err);

//==================================================================================================================================
int
md_waitany(int * const indexCompleted,
           int * const source,
           int * const type,
           int requestVlen,
           M1D_MPI_Request * const requestVector);
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
