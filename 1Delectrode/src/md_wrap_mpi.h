/**
 *  @file md_wrap_mpi.h  Low level mpi wrappers
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef _MD_WRAP_MPI_H
#define _MD_WRAP_MPI_H

#include "m1d_defs.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_MPI
//! MPI object that specifies a request type
#define M1D_MPI_Request MPI_Request
//! MPI object that specifies a communication type
#define M1D_MPI_Comm    MPI_Comm
#else
#define M1D_MPI_Request int
#define M1D_MPI_Comm    int
#endif
//----------------------------------------------------------------------------------------------------------------------------------
namespace M1D
{
//==================================================================================================================================
//! Blocking read to the current processor from a source processor specifying specifically the source processor 
//! and typeid of the message that is to be expected
/*!
 *  Sends bytes bytes of data from the address buf to the destination processor, dest, which is the calling processor of this
 *  function.  The calling processor will hang on this send until the message is received as specified. 
 *
 *  @param[out]              buf                 Location of the data buffer to receive the message
 *  @param[in]               bytes               Number of bytes of data to be read
 *  @param[in]               source              Source Processor. This is a positive number specifying the source processor
 *  @param[in]               type                type ID of the sending process.  This is a positive number
 *                                               specifying the type 
 *
 *  @return                                      The number of bytes read
 */
int md_read_specific(void* const buf, int bytes, const int source, const int type);

//==================================================================================================================================
//! Blocking read to the current processor from a source processor
/*!
 *  Sends bytes bytes  of data from the address buf to the destination processor, dest.
 *  The processor will hang on this send until it is finished.
 *
 *  @param[out]              buf                 Location of the data
 *  @param[in]               bytes               Number of bytes of data to be read
 *  @param[in,out]           source              Source Processor. Either this is a positive number specifying the source processor
 *                                               or it can be equal to -1, in which case MPI_ANY_SOURCE is used.
 *                                               On return, source is equal to the procID of the sending processor.
 *  @param[in,out]           type                type ID of the sending process.  Either this is a positive number
 *                                               specifying the type or it can be equal to -1, in which case MPI_ANY_TAG is used.
 *                                               On return, type is equal to the type ID of the message.
 *
 *  @return                                      The number of bytes read
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

//! Machine dependent wrapped message-reading communication routine for MPI. This call does not block.
/*!
 *
 *  @param[in,out]           buf                 Location of the data
 *  @param[in]               bytes               Number of bytes of data to be expected to be received
 *  @param[in,out]           source              Source Processor. Either this is a positive number specifying the source processor
 *                                               or it can be equal to -1, in which case MPI_ANY_SOURCE is used.
 *                                               On return, source is equal to the procID of the sending processor.
 *  @param[in,out]           type                type ID of the sending process.  Either this is a positive number
 *                                               specifying the type or it can be equal to -1, in which case MPI_ANY_TAG is used.
 *                                               On return, type is equal to the type ID of the message.
 *
 *  @param[out]              request             Pointer to the MPI_Request object that is initialized and that can be referred
 *                                               to later as an identifier for this MPI_Irecv() request
 *
 *  @return                                      returns the error code produced by MPI_Irecv(). 0 indicates no error
 */
int md_wrap_iread(void* const buf, int bytes, int* const source, int* const type, M1D_MPI_Request* const request);

//==================================================================================================================================
//! Machine dependent wrapped message-reading communication routine for MPI specifying specifically the source processor 
//! and typeid of the message that is to be expected. This call does not block.
/*!
 *
 *  @param[in,out]           buf                 Location of the data
 *  @param[in]               bytes               Number of bytes of data to be expected to be received
 *  @param[in]               sourceC             Source Processor. This is a positive number specifying the source processor
 *  @param[in]               typeC               type ID of the sending process.  This is a positive number that must be matched.
 *
 *  @param[out]              request             Pointer to the MPI_Request object that is initialized and that can be referred
 *                                               to later as an identifier for this MPI_Irecv() request
 *
 *  @return                                      returns the error code produced by MPI_Irecv(). 0 indicates no error
 */
int md_wrap_iread_specific(void* const buf, const int bytes, const int sourceC, const int typeC, M1D_MPI_Request* const request);

//==================================================================================================================================
//! Blocking send from the current processor to a destination processor
/*!
 *  Sends bytes  of data from the address buf to the destination processor, dest. 
 *  The processor will hang on this send until it is finished.
 *
 *  @param[in]               buf                 Location of the data to be sent
 *  @param[in]               bytes               Number of bytes of data to be sent. Can be zero.
 *  @param[in]               dest                Destination processor
 *  @param[in]               type                type ID of the sending process.
 *
 *  @return                                      returns an error code from MPI_Send. Should be 0 for a successful write.
 */
int md_wrap_write(void* const buf, int bytes, int dest, int type);

//==================================================================================================================================
//! Blocking wait for the current processor based on a request object. This is a wrapper around MPI_Wait().
/*!
 *  @param[out]               source             Returns the processor ID of the received message
 *  @param[out]               type               Returns the type ID of the received message.
 *  @param[in]                request            Pointer to the request object to use in the MPI_Waite() call
 *
 *  @return                                      returns the number of bytes that was received in the message
 */
int md_wrap_wait(int* const source, int* const type, M1D_MPI_Request* const request);

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
