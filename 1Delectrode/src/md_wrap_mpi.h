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
 *  @param[in]                request            Pointer to the request object to use in the MPI_Wait() call
 *
 *  @return                                      returns the number of bytes that was received in the message
 */
int md_wrap_wait(int* const source, int* const type, M1D_MPI_Request* const request);

//==================================================================================================================================
//! Blocking wait for a specific nonblocking read/write for the current processor based on a request object. This is a wrapper around MPI_Wait().
/*!
 *  If the expected processor ID or Type ID is not what is actually received a print message is sent to stdout, but otherwise
 *  no action is taken.
 *
 *  @param[in]                sourceC            The processor ID of the received message to be expected
 *  @param[in]                typeC              The type ID of the received message to be expected
 *  @param[in]                request            Pointer to the request object to use in the MPI_Wait() call
 *
 *  @return                                      returns the number of bytes that was received in the message
 */
int md_wrap_wait_specific(const int sourceC, const int typeC, M1D_MPI_Request* const request);

//==================================================================================================================================
//! Non-blocking send from the current processor to a destination processor
/*!
 *  Sends bytes  of data from the address buf to the destination processor, dest. 
 *  The processor will return from this send before the send is finished.
 *
 *  @param[in]               buf                 Location of the data to be sent
 *  @param[in]               bytes               Number of bytes of data to be sent. Can be zero.
 *  @param[in]               dest                Destination processor number
 *  @param[in]               type                type ID of the sending process.
 *  @param[in,out]           request             Pointer to the request object to use in the MPI_Wait() call
 *
 *  @return                                      returns an error code from MPI_Isend. Should be 0 for a successful write.
 */
int md_wrap_iwrite(void* const buf, const int bytes, const int dest, const int type, M1D_MPI_Request * const request);

//==================================================================================================================================
//! Get the basic information about procID and number of processors
/*!
 *  @param[out]              proc                Returns the processor ID
 *  @param[out]              nprocs              Returns the number of procs
 *  @param[out]              dim                 Returns the dimensionality of configuration. Returns 1
 *  @param[in]               comm                Reference to the MPI_Comm object
 */
void md_parallel_info(int* const proc, int* const nprocs, int* const dim, M1D_MPI_Comm& comm);

//==================================================================================================================================
//! Machine dependent wrapped request object deletion routine.
/*!
 *
 *  @param[in]                 request           Pointer to an existing request object that will be freed.
 *
 *  @return                                      Returns 0 if successful
 */
int md_wrap_request_cancel(M1D_MPI_Request* const request);

//==================================================================================================================================
//!  Throw an MPI error when an MPI error is encountered during the routines
/*!
 *   The error object, MPI::Exception, is thrown from this routine. when any of the routines encounters an error
 *
 *  @param[in]               rout                Character message describing the location of the error
 *  @param[in]               err                 Integer error code from the MPI layer.
 */
void md_throw(const char* const rout, const int err);

//==================================================================================================================================
//! Blocking wait for the current processor based on multiple request objects. This is a wrapper around MPI_Wait().
/*!
 *  This routine waits on a vector of request objects. It will return after any of the objects are completed.
 *
 *  @param[out]               indexCompleted     Index ID of the request object that is completed.
 *  @param[out]               source             Returns the processor ID of the received message
 *  @param[out]               type               Returns the type ID of the received message.
 *  @param[in]                requestVlen        Length of the request vector to 
 *  @param[in]                requestVector      Pointer to a vector of request objects to use in the MPI_Wait() call
 *
 *  @return                                      Returns the number of bytes that was received in the message
 */
int md_waitany(int* const indexCompleted, int* const source, int* const type, int requestVlen,
               M1D_MPI_Request* const requestVector);
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
