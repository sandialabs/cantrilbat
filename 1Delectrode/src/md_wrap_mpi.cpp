/**
 *  @file  md_wrap_mpi.cpp
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "md_wrap_mpi.h"

#include <cstdlib>
#include <cstdio>
//----------------------------------------------------------------------------------------------------------------------------------
namespace M1D
{

//! Global buffer for reading a single byte
static int gl_rbuf = 3;

//! Global buffer for sending a single byte
static int gl_sbuf = 3;

//! the processor name
static int the_proc_name = -1;

//==================================================================================================================================
int md_read_specific(void* const buf, int bytes, const int sourceC, const int typeC)
{
#ifdef HAVE_MP 
    int err, buffer = 1;
    MPI_Status status;
    int type = typeC;
    int source = sourceC;
    if (type == -1) {
      type = MPI_ANY_TAG;
    }
    if (source == -1) {
      source = MPI_ANY_SOURCE;
    }

    if (bytes == 0) {
        err = MPI_Recv(&gl_rbuf, 1, MPI_BYTE, source, type, MPI_COMM_WORLD, &status);
    } else {
        err = MPI_Recv(buf, bytes, MPI_BYTE, source, type, MPI_COMM_WORLD, &status);
    }
    if (err != 0) {
        md_throw("MPI_Recv() ERROR from md_read", err);
    }
    MPI_Get_count(&status, MPI_BYTE, &buffer);
    source = status.MPI_SOURCE;
    type = status.MPI_TAG;
    if (typeC != MPI_ANY_TAG && typeC != -1) {
        if (type != typeC) {
            md_throw("md_read_specific() type confusion: ", int2str(type) + " " + int2str(typeC), -1);
        }
    }
    if (sourceC != MPI_ANY_SOURCE && sourceC != -1) {
        if (source != sourceC) {
            md_throw("md_read_specific() source confusion: " + int2str(source) + " " + int2str(sourceC), -1);
        }
    }
    if (bytes != 0) {
        bytes = buffer;
    }
#endif
    return bytes;
}
//==================================================================================================================================
int md_read(void* const buf, int bytes, int* const source, int* const type)
{
#ifdef HAVE_MPI
    int err, buffer = 1;
    MPI_Status status;
    if (*type == -1)
        *type = MPI_ANY_TAG;
    if (*source == -1)
        *source = MPI_ANY_SOURCE;

    if (bytes == 0) {
        err = MPI_Recv(&gl_rbuf, 1, MPI_BYTE, *source, *type, MPI_COMM_WORLD, &status);
    } else {
        err = MPI_Recv(buf, bytes, MPI_BYTE, *source, *type, MPI_COMM_WORLD, &status);
    }
    if (err) {
        md_throw("MPI_Recv ERROR from md_read", err);
    }
    MPI_Get_count(&status, MPI_BYTE, &buffer);
    *source = status.MPI_SOURCE;
    *type = status.MPI_TAG;
    if (bytes != 0)
        bytes = buffer;
#else
   *source = 0;
   *type = 0;
#endif
    return bytes;
}
//==================================================================================================================================
int md_write(void* const buf, const int bytes, const int dest, const int type)
{
#ifdef HAVE_MPI
    int err;
    if (bytes == 0) {
        err = MPI_Send(&gl_sbuf, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD);
    } else {
        err = MPI_Send(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD);
    }
    if (err) {
        md_throw("MPI_Send ERROR from md_write", err);
    }
#endif
    return bytes;
}
//==================================================================================================================================
void md_throw(const char* const rout, const int err)
{
#ifdef HAVE_MPI
    if (err == MPI_ERR_BUFFER) {
        fprintf(stderr, "ERROR %s: MPI_ERR_BUFFER\n", rout);
        exit(-1);
        //throw std::exception(err);
        //throw MPI::Exception(err);
    }
#else
    fprintf(stderr, " MPI error occurred when mpi not compiled in");
    exit(-1);
#endif
}
//==================================================================================================================================
int md_wrap_iread(void* const buf, int bytes, int* const source, int* const type, M1D_MPI_Request* const request)
{
#ifdef HAVE_MPI
    int err;
    if (*type == -1) {
        *type = MPI_ANY_TAG;
    }
    if (*source == -1) {
        *source = MPI_ANY_SOURCE;
    }
    if (bytes == 0) {
        err = MPI_Irecv(&gl_rbuf, 1, MPI_BYTE, *source, *type, MPI_COMM_WORLD, request);
    } else {
        err = MPI_Irecv(buf, bytes, MPI_BYTE, *source, *type, MPI_COMM_WORLD, request);
    }
    return err;
#else
    return 0;
#endif
}
//==================================================================================================================================
int md_wrap_iread_specific(void *buf, int bytes, const int sourceC, const int typeC, M1D_MPI_Request *request)
{
#ifdef HAVE_MPI
    int err = 0;
    if (bytes == 0) {
        err = MPI_Irecv(&gl_rbuf, 1, MPI_BYTE, sourceC, typeC, MPI_COMM_WORLD, request);
    } else {
        err = MPI_Irecv(buf, bytes, MPI_BYTE, sourceC, typeC, MPI_COMM_WORLD, request);
    }
    return err;
#else
    return 0;
#endif
}
//==================================================================================================================================
int md_wrap_write(void* const buf, int bytes, int dest, int type)
{
#ifdef HAVE_MPI
    int err = 0;
    if (bytes == 0) {
        err = MPI_Send(&gl_sbuf, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD);
    } else {
        err = MPI_Send(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD);
    }
    return err;
#else
    return 0;
#endif
} 
//==================================================================================================================================
int md_wrap_wait(int* const source, int* const type, M1D_MPI_Request* const request)
{
#ifdef HAVE_MPI
    int count;
    int err;
    MPI_Status status;
    if ((err = MPI_Wait(request, &status))) {
        md_throw("MPI_Wait ERROR from md_wrap_wait", err);
    }
    MPI_Get_count(&status, MPI_BYTE, &count);
    *source = status.MPI_SOURCE;
    *type = status.MPI_TAG;
    /* return the count, which is in bytes */
    return count;
#else
    return 0;
#endif
}
//==================================================================================================================================
int md_wrap_wait_specific(const int sourceC, const int typeC, M1D_MPI_Request* const request)
{
#ifdef HAVE_MPI
    int err;
    int source = sourceC;
    int type = typeC;
    int count;
    MPI_Status status;
    if ((err = MPI_Wait(request, &status))) {
         md_throw("MPI_Wait ERROR from md_wrap_wait_specific()", err);
    }
    MPI_Get_count(&status, MPI_BYTE, &count);
    source = status.MPI_SOURCE;
    type = status.MPI_TAG;
    if (source != sourceC) {
        printf("md_wrap_wait_specific() error, sources differ: %d %d\n", source, sourceC);
    }
    if (type != typeC) {
        printf("md_wrap_wait_specific() error, types differ: %d %d\n", type, typeC);
    }
    /* return the count, which is in bytes */
    return count;
#else
    return 0;
#endif
}
//==================================================================================================================================
int md_wrap_iwrite(void* const buf, const int bytes, const int dest, const int type, M1D_MPI_Request* const request)
{
#ifdef HAVE_MPI
    int err = 0;
    if (bytes == 0) {
        err = MPI_Isend(&gl_sbuf, 1, MPI_BYTE, dest, type, MPI_COMM_WORLD, request);
    } else {
        err = MPI_Isend(buf, bytes, MPI_BYTE, dest, type, MPI_COMM_WORLD, request);
    }
    return err;
#else
    return 0;
#endif
}
//==================================================================================================================================
void md_parallel_info(int* const proc, int* const nprocs, int* const dim, M1D_MPI_Comm& comm)
{
#ifdef HAVE_MPI
    MPI_Comm_size(comm, nprocs);
    MPI_Comm_rank(comm, proc);
    *dim = 1;
    the_proc_name = *proc;
#else
    *nprocs = 1;
    *dim = 1;
    *proc = 0;
#endif
} 
//==================================================================================================================================
int md_wrap_request_cancel(M1D_MPI_Request* const request)
{
#ifdef HAVE_MPI
    int err = 0;
    if (request != nullptr) {
        err = MPI_Cancel(request);
    }
    err = MPI_Request_free(request);
    return err;
#else
    return 0;
#endif
}
//==================================================================================================================================
int md_waitany(int* const indexCompleted, int* const source, int* const type, int requestVlen, M1D_MPI_Request* const requestVector)
{
#ifdef HAVE_MPI
    MPI_Status status;
    int count, err;
    if ((err = MPI_Waitany(requestVlen, requestVector, indexCompleted, &status))) {
         md_throw("MPI_Wait ERROR from md_waitany", err);
    }
    MPI_Get_count(&status, MPI_BYTE, &count);
    *source = status.MPI_SOURCE;
    *type = status.MPI_TAG;
    return count;
#else
    return 0;
#endif
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
