/**
 *  @file m1d_Comm.cpp  Miscellaneous communications routines.
 *
 */

/*
 *  $Id: m1d_Comm.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */
#include <stdarg.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include <sstream>
#include <string>

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_Import.h>

#include "m1d_Comm.h"
#include "m1d_defs.h"

#include <iostream>
#include <iomanip>

using namespace std;

int statprocID = 0;
bool statInPrint0SyncBlock = false;

//#define DEBUG_COMM

namespace m1d
{

const int typeHandOff = M1D_MSG_TYPE_BASE - 3;
const int type0CommEnd = M1D_MSG_TYPE_BASE - 2;
const int type0CommStart = M1D_MSG_TYPE_BASE - 4;
const int typeP0 = M1D_MSG_TYPE_BASE - 1;
const int type0CommAck = M1D_MSG_TYPE_BASE - 5;
const int type0CommAck2 = M1D_MSG_TYPE_BASE - 6;
//=====================================================================================
int
procChoice_Max(double pmax, const Epetra_Comm& cc, double &gmax)
{
  double gb_me = -1.0;
  double me_me = -1.0;
  cc.MaxAll(&pmax, &gmax, 1);
  if (gmax <= (pmax + 1.0E-13 * fabs(pmax))) {
    me_me = cc.MyPID();
  }
  cc.MaxAll(&me_me, &gb_me, 1);
  int gWinner = (int) gb_me;
  return gWinner;
}
//=====================================================================================
int
procChoice_Max_Brcst1Int(double pmax, const Epetra_Comm& cc, double &gmax, int &pint)
{
  double gb_me = -1.0;
  double me_me = -1.0;

  cc.MaxAll(&pmax, &gmax, 1);
  if (gmax <= (pmax + 1.0E-13 * fabs(pmax))) {
    me_me = cc.MyPID();
  }
  cc.MaxAll(&me_me, &gb_me, 1);
  int gWinner = (int) gb_me;
  cc.Broadcast(&pint, 1, gWinner);
  return gWinner;
}
//=====================================================================================
// Returns the processor which has the minimum value of the first parameter
// and broadcasts an int from the winning processor.
/*
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
procChoice_Min_Brcst1Int(const double pmin, const Epetra_Comm& cc, double &gmin, int &pint)
{
  double gb_me = -1.0;
  double me_me = -1;
  cc.MinAll(const_cast<double *>(&pmin), &gmin, 1);
  if (gmin >= (pmin - 1.0E-13 * fabs(pmin))) {
    me_me = cc.MyPID();
  }
  cc.MaxAll(&me_me, &gb_me, 1);
  int gWinner = (int) gb_me;
  cc.Broadcast(&pint, 1, gWinner);
  return gWinner;
}
//=====================================================================================
int
procChoice_Min(const double pmin, const Epetra_Comm& cc, double &gmin)
{
  double gb_me = -1.0;
  double me_me = -1;
  cc.MinAll(const_cast<double *>(&pmin), &gmin, 1);
  if (gmin >= (pmin - 1.0E-13 * fabs(pmin))) {
    me_me = cc.MyPID();
  }
  cc.MaxAll(&me_me, &gb_me, 1);
  int gWinner = (int) gb_me;
  return gWinner;
}
//=====================================================================================

void
print_sync_start(int do_print_line, Epetra_Comm &comm)
{
#ifdef HAVE_MPI
  //int bytes;
  int flag = 1, type;
  M1D_MPI_Request request;

  type = M1D_MSG_TYPE_BASE - 1;

  int procID = comm.MyPID();
  statprocID = procID;
  if (procID) {
    int from = procID - 1;
    M1D::md_wrap_iread((void *) &flag, sizeof(int), &from, &type, &request);
    //bytes = M1D::md_wrap_wait(&from, &type, &request);
    M1D::md_wrap_wait(&from, &type, &request);
    /*
     * TODO: check return status to see if we answered the write request
     */  
  } else {
    if (do_print_line) {
      printf("\n");
      for (flag = 0; flag < 37; flag++)
        printf("#");
      printf(" PRINT_SYNC_START ");
      for (flag = 0; flag < 25; flag++)
        printf("#");
      (void) printf("\n");
    }
  }
#endif
}

//=====================================================================================
void
print_sync_end(int do_print_line, Epetra_Comm &comm)
{

#ifdef HAVE_MPI
  int flag = 1, from, type, to;
  M1D_MPI_Request request, request2;

  type = M1D_MSG_TYPE_BASE - 1;
  int procID = comm.MyPID();
  int nprocs = comm.NumProc();

  if (procID < nprocs - 1)
    to = procID + 1;
  else {
    to = 0;
    if (do_print_line) {
      (void) printf("\n");
      for (flag = 0; flag < 37; flag++)
        (void) printf("#");
      (void) printf(" PRINT_SYNC_END__ ");
      for (flag = 0; flag < 25; flag++)
        (void) printf("#");
      (void) printf("\n\n");
    }
  }
  // Write to the next processor on the list to say that you are done
  M1D::md_wrap_iwrite((void *) &flag, sizeof(int), to, type, &request);

  if (procID == 0) {
    from = nprocs - 1;
    M1D::md_wrap_iread((void *) &flag, sizeof(int), &from, &type, &request2);
    M1D::md_wrap_wait(&from, &type, &request2);
  }

  M1D::md_wrap_request_cancel(&request); /* Free request object from iwrite call */
  /*
   * Do a final sync amongst all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from
   * Proc (Num_Proc-1).
   */

  comm.Barrier();
#else
  if (do_print_line) {
    (void) printf("\n");
    for (int flag = 0; flag < 37; flag++)
      (void) printf("#");
    (void) printf(" PRINT_SYNC_END__ ");
    for (int flag = 0; flag < 25; flag++)
      (void) printf("#");
    (void) printf("\n\n");
  }
#endif
}
//=====================================================================================

void
print0_sync_start(int do_print_line, stream0 &ss, const Epetra_Comm &comm)
{
  FILE *of = ss.outFile();
  int flag = 1;
  M1D_MPI_Request request;
  int procID = comm.MyPID();
  if (statInPrint0SyncBlock) {
    throw m1d_Error("print0_sync_start ERROR proc " + Cantera::int2str(procID), "Already in a print0_sync block");
  }
  statInPrint0SyncBlock = true;

  statprocID = procID;
  //int num_proc = comm.NumProc();
  if (procID) {
    int from = procID - 1;
    // We wait here for the processor before us to send us a message
    // saying we are ready
    M1D::md_wrap_iread_specific((void *) &flag, sizeof(int), from, typeHandOff, &request);
    M1D::md_wrap_wait_specific(from, typeHandOff, &request);
    // Ok we need to write to processor zero that we are in charge
    int to = 0;
    M1D::md_write((void *) &flag, sizeof(int), to, type0CommStart);
  } else {
    if (do_print_line) {
      fprintf(of,"\n");
      for (flag = 0; flag < 37; flag++)
        fprintf(of,"#");
      fprintf(of," PRINT0_SYNC_START ");
      for (flag = 0; flag < 25; flag++)
        fprintf(of,"#");
      (void) fprintf(of,"\n");
    }
  }
}

//=====================================================================================
void
stream0::print0(const char *format, ...)
{
  va_list ap;
  static char buf[1024];
  char rbuf[8];
  if (!statInPrint0SyncBlock) {
     // HKM - will turn this on next
     //   throw m1d_Error("stream0::print0 ERROR proc " + Cantera::int2str(statprocID), "Not in a print0_sync block");
  }
  va_start(ap, format);
  vsnprintf(buf, 1023, format, ap);
  va_end(ap);
  buf[1023] = '\0';
  int typeP0 = M1D_MSG_TYPE_BASE - 1;
  if (statprocID == 0) {
    fprintf(outFILE_, "%s", buf);
  } else {
    M1D::md_write((void * const ) buf, 1024, 0, typeP0);
    nBufsSent++;
    M1D::md_read_specific((void * const ) rbuf, 4, 0, type0CommAck2);
  }
}
//=====================================================================================
int
stream0::fprintf0only(const char *format, ...)
{
  int n = 0;
  if (statprocID == 0) {
      va_list ap;
      static char buf[1024];
      va_start(ap, format);
      n = vsnprintf(buf, 1023, format, ap);
      if (n < 1023) {
        va_end(ap);
        buf[1023] = '\0';
        fprintf(outFILE_, "%s", buf);
      } else {
        int sze = n+1;
        char *np = (char *) malloc(sizeof(char) * (n+1));
        va_start(ap, format);
        n = vsnprintf(np, sze, format, ap);
        va_end(ap);
        fprintf(outFILE_, "%s", np); 
        free(np);
      }
  }
  return n;
}
//=====================================================================================
void stream0::drawline0(int indentSpaces, int num, char lineChar)
{
  std::string buf;
  for (int i = 0; i < indentSpaces; i++) {
    buf += " ";
  }
  for (int i = 0; i < num; i++) {
     buf += lineChar;
  }
  buf += "\n";
  print0("%s", buf.c_str());
}
//=====================================================================================
int
sprint0(const char *ibuf, FILE *ff)
{
  int iSent = 0;
  const char *jbuf = ibuf;
  char buf[1024];
  char rbuf[8];
  int typeP0 = M1D_MSG_TYPE_BASE - 1;
  if (!statInPrint0SyncBlock) {
    // HKM Will turn this on next
    //throw m1d_Error("sprint0 ERROR proc " + Cantera::int2str(statprocID), "Not in a print0_sync block");
  }
  if (statprocID == 0) {
     fprintf(ff, "%s", ibuf);
  } else {
     int ll = strlen(ibuf);
     if (ll == 0) {
       return 0;
     }
     bool weAreNotDone = true;
     do {
       ll = strlen(jbuf);
       if (ll <= 1023) {
         weAreNotDone = false;
        } else if (ll > 1023) {
          ll = 1023;
        }
        for (int j = 0; j < ll; j++) {
          buf[j] = jbuf[j];
        }
        buf[ll] = '\0';
        M1D::md_write((void * const ) buf, 1024, 0, typeP0);
        iSent++;
        M1D::md_read_specific((void * const ) rbuf, 4, 0, type0CommAck2);
        jbuf += ll;
      } while (weAreNotDone);
  }
  return iSent;
}
//====================================================================================
void
ssprint0(stream0 &ss)
{
  int nSent = sprint0((ss.str()).c_str(), ss.outFILE_);
  ss.nBufsSent += nSent;
  ss.str("");
}
//=====================================================================================
void
print0_sync_end(int do_print_line, stream0 &ss, const Epetra_Comm &comm)
{
  char buf[1025];
  char bufalt[1024];
  //char bufShort[25];
  int flag = 1;
  int from, to;
  int fromReturn;
  int typeReturn;
  //int rbytes;
  int indexCompleted;
  M1D_MPI_Request requestV[2];
  bool weAreNotDone = true;
#ifdef DEBUG_COMM
  int nBufs;
#endif
  FILE *of = ss.outFile();

  int type;

  int procID = comm.MyPID();
  int nprocs = comm.NumProc();
  if (!statInPrint0SyncBlock) {
    throw m1d_Error("print0_sync_end ERROR proc " + Cantera::int2str(procID), "Not already in a print0_sync block");
  }

  if (procID < nprocs - 1) {

    to = procID + 1;
  } else {
    to = 0;
  }
  // Write to processor zero to say you are done with the writes to it.
  if (procID != 0) {
    sprintf(bufalt, "%d", ss.nBufsSent);
    M1D::md_write(bufalt, 1024, 0, type0CommEnd);
    ss.nBufsSent = 0;
    int bytes = sizeof(int);
    M1D::md_read_specific((void * const ) &buf, bytes, 0, type0CommAck);

  }
  // Write to the next processor on the list to say that you are done
  if (nprocs > 1) {
    M1D::md_write((void *) &flag, sizeof(int), to, typeHandOff);
  }
  if (procID == 0) {
    for (int rproc = 1; rproc < nprocs; rproc++) {
#ifdef DEBUG_COMM
      nBufs = 0;
#endif
      // limit the receives to only the rproc
      from = rproc;
      // Limit the type to the 0CommStart type.

      // Do a blocking read from the processor to start off with
      int bytes = sizeof(int);
      //rbytes = M1D::md_read_specific((void * const ) buf, bytes, from, type0CommStart);
      M1D::md_read_specific((void * const ) buf, bytes, from, type0CommStart);
      // Set up a message type that will indicate that we are through reading from the
      // current processor.
      type = type0CommEnd;
      M1D::md_wrap_iread((void *) bufalt, 1024, &from, &type, &(requestV[0]));
      weAreNotDone = true;
      do {
        buf[1024] = '\0';
        type = typeP0;
        M1D::md_wrap_iread((void *) buf, 1024, &from, &type, &(requestV[1]));

        //rbytes = M1D::md_waitany(&indexCompleted, &fromReturn, &typeReturn, 2, requestV);
        M1D::md_waitany(&indexCompleted, &fromReturn, &typeReturn, 2, requestV);
        if (fromReturn != from) {
          fprintf(stderr,"we are confused\n");
          exit(-1);
        }
        if (indexCompleted == 0) {
          weAreNotDone = false;
          if (typeReturn != type0CommEnd) {
            printf("we are confused by type0CommEnd\n");
            exit(-1);
          }
#ifdef DEBUG_COMM
          int bufsExpected = std::atoi(bufalt);
          printf("DC: We have received %d bufs from proc %d\n", nBufs, rproc);
          if (bufsExpected != nBufs) {
            printf("DC:        Error bufsExpected = %d\n", bufsExpected);
          }
#endif
          M1D::md_wrap_request_cancel(&requestV[1]);
          M1D::md_write((void *) &flag, sizeof(int), from, type0CommAck);
        } else {
#ifdef DEBUG_COMM
          nBufs++;
#endif
          if (typeReturn != typeP0) {
            printf("we are confused by typeP0\n");
            exit(-1);
          }
          // Ok, let's actually print the message from Processor 0
          fprintf(of, "%s", buf);
          M1D::md_write((void *) &flag, 4, rproc, type0CommAck2);
        }
      } while (weAreNotDone);
    }

    if (nprocs > 1) {
      from = nprocs - 1;
      //M1D::md_read_specific((void *) &flag, sizeof(int), from, typeHandOff);
      M1D::md_wrap_iread_specific((void *) &flag, sizeof(int), from, typeHandOff, &requestV[0]);
      M1D::md_wrap_wait_specific(from, typeHandOff, &requestV[0]);
    }

    if (do_print_line) {
      (void) fprintf(of, "\n");
      for (flag = 0; flag < 37; flag++) fprintf(of,"#");
      (void) fprintf(of, " PRINT0_SYNC_END__ ");
      for (flag = 0; flag < 25; flag++) fprintf(of, "#"); 
      (void) fprintf(of,"\n\n");
    }
  }

  statInPrint0SyncBlock = false;
  /*
   * Do a final sync amongst all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from
   * Proc (Num_Proc-1).
   */
  comm.Barrier();
}
//=====================================================================================
stream0::stream0(FILE *ff) :
  std::ostringstream(), nBufsSent(0),
  outFILE_(ff)
{
}
//=====================================================================================
stream0&
stream0::flush()
{
  ssprint0(*this);
  return *this;
}
//=====================================================================================
void
print0_epMultiVector(const Epetra_MultiVector &v, const char *label, FILE *of)
{
  stream0 os(of);
  const Epetra_Comm &cc = v.Comm();
  int mypid = cc.MyPID();
  print0_sync_start(false, os, cc);
  if (mypid == 0) {
    fprintf(of, "######################################");
    fprintf(of, " Start of Epetra_MultiVector Printout ");
    fprintf(of, "######################################\n");
    if (label) {
      fprintf(of, "############################### ");
      fprintf(of, "%s", label);
      fprintf(of, " ###################################\n");
    }
  }
  Print0_epMultiVector(os, v);
  print0_sync_end(false, os, cc);
  if (mypid == 0) {
    fprintf(of, "######################################");
    fprintf(of, " End of Epetra_MultiVector Printout ");
    fprintf(of, "######################################\n");
  }
}
//=====================================================================================
void
print0_epIntVector(const Epetra_IntVector &v, const char *label, FILE *of)
{
  stream0 os(of);
  const Epetra_Comm &cc = v.Comm();
  int mypid = cc.MyPID();
  print0_sync_start(false, os, cc);
  if (mypid == 0) {
    fprintf(of, "######################################");
    fprintf(of, " Start of Epetra_MultiVector Printout ");
    fprintf(of, "######################################\n");
    if (label) {
      fprintf(of, "############################### ");
      fprintf(of, "%s", label);
      fprintf(of, " ###################################\n");
    }
  }
  Print0_epIntVector(os, v);
  print0_sync_end(false, os, cc);
  if (mypid == 0) {
    fprintf(of, "######################################");
    fprintf(of, " End of Epetra_MultiVector Printout ");
    fprintf(of, "######################################\n");
  }
}
//=====================================================================================
void
Print0_epMultiVector(stream0 &os, const Epetra_MultiVector &v)
{
  const Epetra_BlockMap &vmap = v.Map();
  int MyPID = vmap.Comm().MyPID();

  int NumVectors1 = v.NumVectors();
  int NumMyElements1 = vmap.NumMyElements();
  int MaxElementSize1 = vmap.MaxElementSize();
  int * MyGlobalElements1 = vmap.MyGlobalElements();
  int * FirstPointInElementList1;
  if (MaxElementSize1 != 1)
    FirstPointInElementList1 = vmap.FirstPointInElementList();
  double ** A_Pointers = v.Pointers();
  if (MyPID == 0) {
    os.width(8);
    os << "     MyPID";
    os << "    ";
    os.width(12);
    if (MaxElementSize1 == 1) {
      os << "       GID                   Value";
    } else {
      os << "     GID/Point";
      for (int j = 0; j < NumVectors1; j++) {
        os.width(20);
        os << "Value  ";
      }
    }
    os << '\n';
  }
  for (int i = 0; i < NumMyElements1; i++) {
    for (int ii = 0; ii < vmap.ElementSize(i); ii++) {
      int iii;
      os.width(10);
      os << MyPID;
      os << "    ";
      os.width(10);
      if (MaxElementSize1 == 1) {
        os << MyGlobalElements1[i] << "    ";
        iii = i;
      } else {
        os << MyGlobalElements1[i] << "/" << ii << "    ";
        iii = FirstPointInElementList1[i] + ii;
      }
      for (int j = 0; j < NumVectors1; j++) {
        os.width(20);
        os << A_Pointers[j][iii];
      }
      os << '\n';
    }
  }
  ssprint0(os);
}
//=====================================================================================
void
Print0_epIntVector(stream0 &os, const Epetra_IntVector &v)
{
  const Epetra_BlockMap &vmap = v.Map();
  int MyPID = vmap.Comm().MyPID();


  int NumMyElements1 = vmap.NumMyElements();
  int MaxElementSize1 = vmap.MaxElementSize();
  int * MyGlobalElements1 = vmap.MyGlobalElements();
  int * FirstPointInElementList1;
  if (MaxElementSize1 != 1)
    FirstPointInElementList1 = vmap.FirstPointInElementList();
  //double ** A_Pointers = v.Pointers();
  int *A_Values = v.Values();
  if (MyPID == 0) {
    os.width(8);
    os << "     MyPID";
    os << "    ";
    os.width(12);
    if (MaxElementSize1 == 1) {
      os << "       GID                   Value";
    } else {
      os << "     GID/Point";
      os.width(20);
      os << "Value  ";
    }
    os << '\n';
  }
  for (int i = 0; i < NumMyElements1; i++) {
    for (int ii = 0; ii < vmap.ElementSize(i); ii++) {
      int iii;
      os.width(10);
      os << MyPID;
      os << "    ";
      os.width(10);
      if (MaxElementSize1 == 1) {
        os << MyGlobalElements1[i] << "    ";
        iii = i;
      } else {
        os << MyGlobalElements1[i] << "/" << ii << "    ";
        iii = FirstPointInElementList1[i] + ii;
      }
   
      os.width(20);
      os << A_Values[iii];
      os << '\n';
    }
  }
  ssprint0(os);
}
//=====================================================================================
void
print0_epBlockMap(const Epetra_BlockMap &map, const char *label, FILE *of)
{
  stream0 os(of);
  const Epetra_Comm &cc = map.Comm();
  int mypid = cc.MyPID();
  print0_sync_start(false, os, cc);
  if (mypid == 0) {
    fprintf(of, "######################################");
    fprintf(of, " Start of Epetra_BlockMap Printout ");
    fprintf(of, "######################################\n");
    if (label) {
      fprintf(of, "############################### ");
      fprintf(of, "%s", label);
      fprintf(of, " ###################################\n");
    }
  }
  Print0_epBlockMap(os, map);
  print0_sync_end(false, os, cc);
  if (mypid == 0) {
    fprintf(of, "######################################");
    fprintf(of, " End of Epetra_BlockMap Printout ");
    fprintf(of, "######################################\n");
  }
}
//=====================================================================================
void
Print0_epBlockMap(stream0 &os, const Epetra_BlockMap &vmap)
{
  int MyPID = vmap.Comm().MyPID();

  int * MyGlobalElements1 = vmap.MyGlobalElements();
  int * FirstPointInElementList1 = 0;
  int * ElementSizeList1 = 0;
  if (!vmap.ConstantElementSize()) {
    FirstPointInElementList1 = vmap.FirstPointInElementList();
    ElementSizeList1 = vmap.ElementSizeList();
  }

  //int NumMyElements1 = vmap.NumMyElements();
  //int MaxElementSize1 = vmap.MaxElementSize();

  if (MyPID == 0) {
    os << "\nNumber of Global Elements  = ";
    os << vmap.NumGlobalElements();
    os << endl;
    os << "Number of Global Points    = ";
    os << vmap.NumGlobalPoints();
    os << endl;
    os << "Maximum of all GIDs        = ";
    os << vmap.MaxAllGID();
    os << endl;
    os << "Minimum of all GIDs        = ";
    os << vmap.MinAllGID();
    os << endl;
    os << "Index Base                 = ";
    os << vmap.IndexBase();
    os << endl;
    if (vmap.ConstantElementSize())
      os << "Constant Element Size      = ";
    os << vmap.ElementSize();
    os << endl;
  }
  os << endl;

  os << "Number of Local Elements   = ";
  os << vmap.NumMyElements();
  os << endl;
  os << "Number of Local Points     = ";
  os << vmap.NumMyPoints();
  os << endl;
  os << "Maximum of my GIDs         = ";
  os << vmap.MaxMyGID();
  os << endl;
  os << "Minimum of my GIDs         = ";
  os << vmap.MinMyGID();
  os << endl;
  os << endl;

  os.width(14);
  os << "     MyPID";
  os << "    ";
  os.width(14);
  os << "       Local Index ";
  os << " ";
  os.width(14);
  os << "      Global Index ";
  os << " ";
  if (!vmap.ConstantElementSize()) {
    os.width(14);
    os << " FirstPointInElement ";
    os << " ";
    os.width(14);
    os << "   ElementSize ";
    os << " ";
  }
  os << endl;

  for (int i = 0; i < vmap.NumMyElements(); i++) {
    os.width(14);
    os << MyPID;
    os << "    ";
    os.width(14);
    os << i;
    os << "    ";
    os.width(14);
    os << MyGlobalElements1[i];
    os << "    ";
    if (!vmap.ConstantElementSize()) {
      os.width(14);
      os << FirstPointInElementList1[i];
      os << "    ";
      os.width(14);
      os << ElementSizeList1[i];
      os << "    ";
    }
    os << endl;
  }

  ssprint0(os);
}
//=====================================================================================
void
print0_epVbrMatrix(const Epetra_VbrMatrix &mat, const char *label, FILE *of)
{
  stream0 os(of);
  const Epetra_Comm &cc = mat.Comm();
  int mypid = cc.MyPID();
  print0_sync_start(true, os, cc);
  if (mypid == 0) {
    fprintf(of, "######################################");
    fprintf(of, " Start of Epetra_BlockMap Printout ");
    fprintf(of, "######################################\n");
    if (label) {
      fprintf(of, "############################### ");
      fprintf(of, "%s", label);
      fprintf(of, " ###################################\n");
    }
  }
  Print0_epVbrMatrix(os, mat);
  print0_sync_end(true, os, cc);
  if (mypid == 0) {
    fprintf(of, "######################################");
    fprintf(of, " End of Epetra_BlockMap Printout ");
    fprintf(of, "######################################\n");
  }
}
//=====================================================================================
void
Print0_epVbrMatrix(stream0 &os, const Epetra_VbrMatrix &mat)
{
  const Epetra_BlockMap &rMap = mat.RowMap();
  int MyPID = rMap.Comm().MyPID();
  //int NumProc = rMap.Comm().NumProc();
  if (MyPID == 0) {
    os << "\nNumber of Global Block Rows  = ";
    os << mat.NumGlobalBlockRows();
    os << endl;
    os << "Number of Global Block Cols  = ";
    os << mat.NumGlobalBlockCols();
    os << endl;
    os << "Number of Global Block Diags = ";
    os << mat.NumGlobalBlockDiagonals();
    os << endl;
    os << "Number of Global Blk Entries = ";
    os << mat.NumGlobalBlockEntries();
    os << endl;
    os << "Global Max Num Block Entries = ";
    os << mat.GlobalMaxNumBlockEntries();
    os << endl;
    os << "\nNumber of Global Rows        = ";
    os << mat.NumGlobalRows();
    os << endl;
    os << "Number of Global Cols        = ";
    os << mat.NumGlobalCols();
    os << endl;
    os << "Number of Global Diagonals   = ";
    os << mat.NumGlobalDiagonals();
    os << endl;
    os << "Number of Global Nonzeros    = ";
    os << mat.NumGlobalNonzeros();
    os << endl;
    os << "Global Maximum Num Entries   = ";
    os << mat.GlobalMaxNumNonzeros();
    os << endl;
    if (mat.LowerTriangular())
      os << " ** Matrix is Lower Triangular **";
    os << endl;
    if (mat.UpperTriangular())
      os << " ** Matrix is Upper Triangular **";
    os << endl;
    if (mat.NoDiagonal())
      os << " ** Matrix has no diagonal     **";
    os << endl;
    os << endl;
  }

  os << "\nNumber of My Block Rows  = ";
  os << mat.NumMyBlockRows();
  os << endl;
  os << "Number of My Block Cols  = ";
  os << mat.NumMyBlockCols();
  os << endl;
  os << "Number of My Block Diags = ";
  os << mat.NumMyBlockDiagonals();
  os << endl;
  os << "Number of My Blk Entries = ";
  os << mat.NumMyBlockEntries();
  os << endl;
  os << "My Max Num Block Entries = ";
  os << mat.MaxNumBlockEntries();
  os << endl;
  os << "\nNumber of My Rows        = ";
  os << mat.NumMyRows();
  os << endl;
  os << "Number of My Cols        = ";
  os << mat.NumMyCols();
  os << endl;
  os << "Number of My Diagonals   = ";
  os << mat.NumMyDiagonals();
  os << endl;
  os << "Number of My Nonzeros    = ";
  os << mat.NumMyNonzeros();
  os << endl;
  os << "My Maximum Num Entries   = ";
  os << mat.MaxNumBlockEntries();
  os << endl;
  os << endl;

  os << flush;
  int NumBlockRows1 = mat.NumMyBlockRows();
  int MaxNumBlockEntries1 = mat.MaxNumBlockEntries();
  int * BlockIndices1 = new int[MaxNumBlockEntries1];
  Epetra_SerialDenseMatrix ** Entries1;
  int RowDim1, NumBlockEntries1;
  int i, j;

  os.width(8);
  os << "   Processor ";
  os.width(10);
  os << "   Block Row Index ";
  os.width(10);
  os << "   Block Col Index";
  os.width(10);
  os << "   Values     ";
  os << endl;

  for (i = 0; i < NumBlockRows1; i++) {
    int BlockRow1 = mat.GRID(i); // Get global row number
    mat.ExtractGlobalBlockRowPointers(BlockRow1, MaxNumBlockEntries1, RowDim1, NumBlockEntries1, BlockIndices1,
        Entries1);

    for (j = 0; j < NumBlockEntries1; j++) {
      os.width(8);
      os << MyPID;
      os << "    ";
      os.width(10);
      os << BlockRow1;
      os << "    ";
      os.width(10);
      os << BlockIndices1[j];
      os << "    ";
      os.width(20);

      if (Entries1[j] == 0) {
        os << "Block Entry == NULL" << endl;
        continue;
      }

      Epetra_SerialDenseMatrix entry(View, Entries1[j]->A(), Entries1[j]->LDA(), RowDim1, Entries1[j]->N());
      int nrow = RowDim1;
      int ncol = entry.N();
      for (int i = 0; i < nrow; i++) {
        if (i != 0) {
          os << "                                        ";
        }
        for (int j = 0; j < ncol; j++) {
          os.width(11);
          os << setprecision(4);
          os << scientific;
          os << entry[j][i] << " ";
        }
        os << endl;
      }

    }
    os << endl;
  }
  /*
   * Cleanup section
   */
  delete[] BlockIndices1;
  /*
   * output the stringstring on Processor 0
   */
  ssprint0(os);
}
//=====================================================================================
stream0 &
operator<<(stream0 &oss, const Epetra_MultiVector &v)
{
  Print0_epMultiVector(oss, v);
  return oss;
}
//=====================================================================================
}
//=====================================================================================
