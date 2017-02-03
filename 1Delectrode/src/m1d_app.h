/**
 * @file m1d_app.h Global data is held here
 *
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_APP_H
#define M1D_APP_H

#include <string>
#include <exception>
#include <vector>
#include <ostream>
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//!  Class that contains a single instance per run
/*!
 *  This class contains a listing of all the errors that are produced by the code. 
 *  The class is a singelton. There is only one instance of this class per job, on each processor.
 *  The class keeps track of the number of processors and the processor id of this job.
 *
 */
class Appl
{
protected:
  //! Constructor is protected because it's a singeton
  /*!
   *  @param[in]             nprocs              Number of processors in the job. Defaults to 1
   *  @param[in]             myprocID            Processor ID of the job. Defaults to 0
   */
  Appl(int nprocs = 1, int myprocID = 0);

public:

  //! Return a pointer to the one and only instance of class Application
  /*!
   * If the an Application object has not yet been created it is created
   *
   *  @param[in]             nprocs              Number of processors in the job. Defaults to 1
   *  @param[in]             myprocID            Processor ID of the job. Defaults to 0
   *
   *  @return                                    Returns a pointer to the application instance
   */
  static Appl* Instance(int nprocs = 1, int myprocID = 0);

  //! Destructor
  ~Appl();

  //! Static function that destroys the application class's data
  static void ApplDestroy();

  //! Set an error condition in the Appl without throwing an exception.
  /*!
   *  This routine adds an error message to the end of the stack of errors that accumulates in the Appl class.
   *
   *  @param[in]             r                   location
   *  @param[in]             msg                 Description of the error
   *
   *  @ingroup errorhandling
   */
  void addError(const std::string &r, const std::string &msg);

  //! Discard the last error message
  /*!
   *  This routine eliminates the last exception to be added to the stack of error messages
   *
   *  @ingroup errorhandling
   */
  void popError();

  //! Return the number of errors that have been encountered so far
  /*!
   *  @return                                    Returns the number of errors
   * \ingroup errorhandling
   */
  int nErrors() const;

  //! Prints all of the errors on an ostream 
  /*!
   *  @param[in]             f                   Reference to an ostream that will be used for the printout
   */ 
  void showErrors(std::ostream& f);

  //! Returns the number of processors
  /*!
   *  @return                                    Returns the number of processors
   */
  int numProcs() const;

  //! Return the processor ID
  /*!
   *  @return                                    Returns the processor ID
   */
  int procID() const;

protected:

  //! Vector of error messages
  std::vector<std::string> errorMessage;

  //! Vector of procedure names for the error messages
  std::vector<std::string> errorRoutine;

public:

  //! Number of processors
  int NumProcs_;

  //! Processor ID of this processor
  int MyProcID_;

private:
  //! Single instance of the Appl class
  static Appl* s_app;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
