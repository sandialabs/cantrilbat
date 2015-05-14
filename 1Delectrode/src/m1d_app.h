/**
 * @file m1d_EpectraExtras.cpp
 *
 */

/*
 *  $Id: m1d_app.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_APP_H
#define M1D_APP_H

#include <string>
#include <exception>
#include <vector>
#include <ostream>

namespace m1d
{

class Appl
{

protected:
  //! constructor is protected because it's a singeton
  Appl(int nprocs = 1, int myprocID = 0);

public:
  //! Return a pointer to the one and only instance of class Application
  /*
   * If the an Application object has not yet been created it is created
   */
  static Appl*
  Instance(int nprocs = 1, int myprocID = 0);

  ~Appl();

  //! Static function that destroys the application class's data
  static void
  ApplDestroy();

  //! Set an error condition in the Appl without
  //! throwing an exception.
  /*!
   * This routine adds an error message to the end of the stack
   * of errors that accumulates in the Appl
   * class.
   * @param r   location
   * @param msg  Description of the error
   * @ingroup errorhandling
   */
  void
  addError(const std::string &r, const std::string &msg);

  //! Discard the last error message
  /*!
   *  This routine eliminates
   * the last exception to be added to the stack of error messages
   *
   * @ingroup errorhandling
   */
  void
  popError();

  //! Return the number of errors that have been encountered so far
  /*!
   * \ingroup errorhandling
   */
  int
  nErrors() const;

  void
  showErrors(std::ostream& f);

protected:

  std::vector<std::string> errorMessage;
  std::vector<std::string> errorRoutine;

  int NumProcs_;
  int MyProcID_;

private:
  //! Single instance of the Appl class
  static Appl* s_app;

};

}

#endif
