/**
 * @file m1d_app.cpp
 *
 */

/*
 *  $Id: m1d_app.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_app.h"
#include "m1d_GlobalIndices.h"
#include "m1d_globals.h"

#include <iostream>
#include <fstream>

using namespace std;
namespace m1d
{

//==============================================================================
// Pointer to the single Application instance
Appl* Appl::s_app = 0;

//==============================================================================
// Protected Constructor
Appl::Appl(int nprocs, int myprocID) :
 NumProcs_(nprocs),
 MyProcID_(myprocID)
{
}

//==============================================================================
Appl*
Appl::Instance(int nprocs, int myprocID)
{
  if (Appl::s_app == 0) {
    Appl::s_app = new Appl(nprocs, myprocID);
  }
  return s_app;
}

//==============================================================================
// Static function that destroys the application class's data
void
Appl::ApplDestroy()
{
  if (Appl::s_app != 0) {
    delete Appl::s_app;
    Appl::s_app = 0;
  }
}

//==============================================================================
Appl::~Appl()
{
}

//==============================================================================
void
Appl::addError(const std::string &r, const std::string &msg)
{
  if (NumProcs_ > 1) {

    std::string msge = msg + " (Proc " + intToString(MyProcID_) + ")";
    errorMessage.push_back(msg);
  } else {
    errorMessage.push_back(msg);
  }
  errorRoutine.push_back(r);
}

//==============================================================================
void
Appl::popError()
{
  if (static_cast<int> (errorMessage.size()) > 0) {
    errorRoutine.pop_back();
    errorMessage.pop_back();
  }
}

//==============================================================================
int
Appl::nErrors() const
{
  return static_cast<int> (errorMessage.size());
}

//==============================================================================
void
Appl::showErrors(std::ostream& f)
{
  int i = static_cast<int> (errorMessage.size());
  if (i == 0)
    return;
  f << std::endl << std::endl;
  f << "************************************************" << std::endl;
  f << "    MultiDomain One Dimension Error!            " << std::endl;
  f << "************************************************" << std::endl << endl;
  for (int j = 0; j < i; j++) {
    f << std::endl;
    f << "Procedure: " << errorRoutine[j] << std::endl;
    f << "Error:     " << errorMessage[j] << std::endl;
  }
  f << std::endl << std::endl;
  errorMessage.clear();
  errorRoutine.clear();
}
//==============================================================================

}
