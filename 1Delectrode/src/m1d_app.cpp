/**
 * @file m1d_app.cpp Singleton containing the application error messages and MPI information
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_app.h"
#include "m1d_globals.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Pointer to the single Application instance
Appl* Appl::s_app = 0;
//==================================================================================================================================
Appl::Appl(int nprocs, int myprocID) :
    NumProcs_(nprocs),
    MyProcID_(myprocID)
{
}
//==================================================================================================================================
Appl* Appl::Instance(int nprocs, int myprocID)
{
    if (Appl::s_app == nullptr) {
        Appl::s_app = new Appl(nprocs, myprocID);
    }
    return s_app;
}
//==================================================================================================================================
void Appl::ApplDestroy()
{
    if (Appl::s_app != nullptr) {
        delete Appl::s_app;
        Appl::s_app = nullptr;
    }
}
//==================================================================================================================================
Appl::~Appl()
{
}
//==================================================================================================================================
void Appl::addError(const std::string& r, const std::string& msg)
{
    if (NumProcs_ > 1) {
        std::string msge = msg + " (Proc " + intToString(MyProcID_) + ")";
        errorMessage.push_back(msg);
    } else {
        errorMessage.push_back(msg);
    }
    errorRoutine.push_back(r);
}
//==================================================================================================================================
void Appl::popError()
{
    if (static_cast<int>(errorMessage.size()) > 0) {
        errorRoutine.pop_back();
        errorMessage.pop_back();
    }
}
//==================================================================================================================================
int Appl::nErrors() const
{
    return static_cast<int>(errorMessage.size());
}
//==================================================================================================================================
void Appl::showErrors(std::ostream& f)
{
    int i = static_cast<int>(errorMessage.size());
    if (i == 0) {
        return;
    }
    f << "\n\n************************************************\n";
    f << "    MultiDomain One Dimension Error!            \n";
    f << "************************************************\n";
    for (int j = 0; j < i; j++) {
        f << "\nProcedure: " << errorRoutine[j] << "\n";
        f <<   "Error:     " << errorMessage[j] << "\n";
    }
    f << "\n" << std::endl;
    errorMessage.clear();
    errorRoutine.clear();
}
//==================================================================================================================================
int Appl::numProcs() const
{
    return NumProcs_;
}
//==================================================================================================================================
int Appl::procID() const
{
    return MyProcID_;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

