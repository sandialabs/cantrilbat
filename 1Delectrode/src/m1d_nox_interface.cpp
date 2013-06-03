//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source: /cvsroot/cads/cantera_apps/1Delectrode/src/m1d_nox_interface.cpp,v $
//  $Author: hkmoffa $
//  $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
//  $Revision: 5 $
// ************************************************************************
//@HEADER

using namespace std;

// ----------   Includes   ----------
#include <iostream>
#include "m1d_nox_interface.h"

// ----------   User Defined Includes   ----------
#include "m1d_nox_problem.h"

//-----------------------------------------------------------------------------
Problem_Interface::Problem_Interface(FiniteElementProblem& Problem) :
  problem(Problem)
{
}

Problem_Interface::~Problem_Interface()
{
}

bool
Problem_Interface::computeF(const Epetra_Vector& x,
                            Epetra_Vector& FVec,
                            FillType flag)
{
  // return problem.evaluate(F_ONLY, &x, &FVec, NULL);
  return true;
}

bool
Problem_Interface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  Epetra_RowMatrix* Jacobian = dynamic_cast<Epetra_RowMatrix*> (&Jac);
  if (Jacobian == NULL) {
    cout << "ERROR: Problem_Interface::computeJacobian() - The supplied"
        << "Epetra_Operator is NOT an Epetra_RowMatrix!" << endl;
    throw ;
  }
  return problem.evaluate(MATRIX_ONLY, &x, NULL, Jacobian);
}

bool
Problem_Interface::computePrecMatrix(const Epetra_Vector& x,
                                     Epetra_RowMatrix& M)
{
  Epetra_RowMatrix* precMatrix = dynamic_cast<Epetra_RowMatrix*> (&M);
  if (precMatrix == NULL) {
    cout << "ERROR: Problem_Interface::computePreconditioner() - The supplied"
        << "Epetra_Operator is NOT an Epetra_RowMatrix!" << endl;
    throw ;
  }
  return problem.evaluate(MATRIX_ONLY, &x, NULL, precMatrix);
}
bool
Problem_Interface::computePreconditioner(const Epetra_Vector& x,
                                         Epetra_Operator& M,
                                         Teuchos::ParameterList* precParams)
{
  cout
      << "ERROR: Problem_Interface::preconditionVector() - Use Explicit Jaciban only for this test problem!"
      << endl;
  throw 1;
}
//-----------------------------------------------------------------------------

