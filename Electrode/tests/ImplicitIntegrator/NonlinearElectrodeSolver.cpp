/**
 *
 *  @file NonlinearElectrodeSolver.cpp
 *
 *  Damped Newton solver for 0D and 1D problems
 */

/*
 *  $Date: 2013-01-07 14:15:37 -0700 (Mon, 07 Jan 2013) $
 *  $Revision: 496 $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CANTERA_APP
#define CANTERA_APP
#endif

#include <limits>

#include "zuzax/numerics/SquareMatrix.h"
#include "NonlinearElectrodeSolver.h"

#include "zuzax/base/clockWC.h"
#include "zuzax/kernel/vec_functions.h"
#include <ctime>

#include <cfloat>

extern void print_line(const char *, int);

#include <vector>
#include <cstdio>
#include <cmath>


#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif

using namespace Zuzax;
using namespace std;

namespace Zuzax {


 //====================================================================================================================
  //-----------------------------------------------------------
  //                 Constants
  //-----------------------------------------------------------

  const double DampFactor = 4;
  const int NDAMP = 7;
 //====================================================================================================================


  //===================================================================================================================
  ResidJacElectrodeEval::ResidJacElectrodeEval(Electrode * elec, doublereal atol) :
      ResidJacEval(atol)
  {
    neq_ = 2;
    m_electrodePtr_ = elec;
  }
  //===================================================================================================================
  int ResidJacElectrodeEval::nEquations() const
  {
    return neq_;
  }
  //===================================================================================================================

  void ResidJacElectrodeEval::evalResidNJ(doublereal t, const doublereal deltaT,
			   const doublereal * const y, const doublereal * const ydot,
			   doublereal * const resid,  const ResidEval_Type_Enum evalType,
			   int id_x, doublereal delta_x)
  {

    double x1, x2;
    double eq1, eq2;

    x1 = y[0]-0.5; 
    x2 = y[1]-0.9;

#ifdef SOLVE_SS
    double yold1 = ydot[0];
    double yold2 = ydot[1];
    double deltat = 0.1;

    eq1 = (-101*x1 + 99*x2)/2 - (y[0] - yold1) / deltat; 
    eq2 = (99*x1 - 101*x2)/2 -  (y[1] - yold2) / deltat;

#else

    eq1 = (-1001*x1 + 999*x2)/2 - ydot[0]; 
    eq2 = (999*x1 - 1001*x2)/2 -  ydot[1];

#endif

//     eq1 = (-1001*x1 + 999*x2)/2; 
//     eq2 = (999*x1 - 1001*x2)/2;
  
    resid[0] = eq1;
    resid[1] = eq2; 
//     cout << "resid0 = " << resid[0] << ", resid1 = " << resid[1] << endl;
//     if (fabs(resid[0]) + fabs(resid[1]) < 1.0E-1) {
//       outfile << t << "," << x1 << endl;
//     }
    
    
  }
}
