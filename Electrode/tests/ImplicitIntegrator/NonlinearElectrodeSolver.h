/**
 *  @file NonlinearElectrodeSolver.h
 *    Class that calculates the solution to a nonlinear, dense, set
 *    of equations (see \ref numerics
 *    and class \link Cantera::NonlinearElectrodeSolver NonlinearElectrodeSolver\endlink).
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

#ifndef CT_NONLINEARELECTRODESOLVER_H
#define CT_NONLINEARELECTRODESOLVER_H

#include "cantera/numerics/NonlinearSolver.h"
#include "cantera/numerics/ResidJacEval.h"
#include "Electrode.h"

using namespace Cantera;
using namespace std;

namespace Cantera {

  // I think steady state is the only option I'm gunning for
#define  NSOLN_TYPE_PSEUDO_TIME_DEPENDENT 2
#define  NSOLN_TYPE_TIME_DEPENDENT 1
#define  NSOLN_TYPE_STEADY_STATE   0
  
#define NSOLN_JAC_NUM 1
#define NSOLN_JAC_ANAL 2
  
  class ResidJacElectrodeEval : public ResidJacEval {
    
  public:

    ResidJacElectrodeEval(Electrode * elec, doublereal atol = 1.0e-13);
    
    int nEquations() const;
    
    virtual void evalResidNJ(const doublereal t, const doublereal delta_t,
			     const doublereal * const y,
			     const doublereal * const ydot,
			     doublereal * const resid,
			     const ResidEval_Type_Enum evalType = Base_ResidEval,
			     const int id_x = -1, 
			     const doublereal delta_x = 0.0);

  protected:
    
    //! Electrode object pointer
    Electrode * m_electrodePtr_;
    
  };
  
}

#endif
