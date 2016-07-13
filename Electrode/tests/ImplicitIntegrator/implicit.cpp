/**
 *  @file mixGasTransport.cpp
 *       test problem for mixture transport
 */

//  Example 
//
// Test case for mixture transport in a gas
// The basic idea is to set up a gradient of some kind.
// Then the resulting transport coefficients out.
// Essentially all of the interface routines should be
// exercised and the results dumped out.
//
// A blessed solution test will make sure that the actual
// solution doesn't change as a function of time or 
// further development. 

// perhaps, later, an analytical solution could be added

// An open Rankine cycle

#ifndef SOLVE_SS
#define SOLVE_SS
#endif

#include <string>
#include <map>
#include <iomanip>

#include "cantera/numerics.h"
#include "cantera/numerics/ResidJacEval.h"    // defines class Water
#include "cantera/numerics/NonlinearSolver.h"
#include "Electrode.h"
//#include "NonlinearElectrodeSolver.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
using namespace std;

ofstream outfile;
double deltaT = 0.01;

void printDbl(double val) {
  if (fabs(val) < 5.0E-17) {
    cout << " nil";
  } else {
    cout << val;
  }
}

class ResidJacElectrodeEval : public ResidJacEval {
  
protected:
  
  //! Electrode object pointer
  Electrode * m_electrodePtr_;

public:
  
  ResidJacElectrodeEval(Electrode * elec, doublereal atol = 1.0e-13) :
    ResidJacEval(atol)
  {
    neq_ = 2;
    m_electrodePtr_ = elec;
  }
  
  int nEquations() const
  {
    return neq_;
  }
  
  virtual void evalResidNJ(const doublereal t, const doublereal delta_t,
			   const doublereal * const y,
			   const doublereal * const ydot,
			   doublereal * const resid,
			   const ResidEval_Type_Enum evalType = Base_ResidEval,
			   const int id_x = -1, 
			   const doublereal delta_x = 0.0)
  {
    
    double x1, x2;
    double eq1, eq2;
    
    x1 = y[0]; 
    x2 = y[1];
    
#ifdef SOLVE_SS
    double yold1 = ydot[0];
    double yold2 = ydot[1];
    
    eq1 = (-101*x1 + 99*x2)/2 - (y[0] - yold1) / deltaT; 
    eq2 = (99*x1 - 101*x2)/2 -  (y[1] - yold2) / deltaT;
    
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
      
};



int main(int argc, char** argv) {
  int k;
  string infile = "diamond.xml";
  outfile.open("outfile.dat");
  //string outputName = "solution.csv";
  //File* fp = fopen(outputName.c_str(),"w");

  try {
    int num_newt_its = 0;
    int num_linear_solves = 0;
    int numBacktracks = 0;
    int loglevelInput = 9;
    Electrode * elec = 0;
    ResidJacElectrodeEval r1(elec);
    int neq = r1.nEquations();
    double y_comm[neq];
    double ylow[neq];
    double yhigh[neq];
    double atol[neq];
    double yold_comm[neq];
    double ydot_comm[neq];;
    double time_curr = 0.0;
    SquareMatrix jac(neq);
    double CJ = 1.0/deltaT;
    for (k = 0; k < neq; k++) {
      //y_comm[k] = 0.1;
      ydot_comm[k] = 0.0;
    }
    y_comm[0] = 1.4142;
    y_comm[1] = 0.0;
    atol[0] = 1.0E-12;
    atol[1] = 1.0E-12;
    ylow[0] = 0.0;
    ylow[1] = 0.0;
    yhigh[0] = 100.0;
    yhigh[1] = 2.0;
  
    NonlinearSolver *nls = new NonlinearSolver(&r1);

#ifdef SOLVE_SS
    int solnType =  NSOLN_TYPE_STEADY_STATE ;
#else
    int solnType =  NSOLN_TYPE_TIME_DEPENDENT;
#endif

    outfile << setw(15) << "time";
    for (int i = 0; i < neq; i++) {
      outfile << "," << setw(13) << "y[" << i << "]";
    }
    outfile << endl;
    
    outfile << setw(15) << time_curr;
    for (int i = 0; i < neq; i++) {
      outfile << "," << setw(15) << y_comm[i];
    }
    outfile << endl;
  
    for ( int s = 0; s < (int)(1.0/deltaT); s++ ) {

#ifdef SOLVE_SS
      yold_comm[0] = y_comm[0];
      yold_comm[1] = y_comm[1];
#endif

      nls->setAtol(DATA_PTR(atol));
      nls->setBoundsConstraints(&ylow[0], &yhigh[0]);
      int flag = nls->solve_nonlinear_problem(solnType, y_comm,
					      yold_comm, CJ,
					      time_curr,
					      jac,
					      num_newt_its,
					      num_linear_solves,
					      numBacktracks,
					      loglevelInput);
      if (flag < 0) {
	cout << "Nonlinear solver did not converge" << endl;
	break;
      }

      ydot_comm[0] = (y_comm[0] - yold_comm[0] ) * CJ;
      ydot_comm[1] = (y_comm[1] - yold_comm[1] ) * CJ;

      time_curr += deltaT;
      outfile << setw(15) << time_curr;
      for (int i = 0; i < neq; i++) {
	outfile << "," << setw(15) << y_comm[i];
      }
      outfile << endl;

    }

    double res[10];

    //r1.evalResidNJ(0.0, 0.0, y_comm, ydot_comm, res);
 
    printf("\nSOLN:\n");

    printf(" Unk    Value        residual if GE 1.0E-15\n");
 
    for (k = 0; k < neq; k++) {
   
      double tmp = res[k];
      if (fabs(tmp) < 1.0E-15) {
	tmp = 0.0;
      }
      printf("    %d %13.5E %13.5E %13.5E\n", k, y_comm[k], ydot_comm[k], tmp);
 
    }
 

  }
  catch (CanteraError) {
    showErrors(cout);
  }
  //fclose(fp);
  return 0;
}
/***********************************************************/
