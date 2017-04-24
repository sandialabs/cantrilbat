/**
 *  @file m1d_ProblemResidEval.cpp
 *
 **/

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_FlatBatteryResidEval.h"

//#include "m1d_SurDomain_AnodeCollector.h"
//#include "m1d_SurDomain_CathodeCollector.h"
#include "m1d_SurDomain_FlatFeS2Cathode.h"
#include "m1d_SurDomain_FlatLiSiAnode.h"

#include "m1d_porousLiKCl_dom1D.h"
#include "m1d_porousElectrode_dom1D.h"

#include "m1d_DomainLayout.h"
#include "m1d_ProblemStatementCell.h"

#include "m1d_Comm.h"
#include "m1d_GlobalIndices.h"
#include "m1d_globals.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

extern m1d::ProblemStatementCell PSinput;

namespace m1d
{
  //=====================================================================================================================
  // Default constructor
  /*
   *
   * @param atol   Absolute tolerance for calculation
   */
  FlatBatteryResidEval::FlatBatteryResidEval(double atol) :
      BatteryResidEval(atol)
  {
      // Set the anode and cathode types to a flat domain boundary
      anodeType_ = 1;
      cathodeType_ = 1;
  }
  //=====================================================================================================================
  // Destructor
  /*
   *
   */
  FlatBatteryResidEval::~FlatBatteryResidEval()
  {
  }
  //=====================================================================================================================
  /*!
   *
   * @param r  Object to be copied
   */
  FlatBatteryResidEval::FlatBatteryResidEval(const FlatBatteryResidEval &r) :
      BatteryResidEval(r.m_atolDefault)
  {
    *this = r;
  }
  //=====================================================================================================================
  // Assignment operator
  /*
   *
   * @param r  Object to be copied
   * @return   Returns a copy of the current problem
   */
  FlatBatteryResidEval &
  FlatBatteryResidEval::operator=(const FlatBatteryResidEval &r)
  {
    if (this == &r) {
      return *this;
    }

    BatteryResidEval::operator=(r);

    return *this;
  }

//=====================================================================================================================
  //=====================================================================================================================
  // This function may be used to create output at various points in the
  // execution of an application.
  /*
   *   These functions are not affected by the print controls of the nonlinear solver
   *   and the time stepper.
   *
   *      ievent is a description of the event that caused this
   *      function to be called.
   *
   *      @param ievent  Event that's causing this routine to be called.
   *                     =  0 Initial conditions for a calculation
   *                     =  1 Completion of a successful intermediate time step.
   *                     =  2 Completion of a successful Final time or final calculation.
   *                     =  3 Completion of a successful Intermediate nonlinear step
   *                     = -1 unsuccessful time step that converged, but failed otherwise
   *                     = -2 unsuccessful nonlinear step.
   *
   *      @param time_current      Current time
   *      @param delta_t_n         Current value of delta_t
   *      @param istep             Current step number
   *      @param y_n               Current value of the solution vector
   *      @param ydot_n_ptr        Current value of the time deriv of the solution vector
   */
  void
  FlatBatteryResidEval::user_out(const int ievent,
			     const double time_current,
			     const double delta_t_n,
			     const int istep,
			     const Epetra_Vector_Ghosted &y_n,
			     const Epetra_Vector_Ghosted * const ydot_n_ptr)
  {
    BatteryResidEval::user_out(ievent, time_current, delta_t_n, istep, y_n, ydot_n_ptr);
  }

  //=====================================================================================================================
  // Write the solution to either the screen or to a log file
  /*
   *
   * @param ievent  Type of the event. The following form is used:
   *             0 Initial conditions
   *             1 Completion of a successful intermediate step.
   *             2 Final successful conditions.
   *             3 Intermediate nonlinear step
   *            -1 unsuccessful step
   *
   * @param m_y_n    Current value of the solution vector
   * @param m_ydot_n  Current value of the derivative of the solution vector
   */
  void
  FlatBatteryResidEval::writeSolution(const int ievent,
				  const bool doTimeDependentResid,
				  const double time_current,
				  const double delta_t_n,
				  int istep,
				  const Epetra_Vector_Ghosted &y_n,
				  const Epetra_Vector_Ghosted * const ydot_n_ptr,
				  const Solve_Type_Enum solveType, 
                                  const double delta_t_np1)
  {
    ProblemResidEval::writeSolution(ievent, doTimeDependentResid, time_current, delta_t_n, istep, y_n,
				    ydot_n_ptr, solveType, delta_t_np1);
    if (ievent == 0 || ievent == 1 || ievent == 2) {
      write_IV(ievent, doTimeDependentResid, time_current, delta_t_n, istep, y_n, ydot_n_ptr);
    }
  }
  //=====================================================================================================================
  void
  FlatBatteryResidEval::write_IV(const int ievent,
			     const bool doTimeDependentResid,
			     const double time_current,
			     const double delta_t_n,
			     int istep,
			     const Epetra_Vector_Ghosted &y_n,
			     const Epetra_Vector_Ghosted * const ydot_n_ptr)
  {
    DomainLayout &DL = *DL_ptr_;

    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    //SurDomain_CathodeCollector *c_ptr = dynamic_cast<SurDomain_CathodeCollector *>(d_ptr);
    SurDomain_FlatFeS2Cathode *c_ptr = dynamic_cast<SurDomain_FlatFeS2Cathode *>(d_ptr);
    ZZCantera::Electrode *ElectrodeC_ = c_ptr->ElectrodeC_;

    /*
     * Calculate the rates of production of all species in the Electrode
     */
    double icurr = ElectrodeC_->getNetSurfaceProductionRatesCurrent(0, &c_ptr->electrodeSpeciesProdRates_[0]);


    double c_icurr =  icurr * c_ptr->surfaceArea_;
    double phiCath = c_ptr->phiCathode_;

    SurDomain1D *ad_ptr = DL.SurDomain1D_List[0];
    SurDomain_FlatLiSiAnode *ac_ptr = dynamic_cast<SurDomain_FlatLiSiAnode *>(ad_ptr);

    ZZCantera::Electrode *ElectrodeA_ = ac_ptr->ElectrodeA_;
    double icurrA = ElectrodeA_->getNetSurfaceProductionRatesCurrent(0, &ac_ptr->electrodeSpeciesProdRates_[0]);
   
    double a_icurr = icurrA * ac_ptr->surfaceArea_;
    

    int procID = Comm_ptr->MyPID();

    if (!procID) {
    //looking for cathode capacity and depth of discharge
    BulkDomain1D *cd_ptr = DL.BulkDomain1D_List.back();
    //porousLiKCl_FeS2Cathode_dom1D *cc_ptr = dynamic_cast<porousLiKCl_FeS2Cathode_dom1D *>(cd_ptr);
    double capacityZeroDoD, spec_capacityZeroDoD;
    double dischargedCapacity, spec_dischargedCapacity;
    cd_ptr->reportSolutionParam( "CapacityZeroDoD", &capacityZeroDoD );
    cd_ptr->reportSolutionParam( "DepthOfDischarge", &dischargedCapacity );
    cd_ptr->reportSolutionParam( "SpecificCapacityZeroDoD", &spec_capacityZeroDoD );
    cd_ptr->reportSolutionParam( "SpecificDepthOfDischarge", &spec_dischargedCapacity );

    FILE *fp = 0;
    if (ievent == 0) {
      fp = fopen("timeDep_IV.dat", "w");
      fprintf(fp, "TITLE = \"Time versus Current or Voltage\"\n");
      fprintf(fp, "VARIABLES = \" T [s]\"\n");
      fprintf(fp, "\"Voltage [volts] \"\n");
      fprintf(fp, "\"CathodeCurrent [mA/cm2]\"\n");
      fprintf(fp, "\"Initial Specific Cathode Capacity [mA-hr/g] \"\n");
      fprintf(fp, "\"Discharged Specific Cathode Capacity [mA-hr/g] \"\n");
      fprintf(fp, "\"Initial Cathode Capacity [A-s/m2] \"\n");
      fprintf(fp, "\"Discharged Cathode Capacity [A-s/m2] \"\n");
      fprintf(fp, "\"AnodeCurrent [mA/cm2]\"\n");
    } else {
      fp = fopen("timeDep_IV.dat", "a");
    }
    fprintf(fp, "   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E,   %15.5E \n", 
	    time_current, phiCath, 0.1 * c_icurr, spec_capacityZeroDoD/3.6, spec_dischargedCapacity/3.6, 
	    capacityZeroDoD, dischargedCapacity, 0.1 * a_icurr );
    fflush(fp);
    fclose(fp);
    }
    Comm_ptr->Barrier();
  }
//================================================================================================================================
double
FlatBatteryResidEval::reportCathodeVoltage() const {
    DomainLayout &DL = *DL_ptr_;
    // we want the last surface, but be careful when we go to double tap batteries
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_FlatFeS2Cathode* c_ptr = dynamic_cast<SurDomain_FlatFeS2Cathode*>(d_ptr);
    // might have to update the SurDomain.
    double phi =  c_ptr->phiCathode_;
    
    return phi;
} 
//====================================================================================================================
double
FlatBatteryResidEval::reportCathodeCurrent() const {
    DomainLayout &DL = *DL_ptr_;
  
    SurDomain1D *d_ptr = DL.SurDomain1D_List.back();
    SurDomain_FlatFeS2Cathode* c_ptr = dynamic_cast<SurDomain_FlatFeS2Cathode*>(d_ptr);
    double icurr = c_ptr->icurrCollector_;
    return icurr;
}


//====================================================================================================================
//=====================================================================================================================
} // end of m1d namespace
//=====================================================================================================================
