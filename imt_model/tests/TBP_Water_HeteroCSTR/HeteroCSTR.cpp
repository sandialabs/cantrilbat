/*
 * $Id: AtoB_1.cpp 222 2012-06-26 21:55:32Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */



#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/numerics/ResidEval.h"
#include "cantera/numerics/NonlinearSolver.h"
#include "cantera/numerics/CVodesIntegrator.h"

#include "cantera/kinetics/AqueousKinetics.h"

#include "InterfacialMassTransfer_input.h"
#include "InterfacialMassTransfer.h"
#include "InterfacialMassTransfer_1to1Distrib.h"
#include "imtPSS_NoSurf.h"
#include "imtPSS_NoSurf.h"
#include "ExtraGlobalRxn.h"
#include "RxnMolChange.h"
#include "BlockEntry.h"

#include <stdio.h>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace Cantera;
using namespace mdpUtil;


// a lvl of one prints out the .csv file
int mpequil_debug_print_lvl = 1;
int VCS_Debug_Print_Lvl = 3;

void printUsage() {
    cout << "usage: electrodeCell [-h] [-help_cmdfile] [-d #] [mpequil.inp]"
         <<  endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << "  electrodeCell.inp    : command file" << endl;
    cout << "                     : (if missing, assume mpequil.inp)" 
	 << endl;
    cout << endl;
}

namespace Cantera {

/**
 *  Virtual base class for ODE right-hand-side function evaluators.
 *  Classes derived from FuncEval evaluate the right-hand-side function
 * \f$ \vec{F}(t,\vec{y})\f$ in 
 * \f[
 *  \dot{\vec{y}} = \vec{F}(t,\vec{y}).
 * \f] 
 *  @ingroup odeGroup 
 */
  class HeteroCSTRFunc : public FuncEval {

  public:

    HeteroCSTRFunc() :
      FuncEval(),
      setOnce(0),
      t_init_(0.0),
      t_final_(0.0),
      deltaT(0.0),
      tpA_ptr_(0),
      tpB_ptr_(0),
      tbpKinetics_(0)
    {
      Temp_ = 298.15;
      Pres_ = OneAtm;
    }

    virtual ~HeteroCSTRFunc() 
    {
    }

    void popPointers(ThermoPhase *tpA_ptr, ThermoPhase *tpB_ptr, InterfacialMassTransfer *iface, 
		     AqueousKinetics *tbpKinetics)
    {
      tpA_ptr_ = tpA_ptr;
      tpB_ptr_ = tpB_ptr;
      iface_ = iface;
      tbpKinetics_ = tbpKinetics;
    }

    void setup()
    {
      nSpeciesA_ = tpA_ptr_->nSpecies();

      nSpeciesB_ = tpB_ptr_->nSpecies();
      solnAPhase_ = iface_->solnAPhase_;
      solnBPhase_ = iface_->solnBPhase_;
      neq_ =  nSpeciesA_ +  nSpeciesB_;
      molNumVector_.resize(neq_);
      double vec[100];
      iface_->getMoleNumSpecies(vec);

      for (int k = 0; k < nSpeciesA_; k++) {
	molNumVector_[k] = vec[k];
      }
      netTBP.resize(nSpeciesB_, 0.0);

      int bstart = iface_->getGlobalSpeciesIndex(solnBPhase_);

      for (int k = 0; k < nSpeciesB_; k++) {
	molNumVector_[k + bstart] = vec[k + bstart];
      }

 
      
    }


    /**
     * Evaluate the right-hand-side function. Called by the
     * integrator.
     * @param t time. (input)
     * @param y solution vector. (input)
     * @param ydot rate of change of solution vector. (output)
     * @param p parameter vector
     */
    virtual void eval(double t, double* y, double* ydot, double* p)
    {
      if (t > t_init_) {
	setOnce = 1;
	t_final_ = t;
      } 
      if (t < t_init_) {
	printf("confused t_init = %g, t_final = %g, t = %g\n", t_init_, t_final_, t);
	exit(-1);
      }
      deltaT = t_final_ - t_init_;

      iface_->setSolnA_BoundaryConcentrations(y, Temp_, Pres_);

      int bstart = iface_->getGlobalSpeciesIndex(solnBPhase_);
      iface_->setSolnB_BoundaryConcentrations(&y[bstart], Temp_, Pres_);
   
      // iface_->integrate(deltaT);
    
      doublereal net[30];
      //   iface_->getIntegratedProductionRates(net);

      iface_->getNetProductionRates(net);

      tbpKinetics_->getNetProductionRates(DATA_PTR(netTBP));

      double Bmoles = iface_->phaseMoles(solnBPhase_);
      ThermoPhase *bPhase = & iface_->thermo(solnBPhase_);
      double molV = bPhase->molarVolume();
      double bVol = Bmoles / molV;

      /*
       *  Convert from kmol m-3 s-1   to kmol s-1;
       */
      for (int k = 0; k <  nSpeciesB_; k++) {
	netTBP[k] *= bVol;
      }
 
      
      for (int index = 0; index < neq_; index++) {
	ydot[index] = net[index];

      }



      for (int k = 0; k < nSpeciesB_; k++) {
	ydot[k + bstart] += netTBP[k];
      }

    }

    /**
     * Fill the solution vector with the initial conditions 
     * at initial time t0.
     */
    virtual void getInitialConditions(double t0, size_t leny, double* y)
    {
      for (size_t i = 0; i < leny; i++) {
	y[i] = molNumVector_[i];
      }
    }

    /** 
     * Number of equations. 
     */
    virtual size_t neq()
    {
      return neq_;
    }

    //! Number of parameters.
    virtual size_t nparams() { 
      return 0;
    }

  protected:

    int setOnce;

 
    double t_init_;

    double t_final_;

    double deltaT;

    ThermoPhase *tpA_ptr_;

    ThermoPhase *tpB_ptr_;

    int nSpeciesA_;
    int nSpeciesB_;
    int solnAPhase_;
    int solnBPhase_;
    int neq_;

    std::vector<double> molNumVector_;


    InterfacialMassTransfer *iface_;

    AqueousKinetics *tbpKinetics_;

    std::vector<double> netTBP;

    double Temp_;
    double Pres_;
  private:
  };



  class HeteroCSTRInteg : public CVodesIntegrator {
  
public:
  HeteroCSTRInteg() :
    CVodesIntegrator()
  {
    
  }
  
  
  
  
};

}
//=====================================================================================================


int main(int argc, char **argv)
{
  int ip1 = 0;
  int retn = 0;
  //bool doCathode = false;
  string commandFileNet="";
  bool printInputFormat = false; // print cmdfile.txt format
  bool printedUsage = false; // bool indicated that we have already
  // printed usage

  VCSnonideal::vcs_timing_print_lvl = 0;
  NonlinearSolver::s_TurnOffTiming = true;
  NonlinearSolver::s_print_NumJac = true;
  HeteroCSTRFunc hfunc;
  CVodesIntegrator *  integ = new CVodesIntegrator();

  /*
   * Process the command line arguments
   */ 
  if (argc > 1) {
    string tok;
    for (int j = 1; j < argc; j++) {
      tok = string(argv[j]);
      if (tok[0] == '-') {
	int nopt = static_cast<int>(tok.size());
	for (int n = 1; n < nopt; n++) {
	  if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
	    printInputFormat = true;
	  } else if (tok[n] == 'h') {
	    printUsage();
	    printedUsage = true;
	    exit(1);
	  } else if (tok[n] == 'd') {
	    int lvl = 2;
	    if (j < (argc - 1)) {
	      string tokla = string(argv[j+1]);
	      if (strlen(tokla.c_str()) > 0) {
		lvl = atoi(tokla.c_str());
		n = nopt - 1;
		j += 1;
		if (lvl >= 0 && lvl <= 1000) {
		  if (lvl == 0) ip1 = 0;
		  else          ip1 = lvl; 
		  mpequil_debug_print_lvl = lvl;
		}
	      }  
	    }
	  } else {
	    printUsage();
	    printedUsage = true;
	    exit(1);
	  }
	}
      } else if (commandFileNet == "") {
	commandFileNet = tok;
      } else {
	printUsage();
	printedUsage = true;
	exit(1);
      }
    }
  }


  try {

  
    //  retn = cell_input(commandFileNet);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }
  
  
    Cantera::imtPSS_NoSurf* iface  = new Cantera::imtPSS_NoSurf();
    
    Cantera::IMT_KEY_INPUT *face_input = new IMT_KEY_INPUT();
    
    std::string commandFile = "interface.inp";
  
    /**	
     * Initialize a block input structure for the command file
     */
    BEInput::BlockEntry *cf = new BEInput::BlockEntry("command_file");

    /*	
     * Go get the problem description from the input file
     */
    retn = imt_input(face_input, commandFile, cf);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

 
    retn = iface->model_create(face_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }


    retn = iface->setInitialConditions(face_input);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    ThermoPhase *tpA_ptr = & iface->thermo(iface->solnAPhase_);
    ThermoPhase *tpB_ptr = & iface->thermo(1);

    AqueousKinetics *tbpKinetics = new AqueousKinetics();

    std::vector<ThermoPhase *> aqList;
    aqList.push_back(tpB_ptr);

    XML_Node *xB = iface->volPhaseXMLNode(iface->solnBPhase_);

    importKinetics(*xB, aqList, tbpKinetics);


    hfunc.popPointers(tpA_ptr, tpB_ptr, iface, tbpKinetics);

    hfunc.setup();


    double Temp_ = 298.15;
    double Pres = OneAtm;
    iface->setState_TP(Temp_, Pres);
    integ->initialize(0.0, hfunc);

    integ->setTolerances(1.0E04, 1.0E-16);

    std::vector<double> molNum(30);
    std::vector<double> ydot(30);
    int nT = 500;
    double Tinitial = 0.0;
    double Tfinal = 0.0;
    double Tout = 0.0;
    for (int itimes = 0; itimes < nT; itimes++) {
      iface->resetStartingCondition(Tout);
      Tout = integ->step(10.);

      iface->getMoleNumSpecies(DATA_PTR(molNum));
      Tinitial = Tfinal;


      double *y = integ->solution();

      /*
       *  It seems that the current solution is never in the iface object when returned by cvode
       *  This is a serious issue that indicates that the residual isn't being called at the final value of y.
       *  So, we need to do this manually.
       */
      hfunc.eval(Tout, y, DATA_PTR(ydot), 0);
  

   
      Tfinal = Tout;

      iface->setStateFF_time(Tfinal);

      iface->setFinalFinalStateFromFinal();
      iface->calcIntegratedSourceTerm();

      iface->getMoleNumSpecies(DATA_PTR(molNum));
      doublereal net[30];
      iface->getNetProductionRates(net);
      cout << setw(15) << Tfinal << setw(15) << 0.0 << endl;
      iface->printInterfacialMassTransfer(1, false);
      iface->writeCSVData(2);
    }

    delete cf;
    delete face_input;
    delete iface;
    Cantera::appdelete();

    return retn;

  } catch (CanteraError) {

    showErrors();
    return -1;
  }
} 
