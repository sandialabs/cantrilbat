/*
 * $Id: InterfacialMassTransfer_PseudoSS.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



#include "mdp_allo.h"
#include "cantera/equilibrium.h"

#include "cantera/solvers.h"
#include "cantera/kinetics/ImplicitSurfChem.h"
#include "ReactingSurDomain.h"

#include "PhaseList.h"


//#include "BlockEntryGlobal.h"


#include "InterfacialMassTransfer_PseudoSS.h"

#include "ApplBase_print.h"

using namespace Cantera;
using namespace std;



#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif  

namespace Cantera {
  //======================================================================================================================
  /*
   *  Interface constructor
   *
   *  We initialize the arrays in the structure to the appropriate sizes.
   *  And, we initialize all of the elements of the arrays to defaults.  
   */
  InterfacialMassTransfer_PseudoSS::InterfacialMassTransfer_PseudoSS() :
    InterfacialMassTransfer_Integrator()
  {
  
  }
  //======================================================================================================================
  // Copy Constructor
  /*
   * @param right Object to be copied
   */
  InterfacialMassTransfer_PseudoSS::InterfacialMassTransfer_PseudoSS(const InterfacialMassTransfer_PseudoSS &right) :
    InterfacialMassTransfer_Integrator()
  {
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
  }
  //======================================================================================================================
  // Assignment operator
  /*
   *  @param right object to be copied
   */
  InterfacialMassTransfer_PseudoSS & InterfacialMassTransfer_PseudoSS::operator=(const InterfacialMassTransfer_PseudoSS &right)
  {
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;

    InterfacialMassTransfer_Integrator::operator=(right);
    /*
     * We do a straight assignment operator on all of the
     * data. The vectors are copied.
     */
  
 
    /*
     * Return the reference to the current object
     */
    return *this;
  }
  //======================================================================================================================
  /*
   *  destructor
   *
   * We need to manually free all of the arrays.
   */
  InterfacialMassTransfer_PseudoSS::~InterfacialMassTransfer_PseudoSS() 
  {
     delete isc_prob;
     isc_prob = 0; 
  }
  //======================================================================================================================
  //  Setup the electrode
  /*
   * @param ei    ELECTRODE_KEY_INPUT pointer object
   */
  int 
  InterfacialMassTransfer_PseudoSS::model_create(IMT_KEY_INPUT *ei) {
    
    InterfacialMassTransfer_Integrator::model_create(ei);


    /*
     *  create and initialize the implicit surface chemistry problem.
     *   std::vector<ReactingSurDomain *> RSD_List_
     *      -> For some reason, I can't just send RSD_List_ into ImplicitSurfChem directly- I don't know why
     */
    std::vector<InterfaceKinetics *> IK_List;
    for (size_t i = 0; i < RSD_List_.size(); i++) {
      InterfaceKinetics *ik_ptr = dynamic_cast<InterfaceKinetics *>(RSD_List_[i]);
      IK_List.push_back(ik_ptr);
    }
    isc_prob = new Cantera::ImplicitSurfChem(IK_List);

    neq_isc_prob_ = isc_prob->neq();

    return 0;
  }
  //======================================================================================================================
  //  Set the initial conditions from the input file.
  /*   
   *   (virtual from InterfacialMassTransfer)
   *   (This is a serial virtual function or an overload function)
   *
   *    This is one of the most important routines. It sets up the initial conditions of the interface
   *    from the input file. The interface object itself has been set up from a call to model_create().
   *    After the call to this routine, the interface should be internally ready to be integrated and reacted. 
   *    It takes its input from an IMT_KEY_INPUT object which specifies the setup of the interface
   *    object and the initial state of that object.
   *    
   *    The routine works like an onion initialization. The parent object is initialized before the 
   *    child. This means the child object first calls the parent, before it does its own initializations.
   * 
   * @param ei    IMT_KEY_INPUT pointer object
   *  
   *  @return  Returns zero if successful, and -1 if not successful.
   */
  int InterfacialMassTransfer_PseudoSS::setInitialConditions(IMT_KEY_INPUT *ei)
  {
    return InterfacialMassTransfer_Integrator::setInitialConditions(ei);
  }
  //====================================================================================================================
  // The internal state of the electrode must be kept for the initial and 
  // final times of an integration step.
  /*
   *  This function advances the initial state to the final state that was calculated
   *  in the last integration step. If the initial time is input, then the code doesn't advance
   *  or change anything.
   *
   * @param Tinitial   This is the New initial time. This time is compared against the "old"
   *                   final time, to see if there is any problem.
   */
  void  InterfacialMassTransfer_PseudoSS::resetStartingCondition(double Tinitial) {
    int i;

    if (pendingIntegratedStep_ != 1) {
#ifdef DEBUG_ELECTRODE
      // printf(" InterfacialMassTransfer_PseudoSS::resetStartingCondition WARNING: resetStartingCondition called with no pending integration step\n");
#endif
        return;
    }

    double tbase = MAX(t_init_init_, 1.0E-50);
    if (fabs(Tinitial - t_init_init_) < (1.0E-9 * tbase)) {
      return;
    }

    tbase = MAX(Tinitial, tbase);
    tbase = MAX(tbase, t_final_final_);
    if (fabs(Tinitial - t_final_final_) > (1.0E-9 * tbase)) {
      throw CanteraError("InterfacialMassTransfer_PseudoSS::resetStartingCondition()", "tinit " + fp2str(Tinitial) +" not compat with t_final_final_ "
			 + fp2str(t_final_final_));
    }



    t_init_init_ = Tinitial;

    // reset surface quantities
    for (i = 0; i < m_NumSurPhases; i++) {
      surfaceAreaRS_init_[i] = surfaceAreaRS_final_[i];
    }

    // Reset total species quantities
    for (int k = 0; k < m_NumTotSpecies; k++) {
      spMoles_init_[k] = spMoles_final_[k];
      spMoles_init_init_[k] = spMoles_final_[k];
    }

    mdpUtil::mdp_zero_dbl_1(DATA_PTR(spMoleIntegratedSourceTerm_), m_NumTotSpecies);
    mdpUtil::mdp_zero_dbl_1(DATA_PTR(spMoleIntegratedSourceTermLast_), m_NumTotSpecies);

    // Reset the total phase moles quantities
    for (i = 0; i < m_NumTotPhases; i++) {
      phaseMoles_init_[i] = phaseMoles_final_[i];
      phaseMoles_init_init_[i] = phaseMoles_final_[i];
    }


    /*
     *  Change the initial subcycle time delta here. Note, we should not change it during the integration steps
     *  because we want jacobian calculations to mainly use the same time step history, so that the answers are
     *  comparible irrespective of the time step truncation error.
     */;
    if (deltaTsubcycle_init_next_ < 1.0E299) {
      deltaTsubcycle_init_init_ = deltaTsubcycle_init_next_;
    }
    deltaTsubcycle_init_next_ = 1.0E300;

    pendingIntegratedStep_ = 0;
  }
   // --------------------------- GET MOLE NUMBERS ------------------------------------------------



    // --------------------------- GET THERMO AND VOLTAGE  ------------------------------------------------

 

    // -------------------------  GET VOLUMES -----------------------------------------------------------

    // ---------------------- GET SURFACE AREAS -------------------------------------------------------



    // --------------------------- INTERNAL UPDATE FUNCTIONS --------------------------------------

  //====================================================================================================================
  // Take the state (i.e., the final state) within the InterfacialMassTransfer_Model and push it down
  // to the ThermoPhase Objects
  /*
   *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
   *  objects
   */
  void InterfacialMassTransfer_PseudoSS::updateState()
  {
    InterfacialMassTransfer::updateState();

    /*
     *  Solve the pseudo state problem
     */
    int ifuncOverride = 0;
    
    solvePseudoSteadyStateProblem(ifuncOverride, 0.0);
  } 
  //====================================================================================================================
  //! Carry out the pseudo steady state calculation
  /*!
   *
   */
  void InterfacialMassTransfer_PseudoSS::solvePseudoSteadyStateProblem(int ifuncOverride, doublereal timeScaleOverride)
  {
    isc_prob->solvePseudoSteadyStateProblem(ifuncOverride,timeScaleOverride);

  }
  //====================================================================================================================
  void InterfacialMassTransfer_PseudoSS::printInterfacialMassTransfer_List(int pSrc, bool subTimeStep) {
    printf("\n");
    printf("         PName                   MoleNum      molarVol     Volume       FractVol     Voltage   \n");
    printf("     ============================================================================================\n");
    double egv = TotalVol();
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
      std::string pname = PhaseNames_[iph];
      ThermoPhase &tp = thermo(iph);
      double mv = tp.molarVolume();
      double pv = mv * phaseMoles_final_[iph];
      printf("     ");
      ca_ab::pr_sf_lj(pname, 24, 1);
      printf(" %12.3E", phaseMoles_final_[iph]);
      printf(" %12.3E", mv);
      printf(" %12.3E", pv);
      printf(" %12.3E", pv / egv);
  
      printf("\n");
    }
    printf("     ============================================================================================\n");


  }
  //====================================================================================================================
  // Print conditions of the electrode for the current integration step to stdout
  /*
   *  @param pSrc          Print Source terms that have occurred during the step from the initial_initial
   *                       to the final_final time.
   *                       The default is to print out the source terms
   *  @param  subTimeStep  Print out conditions from the most recent subTimeStep and not the global
   *                       time step. The default is to print out the global values 
   */
  void InterfacialMassTransfer_PseudoSS::printInterfacialMassTransfer(int pSrc, bool subTimeStep) {
    int iph;
    double egv = TotalVol();
    printf("   ==============================================================================================\n");
    if (subTimeStep) {
      printf("      InterfacialMassTransfer_PseudoSS at intermediate-step time final = %g\n", t_final_);
      printf("                   intermediate-step time init  = %g\n", t_init_);
    } else {
      printf("      InterfacialMassTransfer_PseudoSS at time final = %g\n", t_final_final_);
      printf("                   time init  = %g\n", t_init_init_);
    }
    printf("\n");
    printf("                    DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n", 
	   DomainNumber_, CellNumber_, counterNumberIntegrations_);
    printf("   ==============================================================================================\n");
    printf("          Number of surfaces = %d\n", numSurfaces_);
    printf("          Total Volume = %10.3E\n", egv);
    printf("          Temperature = %g\n", Temp_);
    printf("          Pressure = %g\n", Pres_A_Interface_final_);


    printInterfacialMassTransfer_List(pSrc, subTimeStep);

    for (iph = 0; iph < m_NumTotPhases; iph++) {
      printInterfacialMassTransfer_Phase(iph, pSrc, subTimeStep);   
    }
  }
  //===================================================================================================================
 
  void InterfacialMassTransfer_PseudoSS::printInterfacialMassTransfer_Phase(int iph, int pSrc, bool subTimeStep) {
    int isph;
    double *netROP = new double[m_NumTotSpecies];
    ThermoPhase &tp = thermo(iph);
    string pname = tp.id();
    int istart = m_PhaseSpeciesStartIndex[iph];
    int nsp = tp.nSpecies();
    if (printLvl_ <= 1) {
      return;
    }
    printf("     ============================================================================================\n");
    printf("          Phase %d %s \n", iph,pname.c_str() );
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
 

    /*
     * Do specific surface phase printouts
     */
    if (iph >= NumVolPhases_) {
      isph = iph - NumVolPhases_;
      printf("                surface area (final) = %g\n",  surfaceAreaRS_final_[isph]);
      printf("                surface area (init)  = %g\n",  surfaceAreaRS_init_[isph]);
    }
    if (printLvl_ >= 3) {
      printf("\n");
      printf("                Name              MoleFrac_final kMoles_final kMoles_init SrcTermIntegrated(kmol)\n");
      for (int k = 0; k < nsp; k++) {
	string sname = tp.speciesName(k);
	if (pSrc) {
	  if (subTimeStep) {
	    printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
		   spMoles_final_[istart + k], spMoles_init_[istart + k],
		   spMoleIntegratedSourceTermLast_[istart + k]);
	  } else {
	    printf("                %-22s %10.3E %10.3E   %10.3E  %10.3E\n", sname.c_str(), spMf_final_[istart + k],
		   spMoles_final_[istart + k], spMoles_init_init_[istart + k],
		   spMoleIntegratedSourceTerm_[istart + k]);
	  }
	} else {
	  if (subTimeStep) {
	    printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
		   spMoles_final_[istart + k], spMoles_init_[istart + k]);
	  } else {
	    printf("                %-22s %10.3E %10.3E   %10.3E\n", sname.c_str(), spMf_final_[istart + k],
		   spMoles_final_[istart + k], spMoles_init_init_[istart + k]);
	  }
	}
      }
    }
    if (printLvl_ >= 4) {
      if (iph >= NumVolPhases_) {
	const vector<double> &rsSpeciesProductionRates = RSD_List_[isph]->calcNetProductionRates();
	RSD_List_[isph]->getNetRatesOfProgress(netROP);
      
	doublereal * spNetProdPerArea = (doublereal *) spNetProdPerArea_List_.ptrColumn(isph);
	mdpUtil::mdp_zero_dbl_1(spNetProdPerArea, m_NumTotSpecies);
	int nphRS = RSD_List_[isph]->nPhases();
	int kIndexKin = 0;
	for (int kph = 0; kph < nphRS; kph++) {
	  int jph = RSD_List_[isph]->kinOrder[kph];
	  int istart = m_PhaseSpeciesStartIndex[jph];
	  int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	  for (int k = 0; k < nsp; k++) {
	    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
	    kIndexKin++;
	  }
	}
	printf("\n");
	printf("                           spName                  SourceRateLastStep (kmol/m2/s) \n");
	for (int k = 0; k <  m_NumTotSpecies; k++) {
	  string ss = speciesName(k);
	  printf("                           %-22s %10.3E\n", ss.c_str(), spNetProdPerArea[k]);
	}
      }
    }
    printf("     ============================================================================================\n");
    delete [] netROP;

  }


  //====================================================================================================================
  // Write out CSV tabular data on the integrations
  /*
   *  The idea is to print out tabular data about each intermediate run and about each
   *  global run
   *
   *  @param itype Type of the data
   *            - 0      Initialization information
   *            - 1      Normal intermediate information
   *            - 2      Normal global information at the end of a global time step
   *            - -1     Failed intermediate calculation, a failure from a nonlinear solver step
   *            - -2     Failed calculation from the predictor step - not necessarily significant.
   */
  void InterfacialMassTransfer_PseudoSS::writeCSVData(int itype) { 
    if (printLvl_ < 2) return;
    int k;
    static std::string globOutputName;
    static std::string intOutputName;
    static FILE *fpG = 0;
    static FILE *fpI = 0;

    if (itype == 0 || fpI == 0) {
      intOutputName = (Title_ + "_intResults_" + int2str(DomainNumber_) + "_" +
		       int2str(CellNumber_) + ".csv");
       fpI = fopen(intOutputName.c_str(), "w");

       fprintf(fpI, "         Tinit ,        Tfinal ," );

       fprintf(fpI, "      Volts_Soln ,  Volts_InterfacialMassTransfer_PseudoSS ,");

       fprintf(fpI, "      Current ,");

       fprintf(fpI, "  CapDischarged ,");

       for (k = 0; k < m_NumTotSpecies; k++) {
	 string sss = speciesName(k);
	 fprintf(fpI, " MN_%-20.20s,",  sss.c_str());
       }

       for (k = 0; k < m_NumTotSpecies; k++) {
	 string sss = speciesName(k);
	 fprintf(fpI, " SRC_%-20.20s,",  sss.c_str());
       }


       fprintf(fpI, " iType");
       fprintf(fpI, "\n");
       fclose(fpI);
    }

   if (itype == 0 || fpG == 0) {
     globOutputName =( Title_ + "_globalResults_" + int2str(DomainNumber_) + "_" +
		       int2str(CellNumber_) + ".csv");
       fpG = fopen(globOutputName.c_str(), "w");
       fprintf(fpG, "         Tinit ,        Tfinal ," );

       fprintf(fpG, "      Volts_Soln ,  Volts_InterfacialMassTransfer_PseudoSS ,");

       fprintf(fpG, "      Current ,");

       fprintf(fpG, "  CapDischarged ,");

       for (k = 0; k < m_NumTotSpecies; k++) {
	 string sss = speciesName(k);
	 fprintf(fpG, " MN_%-20.20s,",  sss.c_str());
       }

       for (k = 0; k < m_NumTotSpecies; k++) {
	 string sss = speciesName(k);
	 fprintf(fpG, " SRC_%-20.20s,",  sss.c_str());
       }

       fprintf(fpG, " iType");
       fprintf(fpG, "\n");
       fclose(fpG);
    }



   if (itype == 1 || itype == 2 || itype == -1) {
     fpI = fopen(intOutputName.c_str(), "a");
     fprintf(fpI,"  %12.5E ,  %12.5E ,", t_init_, t_final_);
     
 
      for (k = 0; k < m_NumTotSpecies; k++) {
	fprintf(fpI, " %12.5E           ,",  spMoles_final_[k]);
      }

      for (k = 0; k < m_NumTotSpecies; k++) {
	fprintf(fpI, " %12.5E            ,",  spMoleIntegratedSourceTermLast_[k]); 
      }
      
      fprintf(fpI, "   %3d \n", itype);
      fclose(fpI);
    }


    if (itype == 2) {
      fpG = fopen(globOutputName.c_str(), "a");
      fprintf(fpG,"  %12.5E ,  %12.5E ,", t_init_init_, t_final_final_);

 

  
      for (k = 0; k < m_NumTotSpecies; k++) {	
	fprintf(fpG, " %12.5E           ,",  spMoles_final_[k]);
      }

      for (k = 0; k < m_NumTotSpecies; k++) {
	fprintf(fpG, " %12.5E            ,",  spMoleIntegratedSourceTerm_[k]); 
      }
      
      fprintf(fpG, "   %3d \n", itype);
      fclose(fpG);
    }


  }
 

  
  //====================================================================================================================
} // End of namespace Cantera
//======================================================================================================================
