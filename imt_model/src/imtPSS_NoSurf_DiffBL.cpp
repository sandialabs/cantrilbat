/*
 * $Id: InterfacialMassTransfer_1to1Distrib.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



#include "mdp_allo.h"
#include "cantera/equilibrium.h"

#include "cantera/solvers.h"


//#include "PhaseList.h"


#include "imtPSS_NoSurf_DiffBL.h"

#include "ApplBase_print.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif
using namespace std;



#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif  

#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera 
#endif
{
  //======================================================================================================================
  /*
   *  ELECTRODE_INPUT: constructor
   *
   *  We initialize the arrays in the structure to the appropriate sizes.
   *  And, we initialize all of the elements of the arrays to defaults.  
   */
  imtPSS_NoSurf_DiffBL::imtPSS_NoSurf_DiffBL() :
    InterfacialMassTransfer_PseudoSS(),
    cSurfA_(0),
    cSurfB_(0),
    cBoundA_(0),
    cBoundB_(0),
    jKmolfluxA_(0),
    jKmolfluxB_(0), 
    gradX_phaseA_(0),
    gradX_phaseB_(0),
    tranA_(0),
    tranB_(0),
    Velo_A_CVB_init_init_(0.0),
    Velo_A_CVB_init_(0.0),
    Velo_A_CVB_final_(0.0),
    Velo_A_CVB_final_final_(0.0),
    Velo_B_CVB_init_init_(0.0),
    Velo_B_CVB_init_(0.0),
    Velo_B_CVB_final_(0.0),
    Velo_B_CVB_final_final_(0.0),
    grossAVol_(0.0),
    molarVolA_(0.0),
    grossBVol_(0.0),
    molarVolB_(0.0)
  {
    interfaceType_ = PSEUDOSS_NOSURF_IMT;      
  }
  //======================================================================================================================
  // Copy Constructor
  /*
   * @param right Object to be copied
   */
  imtPSS_NoSurf_DiffBL::imtPSS_NoSurf_DiffBL(const imtPSS_NoSurf_DiffBL &right) :
    InterfacialMassTransfer_PseudoSS(),
    cSurfA_(0),
    cSurfB_(0),
    cBoundA_(0),
    cBoundB_(0),
    gradX_phaseA_(0),
    gradX_phaseB_(0),
    tranA_(0),
    tranB_(0),
    Velo_A_CVB_init_init_(0.0),
    Velo_A_CVB_init_(0.0),
    Velo_A_CVB_final_(0.0),
    Velo_A_CVB_final_final_(0.0),
    Velo_B_CVB_init_init_(0.0),
    Velo_B_CVB_init_(0.0),
    Velo_B_CVB_final_(0.0),
    Velo_B_CVB_final_final_(0.0),
    grossAVol_(0.0),
    molarVolA_(0.0),
    grossBVol_(0.0),
    molarVolB_(0.0)
  {
    /*
     * Call the assignment operator.
     */
    *this = operator=(right);
    interfaceType_ = PSEUDOSS_NOSURF_IMT;
  }
  //======================================================================================================================
  // Assignment operator
  /*
   *  @param right object to be copied
   */
  imtPSS_NoSurf_DiffBL & imtPSS_NoSurf_DiffBL::operator=(const imtPSS_NoSurf_DiffBL &right)
  {
    /*ZZ
     * Check for self assignment.
     */
    if (this == &right) return *this;

    InterfacialMassTransfer_PseudoSS::operator=(right);
    /*
     * We do a straight assignment operator on all of the
     * data. The vectors are copied.
     */
    cSurfA_ = right.cSurfA_;
    cSurfB_ = right.cSurfB_;
    cBoundA_ = right.cBoundA_;
    cBoundB_ = right.cBoundB_;
    gradX_phaseA_ = right.gradX_phaseA_;
    gradX_phaseB_ = right.gradX_phaseB_;

    delete tranA_;
    tranA_ = right.tranA_->duplMyselfAsTransport();
    ThermoPhase *tpA = &thermo(solnAPhase_);
    tranA_->setThermo(*tpA);
 

    delete tranB_;
    tranB_ = right.tranB_->duplMyselfAsTransport();
    ThermoPhase *tpB = &thermo(solnBPhase_);
    tranB_->setThermo(*tpB); 
    tranB_->setNDim(1);

    Velo_A_CVB_init_init_ = right.Velo_A_CVB_init_init_;
    Velo_A_CVB_init_ = right.Velo_A_CVB_init_;
    Velo_A_CVB_final_ = right.Velo_A_CVB_final_;
    Velo_A_CVB_final_final_ = right.Velo_A_CVB_final_final_;
    Velo_B_CVB_init_init_ = right.Velo_B_CVB_init_init_;
    Velo_B_CVB_init_ = right.Velo_B_CVB_init_;
    Velo_B_CVB_final_ = right.Velo_B_CVB_final_;
    Velo_B_CVB_final_final_ = right.Velo_B_CVB_final_final_;
    grossAVol_ = right.grossAVol_;
    molarVolA_ = right.molarVolA_;
    grossBVol_ = right.grossBVol_;
    molarVolB_ = right.molarVolB_;
    /*
     * Return the reference to the current object
     */
    return *this;
  }
  //======================================================================================================================
  /*
   *
   *  ELECTRODE_INPUT:destructor
   *
   * We need to manually free all of the arrays.
   */
  imtPSS_NoSurf_DiffBL::~imtPSS_NoSurf_DiffBL() 
  {
  }
  //======================================================================================================================
  // Return the type of interfacial mass transport object 
  /*
   *  Returns the enum type of the object. This is used in the factory routine.
   *
   *  @return Returns an enum type, called   IMT_Types_Enum
   */
  IMT_Types_Enum imtPSS_NoSurf_DiffBL::imtType() const
  {
    return  PSEUDOSS_NOSURF_IMT;
  }
  //======================================================================================================================
  // Set the electrode ID information
  void imtPSS_NoSurf_DiffBL::setID(int domainNum, int cellNum)
  {
 
  }

  //======================================================================================================================
  //  Setup the electrode
  /*
   * @param ei    ELECTRODE_KEY_INPUT pointer object
   */
  int 
  imtPSS_NoSurf_DiffBL::model_create(IMT_KEY_INPUT *ei) {

    /*
     *  For this class, A and B moles are part of the solution vector
     */
    includeABMolesInSolnVector_ = 1;

    InterfacialMassTransfer_PseudoSS::model_create(ei);


    cSurfA_.resize(nSpeciesA_, 0.0);
    cSurfB_.resize(nSpeciesB_, 0.0);
    cBoundA_.resize(nSpeciesA_, 0.0);
    cBoundB_.resize(nSpeciesB_, 0.0);
    jKmolfluxA_.resize(nSpeciesA_, 0.0);
    jKmolfluxB_.resize(nSpeciesB_, 0.0);

    gradX_phaseA_.resize(nSpeciesA_, 0.0);
    gradX_phaseB_.resize(nSpeciesB_, 0.0);

    ThermoPhase *tpA = & thermo(solnAPhase_);

    XML_Node* xmlA = &( volPhaseXMLNode(solnAPhase_) );

    std::string transportModel;
    if (xmlA->hasChild("transport")) {   
      const XML_Node& tranNode = xmlA->child("transport");
      transportModel = tranNode["model"];

      tranA_ = newTransportMgr(transportModel, tpA, 0);
    } else {
      throw CanteraError("  imtPSS_NoSurf_DiffBL::model_create()",
			 "Phase A needs a transport model");
    }
    tranA_->setNDim(1);

    ThermoPhase *tpB = & thermo(solnBPhase_);

    XML_Node* xmlB = &( volPhaseXMLNode(solnBPhase_) );

    std::string transportModelB;
    if (xmlB->hasChild("transport")) {   
      const XML_Node& tranNode = xmlB->child("transport");
      transportModelB = tranNode["model"];

      tranB_ = newTransportMgr(transportModelB, tpB, 0);
    } else {
      throw CanteraError("  imtPSS_NoSurf_DiffBL::model_create()",
			 "Phase B needs a transport model");
    }
    tranB_->setNDim(1);

    /*
     *  set the number of equations
     *       time
     *       v_a  and v_b and v_s
     */
    neq_ = 1 + 2 + nSpeciesA_ + nSpeciesB_;


    return 0;
  }
  //======================================================================================================================
  //!  Set the initial conditions from the input file.
  /*!   
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
  int imtPSS_NoSurf_DiffBL::setInitialConditions(IMT_KEY_INPUT *ei)
  {
    int retn = InterfacialMassTransfer_PseudoSS::setInitialConditions(ei);

  
    AssertThrow(BLThickness_A_final_ > 0.0, "imtPSS_NoSurf_DiffBL::setInitialConditions()");
    AssertThrow(BLThickness_B_final_ > 0.0, "imtPSS_NoSurf_DiffBL::setInitialConditions()");

    for (size_t k = 0; k < nSpeciesA_; k++) {
      spMF_solnA_BC_final_[k] = ei->XmfPhaseA_[k]; 
    }
    int astart = 0;
    for (size_t k = 0; k < nSpeciesA_; k++) {
      spMf_final_[k + astart] = ei->XmfPhaseA_[k]; 
    }
    ThermoPhase *tpA = & thermo(solnAPhase_);
    tpA->setState_TPX(Temp_, Pres_A_Interface_final_,  &spMf_final_[astart]);
    double mv_A = tpA->molarVolume();
    double grossAVol = surfaceAreaRS_final_[0] * BLThickness_A_final_ * 0.5;
    
    for (size_t k = 0; k < nSpeciesA_; k++) {
      spMoles_final_[k + astart] = spMf_final_[k + astart] * grossAVol / mv_A;
    }


    for (size_t k = 0; k < nSpeciesB_; k++) {
      spMF_solnB_BC_final_[k] = ei->XmfPhaseB_[k]; 
    }
    int bstart =  m_PhaseSpeciesStartIndex[solnBPhase_];
    for (size_t k = 0; k < nSpeciesB_; k++) {
      spMf_final_[k + bstart] = ei->XmfPhaseB_[k]; 
    }
    ThermoPhase *tpB = & thermo(solnBPhase_);
    tpB->setState_TPX(Temp_, Pres_B_Interface_final_,  &spMf_final_[bstart]);
    double mv_B = tpB->molarVolume();
    double grossBVol = surfaceAreaRS_final_[0] * BLThickness_B_final_ * 0.5;
    
    for (size_t k = 0; k < nSpeciesB_; k++) {
      spMoles_final_[k + bstart] = spMf_final_[k + bstart] * grossBVol / mv_B;
    }
    updateState();
    setInitStateFromFinal(true);
    return retn;
  }
  // ---------------------------------------------------------------------------------------------
  // ----------------------- SPECIFY AND OBTAIN PROBLEM PARAMETERS -------------------------------
  // ---------------------------------------------------------------------------------------------

  // ------------------------------ OBTAIN STATIC PROBLEM INFORMATION ----------------------------


  // ------------------------------ SPECIFY BASIC THERMO CONDITIONS  ------------------------------------
  
  // ------------------------------ SPECIFY PROBLEM PARAMETERS ------------------------------------

  //================================================================================================  
  // Set the mole numbers and internal state of the solnA phase at the final time of the
  // global step.
  /*
   *  We set the mole numbers of the solnA phase separately from 
   *  the rest of the phases.
   *
   *  We always make sure that mole numbers are positive by clipping. We
   *  always make sure that mole fractions sum to one.
   *
   *  If we are not following mole numbers in the electrode, we set the
   *  total moles to the internal constant, electrolytePseudoMoles_, while
   *  using this vector to set the mole fractions.
   *
   * @param solnAMoleNum vector of mole numbers of the species in the
   *                     electrolyte phase. 
   *                     units = kmol
   *                     size = number of species in the electrolyte phase 
   *
   * @param tempA        Temperature of phase A
   *
   * @param presA        Pressure of phase A
   *
   * @param extFieldAValues External field values of A, the first one is usually the voltage.
   *
   * @param setInitial   Boolean indicating that we should set the initial values of the
   *                     solnA instead of the final values.(default false)
   */
  void imtPSS_NoSurf_DiffBL::setSolnA_BoundaryConcentrations(const double * const solnAMoleNum, 
							     double tempA, double presA,
							     const double * extFieldAValues,
							     bool setInitial)
  {
    ThermoPhase &tpA = thermo(solnAPhase_);
    size_t nsp = tpA.nSpecies();
    AssertTrace(nsp == ( m_PhaseSpeciesStartIndex[solnAPhase_+1] - m_PhaseSpeciesStartIndex[solnAPhase_]));

    double sum = 0.0;
    for (size_t k = 0; k < nsp; k++) {
      sum += solnAMoleNum[k];
    }

    if (setInitial) {
      for (size_t k = 0; k < nsp; k++) {
	spMF_solnA_BC_init_init_[k] =  solnAMoleNum[k] / sum;
      }
      Pres_solnA_BC_init_init_ = presA;
    } else {
      for (size_t k = 0; k < nsp; k++) {
	spMF_solnA_BC_final_final_[k] =  solnAMoleNum[k] / sum;
      }
      Pres_solnA_BC_final_final_ = presA;
    }
  }
  //================================================================================================  
  // Set the mole numbers and internal state of the solnB phase at the final time of the
  // global step.
  /*
   *  We set the mole numbers of the solnA phase separately from 
   *  the rest of the phases.
   *
   *  We always make sure that mole numbers are positive by clipping. We
   *  always make sure that mole fractions sum to one.
   *
   *  If we are not following mole numbers in the electrode, we set the
   *  total moles to the internal constant, electrolytePseudoMoles_, while
   *  using this vector to set the mole fractions.
   *
   * @param solnBMoleNum vector of mole numbers of the species in the
   *                     electrolyte phase. 
   *                     units = kmol
   *                     size = number of species in the electrolyte phase 
   *
   * @param tempB        Temperature of phase B
   *
   * @param presB        Pressure of phase B
   *
   * @param extFieldAValues External field values of B, the first one is usually the voltage.
   *
   * @param setInitial   Boolean indicating that we should set the initial values of the
   *                     solnB instead of the final values.(default false)
   */
  void imtPSS_NoSurf_DiffBL::setSolnB_BoundaryConcentrations(const double * const solnBMoleNum, double tempB, double presB,
							     const double * extFieldBValues, bool setInitial)
  {
    ThermoPhase &tpB = thermo(solnBPhase_);
    size_t nsp = tpB.nSpecies();
    AssertTrace(nsp == ( m_PhaseSpeciesStartIndex[solnBPhase_+1] - m_PhaseSpeciesStartIndex[solnBPhase_]));

    double sum = 0.0;
    for (size_t k = 0; k < nsp; k++) {
      sum += solnBMoleNum[k];
    }

    if (setInitial) {
      for (size_t k = 0; k < nsp; k++) {
	spMF_solnB_BC_init_init_[k] =  solnBMoleNum[k] / sum;
      }
      Pres_solnB_BC_init_init_ = presB;
    } else {
      for (size_t k = 0; k < nsp; k++) {
	spMF_solnB_BC_final_final_[k] =  solnBMoleNum[k] / sum;
      }
      Pres_solnB_BC_final_final_ = presB;
    }
  }
  // ---------------------------------------------------------------------------------------------
  // ----------------------------- CARRY OUT INTEGRATIONS -----------------------------------------
  // ---------------------------------------------------------------------------------------------
  int  imtPSS_NoSurf_DiffBL::predictSoln_F1()
  {
  /*
     *    Here we essentially duplicate the calcResid logic, using a forwards Euler approach and putting 
     *    the results in soln_predict and spMoles_final_[]
     */

    double advecS_A =  - (Velo_S_Interface_final_) *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecB_A  = - Velo_A_Interface_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    /*
     *  source terms for the A phase
     *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
     *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
     */
    size_t astart = 0;
    double moleASum = 0.0;
    for (size_t k = 1; k < nSpeciesA_; k++)  {
     
      /*
       *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      double delta =  (deltaTsubcycleCalc_ * DspMoles_final_[k + astart]
		+	jKmolfluxA_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_);

      if ((Velo_S_Interface_final_ -Velo_A_Interface_final_) < 0) {
	// We gain material from the surrounding boundary region
	delta += (advecS_A - advecB_A) * (cBoundA_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecA is negative)
	delta += (advecS_A - advecB_A) * ( cSurfA_[k]);
      }


      spMoles_final_[astart + k] = spMoles_init_[astart + k] + delta;
      moleASum += spMoles_final_[astart + k];
    }
    spMoles_final_[astart + 0] = grossAVol_ / molarVolA_ - moleASum; 

    /*
     *  source terms for the A phase
     *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
     *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
     */
    int bstart = m_PhaseSpeciesStartIndex[solnBPhase_];
  
    double advecS_B = Velo_S_Interface_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecB_B  =  Velo_B_Interface_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    
    double moleBSum = 0.0;
    for (size_t k = 1; k < nSpeciesB_; k++)  {
     
      /*

       *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      double  delta =  (deltaTsubcycleCalc_ * DspMoles_final_[k + bstart]
		+	jKmolfluxB_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_);

      if ((Velo_S_Interface_final_ -  Velo_B_Interface_final_) > 0) {
	// We gain material from the surrounding boundary region
	delta += (advecS_B - advecB_B) * ( cBoundB_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecB is negative)
	delta += (advecS_B - advecB_B) * ( cSurfB_[k]);
      }
  

      spMoles_final_[bstart + k] = spMoles_init_[bstart + k] + delta;
      moleBSum += spMoles_final_[bstart + k];
    }
    spMoles_final_[bstart + 0] = grossBVol_ / molarVolB_ - moleBSum; 

    /*
     *  Create a failure if we significantly decrease the moles to near zero. This will never pass
     *  the predictor corrector
     */

    soln_predict_[0] =  deltaTsubcycleCalc_;

    double retn = 1;
    for (size_t k = 0; k < nSpeciesA_; k++)  {
      soln_predict_[1 + astart + k] = spMoles_final_[astart + k];
      if (spMoles_init_[astart + k] > 0.0) {
	if (spMoles_final_[astart + k] < 0.33 * spMoles_init_[astart + k]) {
	  retn = -1;
	}
      }
    }

    for (size_t k = 0; k < nSpeciesB_; k++)  {
      soln_predict_[1+ bstart + k] = spMoles_final_[bstart + k];
      if (spMoles_init_[bstart + k] > 0.0) {
	if (spMoles_final_[bstart + k] < 0.33 * spMoles_init_[bstart + k]) {
	  retn = -1;
	}
      }
    }
    return retn;
  }


  //====================================================================================================================
  // Predict the solution
  /*
   *  (virtual from InterfacialMassTransfer_Integrator)
   *
   * Ok at this point we have a time step 
   * and initial conditions consisting of phaseMoles_init_ and spMF_init_.
   * We now calculate predicted solution components from these conditions.
   * Additionally, we evaluate whether any multispecies phases are going to pop into existence,
   * adding in potential seed values, or die, creating a different equation set and time step equation.
   * We set the phaseExistence flags in the kinetics solver to reflect phase pops.
   *
   * @return   Returns the success of the operation
   *                 1  A predicted solution is achieved 
   *                 2  A predicted solution with a multispecies phase pop is acheived
   *                 0  A predicted solution is not achieved, but go ahead anyway
   *                -1  The predictor suggests that the time step be reduced and a retry occur.
   */
  int  imtPSS_NoSurf_DiffBL::predictSoln()
  {
    if ((size_t) neq_ != 1 + 2 + nSpeciesA_ + nSpeciesB_) {
      printf("we shouldn't be here\n");
    }
    deltaTsubcycleCalc_ = deltaTsubcycle_;
    soln_predict_[0] = deltaTsubcycleCalc_; 
    size_t astart = 0;
    size_t bstart;
    /*
     *  Take the unpacked solution and calculate the consistent and full final state
     */
    updateState();

    /*
     *  Query Cantera for all of the rate information at the final state (and the initial state if we are doing higher order)
     */
    extractInfo();

    /*
     *   We take the ROP_inner_[] and ROP_outer_[] rates of progress, combine them with the surface area calculation,
     *   and the stoichiometric coefficients to calculate the DspMoles_final_[], which is production rate for
     *   all species in the electrode.
     */
    updateSpeciesMoleChangeFinal();

    /*
     *  Update the Stefan Velocity internal variables for the current final conditions
     */
    updateVelocities();

#ifdef DEBUG_METHOD_F1
    int rf1 =  predictSoln_F1();
    return rf1;
#endif
    for (int III = 0; III < 2; III++) {
      /*
       *    Here we essentially duplicate the calcResid logic, using a forwards Euler approach and putting 
       *    the results in soln_predict and spMoles_final_[]
       */

      double advecS_A = - Velo_S_Interface_final_ * deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
      double advecB_A = - Velo_A_CVB_final_       * deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
      /*
       *  source terms for the A phase
       *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
       *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
       */
      astart = 0;
      double moleASum = 0.0;
      for (size_t k = 1; k < nSpeciesA_; k++)  {
     
	/*
	 *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
	 *  accumulation term.
	 */
	double delta =  (deltaTsubcycleCalc_ * DspMoles_final_[k + astart]
			 +	jKmolfluxA_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_);
	if (massFluxA_ < 0.0) {
	  // We gain material from the surrounding boundary region
	  delta += (advecS_A - advecB_A) * (cBoundA_[k]);
	} else {
	  // We lose material to the surrounding boundary region (advecA is negative)
	  delta += (advecS_A - advecB_A) * ( cSurfA_[k]);
	}


	spMoles_final_[astart + k] = spMoles_init_[astart + k] + delta;
	moleASum += spMoles_final_[astart + k];
      }
      spMoles_final_[astart + 0] = grossAVol_ / molarVolA_ - moleASum; 

      updatePhaseNumbers(solnAPhase_); 
      updateState();
      extractInfo(); 
      updateSpeciesMoleChangeFinal();
      updateVelocities();


      //advecS_A    = - (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
      double advecCVB_A  = - Velo_A_CVB_final_         * surfaceAreaRS_final_[0];
      double resid_1 = (densA_final_ - densA_init_) /  deltaTsubcycleCalc_ * grossAVol_ - massFluxA_;
#ifdef DEBUG_MODE_NOT
      printf("resid_1_drho_dt = %13.7E\n",  (densA_final_ - densA_init_) /  deltaTsubcycleCalc_ * grossAVol_ );
      printf("          densA_final_ = %16.10E",  densA_final_);
      for (size_t ii = 0; ii < nSpeciesA_; ii++) {
	printf(" %16.10E ", spMf_final_[ii]);
      } 
      printf("\n");
      printf("          densA_init_  = %16.10E",  densA_init_);
      for (size_t ii = 0; ii < nSpeciesA_; ii++) {
	printf(" %16.10E ", spMf_init_[ii]);
      } 
      printf("\n");
      printf("resid_1_massF = %13.7E\n",   massFluxA_);
#endif
      if (massFluxA_ < 0.0) {
	// We gain material from the surrounding boundary region
	//    resid_1 -= (advecS_A - advecCVB_A) * densBoundA_;
	resid_1 -= advecCVB_A * densBoundA_;
	double advecS_A = resid_1 / densBoundA_;

	Velo_S_Interface_final_ = - advecS_A /  surfaceAreaRS_final_[0];

	resid_1 = (densA_final_ - densA_init_) /  deltaTsubcycleCalc_ * grossAVol_ - massFluxA_;
	advecS_A    = - (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
	resid_1 -= (advecS_A - advecCVB_A) * densBoundA_;
#ifdef DEBUG_MODE_NOT
	printf("resid_1 = %g\n", resid_1);
#endif
      } else {
	// We lose material to the surrounding boundary region (advecA is negative)
	//  resid_1 -= (advecS_A - advecCVB_A) * densA_final_;

	resid_1 -= advecCVB_A * densA_final_;
	double advecS_A = resid_1 / densA_final_;

	Velo_S_Interface_final_ = - advecS_A /  surfaceAreaRS_final_[0];

      
	resid_1 = (densA_final_ - densA_init_) /  deltaTsubcycleCalc_ * grossAVol_ - massFluxA_;
	advecS_A    = - (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
	resid_1 -= (advecS_A - advecCVB_A) * densA_final_;
#ifdef DEBUG_MODE_NOT
	printf("resid_1_advec = %13.7E\n", (advecS_A - advecCVB_A) * densA_final_);
	printf("resid_1 = %13.7E\n", resid_1);
#endif
      }
   
      /*
       *  source terms for the B phase
       *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
       *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
       */
      bstart = m_PhaseSpeciesStartIndex[solnBPhase_];
  
      double advecS_B = Velo_S_Interface_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
      double advecB_B = Velo_B_CVB_final_ * deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    
      double moleBSum = 0.0;
      for (size_t k = 1; k < nSpeciesB_; k++)  {
     
	/*
	 *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
	 *  accumulation term.
	 */
	double  delta = (deltaTsubcycleCalc_ * DspMoles_final_[k + bstart]
			 +      jKmolfluxB_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_);

	if (massFluxB_ < 0.0) {
	  // We gain material from the surrounding boundary region
	  delta += (advecS_B - advecB_B) * (cBoundB_[k]);
	} else {
	  // We lose material to the surrounding boundary region (advecB is negative)
	  delta += (advecS_B - advecB_B) * ( cSurfB_[k]);
	}
  

	spMoles_final_[bstart + k] = spMoles_init_[bstart + k] + delta;
	moleBSum += spMoles_final_[bstart + k];
      }
      spMoles_final_[bstart + 0] = grossBVol_ / molarVolB_ - moleBSum;

      updatePhaseNumbers(solnBPhase_);
      updateVelocities();

      advecS_B    =  (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
      //double advecCVB_  = - Velo_A_CVB_final_         * surfaceAreaRS_final_[0];
      resid_1 = (densB_final_ - densB_init_) /  deltaTsubcycleCalc_ * grossBVol_ - massFluxB_;
      if (massFluxB_ < 0.0) {
	// We gain material from the surrounding boundary region
	//    resid_1 -= (advecS_A - advecCVB_A) * densBoundA_;
	resid_1 -= advecS_B * densBoundB_;
	double advecCVB_B = - resid_1 / densBoundB_;

	Velo_B_CVB_final_ = advecCVB_B /  surfaceAreaRS_final_[0];

	resid_1 = (densB_final_ - densB_init_) /  deltaTsubcycleCalc_ * grossBVol_ - massFluxB_;
	advecS_B  = (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
	resid_1 -= (advecS_B - advecB_B) * densBoundB_;
#ifdef DEBUG_MODE_NOT
	printf("resid_B = %g\n", resid_1);
#endif
      } else {
	// We lose material to the surrounding boundary region (advecA is negative)
	//  resid_1 -= (advecS_B - advecCVB_B) * densA_final_;

	resid_1 -= advecS_B * densB_final_;
	double advecCVB_B = - resid_1 / densB_final_;

	Velo_B_CVB_final_ =  advecCVB_B /  surfaceAreaRS_final_[0];

	resid_1 = (densB_final_ - densB_init_) /  deltaTsubcycleCalc_ * grossBVol_ - massFluxB_;
	advecS_B    = Velo_S_Interface_final_ * surfaceAreaRS_final_[0];
	advecCVB_B  = Velo_B_CVB_final_       * surfaceAreaRS_final_[0];
    
    
	resid_1 -= (advecS_B - advecCVB_B) * densB_final_;
#ifdef DEBUG_MODE_NOT
	printf("resid_B = %g\n", resid_1);
#endif
      }
    }

    /*
     *  Create a failure if we significantly decrease the moles to near zero. This will never pass
     *  the predictor corrector
     */

    soln_predict_[0] = deltaTsubcycleCalc_;
    soln_predict_[1] = Velo_S_Interface_final_;
    double retn = 1;
    for (size_t k = 0; k < nSpeciesA_; k++)  {
      soln_predict_[2 + astart + k] = spMoles_final_[astart + k];
      if (spMoles_init_[astart + k] > 0.0) {
	if (spMoles_final_[astart + k] < 0.33 * spMoles_init_[astart + k]) {
	  retn = -1;
	}
      }
    }
    soln_predict_[2 + bstart] = Velo_B_CVB_final_;
    for (size_t k = 0; k < nSpeciesB_; k++)  {
      soln_predict_[3 + bstart + k] = spMoles_final_[bstart + k];
      if (spMoles_init_[bstart + k] > 0.0) {
	if (spMoles_final_[bstart + k] < 0.33 * spMoles_init_[bstart + k]) {
	  retn = -1;
	}
      }
    }

    return retn;
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
  void  imtPSS_NoSurf_DiffBL::resetStartingCondition(double Tinitial) {

    if (pendingIntegratedStep_ != 1) {
#ifdef DEBUG_ELECTRODE
      // printf(" imtPSS_NoSurf_DiffBL::resetStartingCondition WARNING: resetStartingCondition called with no pending integration step\n");
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
      throw CanteraError("imtPSS_NoSurf_DiffBL::resetStartingCondition()", "tinit " + fp2str(Tinitial) +" not compat with t_final_final_ "
			 + fp2str(t_final_final_));
    }



    t_init_init_ = Tinitial;

    setInitStateFromFinal(true);

    // reset surface quantities
    for (size_t i = 0; i < m_NumSurPhases; i++) {
      surfaceAreaRS_init_[i] = surfaceAreaRS_final_[i];
    }

    // Reset total species quantities
    for (size_t k = 0; k < m_NumTotSpecies; k++) {
      spMoles_init_[k] = spMoles_final_[k];
      spMoles_init_init_[k] = spMoles_final_[k];
    }

    mdpUtil::mdp_zero_dbl_1(DATA_PTR(spMoleIntegratedSourceTerm_), m_NumTotSpecies);
    mdpUtil::mdp_zero_dbl_1(DATA_PTR(spMoleIntegratedSourceTermLast_), m_NumTotSpecies);

    // Reset the total phase moles quantities
    for (size_t i = 0; i < m_NumTotPhases; i++) {
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

  //====================================================================================================================
  // Pack the nonlinear solver problem
  /*
   *  formulate the nonlinear solver problem to be solved.
   *     Fields to be filled in
   *             yvalNLS_
   *             ylowNLS_
   *             yhighNLS_
   *             atolNLS_
   *             deltaBoundsMagnitudesNLS_     
   */
  void imtPSS_NoSurf_DiffBL::initialPackSolver_nonlinFunction()
  { 
#ifdef DEBUG_METHOD_F1
    yvalNLS_[0] = deltaTsubcycleCalc_;
    ylowNLS_[0] = 0.0;
    yhighNLS_[0] = 1.0E300;
    atolNLS_[0] =  rtol_IntegratedSrc_global_ * 1.0E-3 * (t_final_final_ - t_init_init_);
    deltaBoundsMagnitudesNLS_[0] =  (t_final_final_ - t_init_init_);
    residAtolNLS_[0] = atolNLS_[0];


    yvalNLS_[1] = Velo_S_Interface_final_;
    ylowNLS_[1] = -1.0E50;
    yhighNLS_[1] = 1.0E50;

    size_t astart = 0;
    double sum = 0.0;
    for (size_t k = 0; k < nSpeciesA_; k++) {
      yvalNLS_[ 2 + k] =   spMoles_final_[k + astart];
      sum +=  yvalNLS_[ 2 + k];
      ylowNLS_[ 2 + k] = 0; 
     
      yhighNLS_[2 + k] = 1.0E300;
    }

    yvalNLS_[2 + nSpeciesA_] = Velo_B_CVB_final_;
    ylowNLS_[2 + nSpeciesA_] = -1.0E50;
    yhighNLS_[2 + nSpeciesA_] = 1.0E50;

    size_t bstart = m_PhaseSpeciesStartIndex[solnBPhase_];
    for (size_t k = 0; k < nSpeciesB_; k++) {
      yvalNLS_[ 3 + nSpeciesA_ + k] = spMoles_final_[k + bstart]; 
      sum +=  yvalNLS_[ 1 + nSpeciesA_ + k];
      ylowNLS_[ 3 + nSpeciesA_ + k] = 0; 
      yhighNLS_[3 + nSpeciesA_ + k] = 1.0E300;
    }


    for (size_t k = 0; k < nSpeciesA_; k++) {    
      atolNLS_[ 2 + k] = 1.0E-13 * sum;
      residAtolNLS_[ 2 + k ] =   1.0E-13 * sum;
      deltaBoundsMagnitudesNLS_[ 2 + k] = MAX( 1.0E-9 * sum, 0.5 * spMoles_final_[k + astart]);
    }

    for (size_t k = 0; k < nSpeciesB_; k++) {
      atolNLS_[ 3 + nSpeciesA_ + k] = 1.0E-13 * sum;
      residAtolNLS_[ 3 + nSpeciesA_ + k ] =   1.0E-13 * sum;
      deltaBoundsMagnitudesNLS_[ 3 + nSpeciesA_ + k] = MAX( 1.0E-9 * sum, 0.5 * spMoles_final_[k + bstart]);
    }

    atolNLS_[1] = 1.0E-5;
    residAtolNLS_[1] = 1.0E-5;
    deltaBoundsMagnitudesNLS_[1] = 1.0;

    atolNLS_[2 + nSpeciesA_] = 1.0E-5;
    residAtolNLS_[2 + nSpeciesA_] = 1.0E-5;
    deltaBoundsMagnitudesNLS_[2 + nSpeciesA_] = 1.0;
#else
    // Time
    yvalNLS_[0] = deltaTsubcycleCalc_;
    ylowNLS_[0] = 0.0;
    yhighNLS_[0] = 1.0E300;
    atolNLS_[0] =  rtol_IntegratedSrc_global_ * 1.0E-3 * (t_final_final_ - t_init_init_);
    deltaBoundsMagnitudesNLS_[0] =  (t_final_final_ - t_init_init_);
    residAtolNLS_[0] = atolNLS_[0];

    // Velo_S_Interface_final > A velocity
    double veloMax = MAX(1.0E-5, fabs(Velo_S_Interface_final_));
    veloMax = MAX(veloMax, fabs(Velo_B_CVB_final_));

    yvalNLS_[1] = Velo_S_Interface_final_;
    ylowNLS_[1] = -1.0E50;
    yhighNLS_[1] = 1.0E50;


    // A species 
    size_t astart = 0;
    double sum = 0.0;
    for (size_t k = 0; k < nSpeciesA_; k++) {
      yvalNLS_[ 2 + k] =   spMoles_final_[k + astart];
      sum +=  yvalNLS_[ 2 + k];
      ylowNLS_[ 2 + k] = 0; 
     
      yhighNLS_[2 + k] = 1.0E300;
    }

    // Velo_B_Interface_final > B velocity
    yvalNLS_[2 + nSpeciesA_] = Velo_B_CVB_final_;
    ylowNLS_[2 + nSpeciesA_] = -1.0E50;
    yhighNLS_[2 + nSpeciesA_] = 1.0E50;

    // B species 
    size_t bstart = m_PhaseSpeciesStartIndex[solnBPhase_];
    for (size_t k = 0; k < nSpeciesB_; k++) {
      yvalNLS_[ 3 + nSpeciesA_ + k] = spMoles_final_[k + bstart]; 
      sum +=  yvalNLS_[ 3 + nSpeciesA_ + k];
      ylowNLS_[ 3 + nSpeciesA_ + k] = 0; 
      yhighNLS_[3 + nSpeciesA_ + k] = 1.0E300;
    }

    if (sum <= 0.0) {
      throw CanteraError(" imtPSS_NoSurf_DiffBL::initialPackSolver_nonlinFunction()",
			 "sum <= 0.0");
    }

    for (size_t k = 0; k < nSpeciesA_; k++) {    
      atolNLS_[ 2 + k] = 1.0E-13 * sum;
      residAtolNLS_[ 2 + k ] =   1.0E-13 * sum;
      deltaBoundsMagnitudesNLS_[ 2 + k] = MAX( 1.0E-9 * sum, 0.5 * spMoles_final_[k + astart]);
    }

    for (size_t k = 0; k < nSpeciesB_; k++) {
      atolNLS_[ 3 + nSpeciesA_ + k] = 1.0E-13 * sum;
      residAtolNLS_[ 3 + nSpeciesA_ + k ] =   1.0E-13 * sum;
      deltaBoundsMagnitudesNLS_[ 3 + nSpeciesA_ + k] = MAX( 1.0E-9 * sum, 0.5 * spMoles_final_[k + bstart]);
    }

    atolNLS_[1] = veloMax;
    residAtolNLS_[1] = veloMax;
    deltaBoundsMagnitudesNLS_[1] = 1.0;

    atolNLS_[2 + nSpeciesA_] = veloMax;
    residAtolNLS_[2 + nSpeciesA_] = veloMax;
    deltaBoundsMagnitudesNLS_[2 + nSpeciesA_] = 1.0;

#endif

  }

  //====================================================================================================================
  // Update the velocities within the object at the final conditions.
  /*
   *  (virtual from InterfacialMassTransfer)
   * 
   *  In order to do this, we need the instantaenous source terms and the current value of the
   *  surface area over the local time step. Therefore, this routine must be called 
   *  after DspMoles_final_ has been calculated. 
   */
  void imtPSS_NoSurf_DiffBL::updateVelocities()
  {
    massFluxA_ = getPhaseAMassSourceTerm();
    ThermoPhase &tpA = thermo(solnAPhase_);
    densA_final_ = tpA.density();
    double area = 0.5 * (surfaceAreaRS_init_[0] + surfaceAreaRS_final_[0]);
    double sVelocA = massFluxA_  / area / densA_final_;

    massFluxB_ = getPhaseBMassSourceTerm();
    ThermoPhase &tpB = thermo(solnBPhase_);
    densB_final_ = tpB.density();
    double sVelocB = massFluxB_  / area / densB_final_;

    Velo_A_Interface_final_ =  Velo_S_Interface_final_ - sVelocA;
    Velo_B_Interface_final_ = sVelocB + Velo_S_Interface_final_;

    Velo_A_CVB_final_ = 0.0;
  }

  //==================================================================================================================
  // calculate the residual
  /*
   *   
   *  In this residual we treat species 0 as the solvent. This gets the total concentration equals
   *  the thermo concentration constraint equation. This constraint equation is needed, because we choose not
   *  to solve for the pressure in this formulation.
   *
   *           index                   name
   *    --------------------------------------------------------------------
   *              0                     deltaT
   *              1                     Velo_S_Interface_final_
   *              2                     N_A_final_[0]
   *              2 + nSpeciesA-1       N_A_final_[nSpeciesA_1]
   *    brstart + 1                     Velo_B_CVB_final_
   *    brstart + 2                     N_B_final_[0]
   *    brstart + 2 + nSpeciesA-1       N_B_final_[nSpeciesB_1]
   *
   *   We define that Velo_A_CVB_final_ = 0.0 as a start of the calculation.
   */
  int imtPSS_NoSurf_DiffBL::calcResid_F1(doublevalue* const resid, const ResidEval_Type_Enum evalType)
  {
    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;

    double advecS_A =  - (Velo_S_Interface_final_) *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecCVB_A  = - Velo_A_CVB_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    /*
     *  source terms for the A phase
     *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
     *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
     */
    size_t astart = 0;
    size_t arstart = 2;
    double moleASum = 0.0;
    for (size_t k = 0; k < nSpeciesA_; k++)  {
      moleASum += spMoles_final_[astart + k];
      /*
       *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      resid[k+arstart] = spMoles_final_[astart + k] - spMoles_init_[astart + k]
	- deltaTsubcycleCalc_ * DspMoles_final_[k + astart];

      /*
       *  Add in the diffusive flux to the surface of species. This sums to zero on a mass-averaged velocity basis.
       */
      resid[k+arstart] -= jKmolfluxA_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_;

      /*
       *  Add in the advection term due to the convection of the control volume due to the Stefan velocity of the
       *  interface
       */
      if ((Velo_S_Interface_final_ -Velo_A_CVB_final_) < 0.0) {
	// We gain material from the surrounding boundary region
	resid[k+arstart] -= (advecS_A - advecCVB_A) * (cBoundA_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecA is negative)
	resid[k+arstart] -= (advecS_A - advecCVB_A) * ( cSurfA_[k]);
      }
    }
    /*
     *  We have to close the system using the following total moles
     *  equation and then apply it to find the difference in velocity of the A phase and S interface at the far boundary
     */
    resid[1] = moleASum - grossAVol_ / molarVolA_;
    
    /*
     *  source terms for the B phase
     *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
     *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
     */
    size_t bstart =  m_PhaseSpeciesStartIndex[solnBPhase_];
    size_t brstart =  3 + nSpeciesA_;
    double moleBSum = 0.0;
    double advecS_B = Velo_S_Interface_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecCVB_B  =  Velo_B_CVB_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    for (size_t k = 0; k < nSpeciesB_; k++)  {
      moleBSum += spMoles_final_[bstart + k];
      /*
       *  The residual involves the mol balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      resid[k +brstart] = spMoles_final_[bstart + k] - spMoles_init_[bstart + k]
	- deltaTsubcycleCalc_ * DspMoles_final_[k + bstart];
      
      /*
       *  Add in the diffusive flux to the surface of species. This sums to zero on a mass-averaged velocity basis.
       */
      resid[k + brstart] -= - jKmolfluxB_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_;

      /*
       *  Add in the advection term due to the convection of the control volume due to the Stefan velocity of the
       *  interface
       */
      if ((Velo_S_Interface_final_ -  Velo_B_CVB_final_) > 0) {
	// We gain material from the surrounding boundary region
	resid[k + brstart] -= (advecS_B - advecCVB_B) * ( cBoundB_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecB is negative)
	resid[k + brstart] -= (advecS_B - advecCVB_B) * (  cSurfB_[k]);
      }

    }
    /*
     *  We assume that the solvent is species 0. Then we have to close the system using the following total moles
     *  equation
     */
    resid[brstart - 1] =  moleBSum -  grossBVol_/ molarVolB_;

    return 1;
  }
  //==================================================================================================================
  // calculate the residual
  /*
   *   
   *  In this residual we treat species 0 as the solvent. This gets the total concentration equals
   *  the thermo concentration constraint equation. This constraint equation is needed, because we choose not
   *  to solve for the pressure in this formulation.
   *
   *           index                   name
   *    --------------------------------------------------------------------
   *              0                     deltaT
   *              1                     Velo_S_Interface_final_
   *              2                     N_A_final_[0]
   *              2 + nSpeciesA-1       N_A_final_[nSpeciesA_1]
   *    brstart + 1                     Velo_B_CVB_final_
   *    brstart + 2                     N_B_final_[0]
   *    brstart + 2 + nSpeciesA-1       N_B_final_[nSpeciesB_1]
   *
   *   We define that Velo_A_CVB_final_ = 0.0 as a start of the calculation.
   */
  int imtPSS_NoSurf_DiffBL::calcResid_F2(doublevalue* const resid, const ResidEval_Type_Enum evalType)
  {
    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;

    double advecS_A =  - (Velo_S_Interface_final_) *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecCVB_A  = - Velo_A_CVB_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    /*
     *  source terms for the A phase
     *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
     *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
     */
    size_t astart = 0;
    size_t arstart = 2;
    double moleASum = 0.0;
    for (size_t k = 0; k < nSpeciesA_; k++)  {
      moleASum += spMoles_final_[astart + k];
      /*
       *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      resid[k+arstart] = spMoles_final_[astart + k] - spMoles_init_[astart + k]
	- deltaTsubcycleCalc_ * DspMoles_final_[k + astart];

      /*
       *  Add in the diffusive flux to the surface of species. This sums to zero on a mass-averaged velocity basis.
       */
      resid[k+arstart] -= jKmolfluxA_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_;

      /*
       *  Add in the advection term due to the convection of the control volume due to the Stefan velocity of the
       *  interface
       */
      if ((Velo_S_Interface_final_ -Velo_A_CVB_final_) < 0.0) {
	// We gain material from the surrounding boundary region
	resid[k+arstart] -= (advecS_A - advecCVB_A) * (cBoundA_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecA is negative)
	resid[k+arstart] -= (advecS_A - advecCVB_A) * ( cSurfA_[k]);
      }
    }
    /*
     *  We have to close the system using the following total moles
     *  equation and then apply it to find the difference in velocity of the A phase and S interface at the far boundary
     */
    resid[1] = moleASum - grossAVol_ / molarVolA_;
    
    /*
     *  source terms for the B phase
     *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
     *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
     */
    size_t bstart =  m_PhaseSpeciesStartIndex[solnBPhase_];
    size_t brstart =  3 + nSpeciesA_;
    double moleBSum = 0.0;
    double advecS_B = Velo_S_Interface_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecCVB_B  =  Velo_B_CVB_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    for (size_t k = 0; k < nSpeciesB_; k++)  {
      moleBSum += spMoles_final_[bstart + k];
      /*
       *  The residual involves the mol balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      resid[k +brstart] = spMoles_final_[bstart + k] - spMoles_init_[bstart + k]
	- deltaTsubcycleCalc_ * DspMoles_final_[k + bstart];
      
      /*
       *  Add in the diffusive flux to the surface of species. This sums to zero on a mass-averaged velocity basis.
       */
      resid[k + brstart] -= - jKmolfluxB_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_;

      /*
       *  Add in the advection term due to the convection of the control volume due to the Stefan velocity of the
       *  interface
       */
      if ((Velo_S_Interface_final_ -  Velo_B_CVB_final_) > 0) {
	// We gain material from the surrounding boundary region
	resid[k + brstart] -= (advecS_B - advecCVB_B) * ( cBoundB_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecB is negative)
	resid[k + brstart] -= (advecS_B - advecCVB_B) * (  cSurfB_[k]);
      }

    }
    /*
     *  We assume that the solvent is species 0. Then we have to close the system using the following total moles
     *  equation
     */
    resid[brstart - 1] =  moleBSum -  grossBVol_/ molarVolB_;

    return 1;
  }
  //==================================================================================================================
  // calculate the residual
  /*
   *   
   *  In this residual we treat species 0 as the solvent. This gets the total concentration equals
   *  the thermo concentration constraint equation. This constraint equation is needed, because we choose not
   *  to solve for the pressure in this formulation.
   *
   *  neq_ = 3 +  nSpeciesA +  nSpeciesB
   *
   *           index                  Unknown_name                               Eqn
   *    --------------------------------------------------------------------------------------------------
   *              0                     deltaT                               time eqn
   *              1                     Velo_S_Interface_final_              Total mass continuity equation
   *              2                     N_A_final_[0]                        volume conservation
   *              2 + nSpeciesA-1       N_A_final_[nSpeciesA_1]              individual mass conservation equation
   *    brstart + 1                     Velo_B_CVB_final_
   *    brstart + 2                     N_B_final_[0]
   *    brstart + 2 + nSpeciesA-1       N_B_final_[nSpeciesB_1]
   *
   *   We define that Velo_A_CVB_final_ = 0.0 as a start of the calculation.
   */
  int imtPSS_NoSurf_DiffBL::calcResid(doublevalue* const resid, const ResidEval_Type_Enum evalType)
  {

#ifdef DEBUG_METHOD_F1
    return calcResid_F1(resid, evalType);
#else 
#ifdef DEBUG_METHOD_F2
    return calcResid_F2(resid, evalType);
#endif
#endif

    resid[0] = deltaTsubcycleCalc_ - deltaTsubcycle_;

    double advecS_A    = - (Velo_S_Interface_final_) *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecCVB_A  = - Velo_A_CVB_final_ * deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    /*
     *  source terms for the A phase
     *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
     *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
     */
    size_t astart = 0;
    size_t arstart = 2;
    double mfSum = spMf_final_[astart];
    double moleASum = spMoles_final_[astart];
    //double mfSum = spMf_final_[astart];
    for (size_t k = 1; k < nSpeciesA_; k++) {
      moleASum += spMoles_final_[astart + k];
      mfSum += spMf_final_[astart + k];
      /*
       *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      resid[k+arstart] = spMoles_final_[astart + k] - spMoles_init_[astart + k]
	- deltaTsubcycleCalc_ * DspMoles_final_[k + astart];

      /*
       *  Add in the diffusive flux to the surface of species. This sums to zero on a mass-averaged velocity basis.
       */
      resid[k+arstart] -= jKmolfluxA_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_;

      /*
       *  Add in the advection term due to the convection of the control volume due to the Stefan velocity of the
       *  interface
       */
      if (massFluxA_ < 0.0) {
	// We gain material from the surrounding boundary region
	resid[k+arstart] -= (advecS_A - advecCVB_A) * (cBoundA_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecA is negative)
	resid[k+arstart] -= (advecS_A - advecCVB_A) * ( cSurfA_[k]);
      }
    }
    /*
     *  Throw out one convection equation for the molar volume equation
     */

    resid[arstart] = moleASum - grossAVol_ / molarVolA_;


    /*
     *  We have to close the system using the following total moles
     *  equation and then apply it to find the difference in velocity of the A phase and S interface at the far boundary
     */ 
    advecS_A    = - (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
    advecCVB_A  = - Velo_A_CVB_final_         * surfaceAreaRS_final_[0];
    resid[1] = (densA_final_ - densA_init_) /  deltaTsubcycleCalc_ * grossAVol_ - massFluxA_;
    if (massFluxA_ < 0.0) {
      // We gain material from the surrounding boundary region
      resid[1] -= (advecS_A - advecCVB_A) * densBoundA_;
    } else {
      // We lose material to the surrounding boundary region (advecA is negative)
      resid[1] -= (advecS_A - advecCVB_A) * densA_final_;
    }
#ifdef DEBUG_MODE_NOT
    printf("resid_1_drho_dt = %13.7E\n",  (densA_final_ - densA_init_) /  deltaTsubcycleCalc_ * grossAVol_ );
    printf("          densA_final_ = %16.10E",  densA_final_);
    for (size_t ii = 0; ii < nSpeciesA_; ii++) {
      printf(" %16.10E ", spMf_final_[ii]);
    }
    printf("\n");
    printf("          densA_init_  = %16.10E",  densA_init_);
    for (size_t ii = 0; ii < nSpeciesA_; ii++) {
      printf(" %16.10E ", spMf_init_[ii]);
    } 
    printf("\n");
    printf("resid_1_massF = %13.7E\n",   massFluxA_);
    if (massFluxA_ < 0.0) {
      printf("resid_1_advec = %13.7E\n", (advecS_A - advecCVB_A) * densBoundA_);
    } else {
      printf("resid_1_advec = %13.7E\n", (advecS_A - advecCVB_A) * densA_final_);
    }
    printf("resid_1 = %13.7E\n", resid[1]);
#endif
    /*
     *  source terms for the B phase
     *   - we stick the origin on the interface. Within that reference frame the Stefan velocity causes a
     *     net advective flux which helps the mass transport in one direction and hinders it in the other direction.
     */
    size_t bstart =  m_PhaseSpeciesStartIndex[solnBPhase_];
    size_t brstart =  3 + nSpeciesA_;
    double moleBSum = spMoles_final_[bstart];
    double advecS_B = Velo_S_Interface_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecCVB_B  =  Velo_B_CVB_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];

    for (size_t k = 1; k < nSpeciesB_; k++)  {
      moleBSum += spMoles_final_[bstart + k];
      /*
       *  The residual involves the mol balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      resid[k +brstart] = spMoles_final_[bstart + k] - spMoles_init_[bstart + k]
	- deltaTsubcycleCalc_ * DspMoles_final_[k + bstart];
      
      /*
       *  Add in the diffusive flux to the surface of species. This sums to zero on a mass-averaged velocity basis.
       */
      resid[k + brstart] -= - jKmolfluxB_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_;

      /*
       *  Add in the advection term due to the convection of the control volume due to the Stefan velocity of the
       *  interface
       */
      if (massFluxB_ < 0.0) {
	// We gain material from the surrounding boundary region
	resid[k + brstart] -= (advecS_B - advecCVB_B) * ( cBoundB_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecB is negative)
	resid[k + brstart] -= (advecS_B - advecCVB_B) * (  cSurfB_[k]);
      }

    }
    /*
     *  We assume that the solvent is species 0. Then we have to close the system using the following total moles
     *  equation
     */
    resid[brstart] =  moleBSum - grossBVol_ / molarVolB_; 


    /*
     *  We have to close the system using the following total moles
     *  equation and then apply it to find the difference in velocity of the A phase and S interface at the far boundary
     */ 
    advecS_B    = (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
    advecCVB_B  =  Velo_B_CVB_final_        * surfaceAreaRS_final_[0];

    resid[brstart - 1] = (densB_final_ - densB_init_) /  deltaTsubcycleCalc_ * grossBVol_ - massFluxB_;
    if (massFluxB_ < 0.0) {
      // We gain material from the surrounding boundary region
      resid[brstart - 1] -= (advecS_B - advecCVB_B) * densBoundB_;
    } else {
      // We lose material to the surrounding boundary region (advecA is negative)
      resid[brstart - 1] -= (advecS_B - advecCVB_B) * densB_final_;
    }

   

    return 1;
  }

  //====================================================================================================================

  void imtPSS_NoSurf_DiffBL::printResid_ResidSatisfaction()
  {
    double deltaMoles, rxnMoles, fluxMoles, advMoles, denom, resid1, normResid;
 
    size_t astart = 0;
    //int arstart = 2;
    double mfSum = 0.0;
    double moleASum = 0.0;
    double advecS_A    = - (Velo_S_Interface_final_) *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecCVB_A  = - Velo_A_CVB_final_ * deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];

    printf("\t\t Residual Satisfaction\n");

   
    printf("\t\t      k  deltaMoles      rxnMole         fluxMoles       advMoles        resid           normResid \n");
    for (size_t k = 0; k < nSpeciesA_; k++) {
      moleASum += spMoles_final_[astart + k];
      mfSum += spMf_final_[astart + k];
      /*
       *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      deltaMoles = spMoles_final_[astart + k] - spMoles_init_[astart + k];
      rxnMoles = deltaTsubcycleCalc_ * DspMoles_final_[k + astart];
      fluxMoles = jKmolfluxA_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_;
  

      /*
       *  Add in the advection term due to the convection of the control volume due to the Stefan velocity of the
       *  interface
       */
      if (massFluxA_ < 0.0) {
	// We gain material from the surrounding boundary region

	advMoles = (advecS_A - advecCVB_A) * (cBoundA_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecA is negative)

	advMoles =(advecS_A - advecCVB_A) * ( cSurfA_[k]);
      }
      denom = MAX( spMoles_final_[astart + k],  spMoles_init_[astart + k]);
      denom = MAX( denom,  fabs(rxnMoles));
      denom = MAX(denom, fabs(fluxMoles));
      denom = MAX(denom, fabs(advMoles));
      resid1 = deltaMoles - rxnMoles - fluxMoles - advMoles;
      if (denom > 0.0) {
	normResid = resid1/denom;
      } else {
	normResid = resid1;
      }
      
      printf("\t\t    %3d %14.7E  %14.7E  %14.7E  %14.7E  %14.7E  %14.7E  \n", static_cast<int>(k), deltaMoles,
             rxnMoles, fluxMoles, advMoles, resid1, normResid);
    }
  
    printf("\t\t      k  Mole Sum        GrossVol\n");
    double volSum  = grossAVol_ / molarVolA_;
    resid1 = moleASum - volSum;
    denom = MAX(fabs(moleASum), fabs(volSum));
    if (denom > 0.0) {
      normResid = resid1/denom;
    } else {
      normResid = resid1;
    }
      
    printf("\t\t    %3d %14.7E  %14.7E                                  %14.7E  %14.7E\n", 0, moleASum, volSum, resid1, normResid);

    /*
     *  We have to close the system using the following total moles
     *  equation and then apply it to find the difference in velocity of the A phase and S interface at the far boundary
     */ 
    advecS_A    = - (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
    advecCVB_A  = - Velo_A_CVB_final_         * surfaceAreaRS_final_[0];
    double drhodt = (densA_final_ - densA_init_) /  deltaTsubcycleCalc_ * grossAVol_ ;
    if (massFluxA_ < 0.0) {
      // We gain material from the surrounding boundary region
      advMoles = (advecS_A - advecCVB_A) * densBoundA_;
    } else {
      // We lose material to the surrounding boundary region (advecA is negative)
      advMoles = (advecS_A - advecCVB_A) * densA_final_;
    }
    resid1 = drhodt - massFluxA_ - advMoles;
 
    denom = MAX( fabs(drhodt), fabs(massFluxA_));
    denom = MAX(denom, fabs(advMoles));
    if (denom > 0.0) {
      normResid = resid1/denom;
    } else {
      normResid = resid1;
    }
    printf("\t\t         drhodt          massFluxA_                      advDens  \n");                 
    printf("\t\t    V_S %14.7E  %14.7E                  %14.7E  %14.7E  %14.7E  \n", drhodt, massFluxA_, advMoles, resid1, normResid);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    int bstart =  m_PhaseSpeciesStartIndex[solnBPhase_];
    double moleBSum = 0;
    double advecS_B = Velo_S_Interface_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    double advecCVB_B  =  Velo_B_CVB_final_ *  deltaTsubcycleCalc_ * surfaceAreaRS_final_[0];
    printf("\t\t      k  deltaMoles      rxnMole         fluxMoles       advMoles        resid           normResid \n");
    for (size_t k = 0; k < nSpeciesB_; k++) {
      moleBSum += spMoles_final_[bstart + k];
      mfSum += spMf_final_[bstart + k];
      /*
       *  The residual involves the mole balance of material in the phase A. Here is the reaction term and the
       *  accumulation term.
       */
      deltaMoles = spMoles_final_[bstart + k] - spMoles_init_[bstart + k];
      rxnMoles = deltaTsubcycleCalc_ * DspMoles_final_[k + bstart];
      fluxMoles = - jKmolfluxB_[k] * surfaceAreaRS_final_[0] * deltaTsubcycleCalc_;
  

      /*
       *  Add in the advection term due to the convection of the control volume due to the Stefan velocity of the
       *  interface
       */
      if (massFluxB_ < 0.0) {
	// We gain material from the surrounding boundary region

	advMoles = (advecS_B - advecCVB_B) * (cBoundB_[k]);
      } else {
	// We lose material to the surrounding boundary region (advecA is negative)

	advMoles =(advecS_B - advecCVB_B) * ( cSurfB_[k]);
      }
      denom = MAX(spMoles_final_[bstart + k],  spMoles_init_[bstart + k]);
      denom = MAX(denom,  fabs(rxnMoles));
      denom = MAX(denom, fabs(fluxMoles));
      denom = MAX(denom, fabs(advMoles));
      resid1 = deltaMoles - rxnMoles - fluxMoles - advMoles;
      if (denom > 0.0) {
	normResid = resid1/denom;
      } else {
	normResid = resid1;
      }
      
      printf("\t\t    %3d %14.7E  %14.7E  %14.7E  %14.7E  %14.7E  %14.7E  \n", static_cast<int>(k), deltaMoles, rxnMoles, 
            fluxMoles, advMoles, resid1, normResid);
    }
  
    printf("\t\t      k  Mole Sum        GrossVol\n");
    volSum  = grossBVol_ / molarVolB_;
    resid1 = moleBSum - volSum;
    denom = MAX(fabs(moleBSum), fabs(volSum));
    if (denom > 0.0) {
      normResid = resid1/denom;
    } else {
      normResid = resid1;
    }
      
    printf("\t\t    %3d %14.7E  %14.7E                                  %14.7E  %14.7E\n", 0, moleBSum, volSum, resid1, normResid);

    /*
     *  We have to close the system using the following total moles
     *  equation and then apply it to find the difference in velocity of the A phase and S interface at the far boundary
     */ 
    advecS_B    = (Velo_S_Interface_final_) * surfaceAreaRS_final_[0];
    advecCVB_B  = Velo_B_CVB_final_         * surfaceAreaRS_final_[0];
    drhodt = (densB_final_ - densB_init_) /  deltaTsubcycleCalc_ * grossBVol_ ;
    if (massFluxB_ < 0.0) {
      // We gain material from the surrounding boundary region
      advMoles = (advecS_B - advecCVB_B) * densBoundB_;
    } else {
      // We lose material to the surrounding boundary region (advecA is negative)
      advMoles = (advecS_B - advecCVB_B) * densB_final_;
    }
    resid1 = drhodt - massFluxB_ - advMoles;
 
    denom = MAX( fabs(drhodt), fabs(massFluxB_));
    denom = MAX(denom, fabs(advMoles));
    if (denom > 0.0) {
      normResid = resid1/denom;
    } else {
      normResid = resid1;
    }
    printf("\t\t         drhodt          massFluxB_                      advDens  \n");                 
    printf("\t\t    V_S %14.7E  %14.7E                  %14.7E  %14.7E  %14.7E  \n", drhodt, massFluxB_, advMoles, resid1, normResid);
  }


  //====================================================================================================================
  //     Unpack the soln vector
  /*
   *  (virtual from InterfacialMassTransfer_Integrator)
   *
   *  This function unpacks the solution vector into  phaseMoles_final_,  spMoles_final_, and spMf_final_[]
   */
  void imtPSS_NoSurf_DiffBL::unpackNonlinSolnVector(const double * const y)
  {
#ifdef DEBUG_METHOD_F1
    deltaTsubcycleCalc_ = y[0];
    t_final_ = t_init_ + deltaTsubcycleCalc_;
    Velo_S_Interface_final_ = y[1];
    size_t astart = 0;
    for (size_t k = 0; k < nSpeciesA_; k++) {
      spMoles_final_[k + astart] = y[2 + k];
    }
    Velo_B_CVB_final_ = y[2 + nSpeciesA_];

    size_t bstart = m_PhaseSpeciesStartIndex[solnBPhase_];
    for (size_t k = 0; k < nSpeciesB_; k++) {
      spMoles_final_[k + bstart] = y[3 + nSpeciesA_ +  k];
    }
#else 
    deltaTsubcycleCalc_ = y[0];
    t_final_ = t_init_ + deltaTsubcycleCalc_;
    Velo_S_Interface_final_ = y[1];
    size_t astart = 0;
    for (size_t k = 0; k < nSpeciesA_; k++) {
      spMoles_final_[k + astart] = y[2 + k];
    }
    Velo_B_CVB_final_ = y[2 + nSpeciesA_];

    size_t bstart = m_PhaseSpeciesStartIndex[solnBPhase_];
    for (size_t k = 0; k < nSpeciesB_; k++) {
      spMoles_final_[k + bstart] = y[3 + nSpeciesA_ +  k];
    }

#endif
  }


  //====================================================================================================================
  
  // ----------------------------------------------------------------------------------------------
  // ----------------------------- GET CONDITIONS OUT --------------------------------------------
  // ----------------------------------------------------------------------------------------------
  //
  //       (unless specified this is always at the final conditions and time
  //
  
  // ----------------------------- GET INSTANTANEOUS SOURCE TERMS --------------------------------
  
  // ---------------------------- GET INTEGRATED SOURCE TERMS -------------------------------------

  // --------------------------- GET MOLE NUMBERS ------------------------------------------------

  // --------------------------- GET THERMO AND VOLTAGE  ------------------------------------------------

  // ------------------------- GET VOLUMES -----------------------------------------------------------

  // ---------------------- GET SURFACE AREAS -------------------------------------------------------

    
  // --------------------------- INTERNAL UPDATE FUNCTIONS --------------------------------------

  //====================================================================================================================
  //  Take the state (i.e., the final state) within the InterfacialMassTransfer_Model and push it down
  //  to the ThermoPhase Objects
  /* 
   *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
   *  objects
   */
  void imtPSS_NoSurf_DiffBL::updateState()
  {
    InterfacialMassTransfer_PseudoSS::updateState();
    double gradT = 0.0;
    ThermoPhase *tpA = & thermo(solnAPhase_);
    int astart = 0;
    for (size_t k = 0; k < nSpeciesA_; k++) {
      gradX_phaseA_[k] = (spMf_final_[k + astart] - spMF_solnA_BC_final_[k]) / BLThickness_A_final_;
    }
  
    tranA_->getSpeciesFluxes(1, &gradT,  nSpeciesA_,  DATA_PTR(gradX_phaseA_),
			     nSpeciesA_, DATA_PTR(jKmolfluxA_));
    for (size_t  k = 0; k < nSpeciesA_; k++) {
      double mwk = tpA->molecularWeight(k);
      jKmolfluxA_[k] /= mwk;
    }

    ThermoPhase *tpB = & thermo(solnBPhase_);
    size_t bstart = m_PhaseSpeciesStartIndex[solnBPhase_];
    for (size_t k = 0; k < nSpeciesB_; k++) {
      gradX_phaseB_[k] = (spMF_solnB_BC_final_[k] - spMf_final_[k + bstart]) / BLThickness_B_final_;
    }
  
    tranB_->getSpeciesFluxes(1, &gradT,  nSpeciesB_,  DATA_PTR(gradX_phaseB_),
			     nSpeciesB_, DATA_PTR(jKmolfluxB_));
    for (size_t k = 0; k < nSpeciesB_; k++) {
      double mwk = tpB->molecularWeight(k);
      jKmolfluxB_[k] /= mwk;
    }

    tpA->setState_TPX(Temp_, Pres_A_Interface_final_, DATA_PTR(spMF_solnA_BC_final_));
    tpA->getConcentrations( DATA_PTR(cBoundA_));
    densBoundA_ = tpA->density();

    tpA->setState_TPX(Temp_, Pres_A_Interface_final_, &spMf_final_[astart]);
    tpA->getConcentrations(DATA_PTR(cSurfA_));
    molarVolA_ = tpA->molarVolume();
    grossAVol_ = 0.5 *  BLThickness_A_final_ * surfaceAreaRS_final_[0];
    densA_final_ = tpA->density();

    tpB->setState_TPX(Temp_, Pres_B_Interface_final_, DATA_PTR(spMF_solnB_BC_final_));
    tpB->getConcentrations( DATA_PTR(cBoundB_));
    densBoundB_ = tpB->density();

    tpB->setState_TPX(Temp_, Pres_B_Interface_final_, &spMf_final_[bstart]);
    tpB->getConcentrations(DATA_PTR(cSurfB_));
    molarVolB_ = tpB->molarVolume();
    grossBVol_ = 0.5 *  BLThickness_B_final_ * surfaceAreaRS_final_[0];
    densB_final_ = tpB->density();
  }  

  // ---------------------------------------------------------------------------------------------
  // -------------------------------  SetState Functions -------------------------------------------------------
  // ---------------------------------------------------------------------------------------------
  void imtPSS_NoSurf_DiffBL::setInitStateFromInitInit(bool setFinal)
  {
    InterfacialMassTransfer_PseudoSS::setInitStateFromInitInit(setFinal);
    Velo_A_CVB_init_ = Velo_A_CVB_init_init_;
    Velo_B_CVB_init_ = Velo_B_CVB_init_init_;
    densA_init_ = densA_init_init_;
    densB_init_ = densB_init_init_;


    if (setFinal) {
      Velo_A_CVB_final_ = Velo_A_CVB_init_init_;
      Velo_B_CVB_final_ = Velo_B_CVB_init_init_;
      densA_final_      = densA_init_init_;
      densB_final_      = densB_init_init_;
    }
  }
  //===================================================================================================================
  void imtPSS_NoSurf_DiffBL::setInitStateFromFinal(bool setInitInit) 
  {
    InterfacialMassTransfer_PseudoSS::setInitStateFromFinal(setInitInit);
    Velo_A_CVB_init_ = Velo_A_CVB_final_;
    Velo_B_CVB_init_ = Velo_B_CVB_final_;
    densA_init_ = densA_final_;
    densB_init_ = densB_final_;
    if (setInitInit) {
      Velo_A_CVB_init_init_ = Velo_A_CVB_final_;
      Velo_B_CVB_init_init_ = Velo_B_CVB_final_;
      densA_init_init_ = densA_final_;
      densB_init_init_ = densB_final_;
    }
  }
  //===================================================================================================================
  void imtPSS_NoSurf_DiffBL::setFinalStateFromInit()
  {
    InterfacialMassTransfer_PseudoSS::setFinalStateFromInit();
    Velo_A_CVB_final_ = Velo_A_CVB_init_;
    Velo_B_CVB_final_ = Velo_B_CVB_init_;
    densA_final_ = densA_init_;
    densB_final_ = densB_init_;
  }
  //===================================================================================================================
  void imtPSS_NoSurf_DiffBL::setFinalFinalStateFromFinal() 
  {
    InterfacialMassTransfer_PseudoSS::setFinalFinalStateFromFinal();
    Velo_A_CVB_final_final_ = Velo_A_CVB_final_;
    Velo_B_CVB_final_final_ = Velo_B_CVB_final_;
  }
  // ---------------------------------------------------------------------------------------------
  // ------------------------------ PRINT ROUTINES -----------------------------------------------
  // ---------------------------------------------------------------------------------------------


  //====================================================================================================================
  void imtPSS_NoSurf_DiffBL::printInterfacialMassTransfer_List(int pSrc, bool subTimeStep) {
    printf("\n");
    printf("         PName                   MoleNum      molarVol     Volume       FractVol     Voltage   \n");
    printf("     ============================================================================================\n");
    double egv = TotalVol();
    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
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
  void imtPSS_NoSurf_DiffBL::printInterfacialMassTransfer(int pSrc, bool subTimeStep) {
    double egv = TotalVol();
    printf("   ==============================================================================================\n");
    if (subTimeStep) {
      printf("      imtPSS_NoSurf_DiffBL at intermediate-step time final = %g\n", t_final_);
      printf("                   intermediate-step time init  = %g\n", t_init_);
    } else {
      printf("      imtPSS_NoSurf_DiffBL at time final = %g\n", t_final_final_);
      printf("                   time init  = %g\n", t_init_init_);
    }
    printf("\n");
    printf("                    DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n", 
	   DomainNumber_, CellNumber_, counterNumberIntegrations_);
    printf("   ==============================================================================================\n");
    printf("          Number of surfaces = %d\n", static_cast<int>(numSurfaces_));

    printf("          Total Volume = %10.3E\n", egv);
    printf("          Temperature = %g\n", Temp_);
    printf("          Pressure = %g\n", Pres_A_Interface_final_);

    if (subTimeStep) {
      printf("          Velocities                     FINAL            INIT\n");
      printf("            Interface_Velocity    %12.5E    %12.5E \n", Velo_S_Interface_final_, Velo_S_Interface_init_);
      printf("            Velocity_A_CVB        %12.5E    %12.5E \n", Velo_A_CVB_final_, Velo_A_CVB_init_);
      printf("            Velocity_B_CVB        %12.5E    %12.5E \n", Velo_B_CVB_final_, Velo_B_CVB_init_);

    } else {
      printf("          Velocities               FINAL_FINAL       INIT_INIT\n");
      printf("            Interface_Velocity    %12.5E    %12.5E \n", Velo_S_Interface_final_final_, Velo_S_Interface_init_init_);
      printf("            Velocity_A_CVB        %12.5E    %12.5E \n", Velo_A_CVB_final_final_, Velo_A_CVB_init_init_);
      printf("            Velocity_B_CVB        %12.5E    %12.5E \n", Velo_B_CVB_final_final_, Velo_B_CVB_init_init_);
    }
    if (!subTimeStep) {
      printf("               Number of Subintegrations in last step = %d\n", counterNumberLastSubIntegrations_);
    } else {
      printf("                ID of Subintegration  step = %d\n", counterNumberLastSubIntegrations_);
    }

    printInterfacialMassTransfer_List(pSrc, subTimeStep);

    for (size_t iph = 0; iph < m_NumTotPhases; iph++) {
      printInterfacialMassTransfer_Phase(iph, pSrc, subTimeStep);   
    }


    initialPackSolver_nonlinFunction();

    if (subTimeStep) {
      std::vector<double> residNLS(neq_);
      evalResidNJ(t_final_,deltaTsubcycleCalc_,
			      &yvalNLS_[0], &ydotNLS_[0], &residNLS[0],  Base_ShowSolution);
      
    }

  }
  //===================================================================================================================
 
  void imtPSS_NoSurf_DiffBL::printInterfacialMassTransfer_Phase(size_t iph, int pSrc, bool subTimeStep) {
    size_t isph;
    double *netROP = new double[m_NumTotSpecies];
    ThermoPhase &tp = thermo(iph);
    std::string pname = tp.id();
    size_t istart = m_PhaseSpeciesStartIndex[iph];
    size_t nsp = tp.nSpecies();
    if (printLvl_ <= 1) {
      return;
    }
    printf("     ============================================================================================\n");
    printf("          Phase %d %s \n", static_cast<int>(iph), pname.c_str() );
    printf("                Total moles = %g\n", phaseMoles_final_[iph]);
 

    /*
     * Do specific surface phase printouts
     */
    if (iph >= m_NumVolPhases) {
      isph = iph - m_NumVolPhases;
      printf("                surface area (final) = %g\n",  surfaceAreaRS_final_[isph]);
      printf("                surface area (init)  = %g\n",  surfaceAreaRS_init_[isph]);
    }
    if (printLvl_ >= 3) {
      printf("\n");
      printf("                Name              MoleFrac_final kMoles_final kMoles_init SrcTermIntegrated(kmol)\n");
      for (size_t k = 0; k < nsp; k++) {
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

    // Residual satisfaction



    if (printLvl_ >= 4) {
      if (iph >= m_NumVolPhases) {
	const std::vector<double> &rsSpeciesProductionRates = RSD_List_[isph]->calcNetProductionRates();
	RSD_List_[isph]->getNetRatesOfProgress(netROP);
      
	doublevalue* spNetProdPerArea = (doublevalue*) spNetProdPerArea_List_.ptrColumn(isph);
	mdpUtil::mdp_zero_dbl_1(spNetProdPerArea, m_NumTotSpecies);
	size_t nphRS = RSD_List_[isph]->nPhases();
	size_t kIndexKin = 0;
	for (size_t kph = 0; kph < nphRS; kph++) {
	  size_t jph = RSD_List_[isph]->kinOrder[kph];
	  size_t istart = m_PhaseSpeciesStartIndex[jph];
	  size_t nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	  for (size_t  k = 0; k < nsp; k++) {
	    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
	    kIndexKin++;
	  }
	}
	printf("\n");
	printf("                           spName                  SourceRateLastStep (kmol/m2/s) \n");
	for (size_t k = 0; k <  m_NumTotSpecies; k++) {
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
  void imtPSS_NoSurf_DiffBL::writeCSVData(int itype) { 
    if (printLvl_ < 2) return;
    size_t k;
    static std::string globOutputName;
    static std::string intOutputName;
    static FILE *fpG = 0;
    static FILE *fpI = 0;

    if (itype == 0 || fpI == 0) {
      intOutputName = (Title_ + "_intResults_" + int2str(DomainNumber_) + "_" +
		       int2str(CellNumber_) + ".csv");
       fpI = fopen(intOutputName.c_str(), "w");

       fprintf(fpI, "         Tinit ,        Tfinal ," );

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
