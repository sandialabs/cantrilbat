/*
 * $Id: InterfacialMassTransfer.cpp 520 2013-01-15 21:28:53Z vebruni $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"


#include "cantera/equilibrium.h"

#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/thermo/FixedChemPotSSTP.h"

#include "cantera/numerics/solveProb.h"
#include "cantera/numerics/BEulerInt.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/SurfPhase.h"

#include "cantera/solvers.h"


#include "PhaseList.h"


#include "BlockEntryGlobal.h"

#include "InterfacialMassTransfer.h"
#include "ReactingSurDomain.h"
#include "importAllCTML.h"
#include "RxnMolChange.h"
#include "ExtraGlobalRxn.h"

#include "InterfacialMassTransfer_input.h"

#include "ApplBase_print.h"

using namespace Cantera;
using namespace std;
using namespace BEInput;
using namespace TKInput;
using namespace mdpUtil;


#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif  

namespace Cantera {


  //======================================================================================================================
  /*
   *  IMT_INPUT: constructor
   *
   *  We initialize the arrays in the structure to the appropriate sizes.
   *  And, we initialize all of the elements of the arrays to defaults.  
   */
  InterfacialMassTransfer::InterfacialMassTransfer() :
    PhaseList(),
    pendingIntegratedStep_(0),
    includeABMolesInSolnVector_(0),
    PhaseIndex_mp_(0),
    interfaceType_(BASE_TYPE_IMT),
    externFieldTimeBehaviorType_(EF_TIMEBEHAVIOR_TFINALCONSTANT),
    t_init_init_(0.0),
    t_final_final_(0.0),
    t_init_(0.0),
    t_final_(0.0),
    deltaTsubcycleMax_(1.0E300),
    deltaTsubcycle_init_init_(1.0E300),
    deltaTsubcycleNext_(1.0E300),
    deltaTsubcycle_init_next_(1.0E300),
    choiceDeltaTsubcycle_init_(0),

    Temp_(298.15),
    Pres_A_Interface_init_init_(1.0E5),
    Pres_A_Interface_init_(1.0E5),
    Pres_A_Interface_final_(1.0E5),
    Pres_A_Interface_final_final_(1.0E5),
    Pres_B_Interface_init_init_(1.0E5),
    Pres_B_Interface_init_(1.0E5),
    Pres_B_Interface_final_(1.0E5),
    Pres_B_Interface_final_final_(1.0E5),
    Velo_S_Interface_init_init_(0.0),
    Velo_S_Interface_init_(0.0),
    Velo_S_Interface_final_(0.0),
    Velo_S_Interface_final_final_(0.0),
    Velo_A_Interface_init_init_(0.0),
    Velo_A_Interface_init_(0.0),
    Velo_A_Interface_final_(0.0),
    Velo_A_Interface_final_final_(0.0),
    Velo_B_Interface_init_init_(0.0),
    Velo_B_Interface_init_(0.0),
    Velo_B_Interface_final_(0.0),
    Velo_B_Interface_final_final_(0.0),
    Velo_ReferenceFrame_(0.0),
    
    BLThickness_A_init_(0.0),
    BLThickness_A_final_(0.0),
    BLThickness_A_init_init_(0.0),
    BLThickness_A_final_final_(0.0),
    BLThickness_B_init_(0.0),
    BLThickness_B_final_(0.0),
    BLThickness_B_init_init_(0.0),
    BLThickness_B_final_final_(0.0),

    spMoles_init_init_(0),
    spMoles_init_(0),
    spMoles_final_(0),
    spMoles_final_final_(0),
    
 
    spMf_init_init_(0),
    spMf_init_(0),
    spMf_final_(0),
    spMoles_dot_(0),
    spMoles_predict_(0),
    phaseMolarVolumes_(0),

    Pres_solnA_BC_init_init_(1.0E5),
    Pres_solnA_BC_init_(1.0E5),
    Pres_solnA_BC_final_(1.0E5),
    Pres_solnA_BC_final_final_(1.0E5),
    Pres_solnB_BC_init_init_(1.0E5),
    Pres_solnB_BC_init_(1.0E5),
    Pres_solnB_BC_final_(1.0E5),
    Pres_solnB_BC_final_final_(1.0E5),
    
    spMF_solnA_BC_init_init_(0),
    spMF_solnA_BC_init_(0),
    spMF_solnA_BC_final_(0),
    spMF_solnA_BC_final_final_(0),

    spMF_solnB_BC_init_init_(0),
    spMF_solnB_BC_init_(0),
    spMF_solnB_BC_final_(0),
    spMF_solnB_BC_final_final_(0),
  
    VolPM_(0),
    spChemPot_(0),
    RSD_List_(0),
    ActiveKineticsSurf_(0),
    phaseMoles_init_(0),
    phaseMoles_init_init_(0),
    phaseMoles_final_(0),
    phaseMoles_dot_(0),
    numSurfaces_(0),
    surfaceAreaRS_init_init_(0),
    surfaceAreaRS_init_(0),
    surfaceAreaRS_final_(0),
    spNetProdPerArea_List_(0, 0),
    spMoleIntegratedSourceTerm_(0),
    spMoleIntegratedSourceTermLast_(0),
    Title_("InterfacialMassTransfer"),
 
    solnAPhase_(-1),
    nSpeciesA_(0),
    solnBPhase_(-1),
    nSpeciesB_(0),
    deltaG_(0),
    molarAtol_(1.0E-16),
    DomainNumber_(0),
    CellNumber_(0),
    counterNumberIntegrations_(0),
    counterNumberSubIntegrations_(0),
    counterNumberLastSubIntegrations_(0),
    printLvl_(4),
    detailedResidPrintFlag_(0),
    enableExtraPrinting_(false)
  {

    //m_pl = new PhaseList();
    //m_rSurDomain = new ReactingSurDomain();
  }

  //======================================================================================================================
  // Copy Constructor
  /*
   * @param right Object to be copied
   */
  InterfacialMassTransfer::InterfacialMassTransfer(const InterfacialMassTransfer &right) :
    PhaseList(),  
    pendingIntegratedStep_(0), 
    includeABMolesInSolnVector_(0),
    prob_type(TP),
    PhaseIndex_mp_(0),
    interfaceType_(BASE_TYPE_IMT),
    externFieldTimeBehaviorType_(EF_TIMEBEHAVIOR_TFINALCONSTANT),
    t_init_init_(0.0),
    t_final_final_(0.0),
    t_init_(0.0),
    t_final_(0.0),
    deltaTsubcycleMax_(1.0E300),
    deltaTsubcycle_init_init_(1.0E300),
    deltaTsubcycleNext_(1.0E300),
    deltaTsubcycle_init_next_(1.0E300),
    choiceDeltaTsubcycle_init_(0),

    Temp_(298.15),
    Pres_A_Interface_init_init_(1.0E5),
    Pres_A_Interface_init_(1.0E5),
    Pres_A_Interface_final_(1.0E5),
    Pres_A_Interface_final_final_(1.0E5),
    Pres_B_Interface_init_init_(1.0E5),
    Pres_B_Interface_init_(1.0E5),
    Pres_B_Interface_final_(1.0E5),
    Pres_B_Interface_final_final_(1.0E5),
    Velo_S_Interface_init_init_(0.0),
    Velo_S_Interface_init_(0.0),
    Velo_S_Interface_final_(0.0),
    Velo_S_Interface_final_final_(0.0),
    Velo_A_Interface_init_init_(0.0),
    Velo_A_Interface_init_(0.0),
    Velo_A_Interface_final_(0.0),
    Velo_A_Interface_final_final_(0.0),
    Velo_B_Interface_init_init_(0.0),
    Velo_B_Interface_init_(0.0),
    Velo_B_Interface_final_(0.0),
    Velo_B_Interface_final_final_(0.0),
    Velo_ReferenceFrame_(0.0),

    BLThickness_A_init_(0.0),
    BLThickness_A_final_(0.0),
    BLThickness_A_init_init_(0.0),
    BLThickness_A_final_final_(0.0),
    BLThickness_B_init_(0.0),
    BLThickness_B_final_(0.0),
    BLThickness_B_init_init_(0.0),
    BLThickness_B_final_final_(0.0),


    spMoles_init_init_(0),
    spMoles_init_(0),
    spMoles_final_(0),
    spMoles_final_final_(0),
    
 
    spMf_init_init_(0),
    spMf_init_(0),
    spMf_final_(0),

    spMoles_dot_(0),
    spMoles_predict_(0),
 

    phaseMolarVolumes_(0),

    Pres_solnA_BC_init_init_(1.0E5),
    Pres_solnA_BC_init_(1.0E5),
    Pres_solnA_BC_final_(1.0E5),
    Pres_solnA_BC_final_final_(1.0E5),
    Pres_solnB_BC_init_init_(1.0E5),
    Pres_solnB_BC_init_(1.0E5),
    Pres_solnB_BC_final_(1.0E5),
    Pres_solnB_BC_final_final_(1.0E5),
    
    spMF_solnA_BC_init_init_(0),
    spMF_solnA_BC_init_(0),
    spMF_solnA_BC_final_(0),
    spMF_solnA_BC_final_final_(0),

    spMF_solnB_BC_init_init_(0),
    spMF_solnB_BC_init_(0),
    spMF_solnB_BC_final_(0),
    spMF_solnB_BC_final_final_(0),

    VolPM_(0),
  
    spChemPot_(0),
    RSD_List_(0),
    ActiveKineticsSurf_(0),
    phaseMoles_init_(0),
    phaseMoles_init_init_(0),
    phaseMoles_final_(0),
    phaseMoles_dot_(0),
    numSurfaces_(0),
    surfaceAreaRS_init_init_(0),
    surfaceAreaRS_init_(0),
    surfaceAreaRS_final_(0),
    spNetProdPerArea_List_(0, 0),
    spMoleIntegratedSourceTerm_(0),
    spMoleIntegratedSourceTermLast_(0),
    Title_("InterfacialMassTransfer"),
    solnAPhase_(-1),
    nSpeciesA_(0),
    solnBPhase_(-1),
    nSpeciesB_(0),
    deltaG_(0),

    molarAtol_(1.0E-16),
    DomainNumber_(0),
    CellNumber_(0),
    counterNumberIntegrations_(0),
    counterNumberSubIntegrations_(0),
    counterNumberLastSubIntegrations_(0),
    printLvl_(4),
    detailedResidPrintFlag_(0),
    enableExtraPrinting_(false)
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
  InterfacialMassTransfer & InterfacialMassTransfer::operator=(const InterfacialMassTransfer &right)
  {
    /*
     * Check for self assignment.
     */
    if (this == &right) return *this;

    PhaseList::operator=(right);
    /*
     * We do a straight assignment operator on all of the
     * data. The vectors are copied.
     */
  
    /*
     * Copy over the ReactingSurDomain list
     */
    for (int i = 0; i < m_NumSurPhases; i++) {
      bool idHit = false;
      //if (right.m_rSurDomain == right.RSD_List_[i]) {
      //	idHit = true;
      //}
      if (RSD_List_[i]) {
	delete RSD_List_[i];
	RSD_List_[i] = 0;
      }
      if (right.RSD_List_[i]) {
	RSD_List_[i] = new ReactingSurDomain(*(right.RSD_List_[i]));
	if (idHit) {
	  // m_rSurDomain = RSD_List_[i];
	}
      }
    }
    pendingIntegratedStep_ = right.pendingIntegratedStep_;
    prob_type = right.prob_type;

    PhaseIndex_mp_ = right.PhaseIndex_mp_;
    interfaceType_ = right.interfaceType_;
    externFieldTimeBehaviorType_ = right.externFieldTimeBehaviorType_;
    t_init_init_ = right.t_init_init_;
    t_final_final_ = right.t_final_final_;
    t_init_ = right.t_init_;
    t_final_ = right.t_final_;
    deltaTsubcycleMax_ = right.deltaTsubcycleMax_;
    deltaTsubcycle_init_init_ = right.deltaTsubcycle_init_init_;
    deltaTsubcycleNext_ = right.deltaTsubcycleNext_;
    deltaTsubcycle_init_next_ = right.deltaTsubcycle_init_next_;
    choiceDeltaTsubcycle_init_ = right.choiceDeltaTsubcycle_init_;

    Temp_                                  = right.Temp_;
    Pres_A_Interface_init_init_            = right.Pres_A_Interface_init_init_;
    Pres_A_Interface_init_                 = right. Pres_A_Interface_init_;
    Pres_A_Interface_final_                = right.Pres_A_Interface_final_;
    Pres_A_Interface_final_final_          = right.Pres_A_Interface_final_final_;
    Pres_B_Interface_init_init_            = right.Pres_B_Interface_init_init_;
    Pres_B_Interface_init_                 = right.Pres_B_Interface_init_;
    Pres_B_Interface_final_                = right.Pres_B_Interface_final_;
    Pres_B_Interface_final_final_          = right.Pres_B_Interface_final_final_;
    Velo_S_Interface_init_init_            = right.Velo_S_Interface_init_init_;
    Velo_S_Interface_init_                 = right.Velo_S_Interface_init_;
    Velo_S_Interface_final_                = right.Velo_S_Interface_final_;
    Velo_S_Interface_final_final_          = right.Velo_S_Interface_final_final_;
    Velo_A_Interface_init_init_            = right.Velo_A_Interface_init_init_;
    Velo_A_Interface_init_                 = right.Velo_A_Interface_init_;
    Velo_A_Interface_final_                = right.Velo_A_Interface_final_;
    Velo_A_Interface_final_final_          = right.Velo_A_Interface_final_final_;
    Velo_B_Interface_init_init_            = right.Velo_B_Interface_init_init_;
    Velo_B_Interface_init_                 = right.Velo_B_Interface_init_;
    Velo_B_Interface_final_                = right.Velo_B_Interface_final_;
    Velo_B_Interface_final_final_          = right.Velo_B_Interface_final_final_;
    Velo_ReferenceFrame_                   = right.Velo_ReferenceFrame_;

    BLThickness_A_init_                    = right.BLThickness_A_init_;
    BLThickness_A_final_                   = right.BLThickness_A_final_;
    BLThickness_A_init_init_               = right.BLThickness_A_init_init_;
    BLThickness_A_final_final_             = right.BLThickness_A_final_final_;
    BLThickness_B_init_                    = right.BLThickness_B_init_;
    BLThickness_B_final_                   = right.BLThickness_B_final_;
    BLThickness_B_init_init_               = right.BLThickness_B_init_init_;
    BLThickness_B_final_final_             = right.BLThickness_B_final_final_;


    spMoles_init_init_                     = right.spMoles_init_init_;
    spMoles_init_                          = right.spMoles_init_;
    spMoles_final_                         = right.spMoles_final_;
    spMoles_final_final_                   = right.spMoles_final_final_;
   
    spMf_init_init_                        = right.spMf_init_init_;
    spMf_init_                             = right.spMf_init_;
    spMf_final_                            = right.spMf_final_;

    spMoles_dot_                           = right.spMoles_dot_;
    spMoles_predict_                       = right.spMoles_predict_;
    phaseMolarVolumes_                     = right.phaseMolarVolumes_;

    Pres_solnA_BC_init_init_               = right.Pres_solnA_BC_init_init_;
    Pres_solnA_BC_init_                    = right.Pres_solnA_BC_init_;
    Pres_solnA_BC_final_                   = right.Pres_solnA_BC_final_;
    Pres_solnA_BC_final_final_             = right.Pres_solnA_BC_final_final_;
    Pres_solnB_BC_init_init_               = right.Pres_solnB_BC_init_init_;
    Pres_solnB_BC_init_                    = right.Pres_solnB_BC_init_;
    Pres_solnB_BC_final_                   = right.Pres_solnB_BC_final_;
    Pres_solnB_BC_final_final_             = right.Pres_solnB_BC_final_final_;

    spMF_solnA_BC_init_init_               = right.spMF_solnA_BC_init_init_;
    spMF_solnA_BC_init_                    = right.spMF_solnA_BC_init_;
    spMF_solnA_BC_final_                   = right.spMF_solnA_BC_final_;
    spMF_solnA_BC_final_final_             = right.spMF_solnA_BC_final_final_;
    spMF_solnB_BC_init_init_               = right.spMF_solnB_BC_init_init_;
    spMF_solnB_BC_init_                    = right.spMF_solnB_BC_init_;
    spMF_solnB_BC_final_                   = right.spMF_solnB_BC_final_;
    spMF_solnB_BC_final_final_             = right.spMF_solnB_BC_final_final_;

    VolPM_ = right.VolPM_;
    ActiveKineticsSurf_ = right.ActiveKineticsSurf_;
    phaseMoles_init_ = right.phaseMoles_init_;
    phaseMoles_init_init_ = right.phaseMoles_init_init_;
    phaseMoles_final_ = right.phaseMoles_final_;
    phaseMoles_dot_ = right.phaseMoles_dot_;

    numSurfaces_ = right.numSurfaces_;
    surfaceAreaRS_init_ = right.surfaceAreaRS_init_;
    surfaceAreaRS_init_init_ = right.surfaceAreaRS_init_init_;
    surfaceAreaRS_final_ = right.surfaceAreaRS_final_;
    spNetProdPerArea_List_ = right.spNetProdPerArea_List_;
    spMoleIntegratedSourceTerm_ = right.spMoleIntegratedSourceTerm_;
    spMoleIntegratedSourceTermLast_ = right.spMoleIntegratedSourceTermLast_;
    Title_ = right.Title_;
  

    solnAPhase_ = right.solnAPhase_;
    nSpeciesA_ = right.nSpeciesA_;
    solnBPhase_ = right.solnBPhase_;
    nSpeciesB_ = right.nSpeciesB_;

    deltaG_ = right.deltaG_;

    molarAtol_ = right.molarAtol_;
 
    DomainNumber_ = right.DomainNumber_;
    CellNumber_ = right.CellNumber_;
    printLvl_ = right.printLvl_;
    detailedResidPrintFlag_ = right.detailedResidPrintFlag_;
    enableExtraPrinting_ = right.enableExtraPrinting_;
    counterNumberIntegrations_ = 0;
    counterNumberSubIntegrations_ = 0;
    counterNumberLastSubIntegrations_ = 0;

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
  InterfacialMassTransfer::~InterfacialMassTransfer() 
  {
  
 
 
    for (int i = 0; i < m_NumSurPhases; i++) {
      if (RSD_List_[i]) {
	delete RSD_List_[i];
      }
    }
    //m_rSurDomain = 0;

 
  }  
  //======================================================================================================================
  // Return the type of interfacial mass transport object 
  /*
   *  Returns the enum type of the object. This is used in the factory routine.
   *
   *  @return Returns an enum type, called   IMT_Types_Enum
   */
  IMT_Types_Enum InterfacialMassTransfer::imtType() const
  {
    return BASE_TYPE_IMT; 
  }
  //======================================================================================================================
  // Set the electrode ID information
  void InterfacialMassTransfer::setID(int domainNum, int cellNum)
  {
    DomainNumber_ = domainNum;
    CellNumber_ = cellNum;
  }
  //======================================================================================================================
  void ErrorModelType(int pos, std::string actual, std::string expected) {
    throw Cantera::CanteraError("InterfacialMassTransfer::electrode_model_create() model id",
				"At pos " + int2str(pos) + ", expected phase " + expected + 
		       " but got phase " + actual);
  }
  void pmatch(std::vector<std::string> &pn, int pos, std::string expected) {
    if (pos <= ((int) pn.size() - 1)) {
      if (pn[pos] != expected) {
	ErrorModelType(pos, pn[pos], expected);
      }
    }
  }
  //======================================================================================================================
  //  Setup the electrode
  /*
   * @param ei    ELECTRODE_KEY_INPUT pointer object
   */
  int 
  InterfacialMassTransfer::model_create(IMT_KEY_INPUT *ei) {

    int i, iph;

    // use the assignment operator to transfer for now
    // May get more sophisticated
    PhaseList::operator=(*(ei->m_pl));

    // resize m_NumTotPhases = total number of phase vectors
    phaseMoles_init_.resize(m_NumTotPhases, 0.0);
    phaseMoles_init_init_.resize(m_NumTotPhases, 0.0);
    phaseMoles_final_.resize(m_NumTotPhases, 0.0);
    phaseMoles_dot_.resize(m_NumTotPhases, 0.0);
    phaseMolarVolumes_.resize(m_NumTotPhases, 0.0);

 
    // resize volume phase vectors

    PhaseIndex_mp_.resize(NumVolPhases_, -1);

    // resize surface phase vectors
    numSurfaces_ = m_NumSurPhases;
    surfaceAreaRS_init_init_.resize(m_NumSurPhases, 0.0);
    surfaceAreaRS_init_.resize(m_NumSurPhases, 0.0);
    surfaceAreaRS_final_.resize(m_NumSurPhases, 0.0);
    RSD_List_.resize(m_NumSurPhases, 0);
    numRxns_.resize(m_NumSurPhases,0);
    ActiveKineticsSurf_.resize(m_NumSurPhases, 0);

    // Resize species vectors
    spMoles_init_.resize(m_NumTotSpecies, 0.0);
    spMoles_final_.resize(m_NumTotSpecies, 0.0);
    spMoles_final_final_.resize(m_NumTotSpecies, 0.0);
    spMoles_init_init_.resize(m_NumTotSpecies, 0.0);
    spMoles_dot_.resize(m_NumTotSpecies, 0.0);
    spMoles_predict_.resize(m_NumTotSpecies, 0.0);
    spMf_final_.resize(m_NumTotSpecies, 0.0);
    spMf_init_.resize(m_NumTotSpecies, 0.0);
    spMf_init_init_.resize(m_NumTotSpecies, 0.0);

    spMF_solnA_BC_init_.resize(m_NumTotSpecies, 0.0);
    spMF_solnA_BC_init_init_.resize(m_NumTotSpecies, 0.0);
    spMF_solnA_BC_final_.resize(m_NumTotSpecies, 0.0);
    spMF_solnA_BC_final_final_.resize(m_NumTotSpecies, 0.0);

    spMF_solnB_BC_init_.resize(m_NumTotSpecies, 0.0);
    spMF_solnB_BC_init_init_.resize(m_NumTotSpecies, 0.0);
    spMF_solnB_BC_final_.resize(m_NumTotSpecies, 0.0);
    spMF_solnB_BC_final_final_.resize(m_NumTotSpecies, 0.0);

    VolPM_.resize(m_NumTotSpecies, 0.0);
    spChemPot_.resize(m_NumTotSpecies, 0.0);

    /*
     * OK, Find the first kinetics object
     */   
    int isfound = -1;
    //int ivfound = -1;
    for (i = 0; i < m_NumSurPhases; i++) {
      if (SurPhaseHasKinetics[i]) {
	isfound = i;
	break;
      }
    }
    if (isfound == -1) {
      for (i = 0; i < NumVolPhases_; i++) {
	if (VolPhaseHasKinetics[i]) {
	  //ivfound = i;
	  break;
	}
      }
    }
 
    /*
     * Assign a reacting surface to each surface in the PhaseList object
     */
    for (i = 0; i < m_NumSurPhases; i++) {
      if (SurPhaseHasKinetics[i]) {
	ReactingSurDomain *rsd = new ReactingSurDomain();
	int ok = rsd->importFromPL(this, -1, i);
	if (!ok) {
	  throw CanteraError("cttables main:",
			     "rSurDomain returned an error"); 
	}
	// We carry a list of pointers to ReactingSurDomain
	RSD_List_[i] = rsd;
	numRxns_[i] = rsd->nReactions();
	ActiveKineticsSurf_[i] = 1;
      }
    }

  

 
    /*
     * Resize the species production vector for the electrode
     *  -> per area and total net
     */
    spMoleIntegratedSourceTerm_.resize(m_NumTotSpecies, 0.0);
    spMoleIntegratedSourceTermLast_.resize(m_NumTotSpecies, 0.0);
    DspMoles_final_.resize(m_NumTotSpecies, 0.0);
    spNetProdPerArea_List_.resize(m_NumTotSpecies, m_NumSurPhases, 0.0);

    /*
     * Load up the temperature and pressure
     */
    Temp_ = ei->Temperature;
    Pres_A_Interface_final_ = ei->PressureA;
    Pres_B_Interface_final_ = ei->PressureA;

  
    /*
     *  Loop Over all phases in the PhaseList, adding these
     *  formally to the InterfacialMassTransfer object.
     */
    int nspecies = 0;
    for (iph = 0; iph < m_NumTotPhases; iph++) {
  
      ThermoPhase *tphase = &(thermo(iph));
      int nSpecies = tphase->nSpecies();
  

      // Find the name of the input block
      string phaseBath = "Bath Specification for Phase "; 
      string phaseNm = tphase->name();
      phaseBath += phaseNm; 

      /*
       * Search the ReactingSurDomain to see if the current phase is in
       * the object
       */
      string pname = tphase->id();
   
  
      int kstart = nspecies;
      nspecies += nSpecies;
    
      /*
       *  We copy the mole numbers and mole fraction information from the bath gas section of the input file
       *  here. Note, the total quantitites of electrode material and electrolyte material may be overwritten
       *  later.
       */
      double sum = 0.0;
      for (int k = 0; k < nSpecies; k++) {
	spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] = ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
	spMoles_init_[m_PhaseSpeciesStartIndex[iph] + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k];
	sum += ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
      }
       
      if (sum > 0.0) {
	for (int k = 0; k < nSpecies; k++) {
	  spMf_final_[kstart + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] / sum;
	  spMf_init_[kstart + k] = spMf_final_[kstart + k];
	}
	tphase->setMoleFractions(&(spMf_final_[kstart]));
      } else {
	tphase->getMoleFractions(&(spMf_final_[kstart]));
	tphase->getMoleFractions(&(spMf_init_[kstart]));
      }
    }

    /* set the absolute tolerance to 1e-5 the total number of moles */
    // HEWSON -- 12/1/10 -- this particular tolerance value 
    // and all magic tolerance values throughout are subject to revision
    double tMoles = 0.0;
    for (int k = 0; k < m_NumTotSpecies; k++) {
      tMoles += spMoles_final_[k];
    }
    molarAtol_ = tMoles * 1.0E-5;

    /*
     * Now that spMoles_final_ is sized appropriately, we can call updatePhaseNumbers()
     */
    for (iph = 0; iph < m_NumTotPhases; iph++) {
      updatePhaseNumbers(iph);
    }

    if (strlen(ei->Title.c_str()) > 0) {
      Title_ = ei->Title;
    }
 
    /*
     * Set up the MultiPhase object. Right now it will contain all of the
     * volume phases only.
     */
    for (iph = 0; iph < NumVolPhases_; iph++) {
      PhaseIndex_mp_[iph] = iph;
      phaseMoles_init_[iph] = phaseMoles_final_[iph];
      phaseMoles_init_init_[iph] = phaseMoles_final_[iph];
    }

    /*
     * Identify the phases which will ID the  left and right sides of the domain
     */

    solnAPhase_ = -1;
    solnBPhase_ = -1;

   
    //   InterfaceKinetics *iK = m_rSurDomain;
    //ReactingSurDomain *rsd = RSD_List_[0];
    //   std::vector<RxnMolChange *> & rmcV = rsd->rmcVector;
  
    for (iph = 0; iph < m_NumTotPhases; iph++) {
      ThermoPhase *tp = & (thermo(iph));
      if (tp->name() == ei->PhaseAName) {
        solnAPhase_ = iph;
	nSpeciesA_ = tp->nSpecies();
      }
      if (tp->name() == ei->PhaseBName) {
        solnBPhase_ = iph;
	nSpeciesB_ = tp->nSpecies();
      }
    }
    if (solnAPhase_ == -1) {
      throw CanteraError("InterfacialMassTransfer::model_create()", "Couldn't find solnA phase");
    }
    if (solnBPhase_ == -1) {
      throw CanteraError("InterfacialMassTransfer::model_create()", "Couldn't find solnB phase");
    }
    if (solnAPhase_ == solnBPhase_) {
      throw CanteraError("InterfacialMassTransfer::model_create()", "Phase A and Phase B are the same");
    }


    for (int i = 0; i < m_NumSurPhases; i++) {
      surfaceAreaRS_final_[i] = ei->SurfaceArea;
      surfaceAreaRS_init_[i] = ei->SurfaceArea;
      surfaceAreaRS_init_init_[i] = ei->SurfaceArea;
    }

 
  
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(phaseMoles_init_), DATA_PTR(phaseMoles_final_), m_NumTotPhases);
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(phaseMoles_init_init_), DATA_PTR(phaseMoles_final_), m_NumTotPhases);
    
 
  

    InterfacialMassTransfer::updateState();


    mdpUtil::mdp_copy_dbl_1(DATA_PTR(spMoles_init_init_), DATA_PTR(spMoles_final_), m_NumTotSpecies);
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(spMoles_final_final_), DATA_PTR(spMoles_final_), m_NumTotSpecies);

    /*
     *  Get boundary layer thicknesses
     */
    BLThickness_A_final_ = ei->BLThickness_A;
    BLThickness_B_final_ = ei->BLThickness_B;

    /*
     *  Resize A
     */
    double volA = ei->SurfaceArea * ei->BLThickness_A;
    double sum = 0.0;
    double curr_volA = 0.0;
    ThermoPhase *tp = VolPhaseList[solnAPhase_];
    double mv = tp->molarVolume();
    double phaseM = phaseMoles_final_[solnAPhase_];
    int numSpecA = tp->nSpecies();
    int kstart = getGlobalSpeciesIndex(solnAPhase_, 0);

    if (phaseM <= 0.0) {
       phaseM = volA / mv;
       for (int k = 0; k < numSpecA; k++) {
        spMoles_final_[kstart + k] = phaseM * spMf_final_[kstart + k];
	sum += spMoles_final_[kstart + k];
      }
    } else {
      curr_volA = mv * phaseM;
      double ratio = volA / curr_volA;
      for (int k = 0; k < numSpecA; k++) {
        spMoles_final_[kstart + k] *= ratio; 
	sum += spMoles_final_[kstart + k];
      }
    }
    phaseMoles_final_[solnAPhase_] = sum;

    /*
     *  Resize B
     */
    double volB = ei->SurfaceArea * ei->BLThickness_B;
    sum = 0.0;
    double curr_volB = 0.0;
    tp = VolPhaseList[solnBPhase_];
    mv = tp->molarVolume();
    phaseM = phaseMoles_final_[solnBPhase_];
    kstart = getGlobalSpeciesIndex(solnBPhase_, 0);
    int numSpecB = tp->nSpecies();
    if (phaseM <= 0.0) {
       phaseM = volB / mv;
       for (int k = 0; k < numSpecB; k++) {
        spMoles_final_[kstart + k] = phaseM * spMf_final_[kstart + k];
	sum += spMoles_final_[kstart + k];
      }
    } else {
      curr_volB = mv * phaseM;
      double ratio = volB / curr_volB;
      for (int k = 0; k < numSpecB; k++) {
        spMoles_final_[kstart + k] *= ratio;
	sum += spMoles_final_[kstart + k];
      }
    }
    phaseMoles_final_[solnBPhase_] = sum;

    for (int i = 0; i < m_NumSurPhases; i++) {
      int iph = NumVolPhases_ + i;
      tp = SurPhaseList[i];
      SurfPhase *sp = dynamic_cast<SurfPhase *>(tp);
      double sd = sp->siteDensity();
      phaseM = phaseMoles_final_[iph];
      phaseMoles_final_[iph] = sd *  ei->SurfaceArea;

      kstart = getGlobalSpeciesIndex(iph, 0);
      int numSpec = tp->nSpecies();
      if (phaseM <= 0.0) {
	phaseM = sd *  ei->SurfaceArea;
	for (int k = 0; k < numSpec; k++) {
	  spMoles_final_[kstart + k] = phaseM * spMf_final_[kstart + k];
	}
      } else { 
	double curr_sur = phaseM / sd;
	 
	double ratio = ei->SurfaceArea / curr_sur;
	for (int k = 0; k < numSpecB; k++) {
	  spMoles_final_[kstart + k] *= ratio;
	}
      }
    }


    InterfacialMassTransfer::updateState();
    InterfacialMassTransfer::setInitStateFromFinal(true);
    setFinalFinalStateFromFinal_Oin();

    return 0;
  }
  //====================================================================================================================
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
  int InterfacialMassTransfer::setInitialConditions(IMT_KEY_INPUT *ei) {

     /*
     * Load up the temperature and pressure
     */
    Temp_ = ei->Temperature;
    Pres_A_Interface_final_ = ei->PressureA;
    /*
     *  Loop Over all phases in the PhaseList, adding these
     *  formally to the InterfacialMassTransfer object.
     */
    int nspecies = 0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
  
      ThermoPhase *tphase = &(thermo(iph));
      int nSpecies = tphase->nSpecies();
  

      // Find the name of the input block
      string phaseBath = "Bath Specification for Phase "; 
      string phaseNm = tphase->name();
      phaseBath += phaseNm; 

      /*
       * Search the ReactingSurDomain to see if the current phase is in
       * the object
       */
      string pname = tphase->id();
   
  
      int kstart = nspecies;
      nspecies += nSpecies;
    
      /*
       *  We copy the mole numbers and mole fraction information from the bath gas section of the input file
       *  here. Note, the total quantitites of electrode material and electrolyte material may be overwritten
       *  later.
       */
      double sum = 0.0;
      for (int k = 0; k < nSpecies; k++) {
	spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] = ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
	spMoles_init_[m_PhaseSpeciesStartIndex[iph] + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k];
	sum += ei->MoleNumber[m_PhaseSpeciesStartIndex[iph] + k];
      }
       
      if (sum > 0.0) {
	for (int k = 0; k < nSpecies; k++) {
	  spMf_final_[kstart + k] = spMoles_final_[m_PhaseSpeciesStartIndex[iph] + k] / sum;
	  spMf_init_[kstart + k] = spMf_final_[kstart + k];
	}
	tphase->setMoleFractions(&(spMf_final_[kstart]));
      } else {
	tphase->getMoleFractions(&(spMf_final_[kstart]));
	tphase->getMoleFractions(&(spMf_init_[kstart]));
      }
    }

    for (int i = 0; i < m_NumSurPhases; i++) {
      surfaceAreaRS_final_[i] = ei->SurfaceArea;
      surfaceAreaRS_init_[i] = ei->SurfaceArea;
      surfaceAreaRS_init_init_[i] = ei->SurfaceArea;
    }

    InterfacialMassTransfer::updateState();
   /*
     *  Resize A
     */
    double volA = ei->SurfaceArea * ei->BLThickness_A;
    double sum = 0.0;
    double curr_volA = 0.0;
    ThermoPhase *tp = VolPhaseList[solnAPhase_];
    double mv = tp->molarVolume();
    double phaseM = phaseMoles_final_[solnAPhase_];
    int numSpecA = tp->nSpecies();
    int kstart = getGlobalSpeciesIndex(solnAPhase_, 0);

    if (phaseM <= 0.0) {
       phaseM = volA / mv;
       for (int k = 0; k < numSpecA; k++) {
        spMoles_final_[kstart + k] = phaseM * spMf_final_[kstart + k];
	sum += spMoles_final_[kstart + k];
      }
    } else {
      curr_volA = mv * phaseM;
      double ratio = volA / curr_volA;
      for (int k = 0; k < numSpecA; k++) {
        spMoles_final_[kstart + k] *= ratio; 
	sum += spMoles_final_[kstart + k];
      }
    }
    phaseMoles_final_[solnAPhase_] = sum;

    /*
     *  Resize B
     */
    double volB = ei->SurfaceArea * ei->BLThickness_B;
    sum = 0.0;
    double curr_volB = 0.0;
    tp = VolPhaseList[solnBPhase_];
    mv = tp->molarVolume();
    phaseM = phaseMoles_final_[solnBPhase_];
    kstart = getGlobalSpeciesIndex(solnBPhase_, 0);
    int numSpecB = tp->nSpecies();
    if (phaseM <= 0.0) {
       phaseM = volB / mv;
       for (int k = 0; k < numSpecB; k++) {
        spMoles_final_[kstart + k] = phaseM * spMf_final_[kstart + k];
	sum += spMoles_final_[kstart + k];
      }
    } else {
      curr_volB = mv * phaseM;
      double ratio = volB / curr_volB;
      for (int k = 0; k < numSpecB; k++) {
        spMoles_final_[kstart + k] *= ratio;
	sum += spMoles_final_[kstart + k];
      }
    }
    phaseMoles_final_[solnBPhase_] = sum;

  


    for (int i = 0; i < m_NumSurPhases; i++) {
      int iph = NumVolPhases_ + i;
      tp = SurPhaseList[i];
      SurfPhase *sp = dynamic_cast<SurfPhase *>(tp);
      double sd = sp->siteDensity();
      phaseM = phaseMoles_final_[iph];
      phaseMoles_final_[iph] = sd *  ei->SurfaceArea;

      kstart = getGlobalSpeciesIndex(iph, 0);
      int numSpec = tp->nSpecies();
      if (phaseM <= 0.0) {
	phaseM = sd *  ei->SurfaceArea;
	for (int k = 0; k < numSpec; k++) {
	  spMoles_final_[kstart + k] = phaseM * spMf_final_[kstart + k];
	}
      } else { 
	double curr_sur = phaseM / sd;
	 
	double ratio = ei->SurfaceArea / curr_sur;
	for (int k = 0; k < numSpecB; k++) {
	  spMoles_final_[kstart + k] *= ratio;
	}
      }
    }


    InterfacialMassTransfer::updateState();
    InterfacialMassTransfer::setInitStateFromFinal(true);
    setFinalFinalStateFromFinal_Oin();
    return 0;
  }
  //====================================================================================================================
  // ---------------------------------------------------------------------------------------------
  // ----------------------- SPECIFY AND OBTAIN PROBLEM PARAMETERS -------------------------------
  // ---------------------------------------------------------------------------------------------
  
  // ------------------------------ OBTAIN STATIC PROBLEM INFORMATION ----------------------------
  
  //====================================================================================================================
  Cantera::ReactingSurDomain * InterfacialMassTransfer::currOuterReactingSurface() {
    for (int isk = 0; isk < m_NumSurPhases; isk++) {
      if (ActiveKineticsSurf_[isk]) {
	return RSD_List_[isk];
      }
    }
    return 0;
  }
  //====================================================================================================================
  Cantera::ReactingSurDomain * InterfacialMassTransfer::reactingSurface(int iSurf) {
    return RSD_List_[iSurf];
  }
  //====================================================================================================================
  // Set the phase existence flag in the electrode kinetics object so that kinetics 
  // are calculated correctly
  /*
   *    Flags are set in the kinetics object to tell the kinetics object which phases
   *    have zero moles.  The zero mole indicator is taken from phaseMoles_final_[]. Therefore,
   *    the final state is queried.
   *    There is a special case. Phases that are in justBornPhase_[] vector are allowed to
   *    be set to exist even if their phase moles are zero.
   *
   * @param doInactive   Boolean indicating whether to do inactive surface reactive surface domains 
   * @param assumeStableSingleSpeciesPhases Assume that single phases are stable. This
   *                         allows their production rates to be calculated
   */
  void InterfacialMassTransfer::setPhaseExistenceForReactingSurfaces(bool doInactive,
						       bool assumeStableSingleSpeciesPhases)
  {
    for (int isk = 0; isk < m_NumSurPhases; isk++) {
      /*
       *  Loop over phases, figuring out which phases have zero moles.
       *  Volume phases exist if the initial or final mole numbers are greater than zero
       *  Surface phases exist if the initial or final surface areas are greater than zero. 
       */
      if (ActiveKineticsSurf_[isk] || doInactive) {
	ReactingSurDomain *rsd = RSD_List_[isk];
	int nph = rsd->nPhases();
	for (int jph = 0; jph < nph; jph++) {
	  int iph = rsd->kinOrder[jph];
	  double mm = phaseMoles_init_[iph];
	  double mmf = phaseMoles_final_[iph];
	  ThermoPhase &tp = thermo(iph);
          int nsp = tp.nSpecies();
	  if (iph >=  NumVolPhases_) {
	    // we are in a surface phase
	    int isur = iph -  NumVolPhases_;
	    double sa_init = surfaceAreaRS_init_[isur];
	    double sa_final = surfaceAreaRS_final_[isur];
	    if (sa_init > 0.0 || sa_final > 0.0) {
	      rsd->setPhaseExistence(jph, true);
	    } else {
	      rsd->setPhaseExistence(jph, false);
	    }
	  } else {
	    if (mm <= 0.0 && mmf <= 0.0) {
	      rsd->setPhaseExistence(jph, false);
	      if (nsp == 1 && assumeStableSingleSpeciesPhases) {
		rsd->setPhaseStability(jph, true);
	      }
	    } else {
	      rsd->setPhaseExistence(jph, true);
	    }
	  }

	}
      }
    }
  }
  //================================================================================================
  // Returns the index of a phase in the ReactionSurfaceDomain object
  // given the index of that phase in the PhaseList object
  /*
   * @param PLph index of the phase in the PhaseList object, which is also the
   *             InterfacialMassTransfer_Model object.
   *
   *  @return  Returns the index of the phase in the current ReactingSurDomain
   *           object. A value of -1 in this slot means that the phase doesn't
   *           participate in the  current ReactingSurDomain object
   */
  int InterfacialMassTransfer::ReactingSurfacePhaseIndex(int isk, int PLph) const {
    ReactingSurDomain *rsd = RSD_List_[isk];
    if (!rsd) {
      throw CanteraError("ReactingSurfacePhaseIndex", "Reacting surface not found");
    }
    return rsd->PLtoKinPhaseIndex_[PLph];
  }
  //====================================================================================================================
  //   Reactant stoichiometric coefficient
  /*
   * Get the reactant stoichiometric coefficient for the kth global species
   * in the ith reaction of the reacting surface domain with index isk.
   */
  double InterfacialMassTransfer::reactantStoichCoeff(const int isk, int kGlobal, int i) {
    ReactingSurDomain *rsd= RSD_List_[isk];
    int krsd = rsd->PLtoKinSpeciesIndex_[kGlobal];
    if (krsd == -1) {
      return 0.0;
    }
    double rst = rsd->reactantStoichCoeff(krsd, i);
    return rst;
  }
  //====================================================================================================================
  //   Reactant stoichiometric coefficient
  /* 
   * Get the reactant stoichiometric coefficient for the kth global species
   * in the ith reaction of the reacting surface domain with index isk.
   */
  double InterfacialMassTransfer::productStoichCoeff(const int isk, int kGlobal, int i) {
    ReactingSurDomain *rsd= RSD_List_[isk];
    int krsd = rsd->PLtoKinSpeciesIndex_[kGlobal];
    if (krsd == -1) {
      return 0.0;
    }
    double rst = rsd->productStoichCoeff(krsd, i);
    return rst;
  }
  //====================================================================================================================
  // Specify the external fields are discretized with respect to the time coordinate
  /*
   *       0   Behavior within the global step is akin to backwards Euler. A step jump is 
   *           assumed to the global values at the end of the global time step even for
   *           intermediate times
   *       1   Behaviow within the global step is treated as a linear function between the 
   *           beginning values and the end values. 
   *
   *  @param externFieldTimeBehaviorType  Parameter describing the behavior.
   */
  void InterfacialMassTransfer::
  specifyExternalFieldTimeBehavior(EF_FieldTimeBehavior_Enum  externFieldTimeBehaviorType)
  {
    externFieldTimeBehaviorType_ = externFieldTimeBehaviorType;
  }
  //====================================================================================================================
  // Report  how the external fields are discretized with respect to the time coordinate
  /*
   *       0   Behavior within the global step is akin to backwards Euler. A step jump is 
   *           assumed to the global values at the end of the global time step even for
   *           intermediate times
   *       1   Behaviow within the global step is treated as a linear function between the 
   *           beginning values and the end values. 
   *
   *  @return  Returns a parameter describing the behavior.
   */
  EF_FieldTimeBehavior_Enum InterfacialMassTransfer::reportExternalFieldTimeBehavior() const
  {
    return externFieldTimeBehaviorType_;
  }
  //=====================================================================================================================

  // ------------------------------ SPECIFY BASIC THERMO CONDITIONS  ------------------------------------

  //====================================================================================================================

  //================================================================================================
  void  InterfacialMassTransfer::setState_TP(doublereal temperature, doublereal presA, doublereal presB)
  {
    if (presB == -1.0) {
      presB = presA;
    }
    Temp_ = temperature;
    Pres_A_Interface_final_ = presA;
    Pres_B_Interface_final_ = presB;
    updateState();
  }
  //================================================================================================
  double InterfacialMassTransfer::temperature() const {
    return Temp_;
  }
  //================================================================================================
  double InterfacialMassTransfer::pressure() const {
    return Pres_A_Interface_final_;
  }
  //================================================================================================  

  // Set the global final_final time
  /*!
   *  When we do this we are setting a pending state flag
   *
   *  @param t_final_final  Final time of the global step
   *  @param setFinal       set the final state as well as the final_final
   *                           Defaults to true.
   */
  void InterfacialMassTransfer::setStateFF_time(double t_final_final, bool setFinal)
  {
    pendingIntegratedStep_ = true;
    t_final_final_ = t_final_final;

    if (setFinal) {
      t_final_ = t_final_final;
    }
  }
  //================================================================================================  

  void InterfacialMassTransfer::setTime(double time) {
    if (pendingIntegratedStep_) {
	throw CanteraError("InterfacialMassTransfer::setTime",
			   "called when there is a pending step");
    }
    t_init_init_ = time;
    t_final_ = t_init_init_;
    t_init_  = t_init_init_;
    t_final_final_ = t_init_init_;
  }
  //================================================================================================

  // ------------------------------ SPECIFY PROBLEM PARAMETERS ------------------------------------

  //====================================================================================================================
  // Set the mole numbers in a single phase
  /*
   *  We set the mole numbers of a single phase separately from 
   *  the rest of the phases.
   *
   *  We always make sure that mole numbers are positive by clipping. We
   *  always make sure that mole fractions sum to one.
   *
   *  If we are not following mole numbers in the electrode, we set the
   *  total moles to the internal constant, electrolytePseudoMoles_, while
   *  using this vector to set the mole fractions.
   *
   * @param iph     Phase id.
   * @param moleNum vector of mole numbers of the species in the
   *                     electrolyte phase. 
   *                     units = kmol
   *                     size = number of species in the electrolyte phase 
   */
  void InterfacialMassTransfer::setPhaseMoleNumbers(int iph, const double * const moleNum) {
    if (pendingIntegratedStep_) {
      if (iph != solnAPhase_) {
	throw CanteraError("InterfacialMassTransfer::setPhaseMoleNumbers",
			   "called when there is a pending step");
      }
    }
    int istart = m_PhaseSpeciesStartIndex[iph];
    int nsp =  m_PhaseSpeciesStartIndex[iph+1] - istart;
    for (int k = 0; k < nsp; k++) {
      spMoles_final_[istart + k] = MAX(moleNum[k], 0.0);
      spMoles_init_[istart + k] = spMoles_final_[istart + k];
      spMoles_init_init_[istart + k] = spMoles_final_[istart + k];
    }
    updatePhaseNumbers(iph);

    mdpUtil::mdp_copy_dbl_1(DATA_PTR(spMoles_init_init_), DATA_PTR(spMoles_final_), m_NumTotSpecies);
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(spMoles_init_), DATA_PTR(spMoles_final_), m_NumTotSpecies);
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(spMoles_final_final_), DATA_PTR(spMoles_final_), m_NumTotSpecies);
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(phaseMoles_init_), DATA_PTR(phaseMoles_final_), m_NumTotPhases);
    mdpUtil::mdp_copy_dbl_1(DATA_PTR(phaseMoles_init_init_), DATA_PTR(phaseMoles_final_), m_NumTotPhases);

  
  }

  //================================================================================================  

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
  void InterfacialMassTransfer::setSolnA_BoundaryConcentrations(const double * const solnAMoleNum, 
								double tempA, double presA,
								const double * extFieldAValues,
								bool setInitial)
  {

    int istart = m_PhaseSpeciesStartIndex[solnAPhase_];
    ThermoPhase &tp = thermo(solnAPhase_);
    int nsp = tp.nSpecies();
    AssertTrace(nsp == ( m_PhaseSpeciesStartIndex[solnAPhase_+1] - m_PhaseSpeciesStartIndex[solnAPhase_]));
    double tmp = 0.0;
    double *spMoles = &spMoles_final_[0];
    double *spMf = &spMf_final_[0];
    double *phaseMoles = &phaseMoles_final_[0];
    if (setInitial) {
      spMoles = &spMoles_init_[0];
      spMf = &spMf_init_[0];
      phaseMoles = &phaseMoles_init_[0];
    }

    for (int k = 0; k < nsp; k++) {
      spMoles[istart + k] = MAX(solnAMoleNum[k], 0.0);
      tmp += spMoles[istart + k];
    }
    phaseMoles[solnAPhase_] = tmp;
 
    if (tmp > 1.0E-200) {
      for (int k = 0; k < nsp; k++) {
	spMf[istart + k] = spMoles[istart + k] / tmp;
      }
    } 
    updateState();
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
  void InterfacialMassTransfer::setSolnB_BoundaryConcentrations(const double * const solnBMoleNum, 
								double tempB, double presB,
								const double * extFieldBValues,
								bool setInitial)
  {

    int istart = m_PhaseSpeciesStartIndex[solnBPhase_];
    ThermoPhase &tp = thermo(solnBPhase_);
    int nsp = tp.nSpecies();
    AssertTrace(nsp == ( m_PhaseSpeciesStartIndex[solnBPhase_+1] - m_PhaseSpeciesStartIndex[solnBPhase_]));
    double tmp = 0.0;
    double *spMoles = &spMoles_final_[0];
    double *spMf = &spMf_final_[0];
    double *phaseMoles = &phaseMoles_final_[0];
    if (setInitial) {
      spMoles = &spMoles_init_[0];
      spMf = &spMf_init_[0];
      phaseMoles = &phaseMoles_init_[0];
    }

    for (int k = 0; k < nsp; k++) {
      spMoles[istart + k] = MAX(solnBMoleNum[k], 0.0);
      tmp += spMoles[istart + k];
    }
    phaseMoles[solnBPhase_] = tmp;
 
    if (tmp > 1.0E-200) {
      for (int k = 0; k < nsp; k++) {
	spMf[istart + k] = spMoles[istart + k] / tmp;
      }
    }
    updateState();
  }
  //================================================================================================
  //! Set the surface areas of surfaces within the model at the final time of the global
  //! step.
  /*!
   *   It is up to the object to determine what to do with this information, if anything.
   *
   *  @param surfaceAreasRS   Vector of surface areas (length equal to number of surfaces)
   *
   *  @param setInitial   Boolean indicating that we should set the initial values
   *                      instead of the final values (default false).
   */
  void  InterfacialMassTransfer::setSurfaceAreas(const double * const surfaceAreaRS, bool setInitial)
  {
    double *surfArea = &surfaceAreaRS_final_[0];
    if (setInitial) {
      surfArea = &surfaceAreaRS_init_[0];
    }
    for (int i = 0; i < numSurfaces_; i++) {
      surfArea[i] =  surfaceAreaRS[i];
    }
  }
  //================================================================================================
  // Update all mole numbers in the object from the mole numbers in the spMoles_final_[] vector
  /*
   *  We use the field spMoles_final_[] to set the field phaseMoles_final_[].
   *
   *  We set the mole numbers of a single phase separately from the rest of the phases.
   *
   *  We do not clip the mole numbers to be positive. We allow negative mole numbers.
   *
   *  We make sure that the mole fractions sum to one.
   *
   *   The following fields in this object are set:
   *
   *            spMf_final_[]
   *            VolPM_[]
   *            spElectroChemPot_[]
   *
   *            phaseMoles_final_[iph] 
   *            phaseVoltages_[iph]
   *           	phaseMolarVolumes_[iph]
   *
   *  If we are not following  the mole numbers in the electrode, we set the
   *  total moles to the internal constant, electrolytePseudoMoles_, while
   *  using this vector to set the mole fractions, using the ThermoPhase object
   *  to get the mole fractions.
   *
   * @param iph     Phase id.
   */
  void InterfacialMassTransfer::updatePhaseNumbers(int iph) {
    int istart = m_PhaseSpeciesStartIndex[iph];
    ThermoPhase &tp = thermo(iph);
    int nsp =  m_PhaseSpeciesStartIndex[iph+1] - istart;
    double tmp = 0.0;
    for (int k = 0; k < nsp; k++) {
      tmp += spMoles_final_[istart + k];
    }
    phaseMoles_final_[iph] = tmp;
    if (tmp > 1.0E-200) {
      for (int k = 0; k < nsp; k++) {
	spMf_final_[istart + k] = spMoles_final_[istart + k] / tmp;
      }
      if (iph == solnAPhase_) {

      }
      // Here we set the state within the phase object
      tp.setState_TPX(Temp_, Pres_A_Interface_final_, &spMf_final_[istart]);

  

    } else {
      // We are here when the mole numbers of the phase are zero. In this case, we still need
      // a valid mole fraction vector. This is kept in the ThermoPhase object.
      ThermoPhase &tp = thermo(iph);
      tp.setState_TP(Temp_, Pres_A_Interface_final_);
      tp.getMoleFractions(&(spMf_final_[istart]));
    }
    tp.getPartialMolarVolumes(&(VolPM_[istart]));
    if (iph >= NumVolPhases_) {
      nsp = tp.nSpecies();
      for (int k = 0; k < nsp; k++) {
	VolPM_[istart + k] = 0.0;
      }
    }
    tp.getChemPotentials(&(spChemPot_[istart]));
    
    if (iph < NumVolPhases_) {
      phaseMolarVolumes_[iph] = tp.molarVolume();
    } else {
      phaseMolarVolumes_[iph] = 0.0;
    }
  } 

  //====================================================================================================================
  // Calculates the change in the surface area of all external and internal interfaces within the electrode
  /*
   *  (virtual)
   *  variables to be potentially altered
   *   surfaceAreaRS_[];
   *   isExternalSurface[]
   *   numExternalInterfacialSurfaces_;
   */
  double InterfacialMassTransfer::calcSurfaceAreaChange(double deltaT) {
    double sa_final = surfaceAreaRS_init_[0];
    return sa_final;
  }

  //====================================================================================================================
  doublereal InterfacialMassTransfer::phaseMoles(int iph) const {
    return phaseMoles_final_[iph];
  } 
  
  //================================================================================================
  doublereal InterfacialMassTransfer::elementMoles(int ie) const {
    return 0.0;
  } 
  //================================================================================================
  // Get the total number of moles in the system
  /*
   * @return returns the total number of moles in the system
   */
  double InterfacialMassTransfer::totalMoles() const
  {
    double sum = 0.0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
      sum+=  phaseMoles_final_[iph];
    }
    return sum;
  }
 //================================================================================================
  double InterfacialMassTransfer::speciesChemPotential(int iGlobalSpIndex) const {
    // DANGER -> Check to see this is updated

    for (int iph = 0; iph < m_NumTotPhases; iph++) {
      ThermoPhase *tphase = &(thermo(iph));
      int kStart = m_PhaseSpeciesStartIndex[iph];
      tphase->getChemPotentials(& (spChemPot_[kStart]));
    }
    double tmp = spChemPot_[iGlobalSpIndex];
  
    return tmp;
  }
  //================================================================================================
  void  InterfacialMassTransfer::getMoleFractions(doublereal* const x) const {
    mdpUtil::mdp_copy_dbl_1(x, &(spMf_final_[0]), m_NumTotSpecies);
  }
  //================================================================================================
  void  InterfacialMassTransfer::getMoleNumSpecies(doublereal* const n) const {
    mdpUtil::mdp_copy_dbl_1(n, &(spMoles_final_[0]), m_NumTotSpecies);
  }
  //================================================================================================
  void  InterfacialMassTransfer::getMoleNumPhases(doublereal* const np) const {
    mdpUtil::mdp_copy_dbl_1(np, &(phaseMoles_final_[0]), m_NumTotPhases);
  }
  //================================================================================================
  double  InterfacialMassTransfer::moleFraction(int globalSpeciesIndex) const {
    return spMf_final_[globalSpeciesIndex];
  }
  //================================================================================================
  double  InterfacialMassTransfer::moleNumSpecies(int globalSpeciesIndex) const {
    return spMoles_final_[globalSpeciesIndex];
  }


  //====================================================================================================================
  double InterfacialMassTransfer::SolidVol() const {
    double vol = 0.0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
      int kStart = m_PhaseSpeciesStartIndex[iph];
      ThermoPhase &tp = thermo(iph);
      int nspPhase = tp.nSpecies();
      if (iph != solnAPhase_) {
	for (int k = 0; k < nspPhase; k++) {
	  vol += spMoles_final_[kStart + k] * VolPM_[kStart + k];
	}
      }
    }
    return vol;
  }
  //====================================================================================================================
  double InterfacialMassTransfer::TotalVol() const {
    double vol = 0.0;
   
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
      double psum = 0.0;
      int kStart = m_PhaseSpeciesStartIndex[iph];
      ThermoPhase &tp = thermo(iph);
      int nspPhase = tp.nSpecies();
      for (int k = 0; k < nspPhase; k++) {
	vol += spMoles_final_[kStart + k] * VolPM_[kStart + k];
	psum += spMoles_final_[kStart + k] * VolPM_[kStart + k];
      }
      double mv = tp.molarVolume();
      double palt = mv * phaseMoles_final_[iph];
      if (palt < -1.0E-15) {
	throw CanteraError(" InterfacialMassTransfer::TotalVol() ",
			   " phase volume is negative " + fp2str(palt));
      }
      if (psum < -1.0E-15) {
	throw CanteraError(" InterfacialMassTransfer::TotalVol() ",
			   " phase volume is negative " + fp2str(psum));
      }
      double denom = palt + psum + 1.0E-9;
      if (tp.eosType() != cLattice) {
	if (fabs((palt - psum) / denom) > 1.0E-4) {
	  throw CanteraError(" InterfacialMassTransfer::TotalVol() ",
			     " internal inconsistency " + fp2str(palt) + " " + fp2str(psum));
	}
      }
    }
    return vol;
  } 
  //====================================================================================================================
  void InterfacialMassTransfer::getPhaseVol(double * const phaseVols) const {
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
      phaseVols[iph] = 0.0;
      int kStart = m_PhaseSpeciesStartIndex[iph];
      ThermoPhase &tp = thermo(iph);
      int nspPhase = tp.nSpecies();
      for (int k = 0; k < nspPhase; k++) {
	phaseVols[iph]  += spMoles_final_[kStart + k] * VolPM_[kStart + k];
      }
    }
  }
  //====================================================================================================================
  void InterfacialMassTransfer::getSurfaceAreas(double *const surfArea) const {
    int sz = m_NumSurPhases;
    if ((int) surfaceAreaRS_final_.size() > m_NumSurPhases) {
      sz = (int) surfaceAreaRS_final_.size();
    }
    for (int i = 0; i < sz; i++) {
      surfArea[i] = surfaceAreaRS_final_[i];
    }
  }
  //====================================================================================================================
  void InterfacialMassTransfer::setSurfaceAreas(const double *const surfArea) {
    int sz = m_NumSurPhases;
    if ((int) surfaceAreaRS_final_.size() > m_NumSurPhases) {
      sz = surfaceAreaRS_final_.size();
    }
    for (int i = 0; i < sz; i++) {
      surfaceAreaRS_final_[i] = surfArea[i];
      surfaceAreaRS_init_[i] = surfArea[i];
      if (surfArea[i] <= 0.0) {
	ActiveKineticsSurf_[i] = 0;
      } else {
	if (RSD_List_[i]) {
	  ActiveKineticsSurf_[i] = 1;
	}
      }
    }
  }
  //====================================================================================================================
  //  Take the state (final) within the InterfacialMassTransfer_Model and push it down
  //  to the ThermoPhase Objects 
  /*
   *  We take the values of spMoles_final_[] and propagate them down to the ThermoPhase
   *  objects in the electrode
   */
  void InterfacialMassTransfer::updateState() {
    int iph;

    /*
     * This may be redundant. However, I want to make sure mole fractions are
     * consistent with final moles.
     */
    for (int i = 0; i < m_NumTotPhases; i++) {
      updatePhaseNumbers(i);
    }

    /*
     * Loop over all phases in the object
     */
    for (iph = 0; iph < m_NumTotPhases; iph++) {
      ThermoPhase *tphase = &(thermo(iph));
      std::string pName = tphase->id();
      int kStart = m_PhaseSpeciesStartIndex[iph];
      double pres =  Pres_A_Interface_final_;
      if (iph == solnBPhase_) {
	pres =  Pres_B_Interface_final_;
      }
      /*
       * A direct call to tphase to set the state is unnecessary,
       * because m_mp will do it anyway.
       */
      tphase->setState_TPX(Temp_, pres, &spMf_final_[kStart]);
      /*
       * Ok, we have set the state. Now upload the Vol and mu's.
       */
      tphase->getPartialMolarVolumes(& (VolPM_[kStart]));
      tphase->getChemPotentials(& (spChemPot_[kStart]));

      if (iph < NumVolPhases_) {
	phaseMolarVolumes_[iph] = tphase->molarVolume();
      }
    }
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
  void InterfacialMassTransfer::updateVelocities()
  {
    double massfluxA = getPhaseAMassSourceTerm();
    ThermoPhase &tpA = thermo(solnAPhase_);
    double densA = tpA.density();
    double area = surfaceAreaRS_final_[0];
    double sVelocA = massfluxA  / area / densA;

    double massfluxB = getPhaseBMassSourceTerm();
     ThermoPhase &tpB = thermo(solnBPhase_);
    double densB = tpB.density();
    double sVelocB = massfluxB  / area / densB;

    Velo_A_Interface_final_ = Velo_ReferenceFrame_;
    Velo_S_Interface_final_ = sVelocA + Velo_A_Interface_final_;
    Velo_B_Interface_final_ = sVelocB + Velo_S_Interface_final_;
  }
  
  //------------------------------------------------------------------------------------
  // ----------------------------- GET CONDITIONS OUT ---------------------------------------------
  // ----------------------------------------------------------------------------------------------
  //
  //       (unless specified this is always at the final conditions and time
  //
  //
  // ----------------------------- GET INSTANTANEOUS SOURCE TERMS --------------------------------
  //====================================================================================================================
  // Get the net production rates of all species in the electrode object
  // at the current conditions
  /*
   *
   *  This routine assumes that the underlying objects have been updated.
   *  It uses the default Backwards Euler integration rule, which doesn't take into account of issues with surfaces going away
   */
  void InterfacialMassTransfer::getNetProductionRates(doublereal* const net) const
  {
    mdpUtil::mdp_zero_dbl_1(net, m_NumTotSpecies);
    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the InterfacialMassTransfer.
     */
    /*
     *  For each Reacting surface
     */ 
    for (int isk = 0; isk < numSurfaces_; isk++) {
      if (ActiveKineticsSurf_[isk]) {
	/*
	 *  Just assume surface area is equal to final value
	 */
	double area = surfaceAreaRS_final_[isk];
	/*
	 *  Get the species production rates for the reacting surface
	 */
	const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
	/*
	 * Number of phases in the reacting species.
	 */
	int nphRS = RSD_List_[isk]->nPhases();

	int kIndexKin = 0;
	/*
	 *  Loop over the phases defined in the Reacting Surface Domain object
	 */
	for (int kph = 0; kph < nphRS; kph++) {
	  /*
	   *  Find the phase Id of the current phase in the InterfacialMassTransfer object
	   */
	  int jph = RSD_List_[isk]->kinOrder[kph];

	  /*
	   *  Find the starting species index within the InterfacialMassTransfer object
	   */
	  int istart = m_PhaseSpeciesStartIndex[jph];
	  int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	  for (int k = 0; k < nsp; k++) {
	    net[istart + k] += rsSpeciesProductionRates[kIndexKin] * area;
	    kIndexKin++;
	  }
	}
      }
    }
  } 
 //====================================================================================================================
  // Get the net production rates of all species in the electrode object
  // at the current conditions
  /*
   *
   *  This routine assumes that the underlying objects have been updated.
   *  It uses the default Backwards Euler integration rule, which doesn't take into account of issues with surfaces going away
   */
  void InterfacialMassTransfer::getNetProductionRatesRSD(const int isk, doublereal* const net) const
  {
    mdpUtil::mdp_zero_dbl_1(net, m_NumTotSpecies);
    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the InterfacialMassTransfer.
     */
    /*
     *  For each Reacting surface
     */
    if (ActiveKineticsSurf_[isk]) {
      /*
       *  Just assume surface area is equal to final value
       */
      double area = surfaceAreaRS_final_[isk];
      /*
       *  Get the species production rates for the reacting surface
       */
      const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
      /*
       * Number of phases in the reacting species.
       */
      int nphRS = RSD_List_[isk]->nPhases();

      int kIndexKin = 0;
      /*
       *  Loop over the phases defined in the Reacting Surface Domain object
       */
      for (int kph = 0; kph < nphRS; kph++) {
	/*
	 *  Find the phase Id of the current phase in the InterfacialMassTransfer object
	 */
	int jph = RSD_List_[isk]->kinOrder[kph];

	/*
	 *  Find the starting species index within the InterfacialMassTransfer object
	 */
	int istart = m_PhaseSpeciesStartIndex[jph];
	int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	for (int k = 0; k < nsp; k++) {
	  net[istart + k] += rsSpeciesProductionRates[kIndexKin] * area;
	  kIndexKin++;
	}
      }
    }
  }
  //====================================================================================================================
  //  Returns the current and the net production rates of the phases in kg/s from a single surface
  /*
   *  Returns the net production rates of all phases from reactions on a single surface
   *  
   *  @param isk Surface ID to get the fluxes from.      
   *  @param phaseMassFlux  Returns the mass fluxes of the phases
   */
  void InterfacialMassTransfer::getPhaseMassFlux(doublereal* const phaseMassFlux) const
  {
    mdpUtil::mdp_zero_dbl_1(phaseMassFlux, m_NumTotPhases);

    for (int isk = 0; isk < numSurfaces_; isk++) {
      if (ActiveKineticsSurf_[isk]) { 
	const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
	double mult =  surfaceAreaRS_final_[isk];
	int nphRS = RSD_List_[isk]->nPhases();
	int kIndexKin = 0;
	for (int kph = 0; kph < nphRS; kph++) {
	  int jph = RSD_List_[isk]->kinOrder[kph];
	  ThermoPhase &tp = thermo(jph);
	  int istart = m_PhaseSpeciesStartIndex[jph];
	  int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	  for (int k = 0; k < nsp; k++) {
	    double net = rsSpeciesProductionRates[kIndexKin] * mult;
	    double mw = tp.molecularWeight(k);
	    phaseMassFlux[jph] += net * mw;
	    kIndexKin++;
	  }
	}
      }
    }
  }
  //================================================================================================
  //  Returns the current and the net production rates of the phases in kmol/s from a single surface
  /*
   *  Returns the net production rates of all phases from reactions on a single surface
   *  
   *  @param isk Surface ID to get the fluxes from.      
   *  @param phaseMassFlux  Returns the mass fluxes of the phases
   */
  void InterfacialMassTransfer::getPhaseMoleFlux(const int isk, doublereal* const phaseMoleFlux) const
  { 
    mdpUtil::mdp_zero_dbl_1(phaseMoleFlux, m_NumTotPhases);
    for (int isk = 0; isk < numSurfaces_; isk++) {
    if (ActiveKineticsSurf_[isk]) {
      const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
      double area =  (surfaceAreaRS_final_[isk]);
      int nphRS = RSD_List_[isk]->nPhases();
      int PLph, RSph;
      int kIndexKin = 0;
      for (RSph = 0; RSph < nphRS; RSph++) {
	PLph = RSD_List_[isk]->kinOrder[RSph];
	int istart = m_PhaseSpeciesStartIndex[PLph];
	int nsp = m_PhaseSpeciesStartIndex[PLph+1] - istart;
	for (int k = 0; k < nsp; k++) {
	  double net = rsSpeciesProductionRates[kIndexKin]* area;
	  phaseMoleFlux[PLph] += net;
	  kIndexKin++;
	}
      }
    }
    }
  }
  //====================================================================================================================
  //!  Returns the computed mass flux of species into the A phase at final conditions
  /*!
   *  This is also the creation rate of phase A at the final conditions.
   *
   *      @return Returns mass flux  of species into phase A (kg s-1)
   */
  double InterfacialMassTransfer::getPhaseAMassSourceTerm() const
  {
    double sum =0.0;
    int kIndex = m_PhaseSpeciesStartIndex[solnAPhase_];
    ThermoPhase &tpA = thermo(solnAPhase_);
    int nsp = m_PhaseSpeciesStartIndex[solnAPhase_+1] - kIndex;
    for (int k = 0; k < nsp; k++) {
      double mw = tpA.molecularWeight(k);
      sum += DspMoles_final_[kIndex] * mw;
      kIndex++;
    }    
    return sum;
  }
  //====================================================================================================================
  //!  Returns the computed mass flux of species into the Bphase at final conditions
  /*!
   *  This is also the creation rate of phase B at the final conditions.
   *
   *      @return Returns mass flux  of species into phase B (kg s-1)
   */
  double InterfacialMassTransfer::getPhaseBMassSourceTerm() const
  {
    double sum =0.0;
    int kIndex = m_PhaseSpeciesStartIndex[solnBPhase_];
    ThermoPhase &tpB = thermo(solnBPhase_);
    int nsp = m_PhaseSpeciesStartIndex[solnBPhase_+1] - kIndex;
    for (int k = 0; k < nsp; k++) {
      double mw = tpB.molecularWeight(k);
      sum += DspMoles_final_[kIndex] * mw;
      kIndex++;
    }    
    return sum;
  }
  //====================================================================================================================
  //! Returns the computed Stefan velocity for phase A in the object
  //! at the current final conditions. 
  /*!
   *  Multiply by the surface area to get the stefan volume production rate at the
   *  current final conditions.
   * 
   *    @return returns stefan velocity created in phase A (m s-1)
   */
  double InterfacialMassTransfer::StefanVelocityPhaseA() const
  {
    double massflux = getPhaseAMassSourceTerm();
    ThermoPhase &tpA = thermo(solnAPhase_);
    double dens = tpA.density();
    double area = surfaceAreaRS_final_[0];
    return massflux / area / dens;
  }
  //====================================================================================================================
  //! Returns the computed Stefan velocity for phase Bin the object
  //! at the current final conditions.
  /*!
   *  Multiply by the surface area to get the stefan volume production rate at the
   *  current final conditions.
   * 
   *    @return returns stefan velocity created in phase B (m s-1)
   */
  double InterfacialMassTransfer::StefanVelocityPhaseB() const
  {
    double massflux = getPhaseBMassSourceTerm();
    ThermoPhase &tpB = thermo(solnBPhase_);
    double dens = tpB.density();
    double area = surfaceAreaRS_final_[0];
    return massflux / area / dens;
  }
  //====================================================================================================================
  void InterfacialMassTransfer::getIntegratedProductionRates(doublereal* const net) const
  {
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
      invDelT = 1.0/ (t_final_final_ - t_init_init_);
    }
    for (int k = 0; k < m_NumTotSpecies; k++) {
      net[k] = invDelT * spMoleIntegratedSourceTerm_[k];
    }
  }
  //================================================================================================
  //   Returns the integrated moles created for each phase in the object
  //  over the current global time step
  /*
   *    @param  phaseMolesTransfered vector of moles transfered (length = number of total 
   *            phases in the object)
   *            units = kmol
   */
  void InterfacialMassTransfer::getIntegratedPhaseMoleSourceTerm(doublereal* const phaseMolesCreated) const
  { 
    if (!pendingIntegratedStep_) {
      throw CanteraError(" InterfacialMassTransfer::getIntegratedPhaseMoleSourceTerm",
			   "no pending integration step");
    }
    double sum = 0.0;
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
      phaseMolesCreated[iph] = 0.0;
      ThermoPhase &tp = thermo(iph);
      std::string pname = tp.id();
      int istart = m_PhaseSpeciesStartIndex[iph]; 
      int nsp = tp.nSpecies();
      for (int ik = 0; ik < nsp; ik++) {
	int k = istart + ik;
	phaseMolesCreated[iph] += spMoleIntegratedSourceTerm_[k];
      }
      sum += fabs(phaseMolesCreated[iph]);
    }
  } 
  //================================================================================================
  //    Returns the integrated mass created for a particular phase in the object
  //    over the current global time step
  /*
   *    @param  phaseMolesTransfered vector of moles transfered (length = number of total 
   *            phases in the object)
   *            units = kmol
   */
 double InterfacialMassTransfer::getIntegratedPhaseMassSourceTerm(int iph) const
  {
    if (!pendingIntegratedStep_) {
      throw CanteraError(" InterfacialMassTransfer::getIntegratedPhaseMoleSourceTerm",
			 "no pending integration step");
    }
    double sum = 0.0;
    ThermoPhase &tp = thermo(iph);
    std::string pname = tp.id();
    int istart = m_PhaseSpeciesStartIndex[iph]; 
    int nsp = tp.nSpecies();
    for (int ik = 0; ik < nsp; ik++) {
      int k = istart + ik;
      sum += spMoleIntegratedSourceTerm_[k] * tp.molecularWeight(ik);
    }
    return sum;
  } 
  //================================================================================================
  // Returns the integrated mass created for phase A in the object
  // over the current global time step
  /*
   *    @return returns mass created in phase A
   */
  double InterfacialMassTransfer::getIntegratedPhaseAMassSourceTerm() const
  {
    return getIntegratedPhaseMassSourceTerm(solnAPhase_);
  }
  //================================================================================================
  // Returns the integrated mass created for phase B in the object
  // over the current global time step
  /*
   *    @return returns mass created in phase B
   */
  double InterfacialMassTransfer::getIntegratedPhaseBMassSourceTerm() const
  {
    return getIntegratedPhaseMassSourceTerm(solnBPhase_);
  }
  //===================================================================================================================
  double  InterfacialMassTransfer::integratedStefanVolumeRatePhaseA() const
  {
    double massflux = getIntegratedPhaseAMassSourceTerm();
    ThermoPhase &tp = thermo(solnAPhase_);
    double dens = tp.density();
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
      invDelT = 1.0/ (t_final_final_ - t_init_init_);
    }
    double area = 1.0;

    double sVelocA = massflux * invDelT / area / dens;
    return sVelocA;
  }
  //===================================================================================================================
  double  InterfacialMassTransfer::integratedStefanVolumeRatePhaseB() const
  {
    double massflux = getIntegratedPhaseBMassSourceTerm();
    ThermoPhase &tp = thermo(solnBPhase_);
    double dens = tp.density();
    double invDelT = 1.0;
    if (t_final_final_ > t_init_init_) {
      invDelT = 1.0/ (t_final_final_ - t_init_init_);
    }
    // double area = 0.5 * (surfaceAreaRS_final_final_[0] + surfaceAreaRS_init_init_[0]);

    double sVelocB = massflux * invDelT / dens;
    return sVelocB;
  }
  //===================================================================================================================
  //  Residual calculation for the solution of the Nonlinear integration problem
  /*
   *  
   * @param t             Time                    (input) 
   * @param delta_t       The current value of the time step (input)
   * @param y             Solution vector (input, do not modify)
   * @param ydot          Rate of change of solution vector. (input, do not modify)
   * @param resid         Value of the residual that is computed (output)
   * @param evalType      Type of the residual being computed (defaults to Base_ResidEval)
   * @param id_x          Index of the variable that is being numerically differenced to find
   *                      the jacobian (defaults to -1, which indicates that no variable is being
   *                      differenced or that the residual doesn't take this issue into account)
   * @param delta_x       Value of the delta used in the numerical differencing
   */
  int InterfacialMassTransfer::integrateResid(const doublereal tfinal, const doublereal deltaTsubcycle,
 		                const doublereal * const y, const doublereal * const ydot,
		                doublereal * const resid,
		                const ResidEval_Type_Enum evalType, const int id_x, const doublereal delta_x)
  { 
  
    //double tinit = tfinal - deltaTsubcycle;
    //bool newStep= false;
    //if (fabs(tinit - t_init_) > 1.0E-14) {
     // newStep = true;
  
    //}


    for (int k = 0; k < m_NumTotSpecies; k++) {
      spMoles_final_[k] = y[k];
    }

    vector<doublereal> phaseMoles_tmp(m_NumTotPhases, 0.0);
    vector<doublereal> spMf_tmp(m_NumTotSpecies, 0.0);
    vector<doublereal> spMoles_tmp(m_NumTotSpecies, 0.0);
    vector<doublereal> Xf_tmp(m_NumTotSpecies, 0.0);

    vector<doublereal> srcTerm(m_NumTotSpecies, 0.0);
    static  vector<int> justBornMultiSpecies;
    
 
    updateState();


    t_final_ = t_init_;
    t_final_ += deltaTsubcycle;
  
    /*
     * Advance the initial state so that it matches the final state from the last time iteration
     */
    spMoles_init_ =  spMoles_final_;
    phaseMoles_init_ =  phaseMoles_final_;

   
    justBornMultiSpecies.clear();
     

    /*
     * We start a predictor corrector damping cycle here
     */

    int bornMultiSpecies = -1;



    /*
     *   Set the internal objects to the correct conditions
     *    -> This will be the final conditions.
     */
    updateState();

    /*
     * Loop over surface phases, filling in the phase existence fields within the 
     * kinetics operator
     */
    for (int isk = 0; isk < m_NumSurPhases; isk++) {
      /*
       *  Loop over phases, figuring out which phases have zero moles.
       *  Volume phases exist if the initial or final mole numbers are greater than zero
       *  Surface phases exist if the initial or final surface areas are greater than zero. 
       */
      if (ActiveKineticsSurf_[isk]) {
	ReactingSurDomain *rsd = RSD_List_[isk];
	int nph = rsd->nPhases();
	for (int jph = 0; jph < nph; jph++) {
	  int iph = rsd->kinOrder[jph];

	  double mm = phaseMoles_init_[iph];
	  double mmf = phaseMoles_final_[iph];
	  if (iph >=  NumVolPhases_) {
	    // we are in a surface phase
	    int isur = iph -  NumVolPhases_;
	    double sa_init = surfaceAreaRS_init_[isur];
	    double sa_final = surfaceAreaRS_final_[isur];
	    if (sa_init > 0.0 || sa_final > 0.0) {
	      rsd->setPhaseExistence(jph, true);
	    } else {
	      rsd->setPhaseExistence(jph, false);
	    }
	  } else {
	    if (mm <= 0.0 && mmf <= 0.0) {
	      rsd->setPhaseExistence(jph, false);
	    } else {
	      rsd->setPhaseExistence(jph, true);
	    }
	  }
	  if (iph == bornMultiSpecies) {
	    rsd->setPhaseExistence(jph, true);
	  }
	  for (int iiph = 0; iiph < (int) justBornMultiSpecies.size(); iiph++) {
	    if (iph == justBornMultiSpecies[iiph]) {
	      rsd->setPhaseExistence(jph, true);
	    }
	  }
	}
      }
    }

    /*
     *  This routine basically translates between species lists for the reacting surface
     *  domain and the InterfacialMassTransfer.
     *  Later, when we have more than one reacting surface domain in the electrode object,
     *  this will do a lot more
     */

    for (int isk = 0; isk < m_NumSurPhases; isk++) {
      // Loop over phases, figuring out which phases have zero moles.

      if (ActiveKineticsSurf_[isk]) {
       
	/*
	 *  For each Reacting surface
	 *
	 *  Get the species production rates for the reacting surface
	 */
	//    m_rSurDomain->getNetProductionRates(&RSSpeciesProductionRates_[0]);
	const vector<double> &rsSpeciesProductionRates = RSD_List_[isk]->calcNetProductionRates();
 
    
	double *spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
	/*
	 *  loop over the phases in the reacting surface
	 *  Get the net production vector
	 */ 
	mdpUtil::mdp_zero_dbl_1(spNetProdPerArea, m_NumTotSpecies);
	int nphRS = RSD_List_[isk]->nPhases();
	int jph, kph;
	int kIndexKin = 0;
	for (kph = 0; kph < nphRS; kph++) {
	  jph = RSD_List_[isk]->kinOrder[kph];
	  int istart = m_PhaseSpeciesStartIndex[jph];
	  int nsp = m_PhaseSpeciesStartIndex[jph+1] - istart;
	  for (int k = 0; k < nsp; k++) {
	    spNetProdPerArea[istart + k] += rsSpeciesProductionRates[kIndexKin];
	    if (rsSpeciesProductionRates[kIndexKin] > 0.0) {
	      if ((phaseMoles_init_[jph] <= 0.0)) {
		bool notFound = true;
		for (int iiph = 0; iiph < (int)justBornMultiSpecies.size(); iiph++) {
		  if (jph == justBornMultiSpecies[iiph]) {
		    notFound = false;
		  }
		}
		if (notFound) {
		  if (nsp > 1) {
		    bornMultiSpecies = jph;
		  } else {
		    justBornMultiSpecies.push_back(jph);
		  }
		}
	      }
	    }
	    kIndexKin++;
	  }
	}
      }
    }
    /*
     *  Find the initial surface area to use
     */
    double sa_init = surfaceAreaRS_init_[0];

    /*
     *  Find the final surface area to use
     */
    double sa_final = calcSurfaceAreaChange(deltaTsubcycle);
    surfaceAreaRS_final_[0] = sa_final;

    /*
     *  Calculate the change in the moles of all of the species
     */
  
    for (int k = 0; k < m_NumTotSpecies; k++) {
      spMoles_tmp[k] = spMoles_init_[k];
      for (int isk = 0; isk < m_NumSurPhases; isk++) { 
	if (ActiveKineticsSurf_[isk]) {
	  sa_init =  surfaceAreaRS_init_[isk];
	  sa_final = surfaceAreaRS_final_[isk];
	  double *spNetProdPerArea = spNetProdPerArea_List_.ptrColumn(isk);
	  spMoles_tmp[k]  += 0.5 * deltaTsubcycle * (sa_init + sa_final) * spNetProdPerArea[k];
	  srcTerm[k] = 0.5 * (sa_init + sa_final) * spNetProdPerArea[k];
	}
      }
    }
      
    for (int k = 0; k < m_NumTotSpecies; k++) {
      resid[k] = (spMoles_final_[k] - spMoles_init_[k]) / deltaTsubcycle - srcTerm[k];
    }

    return 0;
  }

  //===================================================================================================================
  //  Calculate the change in the state of the system when integrating from T_global_initial to T_global_final
  /*
   *  All information is kept internal within this routine. This may be done continuously
   *  and the solution is not updated.
   *
   *  Note the tolerance parameters refers to the nonlinear solves within the calculation
   *  They do not refer to time step parameters.
   *
   *  @param deltaT        DeltaT for the global integration step.
   *  @param rtolResid     Relative tolerance for nonlinear solves within the calculation
   *                       Defaults to 1.0E-3
   *  @param atolResid     Absolute tolerance for nonlinear solves within the calculation
   *                       Defaults to 1.0E-12
   *
   *  @return Returns the number of subcycle steps it took to complete the full step.
   */
  int InterfacialMassTransfer::integrate(double deltaT, double GlobalRtolSrcTerm, double GlobalAtolSrcTerm,
					 int fieldInterpolationType, int subIntegrationType, SubIntegrationHistory * sih)
  {
    throw CanteraError("InterfacialMassTransfer::integrate()",
		       "base class called");
    pendingIntegratedStep_ = 1;
    return 0;
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
  void  InterfacialMassTransfer::resetStartingCondition(double Tinitial) {
    int i;

    if (pendingIntegratedStep_ != 1) {
#ifdef DEBUG_ELECTRODE
      // printf(" InterfacialMassTransfer::resetStartingCondition WARNING: resetStartingCondition called with no pending integration step\n");
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
      throw CanteraError("InterfacialMassTransfer::resetStartingCondition()", "tinit " + fp2str(Tinitial) +" not compat with t_final_final_ "
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
  //====================================================================================================================
  // Set the internal initial intermediate and initial global state from the internal final state
  /*
   *  (non-virtual function)  -> function should onionize in-first.
   *
   *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
   *  routine as well.
   *
   * @param setInitInit   Boolean indicating whether you should set the init_init state as well
   */
  void InterfacialMassTransfer::setInitStateFromFinal_Oin(bool setInitInit)
  {
    int i;
    for (int k = 0; k < m_NumTotSpecies; k++) {
      spMoles_init_[k] = spMoles_final_[k];
      spMf_init_[k] = spMf_final_[k];
    }
    for (int k = 0; k < nSpeciesA_; k++) {
      spMF_solnA_BC_init_[k] = spMF_solnA_BC_final_[k];
    }
    for (int k = 0; k < nSpeciesB_; k++) {
      spMF_solnB_BC_init_[k] = spMF_solnB_BC_final_[k];
    }

    for (i = 0; i < m_NumTotPhases; i++) {
      phaseMoles_init_[i] = phaseMoles_final_[i];
    }
    Pres_A_Interface_init_ =  Pres_A_Interface_final_;
    Pres_B_Interface_init_ =  Pres_B_Interface_final_;
    Pres_solnA_BC_init_ = Pres_solnA_BC_final_;
    Pres_solnB_BC_init_ = Pres_solnB_BC_final_;
    Velo_S_Interface_init_ = Velo_S_Interface_final_;
    Velo_A_Interface_init_ = Velo_A_Interface_final_;
    Velo_B_Interface_init_ = Velo_B_Interface_final_;
    BLThickness_A_init_ = BLThickness_A_final_;
    BLThickness_B_init_ = BLThickness_B_final_;
    for (i = 0; i < m_NumSurPhases; i++) {
      surfaceAreaRS_init_[i] = surfaceAreaRS_final_[i];
    }
    if (setInitInit) {
      for (int k = 0; k < m_NumTotSpecies; k++) {
	spMoles_init_init_[k] = spMoles_final_[k];
	spMf_init_init_[k] = spMf_final_[k];
      }
      for (int k = 0; k < nSpeciesA_; k++) {
       spMF_solnA_BC_init_init_[k] = spMF_solnA_BC_final_[k];
      } 
      for (int k = 0; k < nSpeciesB_; k++) {
        spMF_solnB_BC_init_init_[k] = spMF_solnB_BC_final_[k];
      }
      for (i = 0; i < m_NumTotPhases; i++) {
	phaseMoles_init_init_[i] = phaseMoles_final_[i];
      }
      Pres_A_Interface_init_init_ =  Pres_A_Interface_final_;
      Pres_B_Interface_init_init_ =  Pres_B_Interface_final_;
      Pres_solnA_BC_init_init_ = Pres_solnA_BC_final_;
      Pres_solnB_BC_init_init_ = Pres_solnB_BC_final_;
      Velo_S_Interface_init_init_ = Velo_S_Interface_final_;
      Velo_A_Interface_init_init_ = Velo_A_Interface_final_;
      Velo_B_Interface_init_init_ = Velo_B_Interface_final_;
      BLThickness_A_init_init_ = BLThickness_A_final_;
      BLThickness_B_init_init_ = BLThickness_B_final_;

      for (i = 0; i < m_NumSurPhases; i++) {
	surfaceAreaRS_init_init_[i] = surfaceAreaRS_final_[i];
      }
    } 
  }
  //====================================================================================================================
  // Set the internal initial intermediate and initial global state from the internal final state
  /*
   *  (virtual function) 
   *
   *  Set the intial state and the final_final from the final state. We also can set the init_init state from this
   *  routine as well.
   *
   * @param setInitInit   Boolean indicating whether you should set the init_init state as well
   */
  void InterfacialMassTransfer::setInitStateFromFinal(bool setInitInit)
  {
    setInitStateFromFinal_Oin(setInitInit);
  } 
  //====================================================================================================================
  //! Set the internal final intermediate and from the internal init state
  /*!
   *  (non-virtual function)  -> function should onionize in-first.
   *
   *  Set the final state from the init state. This is commonly called during a failed time step
   */
  void InterfacialMassTransfer::setFinalStateFromInit_Oin()
  {
    int i;
    for (int k = 0; k < m_NumTotSpecies; k++) {
      spMoles_final_[k] = spMoles_init_[k];
      spMf_final_[k] = spMf_init_[k];
      spMF_solnA_BC_final_[k] = spMF_solnA_BC_init_[k];
      spMF_solnB_BC_final_[k] = spMF_solnB_BC_init_[k];
    }
    for (i = 0; i < m_NumTotPhases; i++) {
      phaseMoles_final_[i] = phaseMoles_init_[i]; 
    }
    Pres_A_Interface_final_ = Pres_A_Interface_init_;
    Pres_B_Interface_final_ = Pres_B_Interface_init_;
    Pres_solnA_BC_final_ = Pres_solnA_BC_init_;
    Pres_solnB_BC_final_ = Pres_solnB_BC_init_;
    Velo_S_Interface_final_ = Velo_S_Interface_init_; 
    Velo_A_Interface_final_ = Velo_A_Interface_init_;
    Velo_B_Interface_final_ = Velo_B_Interface_init_;
    BLThickness_A_final_ = BLThickness_A_init_;
    BLThickness_B_final_ = BLThickness_B_init_;

    for (i = 0; i < m_NumSurPhases; i++) {
      surfaceAreaRS_final_[i] = surfaceAreaRS_init_[i];
    }
  }
  //====================================================================================================================
  // Set the internal final intermediate state from the internal init state
  /*
   *  (virtual function from Electrode) 
   *
   *  Set the final state from the init state. This is commonly called during a failed time step
   */
  void InterfacialMassTransfer::setFinalStateFromInit()
  {
    setFinalStateFromInit_Oin();
  }
  //====================================================================================================================
  // Set the internal initial intermediatefrom the internal initial global state
  /*
   *  Set the intial state from the init init state. We also can set the final state from this
   *  routine as well.
   *
   *  The final_final is not touched.
   *
   * @param setFinal   Boolean indicating whether you should set the final as well
   */
  void InterfacialMassTransfer::setInitStateFromInitInit(bool setFinal)
  {
    int i;
    for (int k = 0; k < m_NumTotSpecies; k++) {
      spMoles_init_[k] = spMoles_init_init_[k];
      spMf_init_[k] = spMf_init_init_[k];
      spMF_solnA_BC_init_[k] = spMF_solnA_BC_init_init_[k];
      spMF_solnB_BC_init_[k] = spMF_solnB_BC_init_init_[k];
    }
    for (i = 0; i < m_NumTotPhases; i++) {
      phaseMoles_init_[i] = phaseMoles_init_init_[i]; 
    }
    Pres_A_Interface_init_ = Pres_A_Interface_init_init_;
    Pres_B_Interface_init_ = Pres_B_Interface_init_init_;
    Pres_solnA_BC_init_ = Pres_solnA_BC_init_init_;
    Pres_solnB_BC_init_ = Pres_solnB_BC_init_init_;
    Velo_S_Interface_init_ = Velo_S_Interface_init_init_; 
    Velo_A_Interface_init_ = Velo_A_Interface_init_init_;
    Velo_B_Interface_init_ = Velo_B_Interface_init_init_;
    BLThickness_A_init_ = BLThickness_A_init_init_;
    BLThickness_B_init_ = BLThickness_B_init_init_;
    for (i = 0; i < m_NumSurPhases; i++) {
      surfaceAreaRS_init_[i] = surfaceAreaRS_init_init_[i];
    }
    if (setFinal) {
      for (int k = 0; k < m_NumTotSpecies; k++) {
	spMoles_final_[k] = spMoles_init_init_[k];
	spMf_final_[k] = spMf_init_init_[k];
	spMF_solnA_BC_final_[k] = spMF_solnA_BC_init_init_[k];
	spMF_solnB_BC_final_[k] = spMF_solnB_BC_init_init_[k];
      }
      for (i = 0; i < m_NumTotPhases; i++) {
	phaseMoles_final_[i] = phaseMoles_init_init_[i]; 
      }
      Pres_A_Interface_final_ = Pres_A_Interface_init_init_;
      Pres_B_Interface_final_ = Pres_B_Interface_init_init_;
      Pres_solnA_BC_final_ = Pres_solnA_BC_init_init_;
      Pres_solnB_BC_final_ = Pres_solnB_BC_init_init_;
      Velo_S_Interface_final_ = Velo_S_Interface_init_init_; 
      Velo_A_Interface_final_ = Velo_A_Interface_init_init_;
      Velo_B_Interface_final_ = Velo_B_Interface_init_init_;
      BLThickness_A_final_ = BLThickness_A_init_init_;
      BLThickness_B_final_ = BLThickness_B_init_init_;
      for (i = 0; i < m_NumSurPhases; i++) {
	surfaceAreaRS_final_[i] = surfaceAreaRS_init_init_[i];
      }
    }
  }
  //====================================================================================================================
  // Set the internal final global state from the internal final intermediate state
  /*
   *  (non-virtual function from Electrode) 
   *
   *  Set the final_final state from the final state. This is commonly called at the end of successful base integration.
   *  This is an onionize in function
   */
  void  InterfacialMassTransfer::setFinalFinalStateFromFinal_Oin()
  {
    for (int k = 0; k < m_NumTotSpecies; k++) {
      spMoles_final_final_[k] = spMoles_final_[k];
      spMF_solnA_BC_final_final_[k] = spMF_solnA_BC_final_[k];
      spMF_solnB_BC_final_final_[k] = spMF_solnB_BC_final_[k];
    }
    Pres_A_Interface_final_final_ = Pres_A_Interface_final_;
    Pres_B_Interface_final_final_ = Pres_B_Interface_final_;
    Pres_solnA_BC_final_final_ = Pres_solnA_BC_final_;
    Pres_solnB_BC_final_final_ = Pres_solnB_BC_final_;
    Velo_S_Interface_final_final_ = Velo_S_Interface_final_;
    Velo_A_Interface_final_final_ = Velo_A_Interface_final_;
    Velo_B_Interface_final_final_ = Velo_B_Interface_final_;
    BLThickness_A_final_final_ = BLThickness_A_final_;
    BLThickness_B_final_final_ = BLThickness_B_final_;
  }
  //====================================================================================================================
  // Set the internal final global state from the internal final intermediate state
  /*
   *  (virtual function from Electrode) 
   *
   *  Set the final_final state from the final state. This is commonly called at the end of successful base integration
   */
  void InterfacialMassTransfer::setFinalFinalStateFromFinal()
  {
    setFinalFinalStateFromFinal_Oin();
  }
  //====================================================================================================================
  // Report the integrated source term for the electrode over an interval in time.
  /*
   *  This is the net change in the moles of species defined in the electrode over that
   *  interval of time. The conditions at the end of the interval are used to carry
   *  out the integrations.
   *  
   *  @param spMoleDelta The end result in terms of the change in moles of species in the
   *                     electrode.
   *
   *  @return Tfinal    Final time to integrate to.                     
   */
  double InterfacialMassTransfer::integratedSourceTerm(doublereal* const spMoleDelta) const 
  {
    if (t_final_ == t_init_) {
      throw CanteraError(" InterfacialMassTransfer::integratedSourceTerm()", "tfinal == tinit");
    }
    /*
     *  We may do more here to ensure that the last integration is implicit
     */
    mdpUtil::mdp_copy_dbl_1(spMoleDelta, &(spMoleIntegratedSourceTerm_[0]), m_NumTotSpecies);
    return t_final_;
  }
  //====================================================================================================================
  void InterfacialMassTransfer::calcIntegratedSourceTerm() 
  {
    for (int i = 0; i < m_NumTotSpecies; i++) {
      spMoleIntegratedSourceTerm_[i] = spMoles_final_final_[i] - spMoles_init_init_[i];
    }
  }

  //====================================================================================================================
  // Calculate the integrated source term for the electrode over an interval in time.
  /*
   *  This is the net change in the moles of species defined in the electrode over that
   *  interval of time. The conditions at the beginning of the interval are used to carry
   *  out the integrations. This may be used at the start of the time step, since no
   *  unknown conditions are required.
   *
   *  @param deltaT     time to integrate
   *  @param spMoleDelta The end result in terms of the change in moles of species in the
   *                     electrode.
   *
   *  @return Tfinal    Final time to integrate to.
   */
  /*
  double InterfacialMassTransfer::integrateAndPredictSourceTerm(doublereal deltaT, doublereal* const spMoleDelta) 
  {
    integrate(deltaT);
    mdpUtil::mdp_copy_dbl_1(spMoleDelta, &(spMoleIntegratedSourceTerm_[0]), m_NumTotSpecies);
    return t_final_final_;
  }
  */
  //====================================================================================================================
  void InterfacialMassTransfer::printInterfacialMassTransferPhaseList(int pSrc, bool subTimeStep) 
  {
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
  void InterfacialMassTransfer::printInterfacialMassTransfer(int pSrc, bool subTimeStep) {
    int iph;
    double egv = TotalVol();
    printf("   ==============================================================================================\n");
    if (subTimeStep) {
      printf("      InterfacialMassTransfer at intermediate-step time final = %g\n", t_final_);
      printf("                   intermediate-step time init  = %g\n", t_init_);
    } else {
      printf("      InterfacialMassTransfer at time final = %g\n", t_final_final_);
      printf("                   time init  = %g\n", t_init_init_);
    }
    printf("\n");
    printf("                    DomainNumber = %2d , CellNumber = %2d , IntegrationCounter = %d\n", 
	   DomainNumber_, CellNumber_, counterNumberIntegrations_);
    printf("   ==============================================================================================\n");
    printf("          Number of surfaces = %d\n", numSurfaces_);

    printf("          Total Volume = %10.3E\n", egv);
    printf("          Temperature = %g\n", Temp_);
    printf("          PressureA = %g\n", Pres_A_Interface_final_);
    printf("          PressureB = %g\n", Pres_B_Interface_final_);

    printInterfacialMassTransferPhaseList(pSrc, subTimeStep);

    for (iph = 0; iph < m_NumTotPhases; iph++) {
      printInterfacialMassTransferPhase(iph, pSrc, subTimeStep);   
    }
  }
  //===================================================================================================================
 
  void InterfacialMassTransfer::printInterfacialMassTransferPhase(int iph, int pSrc, bool subTimeStep) {
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
  // Determines the level of printing for each step.
  /*
   *   0 -> absolutely nothing is printed for a single call to integrate.
   *   1 -> One line summary per integrate call
   *   2 -> short description, points of interest: Table of nonlinear solve - one line per iteration
   *   3 -> Table is included -> More printing per nonlinear iteration (default) that occurs during the table
   *   4 -> Summaries of the nonlinear solve iteration as they are occurring -> table no longer printed
   *   5 -> Algorithm information on the nonlinear iterates are printed out
   *   6 -> Additional info on the nonlinear iterates are printed out
   *   7 -> Additional info on the linear solve is printed out.
   *   8 -> Info on a per iterate of the linear solve is printed out.
   */
  void InterfacialMassTransfer::setPrintLevel(int printLvl) 
  {
    printLvl_ = printLvl;
  }
  //====================================================================================================================
  //   Set the deltaT used for the subcycle step
  /*
   *  The default setting for this is infinite. If it's infinite then deltaTSubcycle is set to the
   *  deltaT of the integration step. However, there are many times where a smaller natural delta T is
   *  required to solve the equation system
   */
  void InterfacialMassTransfer::setDeltaTSubcycle(doublereal deltaTsubcycle)
  {
    deltaTsubcycleNext_ = deltaTsubcycle;
    deltaTsubcycle_init_next_ = deltaTsubcycle;
    deltaTsubcycle_init_init_ = deltaTsubcycle;
  } 
  //====================================================================================================================
  //   Set the maximum deltaT used for the subcycle step
  void InterfacialMassTransfer::setDeltaTSubcycleMax(doublereal deltaTsubcycle)
  {
    deltaTsubcycleMax_ = deltaTsubcycle;
  }
  //====================================================================================================================
  //  Calculate the time derivatives of the mole numbers at the current subcycle step
  /*
   *   This may be used as a virtual function
   */
  void InterfacialMassTransfer::calculateTimeDerivatives(doublereal deltaTsubcycle)
  {
    if (deltaTsubcycle <= 0.0) {
      deltaTsubcycle = 1.0;
    }
    double invT = 1.0 /  deltaTsubcycle;
    for (int k = 0; k < m_NumTotSpecies; k++) {
      spMoles_dot_[k] = invT * (spMoles_final_[k] - spMoles_init_[k]);
    }
    for (int iph = 0; iph < m_NumTotPhases; iph++) {
      phaseMoles_dot_[iph] = invT * (phaseMoles_final_[iph] - phaseMoles_init_[iph]);
    }
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
  void InterfacialMassTransfer::writeCSVData(int itype) { 
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
