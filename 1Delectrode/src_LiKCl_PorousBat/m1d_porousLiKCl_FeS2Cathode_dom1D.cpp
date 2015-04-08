/*
 * m1d_porousLiKCl_FeS2Cathode_dom1D.cpp
 *
 */

//  This is a heavyweight base class that provides the function
//evaluation for a single bulk domain.
#include "m1d_porousLiKCl_FeS2Cathode_dom1D.h"
#include "m1d_BDT_porCathode_LiKCl.h"

#include "m1d_NodalVars.h"
#include "m1d_LocalNodeIndices.h"

#include "m1d_GlobalIndices.h"
#include "m1d_exception.h"
#include "m1d_DomainLayout.h"
#include "m1d_Comm.h"
#include "m1d_DomainLayout_LiKCl_PorousBat.h" 

#include "Electrode_InfCapacity.h"
#include "Electrode_SimplePhaseChangeDiffusion.h"
#include "Electrode_Factory.h"


#include "cantera/base/ctml.h"
#include "cantera/transport/Tortuosity.h"

//#include "cantera/base/mdp_allo.h"

//next two lines added for salt precipitation
#include "cantera/thermo.h"
#include "cantera/thermo/MargulesVPSSTP.h"
extern int flagPrecipitation;

#include "stdio.h"
#include "stdlib.h"

using namespace std;
using namespace Cantera;

#define DEBUG_HKM_LI

#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;
//=====================================================================================================================
namespace m1d
{

//=====================================================================================================================
porousLiKCl_FeS2Cathode_dom1D::porousLiKCl_FeS2Cathode_dom1D(BDT_porCathode_LiKCl & bdd) :
    porousElectrode_dom1D(bdd), 
    BDT_ptr_(0),
    nph_(0), nsp_(0),
  concTot_cent_(0.0),
  concTot_cent_old_(0.0), 
  icurrInterfacePerSurfaceArea_Cell_(0), xdelCell_Cell_(0),
  concTot_Cell_(0), concTot_Cell_old_(0),
  electrodeCrossSectionalArea_(0),
  capacityDischargedPA_Cell_(0),
  depthOfDischargePA_Cell_(0),
  capacityLeftPA_Cell_(0),
  capacityPA_Cell_(0),     
  depthOfDischarge_Electrode_(0.0),
  capacityZeroDoD_Electrode_(0.0),
  Fleft_cc_(0.0), Fright_cc_(0.0), Vleft_cc_(0.0),
  Vcent_cc_(0.0), Vright_cc_(0.0), VElectrodeLeft_cc_(0.0), VElectrodeCent_cc_(0.0), VElectrodeRight_cc_(0.0),
  t_final_(0.0),
  t_init_(0.0),
  Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0), spCharge_(0), 
  mfElectrolyte_Soln_Curr_(0),
  mfElectrolyte_Thermo_Curr_(0),
  mfElectrolyte_Soln_Cell_old_(0,0),
  mfElectrolyte_SolnDot_Curr_(0),
  pmVolElectrolyte_Curr_(0),
  phiElectrode_Curr_(0.0), 
  concTotDot_Curr_(0.0),
  conductivityElectrode_(1.0E2),
  gradT_trCurr_(0.0), gradV_trCurr_(0.0), gradVElectrode_trCurr_(0.0), gradX_trCurr_(0), Vdiff_trCurr_(0),
  jFlux_trCurr_(0), icurrElectrode_trCurr_(0.0), 
  nSpeciesElectrode_(0), nSurfsElectrode_(0),
  electrodeSpeciesMoleDelta_Cell_(0),
  icurrInterface_Cell_(0),
  phaseMoleFlux_(0), solnMoleFluxInterface_Cell_(0), icurrElectrode_CBL_(0), icurrElectrode_CBR_(0),
  icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0), deltaV_Cell_(0), Ess_Cell_(0), overpotential_Cell_(0),
  icurrRxn_Cell_(0), LiFlux_Cell_(0), solnTemp(0)
{
  BDT_ptr_ =  dynamic_cast<BDT_porCathode_LiKCl *> (&bdd);

  nsp_ = 3;
  nph_ = 1;
}
//=====================================================================================================================
porousLiKCl_FeS2Cathode_dom1D::porousLiKCl_FeS2Cathode_dom1D(const porousLiKCl_FeS2Cathode_dom1D &r) :
    porousElectrode_dom1D((BDT_porCathode_LiKCl &) r.BDD_), 
    BDT_ptr_(0),
    nph_(0), nsp_(0), 
    concTot_cent_(0.0),
  concTot_cent_old_(0.0),  
  icurrInterfacePerSurfaceArea_Cell_(0), xdelCell_Cell_(0),
  concTot_Cell_(0), 
  concTot_Cell_old_(0),
  electrodeCrossSectionalArea_(0),
  capacityDischargedPA_Cell_(0),
  depthOfDischargePA_Cell_(0),
  capacityLeftPA_Cell_(0),
  capacityPA_Cell_(0),
  depthOfDischarge_Electrode_(0.0),
  capacityZeroDoD_Electrode_(0.0),
  Fleft_cc_(0.0), Fright_cc_(0.0), Vleft_cc_(0.0),
  Vcent_cc_(0.0), Vright_cc_(0.0), VElectrodeLeft_cc_(0.0), VElectrodeCent_cc_(0.0), VElectrodeRight_cc_(0.0),
  t_final_(0.0),
  t_init_(0.0),
  Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0), spCharge_(0), 
  mfElectrolyte_Soln_Curr_(0),
  mfElectrolyte_Thermo_Curr_(0),
  mfElectrolyte_Soln_Cell_old_(0,0),
  mfElectrolyte_SolnDot_Curr_(0),
  pmVolElectrolyte_Curr_(0),
  phiElectrode_Curr_(0.0),
  concTotDot_Curr_(0.0),
  conductivityElectrode_(1.0E2),
  gradT_trCurr_(0.0), gradV_trCurr_(0.0), gradVElectrode_trCurr_(0.0), gradX_trCurr_(0), Vdiff_trCurr_(0),
  jFlux_trCurr_(0), icurrElectrode_trCurr_(0.0), 
  nSpeciesElectrode_(0), nSurfsElectrode_(0),
  electrodeSpeciesMoleDelta_Cell_(0),
  icurrInterface_Cell_(0),
  phaseMoleFlux_(0), solnMoleFluxInterface_Cell_(0), icurrElectrode_CBL_(0), icurrElectrode_CBR_(0),
  icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0), deltaV_Cell_(0), Ess_Cell_(0), overpotential_Cell_(0),
  icurrRxn_Cell_(0), LiFlux_Cell_(0), solnTemp(0)
{
  porousLiKCl_FeS2Cathode_dom1D::operator=(r);
}
//=====================================================================================================================
porousLiKCl_FeS2Cathode_dom1D::~porousLiKCl_FeS2Cathode_dom1D()
{ 
}
//=====================================================================================================================
porousLiKCl_FeS2Cathode_dom1D &
porousLiKCl_FeS2Cathode_dom1D::operator=(const porousLiKCl_FeS2Cathode_dom1D &r)
{
  if (this == &r) {
    return *this;
  }
  // Call the parent assignment operator
  porousElectrode_dom1D::operator=(r);

  BDT_ptr_ = r.BDT_ptr_;
  nph_ = r.nph_;
  nsp_ = r.nsp_;
  concTot_cent_ = r.concTot_cent_;
  concTot_cent_old_ = r.concTot_cent_old_;
  icurrInterfacePerSurfaceArea_Cell_ = r.icurrInterfacePerSurfaceArea_Cell_;
  xdelCell_Cell_ = r.xdelCell_Cell_;
  concTot_Cell_old_ = r.concTot_Cell_old_;
  electrodeCrossSectionalArea_ = r.electrodeCrossSectionalArea_;

  capacityDischargedPA_Cell_ = r.capacityDischargedPA_Cell_;
  depthOfDischargePA_Cell_ = r.depthOfDischargePA_Cell_;
  capacityLeftPA_Cell_ = r.capacityLeftPA_Cell_;
  capacityPA_Cell_ = r.capacityPA_Cell_;
  depthOfDischarge_Electrode_ = r.depthOfDischarge_Electrode_;
  capacityZeroDoD_Electrode_ = r.capacityZeroDoD_Electrode_;
  Fleft_cc_ = r.Fleft_cc_;
  Fright_cc_ = r.Fright_cc_;
  Vleft_cc_ = r.Vleft_cc_;
  Vcent_cc_ = r.Vcent_cc_;
  Vright_cc_ = r.Vright_cc_;
  VElectrodeLeft_cc_ = r.VElectrodeLeft_cc_;
  VElectrodeCent_cc_ = r.VElectrodeCent_cc_;
  VElectrodeRight_cc_ = r.VElectrodeRight_cc_;
  t_final_ = r.t_final_;
  t_init_ = r.t_init_;
  Xleft_cc_ = r.Xleft_cc_;
  Xcent_cc_ = r.Xcent_cc_;
  Xright_cc_ = r.Xright_cc_;
  spCharge_ = r.spCharge_;
  mfElectrolyte_Soln_Curr_ = r.mfElectrolyte_Soln_Curr_;
  mfElectrolyte_Thermo_Curr_ = r.mfElectrolyte_Thermo_Curr_;
  mfElectrolyte_Soln_Cell_old_ = r.mfElectrolyte_Soln_Cell_old_;
  mfElectrolyte_SolnDot_Curr_ = r.mfElectrolyte_SolnDot_Curr_;
  pmVolElectrolyte_Curr_ = r.pmVolElectrolyte_Curr_;
  phiElectrode_Curr_ = r.phiElectrode_Curr_;
  concTot_Curr_ = r.concTot_Curr_;
  concTotDot_Curr_ = r.concTotDot_Curr_;
  conductivityElectrode_ = r.conductivityElectrode_;
  gradT_trCurr_ = r.gradT_trCurr_;
  gradV_trCurr_ = r.gradV_trCurr_;
  gradVElectrode_trCurr_ = r.gradVElectrode_trCurr_;
  gradX_trCurr_ = r.gradX_trCurr_;
  Vdiff_trCurr_ = r.Vdiff_trCurr_;
  jFlux_trCurr_ = r.jFlux_trCurr_;
  icurrElectrode_trCurr_ = r.icurrElectrode_trCurr_;
  nSpeciesElectrode_ = r.nSpeciesElectrode_;
  nSurfsElectrode_ = r.nSurfsElectrode_;

  electrodeSpeciesMoleDelta_Cell_ = r.electrodeSpeciesMoleDelta_Cell_;

  icurrInterface_Cell_ = r.icurrInterface_Cell_;
  phaseMoleFlux_ = r.phaseMoleFlux_;
  solnMoleFluxInterface_Cell_ = r.solnMoleFluxInterface_Cell_;
  icurrElectrode_CBL_ = r.icurrElectrode_CBL_;
  icurrElectrode_CBR_ = r.icurrElectrode_CBR_;
  icurrElectrolyte_CBL_ = icurrElectrolyte_CBL_;
  icurrElectrolyte_CBR_ = r.icurrElectrolyte_CBR_;
  deltaV_Cell_ = r.deltaV_Cell_;
  Ess_Cell_ = r.Ess_Cell_;
  overpotential_Cell_ = r.overpotential_Cell_;
  icurrRxn_Cell_ = r.icurrRxn_Cell_;
  LiFlux_Cell_ = r.LiFlux_Cell_;
  solnTemp = r.solnTemp;

  throw CanteraError("", "not implemented");

  return *this;
}
//=====================================================================================================================
// Prepare all of the indices for fast calculation of the residual
/*
 *  Ok, at this point, we will have figured out the number of equations
 *  to be calculated at each node point. The object NodalVars will have
 *  been fully formed.
 *
 *  We use this to figure out what local node numbers/ cell numbers are
 *  needed and to set up indices for their efficient calling.
 *
 *  Child objects of this one will normally call this routine in a
 *  recursive fashion.
 */
void
porousLiKCl_FeS2Cathode_dom1D::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   */
  porousElectrode_dom1D::domain_prep(li_ptr);
  
  BDT_porCathode_LiKCl* BDD_FeS2_Cathode = dynamic_cast<BDT_porCathode_LiKCl *>(&(BDD_));
  if (!BDD_FeS2_Cathode) {
    throw CanteraError(" porousLiKCl_FeS2Cathode_dom1D::domain_prep()", "bad dynamic cast ");
  }

  Cantera::Electrode *ee = BDD_FeS2_Cathode->Electrode_;
  nSpeciesElectrode_ = ee->nSpecies();
  nSurfsElectrode_ = ee->nSurfaces();

  /*
   * Porous electrode domain prep
   */
  //porosity_Cell_.resize(NumLcCells, 0.5);
  //porosity_Cell_old_.resize(NumLcCells, 0.5);
  //surfaceArea_Cell_.resize(NumLcCells, 1.0E5);
  icurrInterfacePerSurfaceArea_Cell_.resize(NumLcCells, 0.0);
  xdelCell_Cell_.resize(NumLcCells, 0.0);
  concTot_Cell_.resize(NumLcCells, 0.0);
  concTot_Cell_old_.resize(NumLcCells, 0.0);
  capacityDischargedPA_Cell_.resize(NumLcCells, 0.0);
  depthOfDischargePA_Cell_.resize(NumLcCells, 0.0);
  capacityLeftPA_Cell_.resize(NumLcCells, 0.0);
  capacityPA_Cell_.resize(NumLcCells, 0.0);
  icurrInterface_Cell_.resize(NumLcCells, 0.0);
  solnMoleFluxInterface_Cell_.resize(NumLcCells, 0.0);
  Xleft_cc_.resize(3, 0.0);
  Xcent_cc_.resize(3, 0.0);
  Xright_cc_.resize(3, 0.0);
  spCharge_.resize(3, 0.0);
  spCharge_[0] = 1.0;
  spCharge_[1] = 1.0;
  spCharge_[2] = -1.0;

  mfElectrolyte_Soln_Curr_.resize(3, 0.0);
  mfElectrolyte_Thermo_Curr_.resize(3, 0.0);
  mfElectrolyte_Soln_Cell_old_.resize(3, NumLcCells, 0.0);
  mfElectrolyte_SolnDot_Curr_.resize(3, 0.0);
  pmVolElectrolyte_Curr_.resize(3, 0.0);

  gradX_trCurr_.resize(3, 0.0);
  Vdiff_trCurr_.resize(3, 0.0);
  jFlux_trCurr_.resize(3, 0.0);
 
  electrodeSpeciesMoleDelta_Cell_.resize(nSpeciesElectrode_ * NumLcCells, 0.0);
  phaseMoleFlux_.resize(20, 0.0);

  solnTemp.resize(10, 0.0);

  icurrElectrode_CBL_.resize(NumLcCells, 0.0);
  icurrElectrode_CBR_.resize(NumLcCells, 0.0);

  icurrElectrolyte_CBL_.resize(NumLcCells, 0.0);
  icurrElectrolyte_CBR_.resize(NumLcCells, 0.0);

  deltaV_Cell_ .resize(NumLcCells, 0.0);
  Ess_Cell_.resize(nSurfsElectrode_ * NumLcCells, 0.0);
  overpotential_Cell_.resize(nSurfsElectrode_ * NumLcCells, 0.0);
  icurrRxn_Cell_.resize(NumLcCells, 0.0);
  LiFlux_Cell_.resize(NumLcCells, 0.0);

  /*
   *  Set the velocity basis of the transport object. We are using
   *  mole-averaged velocities as the basis.
   */
  trans_->setVelocityBasis(ivb_);

  /*
   *  Set up the Electrode Object at each Cell
   */
  instantiateElectrodeCells();

}
//====================================================================================================================
//  An electrode object must be created and initialized for every cell in the domain
/*
 * Create electrode objects for every cell.  
 * Correct the volume and number of moles of 
 * active material within each of these electrode 
 * objects to correspond to the discretized volume.
 */
void
porousLiKCl_FeS2Cathode_dom1D::instantiateElectrodeCells() 
{ 
  /*
   * Indices and nodes used to determine the volume of each cell.
   */  
  int index_RightLcNode;
  int index_LeftLcNode;
  int index_CentLcNode;
  NodalVars *nodeLeft = 0;
  NodalVars *nodeCent = 0;
  NodalVars *nodeRight = 0;
  
  double xCellBoundaryL; //cell boundary left
  double xCellBoundaryR; //cell boundary right

  for (int iCell = 0; iCell <NumLcCells; iCell++) {
    // Note that we need to create a Factory method to instantiate the desired electrode type.
    // int nSurPhases = PSinput.cathode_input_->m_pl->nSurPhases();

    DomainLayout *DL = BDD_.DL_ptr_;

    DomainLayout_LiKCl_infPorousBat *dlc = dynamic_cast<DomainLayout_LiKCl_infPorousBat *>(DL);
    ProblemStatementCell *psc_ptr = dlc->pscInput_ptr_;
    ELECTRODE_KEY_INPUT *ci = psc_ptr->cathode_input_;

    Cantera::Electrode *ee  = newElectrodeObject(ci->electrodeModelName);
    if (!ee) {
      throw  m1d_Error("porousLiKCl_FeS2Cathode_dom1D::instantiateElectrodeCells()",
		       "Electrode factory method failed");
    }

 
    // Turn off printing from the electrode object
    ee->setPrintLevel(0);
    ee->setID(2, iCell);
    /*
     *  Set the initial subcycle deltaT policy
     */
    ee->choiceDeltaTsubcycle_init_ = 1;

    
    int retn = ee->electrode_model_create(PSinput.cathode_input_);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    retn = ee->setInitialConditions(PSinput.cathode_input_);
    if (retn == -1) {
      printf("exiting with error\n");
      exit(-1);
    }

    /* 
     * ComputTRODE_KEY_INPUT * newElectrodeKeyInputObject(std::string model, Electrode_Factory* f)

     */
    //Get the index for the center node 
    index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    // Get the pointer to the NodalVars object for the center node
    nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    
    /*
     *  ------------------- Get the index for the left node -----------------------------
     *    There may not be a left node if we are on the left boundary. In that case
     *    set the pointer to zero and the index to -1. Hopefully, we will get a segfault on an error.
     */
    index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
    if (index_LeftLcNode < 0) {
      nodeLeft = 0;
    } else {
      // get the node structure for the left node
      nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
    }
    /*
     * ------------------------ Get the indexes for the right node ------------------------------------
     */
    index_RightLcNode = Index_RightLcNode_LCO[iCell];
    if (index_RightLcNode < 0) {
      nodeRight = 0;
    } else {
      // get the node structure for the right node
      nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
    }
    
    /*
     * CALCULATE POSITION AND DELTA_X Variables 
     */
    if (nodeLeft) {
      xCellBoundaryL = 0.5 * (nodeLeft->xNodePos() + nodeCent->xNodePos());
    } else {
      xCellBoundaryL = nodeCent->xNodePos();
    }
    if (nodeRight == 0) {
      xCellBoundaryR = nodeCent->xNodePos();
    } else {
      xCellBoundaryR = 0.5 * (nodeRight->xNodePos() + nodeCent->xNodePos());
    }
    /*
     * Calculate the cell width
     */
    xdelCell_Cell_[iCell] = xCellBoundaryR - xCellBoundaryL;
  

    // Compute total electrode volume 
    double totalElectrodeVolume;
    double electrodeArea = -1.0;
    double porosity = -1.0;

    //see if we know the electrode area
    if (PSinput.cathode_input_->electrodeGrossArea > 0.0 ) {
      electrodeArea = PSinput.cathode_input_->electrodeGrossArea;
    } else if  (PSinput.cathode_input_->electrodeGrossDiameter > 0.0 ) {
      electrodeArea = Pi * 0.25 
	* PSinput.cathode_input_->electrodeGrossDiameter * PSinput.cathode_input_->electrodeGrossDiameter ;
    }

    /*
     * if we know the solid and electrode volume, compute porosity
     * otherwise, hopefully we know the porosity
     */
    if (PSinput.cathode_input_->electrodeGrossThickness > 0.0 && electrodeArea > 0.0) {
      totalElectrodeVolume = electrodeArea * PSinput.cathode_input_->electrodeGrossThickness;
      porosity = 1.0 - ee->SolidVol() / totalElectrodeVolume;
      printf("Cathode porosity is %f with %g m^3 solid volume and %g m^3 electrode volume.\n",porosity, ee->SolidVol(), totalElectrodeVolume );
      if ( porosity <= 0.0 ) {
	throw CanteraError("porousLiKCl_LiSiCathode_dom1D::initialConditions()",
			   "Computed porosity is not positive.");
      }	
      /*
       *  reset the moles in the electrode using the computed porosity 
       * and the electrodeWidth FOR THIS node.
       */
      ee->setElectrodeSizeParams(electrodeArea, xdelCell_Cell_[iCell], porosity);
    } else {
      /*
       * Case where we know porosity and not volume
       */
      porosity = PSinput.cathode_input_->porosity;
      if ( porosity <= 0.0 ) {
	throw CanteraError("porousLiKCl_LiSiCathode_dom1D::initialConditions()",
			   "Input porosity is not positive.");
      }
      if (electrodeArea <= 0.0) {
	electrodeArea = 1.0;
      }
      ee->setElectrodeSizeParams(electrodeArea, xdelCell_Cell_[iCell], porosity);
    }    
    /*
     *  Save the cross-sectional area of the elctrode to use throughout the code. It does not change within
     *  this calculation
     */
    electrodeCrossSectionalArea_ = electrodeArea;
    
    // update porosity as computed from electrode input
    porosity_Cell_[iCell] = porosity;
    // porosity_Cell_[iCell] = 0.64007;

    Electrode_Cell_[iCell] = ee;
  }
}
  //====================================================================================================================
  // Function that gets called at end the start of every time step
  /*
   *  This function provides a hook for a residual that gets called whenever a
   *  time step has been accepted and we are about to move on to the next time step.
   *  The call is made with the current time as the time
   *  that is accepted. The old time may be obtained from t and rdelta_t_accepted.
   *
   *  After this call interrogation of the previous time step's results will not be
   *  valid.
   *
   *   @param  doTimeDependentResid  This is true if we are solving a time dependent
   *                                 problem.
   *   @param  soln_ptr              Solution value at the current time
   *   @param  solnDot_ptr           derivative of the solution at the current time.
   *   @param  solnOld_ptr           Solution value at the old time step, n-1
   *   @param  t                     current time to be accepted, n
   *   @param  t_old                 previous time step value
   */
  void
  porousLiKCl_FeS2Cathode_dom1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector *soln_ptr,
						     const Epetra_Vector *solnDot_ptr, const Epetra_Vector *solnOld_ptr,
						     const double t, const double t_old) 
  {
    const Epetra_Vector &soln = *soln_ptr;
    //const doublereal * solnDot_Cellptr = 0;

    for (int iCell = 0; iCell < NumLcCells; iCell++) {
      cIndex_cc_ = iCell;
      /*
       *  ---------------- Get the index for the center node ---------------------------------
       */
      int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
      /*
       *   Get the pointer to the NodalVars object for the center node
       */
      NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];    
      /*
       *  Index of the first equation in the bulk domain of center node
       */
      int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode] + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];

      if (solnDot_ptr) {
	  //solnDot_Cellptr =  &((*solnDot_ptr)[indexCent_EqnStart_BD]);
      } else {
	  //solnDot_Cellptr = 0;
      }      
      //SetupThermoShop1(&(soln[indexCent_EqnStart_BD]), solnDot_Cellptr, 0);
      SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart_BD]));

      concTot_Cell_old_[iCell] = concTot_Curr_;
      porosity_Cell_old_[iCell] = porosity_Curr_;

      double *mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
      mfElectrolyte_Soln_old[0] = mfElectrolyte_Soln_Curr_[0];
      mfElectrolyte_Soln_old[1] = mfElectrolyte_Soln_Curr_[1];
      mfElectrolyte_Soln_old[2] = mfElectrolyte_Soln_Curr_[2];
      /*
       * Tell the electrode object to accept the current step and prep for the next step. 
       *
       * We might at this point do a final integration to make sure we nailed the conditions of the last step.
       * However, we will hold off at implementing this right now
       */
      Electrode * ee = Electrode_Cell_[iCell];
      ee->resetStartingCondition(t);
      //
      //  this is needed for a proper startup 
      //
      ee->updateState();
      //
      //  this is needed for a proper startup - sync initinit with final
      //
      ee->setInitStateFromFinal(true);
    }
  }  
//=====================================================================================================================
// Revert the domain object's conditions to the conditions at the start of the global time step
/*
 *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init_init.
 *
 *  Virtual from m1d_domain.h
 */
void
porousLiKCl_FeS2Cathode_dom1D::revertToInitialGlobalTime()
{
    for (int iCell = 0; iCell < NumLcCells; iCell++) {   
      /*
       * Tell the electrode object to forget about any recent step and revert to t_init_init conditions
       */
      Electrode* ee = Electrode_Cell_[iCell];
      ee->revertToInitialTime(true);
    }
}
//=====================================================================================================================
// Basic function to calculate the residual for the domain.
/*
 *  This class is used just for volumetric domains with an electrolyte.
 *
 *  All residual terms are written with the following sign convention
 *  based on keeping the time derivative term positive.
 *
 *       res = dcdt - dc2 /dx2 - src = 0
 *
 * @param res  Output vector containing the residual
 * @param doTimeDependentResid  boolean indicating whether the time
 *                         dependent residual is requested
 * @param soln_ptr     solution vector at which the residual should be
 *                     evaluated
 * @param solnDot_ptr  solution dot vector at which the residual should
 *                     be evaluated.
 *  @param t           time
 *  @param rdelta_t    inverse of delta_t
 *
 */
void
porousLiKCl_FeS2Cathode_dom1D::residEval(Epetra_Vector &res,
                                         const bool doTimeDependentResid,
                                         const Epetra_Vector *soln_ptr,
                                         const Epetra_Vector *solnDot_ptr,
                                         const Epetra_Vector *solnOld_ptr,
                                         const double t,
                                         const double rdelta_t,
                                         const ResidEval_Type_Enum residType,
					 const Solve_Type_Enum solveType)
{
  residType_Curr_ = residType;
  int index_RightLcNode;
  int index_LeftLcNode;
  int index_CentLcNode;

  t_final_ = t;
  if (rdelta_t < 1.0E-200) {
    t_init_ = t;
  } else {
    t_init_ = t - 1.0/rdelta_t;
  }

  NodalVars *nodeLeft = 0;
  NodalVars *nodeCent = 0;
  NodalVars *nodeRight = 0;

  double xdelL; // Distance from the center node to the left node
  double xdelR; // Distance from the center node to the right node
  double xdelCell; // cell width - right boundary minus the left boundary.
  double xCellBoundaryL; //cell boundary left
  double xCellBoundaryR; //cell boundary right

  //  Electrolyte mole fluxes - this is c V dot n at the boundaries of the cells
  double fluxFright = 0.;
  double fluxFleft;

  //mole fraction fluxes
  std::vector<double> fluxXright(nsp_, 0.0);
  std::vector<double> fluxXleft(nsp_, 0.0);

  double fluxL = 0.0;
  double fluxR = 0.0;

  const Epetra_Vector &soln = *soln_ptr;
  const doublereal * solnDot_Cellptr = 0;
  /*
   * Index of the first equation at the left node corresponding to the first bulk domain, which is the electrolyte
   */
  int indexLeft_EqnStart_BD;
  /*
   * Index of the first equation at the center node corresponding to the first bulk domain, which is the electrolyte
   */
  int indexCent_EqnStart_BD;
  /*
   * Index of the first equation at the right node corresponding to the first bulk domain, which is the electrolyte
   */
  int indexRight_EqnStart_BD;

  /*
   *   Find the species index for the first species in the electrode object pertaining to the electrolyte
   */
  int sf = Electrode_Cell_[0]->solnPhaseIndex();
  int speciesIndexLip = Electrode_Cell_[0]->getGlobalSpeciesIndex(sf, 0);

  /*
   * offset of the electolyte solution unknowns at the current node
   */
  index_CentLcNode = Index_DiagLcNode_LCO[0];
  nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

  /*
   *  Offsets for the equation unknowns in the residual vector for the electrolyte domain
   */
  int EQ_Current_offset_BD = BDD_.EquationIndexStart_EqName[Current_Conservation];
  int EQ_Current_offset_ED = EQ_Current_offset_BD + 1;
  int EQ_TCont_offset_BD = BDD_.EquationIndexStart_EqName[Continuity];
  int EQ_Species_offset_BD = BDD_.EquationIndexStart_EqName[Species_Conservation];
  int EQ_MFSum_offset_BD = BDD_.EquationIndexStart_EqName[MoleFraction_Summation];
  int EQ_ChargeBal_offset_BD = BDD_.EquationIndexStart_EqName[ChargeNeutrality_Summation];

  /*
   * Offsets for the variable unknowns in the solution vector for the electrolyte domain
   */
  int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
  int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
  int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];

  Fright_cc_ = 0.0;

  /*
   *  -------------------- Special section to do the left boundary -------------------------------
   */

  /*
   * Special section if we own the left node of the domain. If we do
   * it will be cell 0
   */
  if (IOwnLeft) {
    DiffFluxLeftBound_LastResid_NE[0] = fluxL;
  }
  /*
   *  Boolean flag that specifies whether the left flux should be recalculated.
   *  Usually you don't need to calculate the left flux because it can be copied from the right flux
   *  from the previous cell.
   */
  bool doLeftFluxCalc = true;
  /*
   *  ------------------------------ LOOP OVER CELL -------------------------------------------------
   *  Loop over the number of Cells in this domain on this processor
   *  This loop is done from left to right.
   */
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;
    Electrode * Electrode_ptr = Electrode_Cell_[iCell];

#ifdef DEBUG_RESID
    if (counterResBaseCalcs_ > 125 && residType == Base_ResidEval) {
      if (iCell == NumLcCells - 1) {
       // printf("we are here\n");
      }
    }
#endif

    /*
     *  ---------------- Get the index for the center node ---------------------------------
     */
    index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    /*
     *  ------------------- Get the index for the left node -----------------------------
     *    There may not be a left node if we are on the left boundary. In that case
     *    set the pointer to zero and the index to -1. Hopefully, we will get a segfault on an error.
     */
    index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
    if (index_LeftLcNode < 0) {
      /*
       *  We assign node object to zero.
       */
      nodeLeft = 0;
      /*
       *  If there is no left node, we assign the left solution index to the center solution index
       */
      indexLeft_EqnStart_BD = indexCent_EqnStart_BD;
    } else {
      // get the node structure for the left node
      nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
      //index of first equation in the electrolyte of the left node
      indexLeft_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_LeftLcNode]
          + nodeLeft->OffsetIndex_BulkDomainEqnStart_BDN[0];
    }
    /*
     * If we are past the first cell, then we have already done the calculation
     * for this flux at the right cell edge of the previous cell
     */
    if (iCell > 0) {
      doLeftFluxCalc = false;
    }

    /*
     * ------------------------ Get the indexes for the right node ------------------------------------
     */
    index_RightLcNode = Index_RightLcNode_LCO[iCell];
    if (index_RightLcNode < 0) {
      nodeRight = 0;
      /*
       *  If there is no right node, we assign the right solution index to the center solution index
       */
      indexRight_EqnStart_BD = indexCent_EqnStart_BD;
    } else {
      //NodalVars
      nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
      //index of first equation of right node
      indexRight_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_RightLcNode]
          + nodeRight->OffsetIndex_BulkDomainEqnStart_BDN[0];
    }

    /*
     * --------------------------- CALCULATE POSITION AND DELTA_X Variables -----------------------------
     * Calculate the distance between the left and center node points
     */
    if (nodeLeft) {
      xdelL = nodeCent->xNodePos() - nodeLeft->xNodePos();
      xCellBoundaryL = 0.5 * (nodeLeft->xNodePos() + nodeCent->xNodePos());
    } else {
      xdelL = 0.0;
      xCellBoundaryL = nodeCent->xNodePos();
    }
    /*
     * Calculate the distance between the right and center node points
     */
    if (nodeRight == 0) {
      xdelR = 0.0;
      xCellBoundaryR = nodeCent->xNodePos();
    } else {
      xdelR = nodeRight->xNodePos() - nodeCent->xNodePos();
      xCellBoundaryR = 0.5 * (nodeRight->xNodePos() + nodeCent->xNodePos());
    }
    /*
     * Calculate the cell width
     */
    xdelCell = xCellBoundaryR - xCellBoundaryL;
    xdelCell_Cell_[iCell] = xdelCell;

    /*
     * --------------------------- DO PRE-SETUPSHOP RASTER OVER LEFT,CENTER,RIGHT -----------------------------
     * Calculate the distance between the left and center node points
     */
    /*
     * Get current velocity, mole fraction, temperature, potential
     * from the solution
     */

    for (int k = 0; k < nsp_; k++) {
      Xcent_cc_[k] = soln[indexCent_EqnStart_BD + iVar_Species_BD + k];
    }
    Vcent_cc_ = soln[indexCent_EqnStart_BD + iVar_Voltage_BD];
    VElectrodeCent_cc_ = soln[indexCent_EqnStart_BD + iVar_Voltage_BD + 1];


    /*
     * Advance the electrode forward and compute the new porosity 
     */ 
    if (solnDot_ptr) {
      solnDot_Cellptr =  &((*solnDot_ptr)[indexCent_EqnStart_BD]);
    } else {
      solnDot_Cellptr = 0;
    }
    SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart_BD]));
    
    
    if (nodeLeft != 0) {
      /*
       * Find the velocity located at the left cell boundary.
       * The left cell boundary velocity is stored at the previous (left)
       * cell index as per our conventions.
       */
      Fleft_cc_ = soln[indexLeft_EqnStart_BD + iVAR_Vaxial_BD];
      for (int k = 0; k < nsp_; k++) {
        Xleft_cc_[k] = soln[indexLeft_EqnStart_BD + iVar_Species_BD + k];
      }
      Vleft_cc_ = soln[indexLeft_EqnStart_BD + iVar_Voltage_BD];
      VElectrodeLeft_cc_ = soln[indexLeft_EqnStart_BD + iVar_Voltage_BD + 1];
    } else {
      /*
       * We are here when we are at the left most part of the boundary. Then, there is no
       * left node. Or, the left node is actually the center node. The boundary conditions
       * on the fluxes are set by boundary condition routines.
       */
      /*
       *   We set the left flux to zero in this case. It may or may not be overridden by
       *   a boundary condition in the adjoining surface domain
       */
      Fleft_cc_ = 0.0;

      for (int k = 0; k < nsp_; k++) {
        Xleft_cc_[k] = Xcent_cc_[k];
      }
      Vleft_cc_ = Vcent_cc_;
      VElectrodeLeft_cc_ = VElectrodeCent_cc_;
    }
    /*
     * The right velocity is stored at the current cell index.
     */
    Fright_cc_ = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];

    if (nodeRight != 0) {

      for (int k = 0; k < nsp_; k++) {
        Xright_cc_[k] = soln[indexRight_EqnStart_BD + iVar_Species_BD + k];
      }
      Vright_cc_ = soln[indexRight_EqnStart_BD + iVar_Voltage_BD];
      VElectrodeRight_cc_ = soln[indexRight_EqnStart_BD + iVar_Voltage_BD + 1];
    } else {

      for (int k = 0; k < nsp_; k++) {
        Xright_cc_[k] = Xcent_cc_[k];
      }
      Vright_cc_ = Vcent_cc_;
      VElectrodeRight_cc_ = VElectrodeCent_cc_;
    }

    /*
     * ------------------- CALCULATE FLUXES AT THE LEFT BOUNDARY -------------------------------
     *
     */
    if (doLeftFluxCalc) {
      if (nodeLeft == 0) {
        /*
         *  We are here if we are at the left node boundary and we
         *  need a flux condition. The default now is to
         *  set the flux to zero. We could put in a more
         *  sophisticated treatment
         */
        fluxFleft = 0.0;
        icurrElectrolyte_CBL_[iCell] = 0.0;
        icurrElectrode_CBL_[iCell] = 0.0;
        for (int k = 0; k < nsp_; k++) {
          fluxXleft[k] = 0.0;
        }
      } else {

        /*
         *  Establish the environment at the left cell boundary
         */
        SetupThermoShop2(&(soln[indexLeft_EqnStart_BD]), &(soln[indexCent_EqnStart_BD]), 0);

        SetupTranShop(xdelL, 0);

        /*
         * Calculate the flux at the left boundary for each equation
         */
        fluxFleft = Fleft_cc_ * concTot_Curr_;

        /*
         * Calculate the flux of species and the flux of charge
         *   - the flux of charge must agree with the flux of species
         */
        icurrElectrolyte_CBL_[iCell] = 0.0;
        icurrElectrode_CBL_[iCell] = icurrElectrode_trCurr_;
        for (int k = 0; k < nsp_; k++) {
          fluxXleft[k] = jFlux_trCurr_[k];
          icurrElectrolyte_CBL_[iCell] += fluxXleft[k] * spCharge_[k];
          if (Fleft_cc_ > 0.0) {
            fluxXleft[k] += Fleft_cc_ * Xleft_cc_[k] * concTot_Curr_;
          } else {
            fluxXleft[k] += Fleft_cc_ * Xcent_cc_[k] * concTot_Curr_;
          }
        }
        icurrElectrolyte_CBL_[iCell] *= (Cantera::Faraday);
      }
    } else {  // !doLeftFluxCalc
      /*
       * Copy the fluxes from the stored right side
       */
      fluxFleft = fluxFright;
      icurrElectrolyte_CBL_[iCell] = icurrElectrolyte_CBR_[iCell - 1];
      icurrElectrode_CBL_[iCell] = icurrElectrode_CBR_[iCell - 1];
      for (int k = 0; k < nsp_; k++) {
        fluxXleft[k] = fluxXright[k];
      }
    }
    /*
     * ------------------- CALCULATE FLUXES AT THE RIGHT BOUNDARY -------------------------------
     *
     */

    /*
     * Calculate the flux at the right cell boundaries
     */
    if (nodeRight == 0) {
      /*
       *  We are here if we are at the right node boundary and we need a flux
       *  condition. The default now is to set the flux to zero. We could
       *  put in a more sophisticated treatment
       */
      AssertTrace(iCell == NumLcCells-1);
      // fluxFright = 0.0;
      if (solnDot_ptr) {
	solnDot_Cellptr =  &((*solnDot_ptr)[indexCent_EqnStart_BD]);
      } else {
	solnDot_Cellptr = 0;
      }
      SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart_BD]));
      fluxFright = Fright_cc_ * concTot_Curr_;
      // fluxVRight = 0.0;
      icurrElectrolyte_CBR_[iCell] = 0.0;
      //  fluxVElectrodeRight = 0.0;
      icurrElectrode_CBR_[iCell] = 0.0;
      for (int k = 0; k < nsp_; k++) {
        fluxXright[k] = 0.0;
      }
    } else {
      /*
       *  Establish the environment at the right cell boundary
       */
      SetupThermoShop2(&(soln[indexCent_EqnStart_BD]), &(soln[indexRight_EqnStart_BD]), 1);

      SetupTranShop(xdelR, 1);

      /*
       * Calculate the flux at the right boundary for each equation
       * This is equal to
       *       Conc * Vaxial * phi
       */
      fluxFright = Fright_cc_ * concTot_Curr_;

      /*
       * Calculate the flux of species and the flux of charge
       *   - the flux of charge must agree with the flux of species
       */
      icurrElectrolyte_CBR_[iCell] = 0.0;
      icurrElectrode_CBR_[iCell] = icurrElectrode_trCurr_;
      for (int k = 0; k < nsp_; k++) {
        fluxXright[k] = jFlux_trCurr_[k];
        // fluxVRight += fluxXright[k] * spCharge_[k];
        icurrElectrolyte_CBR_[iCell] += fluxXright[k] * spCharge_[k];
        if (Fright_cc_ > 0.0) {
          //fluxXright[k] += Fright_cc_ * Xcent_cc_[k] * concTot_Curr_ * porosity_Curr_;
          // Why is this mfElectrolyte_Thermo_Curr and not Xright_cc?
          fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_;
        } else {
          fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_;
        }
      }
      icurrElectrolyte_CBR_[iCell] *= (Cantera::Faraday);
    }

#ifdef DEBUG_RESID
    if (doTimeDependentResid) {

      //     if (residType == Base_ResidEval) {
      //      printf(" Cell = %d, Totalflux_K+ = %10.3e,  Totalflux_Cl- = %10.3e \n", iCell, fluxXright[1], fluxXright[2]);
      //	  printf("           Old porosity = %10.3e,  New porosity = %10.3e \n", porosity_Cell_old_[iCell], porosity_Cell_[iCell] );
      //       printf("           Vmolal = %10.3e, jd_Li+ = %10.3e  jd_K+ = %10.3e jd_Cl- = %10.3e\n", Fright_cc_,
      //             jFlux_trCurr_[0], jFlux_trCurr_[1], jFlux_trCurr_[2]);
      //      printf("           Vmolal = %10.3e, vd_Li+ = %10.3e  vd_K+ = %10.3e vd_Cl- = %10.3e\n", Fright_cc_,
      //             Vdiff_trCurr_[0], Vdiff_trCurr_[1], Vdiff_trCurr_[2]);
      //   }
    }
#endif

    /*
     * --------------- ADD FLUX TERMS INTO THE RESIDUALS --------------------------------------
     */

    /*
     *  Total continuity equation - fluxFright and fluxFleft represent the total mole
     *                              fluxes coming and going from the cell.
     *                    R =   d rho d t + del dot (rho V) = 0
     */
    res[indexCent_EqnStart_BD + EQ_TCont_offset_BD] += (fluxFright - fluxFleft);

    /*
     * Species continuity Equation - 2
     */
    for (int k = 0; k < nsp_ - 2; k++) {
      res[indexCent_EqnStart_BD + EQ_Species_offset_BD + k] += (fluxXright[k] - fluxXleft[k]);
    }

    /*
     * Mole fraction summation equation
     */

    /*
     * Electroneutrality equation
     */

    /*
     *   Current conservation equation - electrolyte
     */
    // res[indexCent_EqnStart_BD + EQ_Current_offset_BD] += (fluxVRight - fluxVLeft);
    res[indexCent_EqnStart_BD + EQ_Current_offset_BD] += (icurrElectrolyte_CBR_[iCell] - icurrElectrolyte_CBL_[iCell]);

    /*
     *   Current conservation equation - electrode
     */
    res[indexCent_EqnStart_BD + EQ_Current_offset_ED] += (icurrElectrode_CBR_[iCell] - icurrElectrode_CBL_[iCell]);

    /*
     *   ------------------- ADD SOURCE TERMS TO THE CURRENT CELL CENTER --------------------------------------
     */ 
    if (solnDot_ptr) {
      solnDot_Cellptr =  &((*solnDot_ptr)[indexCent_EqnStart_BD]);
    } else {
      solnDot_Cellptr = 0;
    }
    //SetupThermoShop1Old(&(soln[indexCent_EqnStart_BD]), solnDot_Cellptr, 0);
    SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart_BD]));
    getElectrolyte_SolnDot(nodeCent, solnDot_Cellptr);
    /*
     * Species 0 Conservation equation
     *    Source terms for the species production rate of Li+.
     */
    res[indexCent_EqnStart_BD + EQ_Species_offset_BD] -= electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * iCell + speciesIndexLip] * rdelta_t
      / electrodeCrossSectionalArea_;

    /*
     * Mole fraction summation equation
     *       For the DAE problem we specify the time derivative formulation. 
     *       For the regular problem, we specify the mole fraction formulation.
     *       Both are equivalent - > set m_isAlgebraic = 2
     */
    if (solveType == DAESystemInitial_Solve) { 
      res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] = 0.0;
      for (int k = 0; k < nsp_; k++) {
	res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] += mfElectrolyte_SolnDot_Curr_[k]; 
      }
    } else {
      res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] = 1.0;
      for (int k = 0; k < nsp_; k++) {
	res[indexCent_EqnStart_BD + EQ_MFSum_offset_BD] -= Xcent_cc_[k];
      }
    }

    /*
     * Electroneutrality equation
     *       For the DAE problem we specify the time derivative formulation. 
     *       For the regular problem, we specify the mole fraction formulation.
     *       Both are equivalent - > set m_isAlgebraic = 2
     */
    if (solveType == DAESystemInitial_Solve) {
      res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] = 0.0;
      for (int k = 0; k < nsp_; k++) {
	res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] += mfElectrolyte_SolnDot_Curr_[k] * spCharge_[k];
      }
    } else {
      res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] = 0.0;
      for (int k = 0; k < nsp_; k++) {
	res[indexCent_EqnStart_BD + EQ_ChargeBal_offset_BD] += Xcent_cc_[k] * spCharge_[k];
      }
    }

    /*
     *  Total continuity equation for moles in the electrolyte phase of the cell
     *     Add in the molar flux from the electrode into the electrolyte phase
     *     We are assuming for the current problem that the volumes stay constant
     */
    res[indexCent_EqnStart_BD + EQ_TCont_offset_BD] -= solnMoleFluxInterface_Cell_[iCell] * rdelta_t;

    /*
     *   Current conservation equation
     *      These are written as a source term. icurrInterface_Curr_ will be positive for the anode
     *      where current flows into the electrolyte and will be negative for the cathode where current
     *      flows out of the cathode
     */
    res[indexCent_EqnStart_BD + EQ_Current_offset_BD] -= icurrInterface_Cell_[iCell];

    /*
     *   Current conservation equation for the current in the electrode material
     *      These are written as a sink term They will be exactly opposite to the electrolyte current
     *      source terms
     */
    res[indexCent_EqnStart_BD + EQ_Current_offset_ED] += icurrInterface_Cell_[iCell];

    /*
     * Special section if we own the left node of the domain. If we do
     * it will be cell 0. Store currents for later use.
     * These are the correct currents that work for the global balances
     */
    if (IOwnLeft && iCell == 0) {
      if (residType == Base_ShowSolution || residType == Base_ResidEval) {
        icurrElectrolyte_CBL_[iCell] = -icurrInterface_Cell_[iCell]  + icurrElectrolyte_CBR_[iCell];
      }
      if (residType == Base_ShowSolution || residType == Base_ResidEval) {
        DiffFluxLeftBound_LastResid_NE[EQ_Current_offset_BD] = icurrElectrolyte_CBL_[iCell];
        DiffFluxLeftBound_LastResid_NE[EQ_Current_offset_ED] = icurrElectrode_CBL_[iCell];
      }
    }
    /*
     * Special section if we own the right node of the domain. If we do
     * it will be cell (NumLcCells - 1). Store currents for later use.
     * These are the correct currents that work for the global balances
     */
    if (IOwnRight && iCell == (NumLcCells - 1)) {

      if (residType == Base_ShowSolution || residType == Base_ResidEval) {
        icurrElectrode_CBR_[iCell] += -icurrInterface_Cell_[iCell]  + icurrElectrode_CBL_[iCell];
      }
      if (residType == Base_ShowSolution || residType == Base_ResidEval) {
        DiffFluxRightBound_LastResid_NE[EQ_Current_offset_BD] = icurrElectrolyte_CBR_[iCell];
        DiffFluxRightBound_LastResid_NE[EQ_Current_offset_ED] = icurrElectrode_CBR_[iCell];
      }
    }
    if (residType == Base_ShowSolution) {
      deltaV_Cell_[iCell] = Electrode_ptr->potentialDrop();
      for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++ ) { 
	Ess_Cell_[nSurfsElectrode_ * iCell + jSurf] = Electrode_ptr->openCircuitVoltage(jSurf);
	overpotential_Cell_[nSurfsElectrode_ * iCell + jSurf] = Electrode_ptr->overpotential(jSurf);
      }
      icurrRxn_Cell_[iCell] = icurrInterface_Cell_[iCell];
      LiFlux_Cell_[iCell] = jFlux_trCurr_[0];
    }

    /*
     *   ------------------ ADD IN THE TIME DEPENDENT TERMS ----------------------------------------------------
     */
    if (doTimeDependentResid) {

      /*
       *   .................... Calculate quantities needed at the current time
       */
#ifdef DEBUG_HKM_NOT
      if (residType == Base_ResidEval) {
        printf(" Cell = %d, Totalflux_Li+_r = %10.3e,  = %10.3e, Totalflux_Li+_l ", iCell, fluxXright[0], fluxXleft[0]);
      }
#endif
#ifdef DEBUG_HKM_LI
      if (residType == Base_ShowSolution) {
	if (iCell == 7) {
	  printf(" Cell = %d, Totalflux_Li+_R = %- 14.7e, -Totalflux_Li+_L = %- 14.7e", iCell, fluxXright[0], -fluxXleft[0]);
          printf(", -Source/dt = %- 14.7e",
		 -electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * iCell + speciesIndexLip] * rdelta_t
		 / electrodeCrossSectionalArea_);
	}
      }
#endif

      double newStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
      double newStuffSpecies0 = Xcent_cc_[0] * newStuffTC;

      /*
       *   .................... Calculate quantities needed at the previous time step
       */
  
      double *mf_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
      double oldStuffTC = concTot_Cell_old_[iCell] * porosity_Cell_old_[iCell] * xdelCell;
      double oldStuffSpecies0 = mf_old[0] * oldStuffTC;

      double tmp = (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;

#ifdef DEBUG_HKM_NOT
      if (residType == Base_ResidEval) {
        printf(" deltaT term = %10.3e BulkSum = %10.3e\n", tmp, tmp + (fluxXright[0] - fluxXleft[0]));
      }
#endif
#ifdef DEBUG_HKM_LI
      if (residType == Base_ShowSolution) {
	if (iCell == 7) {
	  printf(", d(XCP)/dt = %- 14.7e", tmp);
	  double sumResidualX = tmp +  fluxXright[0] - fluxXleft[0] -
	    electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * iCell + speciesIndexLip] * rdelta_t
	    / electrodeCrossSectionalArea_;
	  printf(" ==  %- 14.7e\n", sumResidualX);
	}
      }
#endif

#ifdef DEBUG_HKM_LI
      if (residType == Base_ShowSolution) {
	if (iCell == 7) {
	  printf(" Cell = %d, Totalflux_R = %- 14.7e, -Totalflux_L = %- 14.7e", iCell, fluxFright, -fluxFleft);
          printf(", -Source/dt = %- 14.7e", -solnMoleFluxInterface_Cell_[iCell] * rdelta_t);
	  double sumResidualTC = fluxFright - fluxFleft - solnMoleFluxInterface_Cell_[iCell] * rdelta_t;
	  printf(" ==  %- 14.7e\n", sumResidualTC);
	}
      }
#endif

	/*
	 *  For this calculation we do the time derivative another way
	 */
      if (solveType == DAESystemInitial_Solve) {

	double tmp2;
	/*
	 * Species continuity Equation - 2
	 */
	for (int k = 0; k < nsp_ - 2; k++) {
	  tmp2 = xdelCell * (concTot_Curr_    * porosity_Curr_ * mfElectrolyte_SolnDot_Curr_[k] +
			     concTotDot_Curr_ * porosity_Curr_ * Xcent_cc_[k] );
	  res[indexCent_EqnStart_BD + EQ_Species_offset_BD + k] += tmp2;
	}	
     
      } else {
        /*
         *  Add in the time term for species 0
         */
	res[indexCent_EqnStart_BD + EQ_Species_offset_BD + 0] += tmp;
      }

      /*
       *   Add in the time term for the total continuity equation
       *         note: the current problem will have this term equalling zero always.
       *               However, we put it in here for the next problem.
       */
      res[indexCent_EqnStart_BD + EQ_TCont_offset_BD] += (newStuffTC - oldStuffTC) * rdelta_t;

    }

  }

  /*
   * Special section to do the right boundary
   */
  /*
   * Special section if we own the left node of the domain. If we do
   * it will be last cell
   */
  if (IOwnRight) {
    DiffFluxRightBound_LastResid_NE[0] = fluxR;
  }

}
//====================================================================================================================
// Utility function to calculate quantities before the main residual routine.
/*
 *  This is used for a loop over nodes. All calculated quantities must be internally storred.
 *
 *  Currently this is called during the residual evalultion of the problem.
 * 
 * @param res  Output vector containing the residual
 * @param doTimeDependentResid  boolean indicating whether the time
 *                         dependent residual is requested
 * @param soln_ptr     solution vector at which the residual should be
 *                     evaluated
 * @param solnDot_ptr  solution dot vector at which the residual should
 *                     be evaluated.
 * @param solnOld_ptr  Pointer to the solution vector at the old time step
 *  @param t           time
 *  @param rdelta_t    inverse of delta_t
 */
void
porousLiKCl_FeS2Cathode_dom1D::residEval_PreCalc(const bool doTimeDependentResid,
						 const Epetra_Vector *soln_ptr,
						 const Epetra_Vector *solnDot_ptr,
						 const Epetra_Vector *solnOld_ptr,
						 const double t,
						 const double rdelta_t,
						 ResidEval_Type_Enum residType,
                                                 const Solve_Type_Enum solveType)
{
  residType_Curr_ = residType;
  const Epetra_Vector &soln = *soln_ptr;
  //const doublereal * solnDot_Cellptr = 0;
  int index_RightLcNode;
  int index_LeftLcNode;
  int index_CentLcNode;

  maxElectrodeSubIntegrationSteps_ = 0;

  t_final_ = t;
  if (rdelta_t < 1.0E-200) {
    t_init_ = t;
  } else {
    t_init_ = t - 1.0/rdelta_t;
  }

  NodalVars *nodeLeft = 0;
  NodalVars *nodeCent = 0;
  NodalVars *nodeRight = 0;

  double xCellBoundaryL; //cell boundary left
  double xCellBoundaryR; //cell boundary right
  

  /*
   * Index of the first equation at the center node corresponding to the first bulk domain, which is the electrolyte
   */
  int indexCent_EqnStart_BD;

  /*
   * Offsets for the variable unknowns in the solution vector for the electrolyte domain
   */
  // int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
  int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
  int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
  

  /*
   *  ------------------------------ LOOP OVER CELL -------------------------------------------------
   *  Loop over the number of Cells in this domain on this processor
   *  This loop is done from left to right.
   */
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
      cIndex_cc_ = iCell;

    
      /*
       *  ---------------- Get the index for the center node ---------------------------------
       */
      index_CentLcNode = Index_DiagLcNode_LCO[iCell];
      /*
       *   Get the pointer to the NodalVars object for the center node
       */
      nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

      /*
       *  Index of the first equation in the bulk domain of center node
       */
      indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
	+ nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
      /*
       *  ------------------- Get the index for the left node -----------------------------
       *    There may not be a left node if we are on the left boundary. In that case
       *    set the pointer to zero and the index to -1. Hopefully, we will get a segfault on an error.
       */
      index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
      if (index_LeftLcNode < 0) {
	/*
	 *  We assign node object to zero.
	 */
	nodeLeft = 0;
	/*
	 *  If there is no left node, we assign the left solution index to the center solution index
	 */
      } else {
	// get the node structure for the left node
	nodeLeft = LI_ptr_->NodalVars_LcNode[index_LeftLcNode];
      }
     
      /*
       * ------------------------ Get the indexes for the right node ------------------------------------
       */
      index_RightLcNode = Index_RightLcNode_LCO[iCell];
      if (index_RightLcNode < 0) {
	nodeRight = 0;
      } else {
	//NodalVars
	nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
      }

      /*
       * --------------------------- CALCULATE POSITION AND DELTA_X Variables -----------------------------
       * Calculate the distance between the left and center node points
       */
      if (nodeLeft) {
	xCellBoundaryL = 0.5 * (nodeLeft->xNodePos() + nodeCent->xNodePos());
      } else {
	xCellBoundaryL = nodeCent->xNodePos();
      }
      /*
       * Calculate the distance between the right node and center node points
       */
      if (nodeRight == 0) {
	xCellBoundaryR = nodeCent->xNodePos();
      } else {
	xCellBoundaryR = 0.5 * (nodeRight->xNodePos() + nodeCent->xNodePos());
      }
      /*
       * Calculate the cell width
       */
  
      xdelCell_Cell_[iCell] = xCellBoundaryR - xCellBoundaryL;


      for (int k = 0; k < nsp_; k++) {
	Xcent_cc_[k] = soln[indexCent_EqnStart_BD + iVar_Species_BD + k];
      }
      Vcent_cc_ = soln[indexCent_EqnStart_BD + iVar_Voltage_BD];
      VElectrodeCent_cc_ = soln[indexCent_EqnStart_BD + iVar_Voltage_BD + 1];

      /*
       * Setup the thermo
       */ 
      //if (solnDot_ptr) {
//	solnDot_Cellptr =  &((*solnDot_ptr)[indexCent_EqnStart_BD]);
  //    } else {
//	solnDot_Cellptr = 0;
 //     } 
      //SetupThermoShop1Old(&(soln[indexCent_EqnStart_BD]), solnDot_Cellptr, 0);
      SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart_BD]));
   
      /*
       *  Calculate the electrode reactions.  Also update porosity.
       */
      int numSubcycles = calcElectrode();
      maxElectrodeSubIntegrationSteps_ = std::max(maxElectrodeSubIntegrationSteps_, numSubcycles);
  }
 
}
//=====================================================================================================================
int
porousLiKCl_FeS2Cathode_dom1D::calcElectrode()
{
  Electrode *Electrode_ptr = Electrode_Cell_[cIndex_cc_];

  double deltaT = t_final_ - t_init_;
  if (deltaT == 0.0) {
    deltaT = 1.0E-8;
  }
  Electrode_ptr->updateState();
  /*
   *  Integrate the subgrid Model
   */
  int numSubcycles = Electrode_ptr->integrate(deltaT);
#ifdef DEBUG_HKM
  if (ProblemResidEval::s_printFlagEnv > 0 && BDD_.Residual_printLvl_ > 8) {
      if (numSubcycles > 15) {
          printf("      cathode::calcElectrode(): problematic integration (%d) : Cell = %d,  counterNumberIntegrations_ = %d\n", numSubcycles,
	             cIndex_cc_, Electrode_ptr->counterNumberIntegrations_);
      }
  }
#endif

  /*
   *  Get the kmol 's produced by the electrode reactions. Units = kmol
   */
  Electrode_ptr->integratedSpeciesSourceTerm(&(electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * cIndex_cc_]));

  /*
   * Download the surface area from the electrode object, This is an
   * extrinsic object with units of m2
   */
  double sa[10];
  Electrode_ptr->getSurfaceAreas(sa);
  double saExternal = 0.0;
  const std::vector<bool> &isExternalSurface =  Electrode_ptr->getExternalSurfaceBooleans();
  for (int ph = 0; ph < nSurfsElectrode_; ph++) {
    if (isExternalSurface[ph]) {
      saExternal += sa[ph];
    }
  }

  int kelectron = Electrode_ptr->kSpecElectron();

  /*
   *  Calculate the total extrinsic current coming from the extrinsic electrode object
   *   -> this is in amps. We make the assumption that the current is constant throught
   *      the time interval, so we divide by the deltaT of the time.
   *      We divide by the cross sectional area of the electrode.
   */
  icurrInterface_Cell_[cIndex_cc_] = Faraday * electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * cIndex_cc_ + kelectron]
    / (deltaT * electrodeCrossSectionalArea_);

  /*
   * Get the phase mole flux from the Electrode object  units kmol/s
   */ 
  Electrode_ptr->getIntegratedPhaseMoleTransfer(&phaseMoleFlux_[0]);
  /*
   *  Calculate and store the electrolye molar flux units = kmol/m2/s
   */
  int sf = Electrode_ptr->solnPhaseIndex();
  solnMoleFluxInterface_Cell_[cIndex_cc_] = phaseMoleFlux_[sf] / electrodeCrossSectionalArea_;


  /*
   * What comes out of the surface electrode object is the total surface area in the electrode.
   * Remember the electrode is an extrinsic object that should be based on a m**2 cross-sectional area basis.
   * To get the surface area per m**2 we need to divide by the cross sectional area of the electrode
   */
  surfaceArea_Cell_[cIndex_cc_] = saExternal / (electrodeCrossSectionalArea_);

  /*
   *  Calculate the current per external surface area of electrode. To do this calculation
   *  we divide the total extrinsic current by the total extrinsic surface area
   */
  icurrInterfacePerSurfaceArea_Cell_[cIndex_cc_] = icurrInterface_Cell_[cIndex_cc_]  * electrodeCrossSectionalArea_ / saExternal;

  /* 
   * Compute new porosity based on new volume
   */
  //  porosity_Cell_[cIndex_cc_] = 1 -  Electrode_ptr->SolidVol() / (xdelCell_Cell_[cIndex_cc_] * electrodeCrossSectionalArea_);
  //  porosity_Cell_[cIndex_cc_] = 0.64007;

  if (residType_Curr_ == Base_ShowSolution) {
    capacityDischargedPA_Cell_[cIndex_cc_] = Electrode_ptr->capacityDischarged() / electrodeCrossSectionalArea_;
    depthOfDischargePA_Cell_[cIndex_cc_] = Electrode_ptr->depthOfDischarge() / electrodeCrossSectionalArea_;
    capacityLeftPA_Cell_[cIndex_cc_] = Electrode_ptr->capacityLeft() / electrodeCrossSectionalArea_;
    capacityPA_Cell_[cIndex_cc_]= Electrode_ptr->capacity() / electrodeCrossSectionalArea_;
  }
  return numSubcycles;
}

//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::SetupThermoShop1(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
{
    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    updateElectrolyte(nv, solnElectrolyte_Curr);
    updateElectrode();
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::SetupThermoShop2(const doublereal * const solnElectrolyte_CurrL,
                                                const doublereal * const solnElectrolyte_CurrR,
                                                int type)
{
  for (int i = 0; i < BDD_.NumEquationsPerNode; i++) {
    solnTemp[i] = 0.5 * (solnElectrolyte_CurrL[i] + solnElectrolyte_CurrR[i]);
  }
  if (type == 0) {
    porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ - 1] + porosity_Cell_[cIndex_cc_]);
  } else {
    porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ + 1] + porosity_Cell_[cIndex_cc_]);
  }
  updateElectrolyte(&solnTemp[0], 0);
  updateElectrode();
}
//=====================================================================================================================
// Function updates the ThermoPhase object for the electrolyte
// given the solution vector
/*
 *   Routine will update the molten salt ThermoPhase object with the current state of the electrolyte
 *
 * @param solnElectrolyte
 */
void
porousLiKCl_FeS2Cathode_dom1D::updateElectrolyte(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
{
    /*
     * Get the temperature: Check to see if the temperature is in the solution vector.
     *   If it is not, then use the reference temperature
     */
    temp_Curr_ = getPointTemperature(nv, solnElectrolyte_Curr);
    /*
     * Get the pressure
     */
    pres_Curr_ = getPointPressure(nv, solnElectrolyte_Curr);

    getMFElectrolyte_soln(nv, solnElectrolyte_Curr);
    getVoltages(nv, solnElectrolyte_Curr);

    ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);
    ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);
    /*
     * Get the partial molar volumes
     */
    ionicLiquid_->getPartialMolarVolumes(DATA_PTR(pmVolElectrolyte_Curr_));

    // Calculate the total concentration of the electrolyte kmol m-3.
    concTot_Curr_ = ionicLiquid_->molarDensity();
}
//=====================================================================================================================
  // Function updates the ThermoPhase object for the electrolyte given the solution vector and solution dot vector
  /*
   *   Routine will update the molten salt ThermoPhase object with the current state of the electrolyte.
   *   It also fills up the _Curr_ fields in this object refering to the Electrolyte properties.
   *
   *  @param solnElectrolyte      Vector of the solution at the current cell and bulk domain
   *  @param solnDotElectrolyte   Vector of the solution dot at the current cell and bulk domain
   */
void
porousLiKCl_FeS2Cathode_dom1D::updateElectrolyte(const doublereal * const solnElectrolyte_Curr,
						 const doublereal * const solnDotElectrolyte_Curr)
{
  /*
   * Get the temperature: Check to see if the temperature is in the solution vector.
   *   If it is not, then use the reference temperature
   */
  temp_Curr_ = TemperatureReference_;
  int iTemp = BDD_.VariableIndexStart_VarName[Temperature];
  if (iTemp >= 0) {
    temp_Curr_ = solnElectrolyte_Curr[iTemp];
  }
  /*
   * Get the pressure
   */
  pres_Curr_ = PressureReference_;
  /*
   * Fill up _Curr_ mf vectors for the electrolyte
   */
  getMFElectrolyte_soln(solnElectrolyte_Curr, solnDotElectrolyte_Curr);
  /*
   * Fill up _Curr_ Voltage values for the electrolyte
   */
  getVoltages(solnElectrolyte_Curr);
  /*
   *  Set the ThermoPhase object with the current solution values
   */
  ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);
  ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);
  /*
   * Get the partial molar volumes
   */
  ionicLiquid_->getPartialMolarVolumes(DATA_PTR(pmVolElectrolyte_Curr_));

  // Calculate the total concentration of the electrolyte kmol m-3.
  concTot_Curr_ = ionicLiquid_->molarDensity();

  // Calculate the time derivative of the concentration of the electrolyte
  concTotDot_Curr_ = 0.0;
  if (solnDotElectrolyte_Curr) {
  for (int k = 0; k < 3; k++) {
    concTotDot_Curr_ += pmVolElectrolyte_Curr_[k] * mfElectrolyte_SolnDot_Curr_[k];
  }
  concTotDot_Curr_ *= (- concTot_Curr_ * concTot_Curr_);
  }
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::updateElectrode()
{ 
  Electrode *Electrode_ptr = Electrode_Cell_[cIndex_cc_];
  /*
   * set the properties in the Electrode object
   *  -> temperature and pressure
   *  -> voltages of the phases
   */
  Electrode_ptr->setState_TP(temp_Curr_, pres_Curr_);
  Electrode_ptr->setVoltages(phiElectrode_Curr_, phiElectrolyte_Curr_);
  /*
   * HKM -> We might should make this a true mole number. Right now we are inputting mole fractions. The
   *        only thing the Electrode object needs from the electrolyte is the thermo ... probably true now ..
   *        maybe true in the future ...
   */
  Electrode_ptr->setElectrolyteMoleNumbers(&(mfElectrolyte_Thermo_Curr_[0]), false);
  /*
   * Set the internal objects within the electrode
   */
  Electrode_ptr->updateState();
}
//=====================================================================================================================
// Retrieves the voltages from the solution vector and puts them into local storage
/*
 * @param solnElectrolyte start of the solution vector at the current node
 */
void
porousLiKCl_FeS2Cathode_dom1D::getVoltages(const double * const solnElectrolyte_Curr)
{
  int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
  phiElectrolyte_Curr_ = solnElectrolyte_Curr[indexVS];
  phiElectrode_Curr_ = solnElectrolyte_Curr[indexVS + 1];
}
//=====================================================================================================================
// Retrieves the voltages from the solution vector and puts them into local storage
/*
 * @param solnElectrolyte start of the solution vector at the current node
 */
void
porousLiKCl_FeS2Cathode_dom1D::getVoltages(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    size_t indexVS = nv->indexBulkDomainVar0(Voltage);
    phiElectrolyte_Curr_ = solnElectrolyte_Curr[indexVS];
    phiElectrode_Curr_ = solnElectrolyte_Curr[indexVS + 1];
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::getMFElectrolyte_soln(const double * const solnElectrolyte_Curr,
						     const double * const solnDotElectrolyte_Curr)
{
  int indexMF = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
  mfElectrolyte_Soln_Curr_[0] = solnElectrolyte_Curr[indexMF];
  mfElectrolyte_Soln_Curr_[1] = solnElectrolyte_Curr[indexMF + 1];
  mfElectrolyte_Soln_Curr_[2] = solnElectrolyte_Curr[indexMF + 2];
  double mf0 = std::max(mfElectrolyte_Soln_Curr_[0], 0.0);
  double mf1 = std::max(mfElectrolyte_Soln_Curr_[1], 0.0);
  double tmp = mf0 + mf1;

  mfElectrolyte_Thermo_Curr_[0] = (mf0) * 0.5 / tmp;
  mfElectrolyte_Thermo_Curr_[1] = (mf1) * 0.5 / tmp;
  mfElectrolyte_Thermo_Curr_[2] = 0.5;

  if (solnDotElectrolyte_Curr) {
    mfElectrolyte_SolnDot_Curr_[0] = solnDotElectrolyte_Curr[indexMF];
    mfElectrolyte_SolnDot_Curr_[1] = solnDotElectrolyte_Curr[indexMF + 1];
    mfElectrolyte_SolnDot_Curr_[2] = solnDotElectrolyte_Curr[indexMF + 2];
  } else {
    mfElectrolyte_SolnDot_Curr_[0] = 0.0;
    mfElectrolyte_SolnDot_Curr_[1] = 0.0;
    mfElectrolyte_SolnDot_Curr_[2] = 0.0;
  }
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    size_t indexMF = nv->indexBulkDomainVar0(MoleFraction_Species);

    mfElectrolyte_Soln_Curr_[0] = solnElectrolyte_Curr[indexMF];
    mfElectrolyte_Soln_Curr_[1] = solnElectrolyte_Curr[indexMF + 1];
    mfElectrolyte_Soln_Curr_[2] = solnElectrolyte_Curr[indexMF + 2];
    double mf0  = std::max(mfElectrolyte_Soln_Curr_[0], 0.0);
    double mf1 = std::max(mfElectrolyte_Soln_Curr_[1], 0.0);
    double tmp = mf0 + mf1;
 
    mfElectrolyte_Thermo_Curr_[0] = (mf0) * 0.5 / tmp;
    mfElectrolyte_Thermo_Curr_[1] = (mf1) * 0.5 / tmp;
    mfElectrolyte_Thermo_Curr_[2] = 0.5;
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::getElectrolyte_SolnDot(const NodalVars* const nv,
			        		      const double * const solnDotElectrolyte_Curr)
{
    size_t indexMF = nv->indexBulkDomainVar0(MoleFraction_Species);
    if (solnDotElectrolyte_Curr) {
	mfElectrolyte_SolnDot_Curr_[0] = solnDotElectrolyte_Curr[indexMF];
	mfElectrolyte_SolnDot_Curr_[1] = solnDotElectrolyte_Curr[indexMF + 1];
	mfElectrolyte_SolnDot_Curr_[2] = solnDotElectrolyte_Curr[indexMF + 2];
    } else {
	mfElectrolyte_SolnDot_Curr_[0] = 0.0;
	mfElectrolyte_SolnDot_Curr_[1] = 0.0;
	mfElectrolyte_SolnDot_Curr_[2] = 0.0;
    }
 
    
    // Calculate the time derivative of the concentration of the electrolyte
    concTotDot_Curr_ = 0.0;
    if (solnDotElectrolyte_Curr) {
	for (int k = 0; k < 3; k++) {
	    concTotDot_Curr_ += pmVolElectrolyte_Curr_[k] * mfElectrolyte_SolnDot_Curr_[k];
	}
	concTotDot_Curr_ *= (- concTot_Curr_ * concTot_Curr_);
    }
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::SetupTranShop(const double xdel, const int type)
{

  /*
   * Determine diffusion velocities
   */

  //set gradients
  gradT_trCurr_ = 0.0;

  if (type == 0) {
    // Left boundary
    gradV_trCurr_ = (Vcent_cc_ - Vleft_cc_) / xdel;
    gradVElectrode_trCurr_ = (VElectrodeCent_cc_ - VElectrodeLeft_cc_) / xdel;

  } else {
    // Right boundary
    gradV_trCurr_ = (Vright_cc_ - Vcent_cc_) / xdel;
    gradVElectrode_trCurr_ = (VElectrodeRight_cc_ - VElectrodeCent_cc_) / xdel;
  }

  if (type == 0) {
    for (int k = 0; k < nsp_; k++) {
      gradX_trCurr_[k] = (Xcent_cc_[k] - Xleft_cc_[k]) / xdel;
    }
  } else {
    for (int k = 0; k < nsp_; k++) {
      gradX_trCurr_[k] = (Xright_cc_[k] - Xcent_cc_[k]) / xdel;
    }
  }

  trans_->getSpeciesVdiffES(1, &gradT_trCurr_, nsp_, &gradX_trCurr_[0], nsp_, &gradV_trCurr_, &Vdiff_trCurr_[0]);

  //  Correct species diffusivities according to the tortuosity
  double bruggemannExp = 2.5;
  Tortuosity tort(bruggemannExp );
  for (int k = 0; k < nsp_; k++) {
    Vdiff_trCurr_[k] *= tort.tortuosityFactor( porosity_Curr_ );
  }
  //Convert from diffusion velocity to diffusion flux
  for (int k = 0; k < nsp_; k++) {
    jFlux_trCurr_[k] = mfElectrolyte_Soln_Curr_[k] * concTot_Curr_ * Vdiff_trCurr_[k];
  }

  double volFSolid = 1.0 - porosity_Curr_;
  icurrElectrode_trCurr_ = -conductivityElectrode_ * pow(volFSolid, 1.5) * gradVElectrode_trCurr_;
}
//=====================================================================================================================
// saving the solution on the domain in an xml node.
/*
 *
 * @param oNode                Reference to the XML_Node
 * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
 * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
 * @param t                    time
 *
 * @param duplicateOnAllProcs  If this is true, all processors will include
 *                             the same XML_Node information as proc 0. If
 *                             false, the xml_node info will only exist on proc 0.
 */
void
porousLiKCl_FeS2Cathode_dom1D::saveDomain(Cantera::XML_Node& oNode,
                                          const Epetra_Vector *soln_GLALL_ptr,
                                          const Epetra_Vector *solnDot_GLALL_ptr,
                                          const double t,
                                          bool duplicateOnAllProcs)
{
  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;

  // Add a child for this domain
  Cantera::XML_Node& bdom = oNode.addChild("domain");

  // Number of equations per node
  int numEquationsPerNode = BDD_.NumEquationsPerNode;

  // Vector containing the variable names as they appear in the solution vector
  std::vector<VarType> &variableNameList = BDD_.VariableNameList;

  //! First global node of this bulk domain
  int firstGbNode = BDD_.FirstGbNode;

  //! Last Global node of this bulk domain
  int lastGbNode = BDD_.LastGbNode;
  int numNodes = lastGbNode - firstGbNode + 1;

  bdom.addAttribute("id", id());
  bdom.addAttribute("points", numNodes);
  bdom.addAttribute("type", "bulk");
  bdom.addAttribute("numVariables", numEquationsPerNode);

  // Dump out the coordinates
  Cantera::XML_Node& gv = bdom.addChild("grid_data");

  std::vector<double> varContig(numNodes);

  int i = 0;
  for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
    NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
    varContig[i] = nv->x0NodePos();
  }
  ctml::addNamedFloatArray(gv, "X0", varContig.size(), &(varContig[0]), "m", "length");

  for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
    VarType vt = variableNameList[iVar];
    i = 0;
    std::string nmm = vt.VariableName(200);
    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
      int istart = nv->EqnStart_GbEqnIndex;
      varContig[i] = (*soln_GLALL_ptr)[istart + ibulk + iVar];
    }
    ctml::addNamedFloatArray(gv, nmm, varContig.size(), &(varContig[0]), "kmol/m3", "concentration");

  }
}
//=====================================================================================================================
static void
drawline(int sp, int ll)
{
  for (int i = 0; i < sp; i++) {
    Cantera::writelog(" ");
  }
  for (int i = 0; i < ll; i++) {
    Cantera::writelog("-");
  }
  Cantera::writelog("\n");
}
//=====================================================================================================================
static void
drawline0(int sp, int ll)
{
  for (int i = 0; i < sp; i++) {
    sprint0(" ");
  }
  for (int i = 0; i < ll; i++) {
    sprint0("-");
  }
  sprint0("\n");
}




//=====================================================================================================================
// Method for writing the header for the surface domain to a tecplot file.
void
porousLiKCl_FeS2Cathode_dom1D::writeSolutionTecplotHeader()
{
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid ; //only proc 0 should write

  if (doWrite) {
    
    //open tecplot file
    FILE* ofp;
    string sss = id();
    char filename[20];
    sprintf(filename,"%s%s",sss.c_str(),".dat");
    ofp = fopen( filename, "w");

    //write title and variable list
    fprintf( ofp, "TITLE = \"Solution on Domain %s\"\n",sss.c_str() );

    // Number of equations per node
    int numVar = BDD_.NumEquationsPerNode;
    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType> &variableNameList = BDD_.VariableNameList;
    //! First global node of this bulk domain

    fprintf(ofp, "VARIABLES = ");
    fprintf(ofp, "\"x [m]\"  \n" );

    for (int k = 0; k < numVar; k++) {
      VarType &vt = variableNameList[k];
      string name = vt.VariableName(15);
      fprintf(ofp, "\"%s\" \t", name.c_str() );
    }
    fprintf(ofp, "\n" );

    /////////////////////////////////////////
    //// end BulkDomain1D section
    //// start application specific output
    /////////////////////////////////////////


    // print potentials within each control volume
    fprintf(ofp, "\"Potential Drop Electrode-Electrolyte [V]\" \t" );
    for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++ ) {
      fprintf(ofp, "\"Surf_%d Equilibrium Potential Drop Electrode-Electrolyte [V]\" \t", jSurf);
      fprintf(ofp, "\"Overpotential_%d [V]\" \t", jSurf );
    }

    //  Print depth of discharge for each control volume
    fprintf(ofp, "\"Capacity discharged [mA-hr/g]\" \t" );
    fprintf(ofp, "\"Depth of discharge [-]\" \t" );
    fprintf(ofp, "\"Capacity remaining [mA-hr/g]\" \t" );
    fprintf(ofp, "\"Initial capacity [mA-hr/g]\" \t" );

    // print porosity, specific surface area, thickness for each control volume
    fprintf(ofp, "\"Porosity []\" \t" );
    fprintf(ofp, "\"Surface Area [] (area per m2)\" \t" );
    fprintf(ofp, "\"Specific current (per particle area) [A/m^2]\" \t" );
    fprintf(ofp, "\"Control volume thickness [m]\" \t" );
    fprintf( ofp, "\n" );

    //set up Electrode objects so that we can plot moles of active materials
    Electrode *ee0 = Electrode_Cell_[0];
    std::vector<double> spmoles(nSpeciesElectrode_, 0.0);
    int solnPhase = ee0->solnPhaseIndex();
    int metalPhase =  ee0->metalPhaseIndex();
    int numVolPhasesE =  ee0->nVolPhases();

    // print mole numbers of active materials
    for (int vph = 0; vph < numVolPhasesE; vph++) {
      ThermoPhase *tp = &(ee0->volPhase(vph));
      int iph = ee0->getGlobalPhaseIndex(tp);
      if (iph == metalPhase || iph == solnPhase) {
	
      } else {
	int nspPhase = tp->nSpecies();
	int kStart =  ee0->getGlobalSpeciesIndex(iph, 0);
	for (int k = 0; k < nspPhase; k++) {
	  string sss = ee0->speciesName(kStart + k);
	  fprintf(ofp, "\"Moles %s [mol/m^2] (per electrode area)\" \t", sss.c_str() );
	}
      }
    }
    fprintf( ofp, "\n" );

    fprintf(ofp, "\n");
    fclose(ofp);
  }

}
//=====================================================================================================================
// Method for writing the solution on the surface domain to a tecplot file.
/*
 *
 * @param soln__GLALL_ptr      Pointer to the Global-All solution vector
 * @param solnDot_ptr          Pointer to the time derivative of the Global-All solution vector
 * @param t                    time
 *
 * @param duplicateOnAllProcs  If this is true, all processors will include
 *                             the same log information as proc 0. If
 *                             false, the loginfo will only exist on proc 0.
 */
void
porousLiKCl_FeS2Cathode_dom1D::writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr,
						    const Epetra_Vector *solnDot_GlAll_ptr,
						    const double t)
{
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid ; //only proc 0 should write

  if (doWrite) {
    
    // get the NodeVars object pertaining to this global node
    GlobalIndices *gi = LI_ptr_->GI_ptr_;
    // Number of equations per node
    int numVar = BDD_.NumEquationsPerNode;
    //! First global node of this bulk domain
    int firstGbNode = BDD_.FirstGbNode;
    //! Last Global node of this bulk domain
    int lastGbNode = BDD_.LastGbNode;
    int numNodes = lastGbNode - firstGbNode + 1;
    
    
    //open tecplot file
    FILE* ofp;
    string sss = id();
    char filename[20];
    sprintf(filename,"%s%s",sss.c_str(),".dat");
    ofp = fopen( filename, "a");
    
    fprintf(ofp, "ZONE T = \"t = %g [s]\" I = %d SOLUTIONTIME = %19.13E\n", t, numNodes, t);

    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];

      //x-position
      fprintf(ofp, "%g \t", nv->xNodePos() );

      for (int iVar = 0; iVar < numVar; iVar++) {

	//other variables
	int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
	int istart = nv->EqnStart_GbEqnIndex;
	fprintf(ofp, "%g \t", (*soln_GlAll_ptr)[istart + ibulk + iVar] );
      }
      fprintf(ofp, "\n");

    /////////////////////////////////////////
    //// end BulkDomain1D section
    /////////////////////////////////////////

      int iCell = iGbNode - firstGbNode;

      fprintf(ofp, "%g \t", icurrRxn_Cell_[iCell]);
      fprintf(ofp, "%g \t", LiFlux_Cell_[iCell]);

      // current flux -- average of left and right fluxes
      fprintf(ofp, "%g \t", 0.5 * ( icurrElectrolyte_CBL_[iCell] + icurrElectrolyte_CBR_[iCell] ) );
      fprintf(ofp, "%g \t", 0.5 * ( icurrElectrode_CBL_[iCell] + icurrElectrode_CBR_[iCell] ) );

      fprintf(ofp, "\n");

      // surface reaction data
      fprintf(ofp, "%g \t", deltaV_Cell_[iCell]);
      for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++ ) {
	fprintf(ofp, "%g \t", Ess_Cell_[nSurfsElectrode_ * iCell + jSurf]);
	fprintf(ofp, "%g \t", overpotential_Cell_[nSurfsElectrode_ * iCell + jSurf]);
      }
      fprintf(ofp, "\n");

      //access Electrode for each CV
      Electrode *ee = Electrode_Cell_[iCell];
      //first we need the FeS2 mass per CV to compare with capacity variables
      double MW_FeS2 = 119.98 * 1e3; //molecular mass in g/kmol
      const  Elements *elem = ee->getGlobalElements();
      int iFe = elem->elementIndex("Fe");
      double mass_FeS2_CV = MW_FeS2 * ee->elementMoles( iFe);
      double ampSecTomAmpHr = 3.6;  //convert from [Amp-sec] to [mAmp-hr]

      // capacity of individual control volumes
      // note that the variables like capacityDischarged_Cell_[iCell] are stored in units of A*s/area
      fprintf(ofp, "%g \t", capacityDischargedPA_Cell_[iCell] * electrodeCrossSectionalArea_ / mass_FeS2_CV / ampSecTomAmpHr );
      fprintf(ofp, "%g \t", depthOfDischargePA_Cell_[iCell]);
      fprintf(ofp, "%g \t", capacityLeftPA_Cell_[iCell] * electrodeCrossSectionalArea_ / mass_FeS2_CV / ampSecTomAmpHr );
      fprintf(ofp, "%g \t", capacityPA_Cell_[iCell] * electrodeCrossSectionalArea_ / mass_FeS2_CV / ampSecTomAmpHr );

      // print porosity, surface area, current source per surface area, thickness for each control volume
      fprintf(ofp, "%g \t", porosity_Cell_[iCell]);
      fprintf(ofp, "%g \t", surfaceArea_Cell_[iCell]);
      fprintf(ofp, "%g \t", icurrInterfacePerSurfaceArea_Cell_[iCell]);
      fprintf(ofp, "%g \t", xdelCell_Cell_[iCell]);
      fprintf( ofp, "\n" );

      // print the properties of the electrode particles
      std::vector<double> spmoles(nSpeciesElectrode_, 0.0);
      int solnPhase = ee->solnPhaseIndex();
      int metalPhase =  ee->metalPhaseIndex();
      int numVolPhasesE =  ee->nVolPhases();

      // print mole numbers of active materials
      ee->getMoleNumSpecies(DATA_PTR(spmoles));
      for (int vph = 0; vph < numVolPhasesE; vph++) {
	ThermoPhase *tp = &(ee->volPhase(vph));
	int iph = ee->getGlobalPhaseIndex(tp);
	if (iph == metalPhase || iph == solnPhase) {
	  
	} else {
	  int nspPhase = tp->nSpecies();
	  int kStart =  ee->getGlobalSpeciesIndex(iph, 0);
	  for (int k = 0; k < nspPhase; k++) {
	    fprintf(ofp, "%g \t", spmoles[kStart + k] / electrodeCrossSectionalArea_);
	  }
	}
      }
      fprintf( ofp, "\n" );

    }
    fclose(ofp);

#undef DAKOTAOUT
#ifdef  DAKOTAOUT
    double firstOutputTime = 0.1;
    static double outputTime = 0.1;
    if ( t > outputTime ) {
      std::ofstream dfp;
      if ( outputTime > firstOutputTime )
	dfp.open( "results.out", std::ios_base::app );
      else
	dfp.open( "results.out", std::ios_base::out );
      int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
      double potentialE_left = (*soln_GlAll_ptr)[indexVS];
      double lithium_left = (*soln_GlAll_ptr)[1];
      
      NodalVars *nv = gi->NodalVars_GbNode[lastGbNode-1];
      int istart = nv->EqnStart_GbEqnIndex;
      double potentialE_right = (*soln_GlAll_ptr)[istart+indexVS];
      double lithium_right = (*soln_GlAll_ptr)[istart+1];
      double ohmicLoss =  potentialE_left - potentialE_right;
      double moleFracDelta = lithium_left - lithium_right;

      // uncomment next line to get ohmic loss in electrolyte
      //dfp << ohmicLoss << "  t" << outputTime << std::endl;

      // uncomment next line to get change in litium mole fraction
      //dfp << moleFracDelta << "  t" << outputTime << std::endl;
      
      double Li_eutectic = 0.585/2.0;
      double Rgas_local = 8.3144;
      double Faraday_local = 9.6485e4;
      double concentrationPotential = Rgas_local * TemperatureReference_ / Faraday_local
	* log( lithium_left / lithium_right ) ;

      //uncomment following line to get concentration potential loss
      dfp << concentrationPotential << "  t" << outputTime << std::endl;

      //uncomment following line to get combined ohmic and concentration losses
      outputTime *= 10.0;
    }
      
#endif
#undef DAKOTAOUT
  }
}



//=====================================================================================================================
// Base class for writing the solution on the domain to a logfile.
/*
 *
 * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
 * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
 * @param soln_ptr             Pointer to the solution vector
 * @param solnDot_ptr          Pointer to the time derivative of the solution vector
 * @param solnOld_ptr          Pointer to the solution vector at the old time step
 * @param residInternal _ptr   Pointer to the current value of the residual just calculated
 *                             by a special call to the residEval()
 * @param t                    time
 * @param rdelta_t             The inverse of the value of delta_t
 * @param indentSpaces         Indentation that all output should have as a starter
 * @param duplicateOnAllProcs  If this is true, all processors will include
 *                             the same log information as proc 0. If
 *                             false, the loginfo will only exist on proc 0.
 */
void
porousLiKCl_FeS2Cathode_dom1D::showSolution(const Epetra_Vector *soln_GlAll_ptr,
                                            const Epetra_Vector *solnDot_GlAll_ptr,
                                            const Epetra_Vector *soln_ptr,
                                            const Epetra_Vector *solnDot_ptr,
                                            const Epetra_Vector *solnOld_ptr,
                                            const Epetra_Vector_Owned *residInternal_ptr,
                                            const double t,
                                            const double rdelta_t,
                                            int indentSpaces,
                                            bool duplicateOnAllProcs)
{

  // nn is the number of block rows in the printout
  int nn = NumDomainEqns / 5;
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool do0Write = (!mypid || duplicateOnAllProcs);
  //bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;

  std::vector<VarType> &variableNameList = BDD_.VariableNameList;
  int iBlock;
  int iGbNode;
  int n;
  int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;

  /*
   * Do some things that should be done after each converged solution (time step)
   * We should put these into another method called from ProblemResidEval::user_out()
   * or similar that is guaranteed to be called after converges solutions.
   */  
  if (mypid == 0) {
    depthOfDischarge_Electrode_ = 0.0;
    capacityZeroDoD_Electrode_  = 0.0;
    for (int iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
      int iCell = iGbNode - BDD_.FirstGbNode;
      depthOfDischarge_Electrode_ += depthOfDischargePA_Cell_[iCell];
      capacityZeroDoD_Electrode_ += capacityPA_Cell_[iCell];
    }
  }
  /*
   * Finished with post converged solution
   */


  std::string indent = "";
  for (int i = 0; i < indentSpaces; i++) {
    indent += " ";
  }
  const char *ind1 = indent.c_str();
  char ind[120];
  strcpy(ind, ind1);
  ind[118] = '\0';
  doublereal v;
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  // Number of points in each vector
  string sss = id();
  stream0 ss;

  if (do0Write) {
    drawline(indentSpaces, 80);
    ss.print0("%s  Solution on Bulk Domain %12s : Number of variables = %d\n", ind, sss.c_str(), NumDomainEqns);
    ss.print0("%s                                         : Number of Nodes = %d\n", ind, nPoints);
    ss.print0("%s                                         : Beginning pos %g\n", ind, BDD_.Xpos_start);
    ss.print0("%s                                         : Ending    pos %g\n", ind, BDD_.Xpos_end);
  }
  if (do0Write) {
    for (iBlock = 0; iBlock < nn; iBlock++) {
      drawline(indentSpaces, 80);
      ss.print0("%s        z   ", ind);
      for (n = 0; n < 5; n++) {
        int ivar = iBlock * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        ss.print0(" %15s", name.c_str());
      }
      ss.print0("\n");
      drawline(indentSpaces, 80);

      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
        NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
        doublereal x = nv->xNodePos();
        ss.print0("\n%s    %-10.4E ", ind, x);
        int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
        int istart = nv->EqnStart_GbEqnIndex;
        for (n = 0; n < 5; n++) {
          int ivar = iBlock * 5 + n;
          VarType vt = variableNameList[ivar];
          v = (*soln_GlAll_ptr)[istart + ibulk + iBlock * 5 + n];
          ss.print0(" %-10.4E ", v);
        }
      }
      ss.print0("\n");
    }

    int nrem = NumDomainEqns - 5 * nn;
    if (nrem > 0) {
      drawline0(indentSpaces, 80);
      ss.print0("%s        z   ", ind);
      for (n = 0; n < nrem; n++) {
        int ivar = nn * 5 + n;
        VarType vt = variableNameList[ivar];
        string name = vt.VariableName(15);
        ss.print0(" %15s", name.c_str());
      }
      ss.print0("\n");
      drawline(indentSpaces, 80);

      for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
        NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
        doublereal x = nv->xNodePos();
        ss.print0("%s    %-10.4E ", ind, x);
        int ibulk = nv->OffsetIndex_BulkDomainEqnStart_BDN[0];
        int istart = nv->EqnStart_GbEqnIndex;
        for (n = 0; n < nrem; n++) {
          int ivar = iBlock * 5 + n;
          VarType vt = variableNameList[ivar];
          v = (*soln_GlAll_ptr)[istart + ibulk + nn * 5 + n];
          ss.print0(" %-10.4E ", v);
        }
        ss.print0("\n");
      }
    }
    drawline(indentSpaces, 80);
  }

  int nLL = 66 + (nSurfsElectrode_ - 1) * 24;
  if (do0Write) {
    // -----------------------------------------------------------------------------------------------------------------
    // --             PRINT REACTION RATES within the cell --
    // -----------------------------------------------------------------------------------------------------------------
    ss.print0("\n");
    drawline0(indentSpaces, nLL);
    ss.print0("%s        z       Delta_V     ", ind);
    for (int isurf = 0; isurf < nSurfsElectrode_; isurf++) {
       ss.print0("   Ess[%1d]   Overpot[%1d]  ", isurf, isurf);
    }
    ss.print0("icurrRxnCell");
    ss.print0("\n");
    drawline0(indentSpaces, nLL);
  }
  doublereal x;
  int iCell;
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      ss.print0("%11.4E ", deltaV_Cell_[iCell]);
      for (int isurf = 0; isurf < nSurfsElectrode_; isurf++) {
        ss.print0("%11.4E ", Ess_Cell_[nSurfsElectrode_ * (iCell) + isurf]);
        ss.print0("%11.4E ", overpotential_Cell_[nSurfsElectrode_ * (iCell) + isurf]);
      }
      ss.print0("%11.4E ", icurrRxn_Cell_[iCell]);
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }
  drawline0(indentSpaces, nLL);


  if (do0Write) {
    // -----------------------------------------------------------------------------------------------------------------
    // --             PRINT VOLUME DETAILS --
    // -----------------------------------------------------------------------------------------------------------------
    ss.print0("\n");
    drawline0(indentSpaces, 80);
    ss.print0("%s        z      Porosity", ind);
    ss.print0("\n");
    drawline0(indentSpaces, 80);
  }
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      ss.print0("%11.4E ", porosity_Cell_[iCell]);
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }

  if (do0Write) {
    // -----------------------------------------------------------------------------------------------------------------
    // --             PRINT SURFACE REACTION DETAILS ABOUT EACH CELL --
    // -----------------------------------------------------------------------------------------------------------------
    ss.print0("\n");
    drawline0(indentSpaces, 80);
    ss.print0("%s        z     SurfaceArea  SurfAreaDens  currPerSA  DeltaZ", ind);
    ss.print0("\n");
    drawline0(indentSpaces, 80);
  }
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      ss.print0("%11.4E ", surfaceArea_Cell_[iCell]);
      ss.print0("%11.4E ", surfaceArea_Cell_[iCell]/xdelCell_Cell_[iCell]);
      ss.print0("%11.4E ", icurrInterfacePerSurfaceArea_Cell_[iCell]);
      ss.print0("%11.4E ", xdelCell_Cell_[iCell]);
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }

  if (do0Write) {
    // -----------------------------------------------------------------------------------------------------------------
    // --             PRINT DEPTH OF DISCHARGE VALUES FOR EACH CELL --
    // -----------------------------------------------------------------------------------------------------------------
    ss.print0("\n");
    drawline0(indentSpaces, 80);
    ss.print0("%s        z    capDischarged DepthDischge  capLeft  capZeroDOD  ", ind);
    ss.print0("\n");
    drawline0(indentSpaces, 80);
  }
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      ss.print0("%11.4E ", capacityDischargedPA_Cell_[iCell]);
      ss.print0("%11.4E ", depthOfDischargePA_Cell_[iCell]);
      ss.print0("%11.4E ", capacityLeftPA_Cell_[iCell]);
      ss.print0("%11.4E ", capacityPA_Cell_[iCell]);
  
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }
   
  Electrode *ee0 = Electrode_Cell_[0];
  std::vector<double> spmoles(nSpeciesElectrode_, 0.0);
  int solnPhase = ee0->solnPhaseIndex();
  int metalPhase =  ee0->metalPhaseIndex();
  int numVolPhasesE =  ee0->nVolPhases();
  if (do0Write) {
    // -----------------------------------------------------------------------------------------------------------------
    // --             PRINT ELECTRODE SOLID MOLE NUMBERS --
    // -----------------------------------------------------------------------------------------------------------------
    ss.print0("\n");
    drawline0(indentSpaces, 120);
    ss.print0("%s        z       ", ind);

    for (int vph = 0; vph < numVolPhasesE; vph++) {
      ThermoPhase *tp = &(ee0->volPhase(vph));
      int iph = ee0->getGlobalPhaseIndex(tp);
      if (iph == metalPhase || iph == solnPhase) {
	
      } else {
	int nspPhase = tp->nSpecies();
	int kStart =  ee0->getGlobalSpeciesIndex(iph, 0);
	for (int k = 0; k < nspPhase; k++) {
	  string sss = ee0->speciesName(kStart + k);
	  ss.print0("%-15.15s ", sss.c_str());
	}
      }
    }
    ss.print0("\n");
    drawline0(indentSpaces, 120);
  }
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      iCell = iGbNode - BDD_.FirstGbNode;
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      x = nv->xNodePos();
      ss.print0("%s    %-10.4E ", ind, x);
      Electrode *ee = Electrode_Cell_[iCell];
      ee->getMoleNumSpecies(DATA_PTR(spmoles));
      for (int vph = 0; vph < numVolPhasesE; vph++) {
	ThermoPhase *tp = &(ee->volPhase(vph));
	int iph = ee->getGlobalPhaseIndex(tp);
	if (iph == metalPhase || iph == solnPhase) {

	} else {
	  int nspPhase = tp->nSpecies();
	  int kStart =  ee->getGlobalSpeciesIndex(iph, 0);
	  for (int k = 0; k < nspPhase; k++) {
	    ss.print0("% -15.6E ", spmoles[kStart + k] / electrodeCrossSectionalArea_);
	  }
	  
	}
      }
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }
  

  // -----------------------------------------------------------------------------------------------------------------
  // --             PRINT ELECTRODE MODEL TABLES
  // -----------------------------------------------------------------------------------------------------------------
  int nETables = ee0->getNumPrintTables();
  std::vector<std::string> colHeaders;
  std::vector<double> colValues;
  int ncols = 0;
  nLL = 66;
  for (int iETable = 0; iETable < nETables; iETable++) {
    if (do0Write) {
      ee0->getPrintTable(iETable, colHeaders, colValues);
      ncols = colHeaders.size();
      nLL = 40 + ncols * 18;
      ss.print0("\n");
      drawline0(indentSpaces, nLL);
      ss.print0("%s        z       ", ind);
   
      for (int icol = 0; icol < ncols; icol++) {
	ss.print0("%-15.15s ", colHeaders[icol].c_str());
      }
      ss.print0("\n");
      drawline0(indentSpaces, nLL);
  
    }
    for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
      print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
      if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
	iCell = iGbNode - BDD_.FirstGbNode;
	NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	x = nv->xNodePos();
	ss.print0("%s    %-10.4E ", ind, x);
	Electrode *ee = Electrode_Cell_[iCell];
	ee->getPrintTable(iETable, colHeaders, colValues);
	for (int icol = 0; icol < ncols; icol++) {
	  ss.print0("% -15.6E ", colValues[icol]);
	}
	ss.print0("\n");
      }

      print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    }
    drawline0(indentSpaces, nLL);
  }

  if (do0Write) {
    // ----------------------------------------------------
    // --             PRINT FLUXES AT THE CELL BOUNDARIES --
    // ----------------------------------------------------   
    ss.print0("\n");
    drawline0(indentSpaces, 80);
    ss.print0("%s    CellBound    z      IcurrElectrolyte  IcurrElectrode ", ind);
    ss.print0("\n");
    drawline(indentSpaces, 80);
  }
  NodalVars *nvl;
  NodalVars *nvr;
  for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
      ss.print0("%s    ", ind);
      if (iGbNode == BDD_.FirstGbNode) {

        iCell = 0;
        nvr = gi->NodalVars_GbNode[BDD_.FirstGbNode];
        x = nvr->xNodePos();
        ss.print0("Lft-0     %11.4E ", x);
        ss.print0("%11.4E ", icurrElectrolyte_CBL_[0]);
        ss.print0("%11.4E", icurrElectrode_CBL_[0]);
	ss.print0("\n");
	ss.print0("%s    ", ind);
	nvl = gi->NodalVars_GbNode[iGbNode];
        nvr = gi->NodalVars_GbNode[iGbNode + 1];
        x = 0.5 * (nvl->xNodePos() + nvr->xNodePos());
        ss.print0("%3d-%-3d   %11.4E ", iCell, iCell + 1, x);
        ss.print0("%11.4E ", icurrElectrolyte_CBR_[iCell]);
        ss.print0("%11.4E ", icurrElectrode_CBR_[iCell]);

      } else if (iGbNode == BDD_.LastGbNode) {
        iCell = BDD_.LastGbNode - BDD_.FirstGbNode;
        nvr = gi->NodalVars_GbNode[BDD_.LastGbNode];
        x = nvr->xNodePos();
        ss.print0("%3d-Rgt   %11.4E ", iCell, x);
        ss.print0("%11.4E ", icurrElectrolyte_CBR_[iCell]);
        ss.print0("%11.4E ", icurrElectrode_CBR_[iCell]);
      } else {
        iCell = iGbNode - BDD_.FirstGbNode;
        nvl = gi->NodalVars_GbNode[iGbNode];
        nvr = gi->NodalVars_GbNode[iGbNode + 1];
        x = 0.5 * (nvl->xNodePos() + nvr->xNodePos());
        ss.print0("%3d-%-3d   %11.4E ", iCell, iCell + 1, x);
        ss.print0("%11.4E ", icurrElectrolyte_CBR_[iCell]);
        ss.print0("%11.4E ", icurrElectrode_CBR_[iCell]);
      }
      ss.print0("\n");
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
  }

}
//=======================================================================================================================
// Set the underlying state of the system from the solution vector
/*
 *   Note this is an important routine for the speed of the solution.
 *   It would be great if we could supply just exactly what is changing here.
 *   This routine is always called at the beginning of the residual evaluation process.
 *
 *   This is a natural place to put any precalculations of nodal quantities that
 *   may be needed by the residual before its calculation.
 *
 *   Also, this routine is called with delta_t = 0. This implies that a step isn't being taken. However, the
 *   the initial conditions must be propagated.
 *
 *   Note, in general t may not be equal to t_old + delta_t. If this is the case, then the solution is
 *   interpolated across the time interval and then the solution applied.
 *
 *   If doTimeDependentResid then delta_t > 0. 
 *   If !doTimeDependentResid then usually delta_t = 0 but not necessarily
 *
 * @param doTimeDependentResid
 * @param soln
 * @param solnDot
 * @param t
 * @param delta_t delta t. If zero then delta_t equals 0.
 * @param t_old
 */
void
porousLiKCl_FeS2Cathode_dom1D::setStateFromSolution(const bool doTimeDependentResid, 
                                                    const Epetra_Vector_Ghosted *soln_ptr, 
						    const Epetra_Vector_Ghosted *solnDot_ptr,
						    const double t, const double delta_t, const double t_old)
{
    int indexCent_EqnStart_BD;
    bool doAll = true;
    //bool doInit = false;
    //bool doFinal = true; // This is always true as you can only set the final Electrode object
    //size_t  indexCent_EqnStart;
    if (doTimeDependentResid) {
	if (delta_t >= 0.0) {
            doAll = false;
	    if (t <= t_old) {
		//doInit = true;
	    }
	}
    }
    const Epetra_Vector_Ghosted& soln = *soln_ptr;

    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;
	//cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];


	// Find the index of the center node
	int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
	// pointer to the NodalVars object for the center node
	NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
	// Index of the first equation in the bulk domain of center node
	indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
				+ nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];

	/*
         *  ---------------- Get the index for the center node ---------------------------------
         *   Get the pointer to the NodalVars object for the center node
         *   Index of the first equation in the bulk domain of center node
         */
	const double *solnCentStart = &(soln[indexCent_EqnStart_BD]);
    
	

        //
        // update the electrolyte quantities 
        //   calculates       temp_Curr_
	//                    pres_Curr_
        //                    mfElectrolyte_Thermo_Curr_[]
	//                    phiElectrolyte_Curr_
        //
        if (solnDot_ptr) {
            const Epetra_Vector_Ghosted& solnDot = *solnDot_ptr;
            //updateElectrolyte(solnCentStart, &(solnDot[indexCent_EqnStart_BD]));
            updateElectrolyte(nodeCent, solnCentStart);
            getElectrolyte_SolnDot(nodeCent, &(solnDot[indexCent_EqnStart_BD]));
        } else {
           // updateElectrolyte(solnCentStart, 0);
            updateElectrolyte(nodeCent, solnCentStart);
        }
	//
	//                    phiElectrode_Curr_
	//
        getVoltages(solnCentStart);
	//
        // Electrode object will be updated and we will compute the porosity
	//
        Electrode* ee = Electrode_Cell_[iCell];
	//
	//  Set the temperature and pressure and voltages in the final_ state
	//
	ee->setState_TP(temp_Curr_, pres_Curr_);
        ee->setVoltages(phiElectrode_Curr_, phiElectrolyte_Curr_);
	//
	//  Set the electrolyte mole fractions
	//
	ee->setElectrolyteMoleNumbers(&(mfElectrolyte_Thermo_Curr_[0]), false);
	//
	// update the internals
	//
	ee->updateState();
	//
	// 
	//
	if (doAll) {
	    ee->setInitStateFromFinal(true);
	    ee->setFinalFinalStateFromFinal();
	}
    }
}
//====================================================================================================
// Generate the initial conditions
/*
 *   The basic algorithm is to loop over the volume domains.
 *   Then, we loop over the surface domains
 *
 * @param doTimeDependentResid    Boolean indicating whether we should
 *                                formulate the time dependent residual
 * @param soln                    Solution vector. This is the input to
 *                                the residual calculation.
 * @param solnDot                 Solution vector. This is the input to
 *                                the residual calculation.
 * @param t                       Time
 * @param delta_t                 delta_t for the initial time step
 */
void
porousLiKCl_FeS2Cathode_dom1D::initialConditions(const bool doTimeDependentResid,
                                                 Epetra_Vector *soln_ptr,
                                                 Epetra_Vector *solnDot,
                                                 const double t,
                                                 const double delta_t)
{
  Epetra_Vector &soln = *soln_ptr;
  
  int index_CentLcNode;
  NodalVars *nodeCent = 0;
  int indexCent_EqnStart_BD;
  
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;
    
    index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    // pointer to the NodalVars object for the center node
    nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
    // Index of the first equation in the bulk domain of center node
    indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
      + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 0.0;

    //soln[indexCent_EqnStart_BD + iVar_Species_BD + 0] = 0.3;
    //soln[indexCent_EqnStart_BD + iVar_Species_BD + 1] = 0.2;
    // HACK from HEWSON
    // WE ARE USING INPUT FROM THE SEPARATOR MOLE FRACTIONS HERE
    /*
     * Get initial mole fractions from PSinput
     */
    int iLiCl = PSinput.PhaseList_->globalSpeciesIndex ("LiCl(L)");
    int iKCl = PSinput.PhaseList_->globalSpeciesIndex ("KCl(L)");
    int iLi = PSinput.PhaseList_->globalSpeciesIndex ("Li+");
    int iK = PSinput.PhaseList_->globalSpeciesIndex ("K+");
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 0] =
      PSinput.electrolyteMoleFracs_[iLiCl] / 2.0 +
      PSinput.electrolyteMoleFracs_[iLi];
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 1] =
      PSinput.electrolyteMoleFracs_[iKCl] / 2.0 +
      PSinput.electrolyteMoleFracs_[iK];
    soln[indexCent_EqnStart_BD + iVar_Species_BD + 2] = 0.5;

    soln[indexCent_EqnStart_BD + iVar_Voltage_BD] = -0.07;

    //double icurr = PSinput.icurrDischargeSpecified_;

    // For now, pick a voltage near OCV for LiSi-FeS2 cell
    // This is mainly needed for the ee->setVoltages() call.
    // adjust below with specific BC.
    soln[indexCent_EqnStart_BD + iVar_Voltage_ED] = 1.9;

    // Electrode object will beupdated and we will compute the porosity
    Electrode * ee = Electrode_Cell_[iCell];

    double porosity = 1.0 - ee->SolidVol() / (xdelCell_Cell_[iCell] * electrodeCrossSectionalArea_);
    // update porosity as computed from electrode input
    porosity_Cell_[iCell] = porosity;
    //    porosity_Cell_[iCell] = 0.64007;

    // Let's compute a characteristic electrical conductivity based on these initial conditions.
    //We only need to do this for one control volume
    if (iCell == 0 ) {
      //Fill thermodynamic state
	SetupThermoShop1( nodeCent, &(soln[indexCent_EqnStart_BD]));
      //This computes conductivity based on unit potential gradient and no species gradients
      electrolyteConduct_ = trans_->getElectricConduct();
      fprintf(stderr, "The cathode conductivity is %g S/m before porosity\n",electrolyteConduct_ );
      //  Correct species diffusivities according to the tortuosity
      // Porosity is known from porosity_Cell_ set in domain_prep
      double bruggemannExp = 2.5;
      Tortuosity tort(bruggemannExp );
      electrolyteConduct_ *= tort.tortuosityFactor( porosity_Cell_[iCell] );
      fprintf(stderr, "and %g S/m after porosity\n",electrolyteConduct_ );
    }  //end if (iCell==0)
    

    // update the electrode object--will be used for better estimate of initial voltage
    ee->setState_TP(temp_Curr_, pres_Curr_);
    ee->setVoltages( soln[indexCent_EqnStart_BD + iVar_Voltage_ED], soln[indexCent_EqnStart_BD + iVar_Voltage_BD] );
    ee->setElectrolyteMoleNumbers( &(soln[indexCent_EqnStart_BD + iVar_Species_BD ]), true );
    ee->updateState();
    ee->setInitStateFromFinal(true);
    ee->setFinalFinalStateFromFinal();
    

    //better estimate for initial voltage
    //first we need to determine what BC is specified
    switch (PSinput.cathodeBCType_)
      {
      case 0:
	//update the electrolyte voltage with the open circuit voltage    
	soln[indexCent_EqnStart_BD + iVar_Voltage_ED] = PSinput.CathodeVoltageSpecified_;
	soln[indexCent_EqnStart_BD + iVar_Voltage_BD] = soln[indexCent_EqnStart_BD + iVar_Voltage_ED] - ee->openCircuitVoltage(0);
	break;
      case 2:
	soln[indexCent_EqnStart_BD + iVar_Voltage_ED] = PSinput.CathodeVoltageSpecified_ * PSinput.TimeDepFunction_(0.0);
	soln[indexCent_EqnStart_BD + iVar_Voltage_BD] = soln[indexCent_EqnStart_BD + iVar_Voltage_ED] - ee->openCircuitVoltage(0);
	break;
      case 4:
      case 6:
      case 8:
	soln[indexCent_EqnStart_BD + iVar_Voltage_ED] = PSinput.CathodeVoltageSpecified_ * PSinput.BC_TimeDep_->value(0.0);
	soln[indexCent_EqnStart_BD + iVar_Voltage_BD] = soln[indexCent_EqnStart_BD + iVar_Voltage_ED] - ee->openCircuitVoltage(0);
	break;

      }

  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void porousLiKCl_FeS2Cathode_dom1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln,
						  Epetra_Vector_Ghosted & atolVector,
						  const Epetra_Vector_Ghosted * const atolV)
{
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

    /*
     *    Get the index for the center node
     */
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];

    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled -> so this is a characteristic value
     */
    double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
    atolVector[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = std::max(1.0E-4, 1.0E-1 * vax);

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E-2;
    if (val < 1.0E-12) {
      val = 1.0E-12;
    }
    for (int k = 0; k < nsp_; k++) {
      atolVector[indexCent_EqnStart_BD + iVar_Species_BD + k] = val;
    }

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.-> so this is a characteristic value
     *         1 kcal gmol-1 = 0.05 volts
     */
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05;
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.05;
  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 * @param atolVector Reference for the atol vector to fill up
 */
void porousLiKCl_FeS2Cathode_dom1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln,
							  const Epetra_Vector_Ghosted & solnDot,
							  Epetra_Vector_Ghosted & atolVector,
							  const Epetra_Vector_Ghosted * const atolV)
{
 
  for (int iCell = 0; iCell < NumLcCells; iCell++) {
    cIndex_cc_ = iCell;

    /*
     *  ---------------- Get the index for the center node ---------------------------------
     */
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled.
     */
    double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
    atolVector[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = std::max(1.0E-4, 1.0E-1 * vax);

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E-2;
    if (val < 1.0E-6) {
      val = 1.0E-6;
    }
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 0] = val;
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 1] = val;
    atolVector[indexCent_EqnStart_BD + iVar_Species_BD + 2] = val;

    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.
     */
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05;
    atolVector[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.05;
  }  
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::setAtolDeltaDamping(double atolDefault, double relcoeff,
						   const Epetra_Vector_Ghosted & soln,
						   Epetra_Vector_Ghosted & atolDeltaDamping,
						   const Epetra_Vector_Ghosted * const atolV)
{

  for (int iCell = 0; iCell < NumLcCells; iCell++) {

    cIndex_cc_ = iCell;

    /*
     *  ---------------- Get the index for the center node ---------------------------------
     */
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled.
     */
    if (iVAR_Vaxial_BD >= 0) {
      //   double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
      atolDeltaDamping[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 1.0E-4 * relcoeff;
    }

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E3;
    if (val < 1.0E-4) {
      val = 1.0E-4;
    }
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 0] = val * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 1] = val * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 2] = val * relcoeff;
    
    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.
     */
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05 * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.05 * relcoeff;
  }
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff,
							   const Epetra_Vector_Ghosted & soln,
							   const Epetra_Vector_Ghosted & solnDot,
							   Epetra_Vector_Ghosted & atolDeltaDamping,
							   const Epetra_Vector_Ghosted * const atolV)
{

  for (int iCell = 0; iCell < NumLcCells; iCell++) {

    cIndex_cc_ = iCell;

    /*
     *  ---------------- Get the index for the center node ---------------------------------
     */
    int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    /*
     *   Get the pointer to the NodalVars object for the center node
     */
    NodalVars * nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Index of the first equation in the bulk domain of center node
     */
    int indexCent_EqnStart_BD = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode]
        + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
    int iVAR_Vaxial_BD = BDD_.VariableIndexStart_VarName[Velocity_Axial];
    int iVar_Species_BD = BDD_.VariableIndexStart_VarName[MoleFraction_Species];
    int iVar_Voltage_BD = BDD_.VariableIndexStart_VarName[Voltage];
    int iVar_Voltage_ED = iVar_Voltage_BD + 1;

    /*
     * Set the atol value for the axial velocity
     *   arithmetically scaled.
     */
    if (iVAR_Vaxial_BD >= 0) {
      //   double vax = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];
      atolDeltaDamping[indexCent_EqnStart_BD + iVAR_Vaxial_BD] = 1.0E-2 * relcoeff;
    }

    /*
     * Set atol values for the species mole fractions
     */
    double val = atolDefault * 1.0E7;
    if (val < 1.0E-1) {
      val = 1.0E-1;
    }
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 0] = val * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 1] = val * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Species_BD + 2] = val * relcoeff;
    
    /*
     * Set the atol value for the electrolyte voltage
     *      arithmetically scaled.
     */
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_BD] = 0.05 * relcoeff;
    atolDeltaDamping[indexCent_EqnStart_BD + iVar_Voltage_ED] = 0.05 * relcoeff;
  }
}

//=====================================================================================================================
/**
 * Method to check for precipitation of the salts.
 * Returns index of offending cation or -1 if no precipitation
 */
int 
porousLiKCl_FeS2Cathode_dom1D::checkPrecipitation() {

  // molten salt phase
  string id_salt = "LiKCl_Margules";
  int iph = (PSinput.PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_FeS2Cathode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  ThermoPhase* tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  
  MargulesVPSSTP *salt ;
  salt = dynamic_cast<Cantera::MargulesVPSSTP *>( tmpPhase->duplMyselfAsThermoPhase() );
  
  int iKCl_l = salt->speciesIndex("KCl(L)");
  int iLiCl_l = salt->speciesIndex("LiCl(L)");
  int iK_ion = ionicLiquid_->speciesIndex("K+");
  int iLi_ion = ionicLiquid_->speciesIndex("Li+");

  //solid LiCl phase
  id_salt = "LiCl(S)";
  iph = (PSinput.PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_FeS2Cathode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  Cantera::ThermoPhase *LiCl_solid = tmpPhase->duplMyselfAsThermoPhase() ;

  //solid KCl phase
  id_salt = "KCl(S)";
  iph = (PSinput.PhaseList_)->globalPhaseIndex(id_salt);
  if (iph < 0) {
    throw CanteraError("porousLiKCl_FeS2Cathode_dom1D::checkPrecipitation()",
                       "Can't find the phase in the phase list: " + id_salt);
  }
  tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
  Cantera::ThermoPhase *KCl_solid = tmpPhase->duplMyselfAsThermoPhase() ;



  //set current states
  //mole fraction of the electrolyte ions are held in
  // mfElectrolyte_Thermo_Curr_[3]
  double x[2];
  x[iKCl_l] = 2.0 * mfElectrolyte_Thermo_Curr_[iK_ion];
  x[iLiCl_l] = 2.0 * mfElectrolyte_Thermo_Curr_[iLi_ion];
  salt->setState_TPX(temp_Curr_, pres_Curr_, x);
  
  /*  
  string f_licl = "LiCl_solid.xml";
  string id = "LiCl(S)";
  Cantera::ThermoPhase *LiCl_solid = Cantera::newPhase(f_licl, id);
  */
  LiCl_solid->setState_TP(temp_Curr_, pres_Curr_);
  
  /*  
  string f_kcl = "KCl_solid.xml";
  id = "KCl(S)";
  Cantera::ThermoPhase *KCl_solid = Cantera::newPhase(f_kcl, id);
  */
  KCl_solid->setState_TP(temp_Curr_, pres_Curr_);
  
  //get chemical potentials
  double mu[3];
  //molten salt
  salt->getChemPotentials(mu);
  double mu_LiCl_liq = mu[iLiCl_l];
  double mu_KCl_liq = mu[iKCl_l];
  //solid LiCl
  LiCl_solid->getChemPotentials(mu);
  double mu_LiCl_solid = mu[0];
  //solid KCl
  KCl_solid->getChemPotentials(mu);
  double mu_KCl_solid = mu[0];
  
  if ( mu_KCl_solid < mu_KCl_liq ) {
    std::cout << "KCl is precipitating at mole fraction "
	      << 2.0 * mfElectrolyte_Thermo_Curr_[iK_ion] << std::endl;
    return iK_ion;
  }
  if ( mu_LiCl_solid < mu_LiCl_liq ) {
    std::cout << "LiCl is precipitating at mole fraction "
	      << 2.0 * mfElectrolyte_Thermo_Curr_[iLi_ion] << std::endl;
    return iLi_ion;
  }
  
  //no precipitation
  return -1;
  
}
//=====================================================================================================================
// returns value for depthOfDischarge_Electrode_ and capacityZeroDoD_Electrode_
void
porousLiKCl_FeS2Cathode_dom1D::getCapacity( double *initialCap, double *dischargedCap )
{
  depthOfDischarge_Electrode_ = 0.0;
  capacityZeroDoD_Electrode_  = 0.0;
  for (int iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
    int iCell = iGbNode - BDD_.FirstGbNode;
    depthOfDischarge_Electrode_ += depthOfDischargePA_Cell_[iCell];
    capacityZeroDoD_Electrode_ += capacityPA_Cell_[iCell];
  }
  *initialCap = capacityZeroDoD_Electrode_;
  *dischargedCap = depthOfDischarge_Electrode_;
}
//=====================================================================================================================
double
porousLiKCl_FeS2Cathode_dom1D::getInitialCathodeMass()
{
  //something to find the mass of the cathode originally
  BDT_porCathode_LiKCl * BDD_FeS2_Cathode = dynamic_cast<BDT_porCathode_LiKCl *>(&(BDD_));
  if (!BDD_FeS2_Cathode) {
    throw CanteraError(" porousLiKCl_FeS2Cathode_dom1D::getInitialCathodeMass()", "bad dynamic cast ");
  }

  Cantera::Electrode *ee = BDD_FeS2_Cathode->Electrode_;

  double molesFe = ee->elementSolidMoles("Fe");
  double MW_FeS2 = 119.98;
  double mass_FeS2 = 0.0;
  if (molesFe > 0.0) {
      mass_FeS2 = molesFe * MW_FeS2 * 1.0e3; //mass in grams
  }
  //  fprintf(stderr,"got mass %g\n",mass_FeS2);
  return mass_FeS2;
}  
//=====================================================================================================================
// Get parameters specified by text strings
int 
porousLiKCl_FeS2Cathode_dom1D::getSolutionParam(std::string paramName, double * const paramVal)
{
  if (paramName == "CapacityZeroDoD") {
    //make sure this is updated
    paramVal[0] = capacityZeroDoD_Electrode_ ;
  }
  else if (paramName == "DepthOfDischarge") {
    //make sure this is updated
    paramVal[0] = depthOfDischarge_Electrode_ ;
  }
  else if (paramName == "SpecificCapacityZeroDoD") {
    //make sure this is updated
    //this is per gram of FeS2
    paramVal[0] = capacityZeroDoD_Electrode_ * electrodeCrossSectionalArea_ / getInitialCathodeMass();
  }
  else if (paramName == "SpecificDepthOfDischarge") {
    //make sure this is updated
    //this is per gram of FeS2
    paramVal[0] = depthOfDischarge_Electrode_ * electrodeCrossSectionalArea_ / getInitialCathodeMass();
  }
  else {
    paramVal[0] = 0.0;
    return -1 ;
  }
  return 1;
}
//=====================================================================================================================
void
porousLiKCl_FeS2Cathode_dom1D::err(const char *msg)
{
  printf("porousLiKCl_FeS2Cathode_dom1D: function not implemented: %s\n", msg);
  exit(-1);
}
//=====================================================================================================================
//=====================================================================================================================
double porousLiKCl_FeS2Cathode_dom1D::capacityPA(int platNum) const
{
    double totalCapacity = 0.0;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
        totalCapacity += ee->capacity(platNum);
        capacityPA_Cell_[iCell]= ee->capacity() / crossSectionalArea_;
    }
    totalCapacity /= crossSectionalArea_;
    return totalCapacity;
}
//=====================================================================================================================
double porousLiKCl_FeS2Cathode_dom1D::capacityDischargedPA(int platNum) const
{
    double totalCapacity = 0.0;
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
        totalCapacity += ee->capacityDischarged(platNum);
        capacityDischargedPA_Cell_[iCell] = ee->capacityDischarged() / crossSectionalArea_;
    }
    totalCapacity /= crossSectionalArea_;
    return totalCapacity;
}
//=====================================================================================================================
double porousLiKCl_FeS2Cathode_dom1D::capacityLeftPA(int platNum, double voltsMax, double voltsMin) const
{
    double totalCapacity = 0.0;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
        totalCapacity += ee->capacityLeft(platNum, voltsMax, voltsMin);
        capacityLeftPA_Cell_[iCell] = ee->capacityLeft() / crossSectionalArea_;
    }
    totalCapacity /= crossSectionalArea_;
    return totalCapacity;
}
//=====================================================================================================================
double porousLiKCl_FeS2Cathode_dom1D::depthOfDischargePA(int platNum) const
{
    double dod = 0.0;
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
        dod += ee->depthOfDischarge(platNum);
        depthOfDischargePA_Cell_[iCell] = ee->depthOfDischarge() / crossSectionalArea_;
    }
    dod /= crossSectionalArea_;
    return dod;
}
//=====================================================================================================================
double porousLiKCl_FeS2Cathode_dom1D::depthOfDischargeStartingPA(int platNum) const
{
    double dodStarting = 0.0;
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; iCell++) {
        dodStarting += Electrode_Cell_[iCell]->depthOfDischargeStarting(platNum);
    }
    return dodStarting / crossSectionalArea_;
}
//=====================================================================================================================
// Reset the counters that keep track of the amount of discharge to date
void porousLiKCl_FeS2Cathode_dom1D::resetCapacityDischargedToDate()
{
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
        ee->resetCapacityDischargedToDate();
    }
}
//=====================================================================================================================
//! Return a value for the open circuit potential without doing a formally correct calculation
/*!
 *  Currently this is defined as the open circuit potential on the outside electrode.
 *
 *   @return return the open circuit potential 
 */
double porousLiKCl_FeS2Cathode_dom1D::openCircuitPotentialQuick() const
{
    Electrode* ee = Electrode_Cell_[NumLcCells-1];
    double ocv = ee->openCircuitVoltage(nSurfsElectrode_ - 1);
    return ocv;
}
//=====================================================================================================================

} //namespace m1d
//=====================================================================================================================
