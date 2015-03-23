/**
 * @file m1d_porousLiIon_Cathode_dom1D.cpp
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_porousLiIon_Cathode_dom1D.h"

#include "LiIon_PorousBat.h"
#include "m1d_BDT_porCathode_LiIon.h"
#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_DomainLayout_LiIon_PorousBat.h"
#include "m1d_ProblemStatementCell.h"
#include "Electrode.h"
#include "Electrode_Factory.h"

#include "cantera/transport/Tortuosity.h"

using namespace std;

//
//  Necessary expediency until we model dUdT correctly and fully, matching to experiment
//
#define DELTASHEAT_ZERO 1
//=====================================================================================================================
namespace m1d
{

//=====================================================================================================================
porousLiIon_Cathode_dom1D::porousLiIon_Cathode_dom1D(BDT_porCathode_LiIon& bdd) :
    porousElectrode_dom1D(bdd),
    ionicLiquid_(0), 
    trans_(0), nph_(0), nsp_(0),
    concTot_cent_(0.0),
    concTot_cent_old_(0.0),
    icurrInterfacePerSurfaceArea_Cell_(0), xdelCell_Cell_(0),
    concTot_Cell_(0), concTot_Cell_old_(0),
    capacityDischargedPA_Cell_(0),
    depthOfDischargePA_Cell_(0),
    capacityLeftPA_Cell_(0),
    capacityPA_Cell_(0),
    Fleft_cc_(0.0), Fright_cc_(0.0), Vleft_cc_(0.0),
    Vcent_cc_(0.0), Vright_cc_(0.0), VElectrodeLeft_cc_(0.0), VElectrodeCent_cc_(0.0), VElectrodeRight_cc_(0.0),
    t_final_(0.0),
    t_init_(0.0),
    Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0), spCharge_(0), mfElectrolyte_Soln_Curr_(0),
    mfElectrolyte_Thermo_Curr_(0), 
    phiElectrode_Curr_(0.0), conductivityElectrode_(1.0E2),
    gradT_trCurr_(0.0), gradV_trCurr_(0.0), gradVElectrode_trCurr_(0.0), gradX_trCurr_(0), Vdiff_trCurr_(0),
    jFlux_trCurr_(0), icurrElectrode_trCurr_(0.0),
    nSpeciesElectrode_(0), nSurfsElectrode_(0),
    electrodeSpeciesMoleDelta_Cell_(0),
    icurrInterface_Cell_(0),
    phaseMoleTransfer_(0), solnMoleFluxInterface_Cell_(0), icurrElectrode_CBL_(0), icurrElectrode_CBR_(0),
    icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0), deltaV_Cell_(0), 
    Ess_Surf_Cell_(0), 
    overpotential_Surf_Cell_(0),
    icurrRxn_Cell_(0), LiFlux_Cell_(0),
    iECDMC_(-1),
    iLip_(-1),
    iPF6m_(-1),
    solnTemp(0)
{
    BDT_porCathode_LiIon* fa = dynamic_cast<BDT_porCathode_LiIon*>(&bdd);
    if (!fa) {
        throw m1d_Error("confused", "confused");
    }
    /*
     * This is a shallow pointer copy. The BDT object owns the ionicLiquid_ object
     */
    ionicLiquid_ = fa->ionicLiquid_;
    /*
     *  This is a shallow pointer copy. The BDT object owns the transport object
     */
    trans_ = fa->trans_;
    /*
     *  This is a shallow pointer copy. The BDT object owns the Electrode object
     */
    // Electrode_ = fa->Electrode_;
    nsp_ = 3;
    nph_ = 1;

    iECDMC_ = ionicLiquid_->speciesIndex("ECDMC");
    if (iECDMC_ < 0) {
        throw CanteraError("confused", "confused");
    }
    iLip_ = ionicLiquid_->speciesIndex("Li+");
    if (iLip_ < 0) {
        throw CanteraError("confused", "confused");
    }
    iPF6m_ = ionicLiquid_->speciesIndex("PF6-");
    if (iPF6m_ < 0) {
        throw CanteraError("confused", "confused");
    }

    conductivityElectrode_ = PSinput.conductivityCathode_;

}
//=====================================================================================================================
porousLiIon_Cathode_dom1D::porousLiIon_Cathode_dom1D(const porousLiIon_Cathode_dom1D& r) :
    porousElectrode_dom1D((m1d::BDD_porousElectrode&)   r.BDD_),
    ionicLiquid_(0), 
    trans_(0), nph_(0), nsp_(0),
    concTot_cent_(0.0),
    concTot_cent_old_(0.0),
    icurrInterfacePerSurfaceArea_Cell_(0), xdelCell_Cell_(0),
    concTot_Cell_(0),
    concTot_Cell_old_(0),
    capacityDischargedPA_Cell_(0),
    depthOfDischargePA_Cell_(0),
    capacityLeftPA_Cell_(0),
    capacityPA_Cell_(0),
    Fleft_cc_(0.0), Fright_cc_(0.0), Vleft_cc_(0.0),
    Vcent_cc_(0.0), Vright_cc_(0.0), VElectrodeLeft_cc_(0.0), VElectrodeCent_cc_(0.0), VElectrodeRight_cc_(0.0),
    t_final_(0.0),
    t_init_(0.0),
    Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0), spCharge_(0), mfElectrolyte_Soln_Curr_(0),
    mfElectrolyte_Thermo_Curr_(0),
    phiElectrode_Curr_(0.0), conductivityElectrode_(1.0E2),
    gradT_trCurr_(0.0), gradV_trCurr_(0.0), gradVElectrode_trCurr_(0.0), gradX_trCurr_(0), Vdiff_trCurr_(0),
    jFlux_trCurr_(0), icurrElectrode_trCurr_(0.0),
    nSpeciesElectrode_(0), nSurfsElectrode_(0),
    electrodeSpeciesMoleDelta_Cell_(0),
    icurrInterface_Cell_(0),
    phaseMoleTransfer_(0), solnMoleFluxInterface_Cell_(0), icurrElectrode_CBL_(0), icurrElectrode_CBR_(0),
    icurrElectrolyte_CBL_(0), icurrElectrolyte_CBR_(0), deltaV_Cell_(0), 
    Ess_Surf_Cell_(0), 
    overpotential_Surf_Cell_(0),
    icurrRxn_Cell_(0), LiFlux_Cell_(0),
    iECDMC_(-1),
    iLip_(-1),
    iPF6m_(-1),
    solnTemp(0)
{
    porousLiIon_Cathode_dom1D::operator=(r);

    conductivityElectrode_ = PSinput.conductivityCathode_;
}
//=====================================================================================================================
porousLiIon_Cathode_dom1D::~porousLiIon_Cathode_dom1D()
{
}
//=====================================================================================================================
porousLiIon_Cathode_dom1D&
porousLiIon_Cathode_dom1D::operator=(const porousLiIon_Cathode_dom1D& r)
{
    if (this == &r) {
        return *this;
    }
    // Call the parent assignment operator
    porousElectrode_dom1D::operator=(r);

    ionicLiquid_ = r.ionicLiquid_;
    trans_ = r.trans_;
    //Electrode_ = r.Electrode_;

    nph_ = r.nph_;
    nsp_ = r.nsp_;
    concTot_cent_ = r.concTot_cent_;
    concTot_cent_old_ = r.concTot_cent_old_;
    icurrInterfacePerSurfaceArea_Cell_ = r.icurrInterfacePerSurfaceArea_Cell_;
    xdelCell_Cell_ = r.xdelCell_Cell_;
    concTot_Cell_old_ = r.concTot_Cell_old_;

    capacityDischargedPA_Cell_ = r.capacityDischargedPA_Cell_;
    depthOfDischargePA_Cell_ = r.depthOfDischargePA_Cell_;
    capacityLeftPA_Cell_ = r.capacityLeftPA_Cell_;
    capacityPA_Cell_ = r.capacityPA_Cell_;
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
    mfElectrolyte_Soln_Cell_old_          = r.mfElectrolyte_Soln_Cell_old_;
    phiElectrode_Curr_ = r.phiElectrode_Curr_;
    concTot_Curr_ = r.concTot_Curr_;
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
    phaseMoleTransfer_ = r.phaseMoleTransfer_;
    solnMoleFluxInterface_Cell_ = r.solnMoleFluxInterface_Cell_;
    icurrElectrode_CBL_ = r.icurrElectrode_CBL_;
    icurrElectrode_CBR_ = r.icurrElectrode_CBR_;
    icurrElectrolyte_CBL_ = icurrElectrolyte_CBL_;
    icurrElectrolyte_CBR_ = r.icurrElectrolyte_CBR_;
    deltaV_Cell_ = r.deltaV_Cell_;
    Ess_Surf_Cell_ = r.Ess_Surf_Cell_;
    overpotential_Surf_Cell_ = r.overpotential_Surf_Cell_;
    icurrRxn_Cell_ = r.icurrRxn_Cell_;
    LiFlux_Cell_ = r.LiFlux_Cell_;
    iECDMC_ = r.iECDMC_;
    iLip_ = r.iLip_;
    iPF6m_ = r.iPF6m_;
    maxElectrodeSubIntegrationSteps_ = r.maxElectrodeSubIntegrationSteps_;
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
porousLiIon_Cathode_dom1D::domain_prep(LocalNodeIndices* li_ptr)
{
    /*
     * First call the parent domain prep to get the node information
     */
    porousElectrode_dom1D::domain_prep(li_ptr);

    BDT_porCathode_LiIon* BDD_FeS2_Cathode = dynamic_cast<BDT_porCathode_LiIon*>(&(BDD_));
    if (!BDD_FeS2_Cathode) {
        throw CanteraError(" porousLiIon_Cathode_dom1D::domain_prep()", "bad dynamic cast ");
    }

    /*
     *  In the bulk domain description object, we have a prototypical Electrode object already instantiated
     *  for the cathode. We will use this to figure out the number of Electrode species and the number of
     *  Electrode surfaces.
     */
    Cantera::Electrode* ee = BDD_FeS2_Cathode->Electrode_;
    nSpeciesElectrode_ = ee->nSpecies();
    nSurfsElectrode_ = ee->nSurfaces();

    for (int i = 0; i < NumLcCells; i++) {
        porosity_Cell_[i] = 0.5;
        porosity_Cell_old_[i] = 0.5;
    }
    /*
     * Porous electrode domain prep
     */
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
    Xleft_cc_.resize(nsp_, 0.0);
    Xcent_cc_.resize(nsp_, 0.0);
    Xright_cc_.resize(nsp_, 0.0);


    spCharge_.resize(nsp_, 0.0);
    for (int k = 0; k < nsp_; k++) {
        spCharge_[k] = ionicLiquid_->charge(k);
    }

    mfElectrolyte_Soln_Curr_.resize(nsp_, 0.0);
    mfElectrolyte_Thermo_Curr_.resize(nsp_, 0.0);

    gradX_trCurr_.resize(3, 0.0);
    Vdiff_trCurr_.resize(3, 0.0);
    jFlux_trCurr_.resize(3, 0.0);

    electrodeSpeciesMoleDelta_Cell_.resize(nSpeciesElectrode_ * NumLcCells, 0.0);
    phaseMoleTransfer_.resize(20, 0.0);

    solnTemp.resize(10, 0.0);

    icurrElectrode_CBL_.resize(NumLcCells, 0.0);
    icurrElectrode_CBR_.resize(NumLcCells, 0.0);

    icurrElectrolyte_CBL_.resize(NumLcCells, 0.0);
    icurrElectrolyte_CBR_.resize(NumLcCells, 0.0);

    deltaV_Cell_ .resize(NumLcCells, 0.0);
    Ess_Surf_Cell_.resize(nSurfsElectrode_ * NumLcCells, 0.0);
    overpotential_Surf_Cell_.resize(nSurfsElectrode_ * NumLcCells, 0.0);
    icurrRxn_Cell_.resize(NumLcCells, 0.0);
    LiFlux_Cell_.resize(NumLcCells, 0.0);

    mfElectrolyte_Soln_Cell_old_.resize(nsp_, NumLcCells, 0.0);

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
porousLiIon_Cathode_dom1D::instantiateElectrodeCells()
{
    /*
     * Indices and nodes used to determine the volume of each cell.
     */
    int index_RightLcNode;
    int index_LeftLcNode;
    int index_CentLcNode;
    NodalVars* nodeLeft = 0;
    NodalVars* nodeCent = 0;
    NodalVars* nodeRight = 0;

    double xCellBoundaryL; //cell boundary left
    double xCellBoundaryR; //cell boundary right

    for (int iCell = 0; iCell <NumLcCells; iCell++) {
        // Note that we need to create a Factory method to instantiate the desired electrode type.
        // int nSurPhases = PSinput.cathode_input_->m_pl->nSurPhases();

        DomainLayout* DL = BDD_.DL_ptr_;

        DomainLayout_LiIon_PorousBat* dlc = dynamic_cast<DomainLayout_LiIon_PorousBat*>(DL);
        ProblemStatementCell* psc_ptr = dlc->pscInput_ptr_;
        ELECTRODE_KEY_INPUT* ci = psc_ptr->cathode_input_;

        Cantera::Electrode* ee  = newElectrodeObject(ci->electrodeModelName);
        if (!ee) {
            throw  m1d_Error("porousLiIon_Cathode_dom1D::instantiateElectrodeCells()",
                             "Electrode factory method failed");
        }
        /*
         *  Do thermal property calculations within the electrode when we are doing heat
         *  transfer calculations
         */
        if (PSinput.doHeatSourceTracking_) {
         ee->doThermalPropertyCalculations_ = true;
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
         *  turn on the saving of states in xml files. We need this for restarts.
         */
        ee->electrode_stateSave_create();

        /*
         * Compute cell volume
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
        double electrodeGrossArea = -1.0;
        double porosity = -1.0;
        //
        //  see if we know the electrode area
        //
        if (PSinput.cathode_input_->electrodeGrossArea > 0.0) {
            electrodeGrossArea = PSinput.cathode_input_->electrodeGrossArea;
        } else if (PSinput.cathode_input_->electrodeGrossDiameter > 0.0) {
            electrodeGrossArea = Pi * 0.25 * PSinput.cathode_input_->electrodeGrossDiameter *
                                 PSinput.cathode_input_->electrodeGrossDiameter ;
        }

       /*
         * if we know the solid and electrode volume, compute porosity
         * otherwise, hopefully we know the porosity
         */
	if (PSinput.cathode_input_->porosity > 0.0) {
	    porosity = PSinput.cathode_input_->porosity;
	} else {
            if (PSinput.cathode_input_->electrodeGrossThickness > 0.0 && electrodeGrossArea > 0.0) {
               totalElectrodeVolume = electrodeGrossArea * PSinput.cathode_input_->electrodeGrossThickness;
               porosity = 1.0 - ee->SolidVol() / totalElectrodeVolume;
               printf("Cathode porosity is %f with %g m^3 solid volume and %g m^3 electrode volume.\n",porosity, ee->SolidVol(),
                      totalElectrodeVolume);
               if (porosity <= 0.0) {
                   throw CanteraError("porousLiIon_Cathode_dom1D::instantiateElectrodeCells()",
                                      "Computed porosity is not positive.");
               }
            } else {
               throw m1d_Error("porousLiIon_Cathode_dom1D::instantiateElectrodeCells()",
                               "Unknown Porosity");
            }
        } 

        /*
         *  Save the cross-sectional area of the electrode to use throughout the code. It does not change within
         *  this calculation
         */
        if (electrodeGrossArea > 0.0) {
	    if (!doubleEqual(crossSectionalArea_, electrodeGrossArea, 1.0E-30, 9)) {
		throw m1d_Error("porousLiIon_Cathode_dom1D::initialConditions()",
				"crossSectional Area, " + fp2str(crossSectionalArea_) + 
				", differs from cross sectional area input from cathode file, " + fp2str(electrodeGrossArea));
	    }
	    electrodeGrossArea = crossSectionalArea_;
        } else {
	    electrodeGrossArea = crossSectionalArea_;
	}

	/*
	 *  reset the moles in the electrode using the computed porosity
	 * and the electrodeWidth FOR THIS node.
	 */
	ee->setElectrodeSizeParams(electrodeGrossArea, xdelCell_Cell_[iCell], porosity);

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
porousLiIon_Cathode_dom1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                                               const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr,
                                               const double t, const double t_old)
{
    const Epetra_Vector& soln = *soln_ptr;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;
        /*
         *  ---------------- Get the index for the center node ---------------------------------
         */
        int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
        /*
         *   Get the pointer to the NodalVars object for the center node
         */
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];

        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

        concTot_Cell_old_[iCell] = concTot_Curr_;
        porosity_Cell_old_[iCell] = porosity_Curr_;

        double* mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
        mfElectrolyte_Soln_old[0] = mfElectrolyte_Soln_Curr_[0];
        mfElectrolyte_Soln_old[1] = mfElectrolyte_Soln_Curr_[1];
        mfElectrolyte_Soln_old[2] = mfElectrolyte_Soln_Curr_[2];
        /*
         * Tell the electrode object to accept the current step and prep for the next step.
         *
         * We might at this point do a final integration to make sure we nailed the conditions of the last step.
         * However, we will hold off at implementing this right now
         */
        Electrode* ee = Electrode_Cell_[iCell];
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
porousLiIon_Cathode_dom1D::revertToInitialGlobalTime()
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
 *  All residual terms are written on a per area basis
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
porousLiIon_Cathode_dom1D::residEval(Epetra_Vector& res,
                                     const bool doTimeDependentResid,
                                     const Epetra_Vector* soln_ptr,
                                     const Epetra_Vector* solnDot_ptr,
                                     const Epetra_Vector* solnOld_ptr,
                                     const double t,
                                     const double rdelta_t,
                                     const  ResidEval_Type_Enum residType,
                                     const Solve_Type_Enum solveType)
{
    static int tmpsSetup = 0;
    if (!tmpsSetup) {
	residSetupTmps();
	tmpsSetup = 1;
    }
    residType_Curr_ = residType;

    int index_CentLcNode;

    t_final_ = t;
    if (rdelta_t < 1.0E-200) {
        t_init_ = t;
    } else {
        t_init_ = t - 1.0/rdelta_t;
    }

    NodalVars* nodeLeft = 0;
    NodalVars* nodeCent = 0;
    NodalVars* nodeRight = 0;

    double xdelL; // Distance from the center node to the left node
    double xdelR; // Distance from the center node to the right node
    double xdelCell; // cell width - right boundary minus the left boundary.

    //  Electrolyte total molar fluxe on Right side of cell - this is c V dot n at the boundaries of the cells - kmol s-1 m-2
    double fluxFright = 0.;
    //  Electrolyte total molar fluxe on Left side of cell - this is c V dot n at the boundaries of the cells - kmol s-1 m-2
    double fluxFleft;

    //mole fraction fluxes
    std::vector<double> fluxXright(nsp_, 0.0);
    std::vector<double> fluxXleft(nsp_, 0.0);

    double fluxL = 0.0;
    double fluxR = 0.0;

    const Epetra_Vector& soln = *soln_ptr;
    //  const Epetra_Vector &solnOld = *solnOld_ptr;

    /*
     * Index of the first equation at the left node corresponding to the first bulk domain, which is the electrolyte
     */
    int indexLeft_EqnStart;
    /*
     * Index of the first equation at the center node corresponding to the first bulk domain, which is the electrolyte
     */
    int indexCent_EqnStart;
    /*
     * Index of the first equation at the right node corresponding to the first bulk domain, which is the electrolyte
     */
    int indexRight_EqnStart;

    /*
     *   Find the species index for the first species in the electrode object pertaining to the electrolyte
     */
    int sf = Electrode_Cell_[0]->solnPhaseIndex();
    int indexStartEOelectrolyte = Electrode_Cell_[0]->getGlobalSpeciesIndex(sf, 0);
    /*
     * offset of the electolyte solution unknowns at the current node
     */
    index_CentLcNode = Index_DiagLcNode_LCO[0];
    nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

    /*
     *  Offsets for the equation unknowns in the residual vector for the electrolyte domain
     */

    /*
     * Offsets for the variable unknowns in the solution vector for the electrolyte domain
     */
  
  

    Fright_cc_ = 0.0;

    /*
     *  -------------------- Special section to do the left boundary -------------------------------
     */

    /*
     * Special section if we own the left node of the domain. If we do
     * it will be cell 0
     */
    if (IOwnLeft) {
	cellTmps& cTmps          = cellTmpsVect_Cell_[0];
        NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
        DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity] = fluxL;

        //DiffFluxLeftBound_LastResid_NE[EQ_TCont_offset] = fluxL;
        indexCent_EqnStart =
            LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode] + nodeCent->OffsetIndex_BulkDomainEqnStart_BDN[0];
        for (int k = 0; k < nsp_; k++) {
            DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset  + k] = 0.0;
        }
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
        Electrode* Electrode_ptr = Electrode_Cell_[iCell];

	cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
	NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

#ifdef DEBUG_RESID
        if (counterResBaseCalcs_ > 125 && residType == Base_ResidEval) {
            if (iCell == NumLcCells - 1) {
                // printf("we are here\n");
            }
        }
#endif
        /*
         *  Placeholder for debugging
         */
        if (IOwnLeft && iCell == 0) {
            if (residType == Base_ShowSolution) {
                cIndex_cc_ = iCell;
            }
        }
	/*
         *  ---------------- Get the index for the center node ---------------------------------
         *   Get the pointer to the NodalVars object for the center node
         *   Index of the first equation in the bulk domain of center node
         */
        nodeCent = cTmps.nvCent_;
        indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
        /*
         *  ------------------- Get the index for the left node -----------------------------
         *    There may not be a left node if we are on the left boundary. In that case
         *    set the pointer to zero and the index to -1. Hopefully, we will get a segfault on an error.
         */
	nodeLeft = cTmps.nvLeft_;
	indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;
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
        nodeRight = cTmps.nvRight_;
        indexRight_EqnStart = nodeTmpsRight.index_EqnStart;
        /*
         * --------------------------- CALCULATE POSITION AND DELTA_X Variables -----------------------------
         * Calculate the distance between the left and center node points
         */
        /*
         * --------------------------- CALCULATE POSITION AND DELTA_X Variables -----------------------------
         * Calculate the distance between the left and center node points
         */
        xdelL = cTmps.xdelL_;
        /*
         * Calculate the distance between the right and center node points
         */
        xdelR = cTmps.xdelR_;
        /*
         * Calculate the cell width
         */
        xdelCell = cTmps.xdelCell_;
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
            Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
        }
        Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];
	// HKM fix up - technically correct
        VElectrodeCent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage + 1];


        /*
         * Advance the electrode forward and compute the new porosity
         */

        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));


        if (nodeLeft != 0) {
            /*
             * Find the velocity located at the left cell boundary.
             * The left cell boundary velocity is stored at the previous (left)
             * cell index as per our conventions.
             */
            Fleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Velocity_Axial];
            for (int k = 0; k < nsp_; k++) {
                Xleft_cc_[k] = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_MoleFraction_Species + k];
            }
            Vleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Voltage];
            VElectrodeLeft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Voltage + 1];
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
        Fright_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Velocity_Axial];

        if (nodeRight != 0) {

            for (int k = 0; k < nsp_; k++) {
                Xright_cc_[k] = soln[indexRight_EqnStart + nodeTmpsRight.Offset_MoleFraction_Species + k];
            }
            Vright_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage];
            VElectrodeRight_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage + 1];
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
                SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent, &(soln[indexCent_EqnStart]), 0);

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
                        fluxXleft[k] += Fleft_cc_ * Xleft_cc_[k];
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
            SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
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
            SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart]), nodeRight, &(soln[indexRight_EqnStart]), 1);

            SetupTranShop(xdelR, 1);

            /*
             * Calculate the molar flux at the right boundary for each equation
             * This is equal to
             *       Conc * VMolaraxial * phi
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
#ifdef DEBUG_RESID
        double residBefore = 0.0;
        if (IOwnLeft && iCell == 0) {
            if (residType == Base_ShowSolution) {
                residBefore =  res[indexCent_EqnStart + nodeTmpsCenter.RO_Electrolyte_Continuity ];
                double tmp3 = solnMoleFluxInterface_Cell_[iCell];
                double sum5 = residBefore + fluxFright - tmp3;
                printf(" Cell = %d, ResBefore = -fluxleft = %10.3e, FluxRight = %10.3e, Prod = %10.3e total = %10.3e \n", iCell,
                       residBefore, fluxFright, tmp3, sum5);

            }
        }
#endif
        /*
         *  Total continuity equation - fluxFright and fluxFleft represent the total mole
         *                              fluxes coming and going from the cell.
         *                    R =   d rho d t + del dot (rho V) = 0
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Electrolyte_Continuity] += (fluxFright - fluxFleft);

        /*
         * Species continuity Equation
         */
        for (int k = 0; k < nsp_; k++) {
            if (k != iECDMC_ && k != iPF6m_) {
                res[indexCent_EqnStart + nodeTmpsCenter.RO_Species_Eqn_Offset + k] += (fluxXright[k] - fluxXleft[k]);
            }
        }

        /*
         * Electroneutrality equation
         */

        /*
         *   Current conservation equation - electrolyte
         */
        // res[indexCent_EqnStart + EQ_Current_offset_BD] += (fluxVRight - fluxVLeft);
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation] += (icurrElectrolyte_CBR_[iCell] - icurrElectrolyte_CBL_[iCell]);

        /*
         *   Current conservation equation - electrode
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation + 1] += (icurrElectrode_CBR_[iCell] - icurrElectrode_CBL_[iCell]);

        /*
         *   ------------------- ADD SOURCE TERMS TO THE CURRENT CELL CENTER --------------------------------------
         */
        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

        /*
         * Species Conservation equations
         *    Source terms for the species production rate of Li+.
         */
        for (int k = 0; k < nsp_; k++) {
            if (k != iECDMC_ && k != iPF6m_) {
                res[indexCent_EqnStart + nodeTmpsCenter.RO_Species_Eqn_Offset  + k] -=
                    electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * iCell + indexStartEOelectrolyte + k] * rdelta_t /
                    crossSectionalArea_;
            }
        }

        /*
         * Mole fraction summation equation
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_MFSum_offset] = 1.0;
        for (int k = 0; k < nsp_; k++) {
            res[indexCent_EqnStart + nodeTmpsCenter.RO_MFSum_offset] -= Xcent_cc_[k];
        }

        /*
         * Electroneutrality equation
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_ChargeBal_offset] = 0.0;
        for (int k = 0; k < nsp_; k++) {
            res[indexCent_EqnStart + nodeTmpsCenter.RO_ChargeBal_offset] += Xcent_cc_[k] * spCharge_[k];
        }

        /*
         *  Total continuity equation for moles in the electrolyte phase of the cell
         *     Add in the molar flux from the electrode into the electrolyte phase
         *     We are assuming for the current problem that the volumes stay constant
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Electrolyte_Continuity] -= solnMoleFluxInterface_Cell_[iCell];

        /*
         *   Current conservation equation
         *      These are written as a source term. icurrInterface_Curr_ will be positive for the anode
         *      where current flows into the electrolyte and will be negative for the cathode where current
         *      flows out of the cathode
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation] -= icurrInterface_Cell_[iCell];

        /*
         *   Current conservation equation for the current in the electrode material
         *      These are written as a sink term They will be exactly opposite to the electrolyte current
         *      source terms
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation + 1] += icurrInterface_Cell_[iCell];

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
                DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Current_Conservation] = icurrElectrolyte_CBL_[iCell];
                DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Current_Conservation + 1] = 0.0;
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
                DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Current_Conservation] = icurrElectrolyte_CBR_[iCell];
                DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Current_Conservation + 1] = icurrElectrode_CBR_[iCell];
            }
        }
        if (residType == Base_ShowSolution) {
            deltaV_Cell_[iCell] = Electrode_ptr->potentialDrop();
            for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++) {
                Ess_Surf_Cell_[nSurfsElectrode_ * iCell + jSurf]           = Electrode_ptr->openCircuitVoltage(jSurf);
                overpotential_Surf_Cell_[nSurfsElectrode_ * iCell + jSurf] = Electrode_ptr->overpotential(jSurf);
            }
            icurrRxn_Cell_[iCell] = icurrInterface_Cell_[iCell];
            LiFlux_Cell_[iCell] = jFlux_trCurr_[iLip_];
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
                printf(" Cell = %d, Totalflux_Li+_r = %10.3e,  = %10.3e, Totalflux_Li+_l ", iCell, fluxXright[iLip_], fluxXleft[iLip_]);
            }
#endif
#ifdef DEBUG_HKM_LI
            if (residType == Base_ShowSolution) {
                if (iCell == 7) {
                    printf(" Cell = %d, Totalflux_Li+_R = %- 14.7e, -Totalflux_Li+_L = %- 14.7e", iCell, fluxXright[iLip_],
                           -fluxXleft[iLip_]);
                    printf(", -Source/dt = %- 14.7e",
                           -electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * iCell + indexStartEOelectrolyte + 1]
                           * rdelta_t / crossSectionalArea_);
                }
            }
#endif

            double newStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
            double newStuffSpecies0 = Xcent_cc_[iLip_] * newStuffTC;

            /*
             *   .................... Calculate quantities needed at the previous time step
             */

            double* mf_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
            double oldStuffTC = concTot_Cell_old_[iCell] * porosity_Cell_old_[iCell] * xdelCell;
            double oldStuffSpecies0 = mf_old[iLip_] * oldStuffTC;

            double tmp = (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;

#ifdef DEBUG_HKM_NOT
            if (residType == Base_ResidEval) {
                printf(" deltaT term = %10.3e BulkSum = %10.3e\n", tmp, tmp + (fluxXright[iLip_] - fluxXleft[iLip_]));
            }
#endif
#ifdef DEBUG_HKM_LI
            if (residType == Base_ShowSolution) {
                if (iCell == 7) {
                    printf(", d(XCP)/dt = %- 14.7e", tmp);
                    double sumResidualX = tmp +  fluxXright[iLip_] - fluxXleft[iLip_] -
                                          electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * iCell + indexStartEOelectrolyte + 1]
                                          * rdelta_t / crossSectionalArea_;
                    printf(" ==  %- 14.7e\n", sumResidualX);
                }
            }
#endif

#ifdef DEBUG_HKM_LI
            if (residType == Base_ShowSolution) {
                if (iCell == 7) {
                    printf(" Cell = %d, Totalflux_R = %- 14.7e, -Totalflux_L = %- 14.7e", iCell, fluxFright, -fluxFleft);
                    printf(", -Source/dt = %- 14.7e", -solnMoleFluxInterface_Cell_[iCell]);
                    double sumResidualTC = fluxFright - fluxFleft - solnMoleFluxInterface_Cell_[iCell];
                    printf(" ==  %- 14.7e\n", sumResidualTC);
                }
            }
#endif


            /*
             *  Add in the time term for species 0
             */
            for (int k = 0; k < nsp_; k++) {
                if (k != iECDMC_ && k != iPF6m_) {
                    res[indexCent_EqnStart + nodeTmpsCenter.RO_Species_Eqn_Offset  + k] += tmp;
                }
            }
            /*
             *   Add in the time term for the total electrolyte mole conservation equation
             *         note: the current problem will have this term equalling zero always.
             *               However, we put it in here for the next problem.
             *   Also we add in the source term for that equation here too.
             */
            res[indexCent_EqnStart + nodeTmpsCenter.RO_Electrolyte_Continuity] += (newStuffTC - oldStuffTC) * rdelta_t
                                  -solnMoleFluxInterface_Cell_[iCell];
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
porousLiIon_Cathode_dom1D::residEval_PreCalc(const bool doTimeDependentResid,
                                             const Epetra_Vector* soln_ptr,
                                             const Epetra_Vector* solnDot_ptr,
                                             const Epetra_Vector* solnOld_ptr,
                                             const double t,
                                             const double rdelta_t,
                                             const ResidEval_Type_Enum residType,
                                             const Solve_Type_Enum solveType)
{
    static int tmpsSetup = 0;
    if (!tmpsSetup) {
	residSetupTmps();
	tmpsSetup = 1;
    }
 
    porousFlow_dom1D::residEval_PreCalc(doTimeDependentResid, soln_ptr, solnDot_ptr,
                      solnOld_ptr, t, rdelta_t, residType, solveType);

    residType_Curr_ = residType;
    const Epetra_Vector& soln = *soln_ptr;

    maxElectrodeSubIntegrationSteps_ = 0;

    t_final_ = t;
    if (rdelta_t < 1.0E-200) {
        t_init_ = t;
    } else {
        // Needed for a correct capacity calculation 
        if (solveType == TimeDependentInitial) {
            t_init_ = t;
        } else {
            t_init_ = t - 1.0/rdelta_t;
        }
    }
  
    /*
     * Index of the first equation at the center node corresponding to the first bulk domain, which is the electrolyte
     */
    int indexCent_EqnStart;

    /*
     *  ------------------------------ LOOP OVER CELL -------------------------------------------------
     *  Loop over the number of Cells in this domain on this processor
     *  This loop is done from left to right.
     */
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

        cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
        NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	NodalVars* nodeCent = cTmps.nvCent_;

        /*
         *   Get the pointer to the NodalVars object for the center node
         */
        indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;

        /*
         * --------------------------- CALCULATE POSITION AND DELTA_X Variables -----------------------------
         * Calculate the distance between the left and center node points
         */
        /*
         * Calculate the cell width
         */
        xdelCell_Cell_[iCell] = cTmps.xCellBoundaryR_ - cTmps.xCellBoundaryL_;

	/*
	 *   Calculate the local quantities at the cell center. These will be used for electrode calculations
	 */
        for (int k = 0; k < nsp_; k++) {
            Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
        }
        Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];
        VElectrodeCent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage + 1];

        /*
         * Setup the thermo at the cell center.
         */
        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

        /*
         *  Calculate the electrode reactions.  Also update porosity.
         */
        int numSubcycles = calcElectrode();
        maxElectrodeSubIntegrationSteps_ = std::max(maxElectrodeSubIntegrationSteps_, numSubcycles);
    }
}
//=====================================================================================================================
int
porousLiIon_Cathode_dom1D::calcElectrode()
{
    Electrode* Electrode_ptr = Electrode_Cell_[cIndex_cc_];
    bool doInstanteous = false;
    int numSubcycles = 0;
    double deltaT = t_final_ - t_init_;
    if (deltaT == 0.0) {
	doInstanteous = true;
    }
    Electrode_ptr->updateState();
    //
    //  Set the init state and the init_init state
    //
    // Electrode_ptr->setInitStateFromFinal(true);

    if (doInstanteous) {
        //
        //  Write into regular array using different units
        // 
        Electrode_ptr->speciesProductionRates(&electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * cIndex_cc_]);
        int kelectron = Electrode_ptr->kSpecElectron();

        icurrInterface_Cell_[cIndex_cc_] = Faraday *
                                           electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * cIndex_cc_ + kelectron]
                                           / (crossSectionalArea_);

        Electrode_ptr->getPhaseProductionRates(&electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * cIndex_cc_],
                                               DATA_PTR(phaseMoleTransfer_));
        int sf = Electrode_ptr->solnPhaseIndex();

        solnMoleFluxInterface_Cell_[cIndex_cc_] = phaseMoleTransfer_[sf] / crossSectionalArea_;

    } else {
        //
        //     Carry out the integration from t_init_init   to t_final_final
        //
        numSubcycles = Electrode_ptr->integrate(deltaT);
#ifdef DEBUG_HKM
        if (ProblemResidEval::s_printFlagEnv > 0 && BDD_.Residual_printLvl_ > 8) {
            if (numSubcycles > 15) {
                printf("      cathode::calcElectrode(): problematic integration ( %d) : Cell = %d, counterNumberIntegrations_ = %d\n",
                       numSubcycles, cIndex_cc_, Electrode_ptr->counterNumberIntegrations_);
            }
        }
#endif
        /*
         *  Get the kmol 's produced by the electrode reactions. Units = kmol
         */
        Electrode_ptr->integratedSpeciesSourceTerm(&electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * cIndex_cc_]);


        int kelectron = Electrode_ptr->kSpecElectron();

        /*
         *  Calculate the total electron current coming from the extrinsic electrode object
         *   -> this is in amps. We make the assumption that the current is constant throught
         *      the time interval, so we divide by the deltaT of the time.
         *      We divide by the cross sectional area of the electrode.
         */
        icurrInterface_Cell_[cIndex_cc_] = Faraday * electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_  * cIndex_cc_ +
                                                                                     kelectron]
                                           / (deltaT * crossSectionalArea_);

#ifdef DEBUG_HKM
        if (residType_Curr_ == Base_ShowSolution) {
            if (cIndex_cc_ == 9) {
            }
        }
#endif
       /*
         * Get the phase mole flux from the Electrode object  units - kmol
         */
        Electrode_ptr->getIntegratedPhaseMoleTransfer(DATA_PTR(phaseMoleTransfer_));
        /*
         *  Calculate and store the electrolye molar flux units = kmol/m2/s
         */
        int sf = Electrode_ptr->solnPhaseIndex();

        solnMoleFluxInterface_Cell_[cIndex_cc_] = phaseMoleTransfer_[sf] / (crossSectionalArea_ * deltaT);
    }

    /*
     * Download the surface area from the electrode object, This is an
     * extrinsic object with units of m2
     */
    double sa[10];
    Electrode_ptr->getSurfaceAreas(sa);
    double saExternal = 0.0;
    const std::vector<bool>& isExternalSurface =  Electrode_ptr->getExternalSurfaceBooleans();
    for (int ph = 0; ph < nSurfsElectrode_; ph++) {
        if (isExternalSurface[ph]) {
            saExternal += sa[ph];
        }
    }
    /*
     * What comes out of the surface electrode object is the total surface area in the electrode.
     * Remember the electrode is an extrinsic object that should be based on a m**2 area basis.
     * To get the surface area per cross-section we need to divide by the cross sectional area of the electrode
     */
    surfaceArea_Cell_[cIndex_cc_] = saExternal / crossSectionalArea_;
    /*
     *  Calculate the current per external surface area of electrode. To do this calculation
     *  we divide the total extrinsic current by the total extrinsic surface area
     */
    icurrInterfacePerSurfaceArea_Cell_[cIndex_cc_] = icurrInterface_Cell_[cIndex_cc_]  * crossSectionalArea_ /
                                                     saExternal;

    /*
     * Compute new porosity based on new volume
     */
    //  porosity_Cell_[cIndex_cc_] = 1 -  Electrode_ptr->SolidVol() / (xdelCell_Cell_[cIndex_cc_] * crossSectionalArea_);
    //  porosity_Cell_[cIndex_cc_] = 0.64007;

    if (residType_Curr_ == Base_ShowSolution) {
        capacityDischargedPA_Cell_[cIndex_cc_] = Electrode_ptr->capacityDischarged() / crossSectionalArea_;
        depthOfDischargePA_Cell_[cIndex_cc_] = Electrode_ptr->depthOfDischarge() / crossSectionalArea_;
        capacityLeftPA_Cell_[cIndex_cc_] = Electrode_ptr->capacityLeft() / crossSectionalArea_;
        capacityPA_Cell_[cIndex_cc_]= Electrode_ptr->capacity() / crossSectionalArea_;
    }

    return numSubcycles;
}
//====================================================================================================================
void
porousLiIon_Cathode_dom1D::eval_PostSoln(
    const bool doTimeDependentResid,
    const Epetra_Vector *soln_ptr,
    const Epetra_Vector *solnDot_ptr,
    const Epetra_Vector *solnOld_ptr,
    const double t,
    const double rdelta_t)
{
    NodalVars* nodeCent = 0;
    NodalVars* nodeLeft = 0;
    NodalVars* nodeRight = 0;
    int indexCent_EqnStart, indexLeft_EqnStart, indexRight_EqnStart;
    const Epetra_Vector& soln = *soln_ptr;
    double xdelL; // Distance from the center node to the left node
    double xdelR; // Distance from the center node to the right node

    double deltaT = 0.0;
    if (rdelta_t > 0.0) {
        deltaT = 1.0 / rdelta_t;
    }

    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

        cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
        NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
        NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
        NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

	qSource_Cell_curr_[iCell] = 0.0;
	jouleHeat_lyte_Cell_curr_[iCell] = 0.0;
	jouleHeat_solid_Cell_curr_[iCell] = 0.0;

	Electrode* Electrode_ptr = Electrode_Cell_[iCell];
	/*
         *  ---------------- Get the index for the center node ---------------------------------
         *   Get the pointer to the NodalVars object for the center node
	 *   Index of the first equation in the bulk domain of center node
         */
	nodeCent = cTmps.nvCent_;
	indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;

	for (int k = 0; k < nsp_; k++) {
            Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
        }
        Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];
	// HKM fix up - technically correct
        VElectrodeCent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage + 1];

	if (iCell == 0) {
            potentialAnodic_ = Vcent_cc_;
        }
        if (iCell == NumLcCells-1) {
            potentialCathodic_ = VElectrodeCent_cc_;
        }

        /*
         *  ------------------- Get the index for the left node -----------------------------
         *    There may not be a left node if we are on the left boundary. In that case
         *    set the pointer to zero and the index to -1.
	 *    The solution index is set to the center solution index in that case as well.
         */
	nodeLeft = cTmps.nvLeft_;
	indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;

	for (int k = 0; k < nsp_; k++) {
	    Xleft_cc_[k] = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_MoleFraction_Species + k];
	}
	Vleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Voltage];
	VElectrodeLeft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Voltage + 1];

        /*
         * If we are past the first cell, then we have already done the calculation
         * for this flux at the right cell edge of the previous cell
         */
    
	xdelL = cTmps.xdelL_;
        /*
         * Calculate the distance between the right and center node points
         */
        xdelR = cTmps.xdelR_;
        /*
         * Calculate the cell width
         */
      

        /*
         * ------------------------ Get the indexes for the right node ------------------------------------
         */
	nodeRight = cTmps.nvRight_;
	indexRight_EqnStart = nodeTmpsRight.index_EqnStart;
	for (int k = 0; k < nsp_; k++) {
	    Xright_cc_[k] = soln[indexRight_EqnStart + nodeTmpsRight.Offset_MoleFraction_Species + k];
	}
	Vright_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage];
	VElectrodeRight_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage + 1];


	if (nodeLeft != 0) {
	    /*
	     *  Establish the environment at the left cell boundary
	     */
	    SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent, &(soln[indexCent_EqnStart]), 0);
	    
	    SetupTranShop(xdelL, 0);
	    /*
	     * Calculate the flux at the left boundary for each equation
	     */
	    gradV_trCurr_ = (Vcent_cc_ - Vleft_cc_) / xdelL;

	    /*
	     * Calculate the flux of species and the flux of charge
	     *   - the flux of charge must agree with the flux of species
	     */
	    icurrElectrolyte_CBL_[iCell] = 0.0;
	    for (int k = 0; k < nsp_; k++) {
		icurrElectrolyte_CBL_[iCell] += jFlux_trCurr_[k] * spCharge_[k];
	    }
	    icurrElectrolyte_CBL_[iCell] *= (Cantera::Faraday);
	    /*
	     *  Joule heating term in electrolyte
	     */
	    qSource_Cell_curr_[iCell]       += - gradV_trCurr_ * icurrElectrolyte_CBL_[iCell] * xdelL * 0.5 * deltaT;
            jouleHeat_lyte_Cell_curr_[iCell]+= - gradV_trCurr_ * icurrElectrolyte_CBL_[iCell] * xdelL * 0.5 * deltaT;
	    /*
	     *  Joule heating term in electrode
	     */
	    qSource_Cell_curr_[iCell]        += - gradVElectrode_trCurr_ * icurrElectrode_trCurr_ * xdelL * 0.5 * deltaT;
	    jouleHeat_solid_Cell_curr_[iCell]+= - gradVElectrode_trCurr_ * icurrElectrode_trCurr_ * xdelL * 0.5 * deltaT;
	}

	if (nodeRight != 0) {
	    /*
             *  Establish the environment at the right cell boundary
             */
            SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart]), nodeRight, &(soln[indexRight_EqnStart]), 1);

            SetupTranShop(xdelR, 1);

            /*
             * Calculate the flux at the right boundary for each equation
             * This is equal to
             *       Conc * Vaxial * phi
             */
	    gradV_trCurr_ = (Vright_cc_ - Vcent_cc_) / xdelR;

            /*
             * Calculate the flux of species and the flux of charge
             *   - the flux of charge must agree with the flux of species
             */
            icurrElectrolyte_CBR_[iCell] = 0.0;
            for (int k = 0; k < nsp_; k++) {
                icurrElectrolyte_CBR_[iCell] += jFlux_trCurr_[k]* spCharge_[k];
            }
            icurrElectrolyte_CBR_[iCell] *= (Cantera::Faraday);
	    /*
	     *  Joule heating term in electrolyte
	     */
	    qSource_Cell_curr_[iCell]       += - gradV_trCurr_ * icurrElectrolyte_CBR_[iCell] * xdelR * 0.5 * deltaT;
	    jouleHeat_lyte_Cell_curr_[iCell]+= - gradV_trCurr_ * icurrElectrolyte_CBR_[iCell] * xdelR * 0.5 * deltaT;
	    /*
	     *  Joule heating term in electrode
	     */
	    qSource_Cell_curr_[iCell        ]  += - gradVElectrode_trCurr_ * icurrElectrode_trCurr_ * xdelR * 0.5 * deltaT;
	    jouleHeat_solid_Cell_curr_[iCell ] += - gradVElectrode_trCurr_ * icurrElectrode_trCurr_ * xdelR * 0.5 * deltaT;
	}
        //
	// Add in the electrode contribution
        //
#ifdef DELTASHEAT_ZERO
	deltaSHeat_Cell_curr_[iCell]= 0.0;
	overPotentialHeat_Cell_curr_[iCell] = Electrode_ptr->getIntegratedThermalEnergySourceTerm_overpotential() / crossSectionalArea_;
	electrodeHeat_Cell_curr_[iCell] = overPotentialHeat_Cell_curr_[iCell];
#else
        electrodeHeat_Cell_curr_[iCell] = Electrode_ptr->getIntegratedThermalEnergySourceTerm() /  crossSectionalArea_;
        overPotentialHeat_Cell_curr_[iCell] = Electrode_ptr->getIntegratedThermalEnergySourceTerm_overpotential() / crossSectionalArea_;
        deltaSHeat_Cell_curr_[iCell]= Electrode_ptr->getIntegratedThermalEnergySourceTerm_reversibleEntropy()/ crossSectionalArea_;
#endif
	qSource_Cell_curr_[iCell] += electrodeHeat_Cell_curr_[iCell];
	qSource_Cell_accumul_[iCell] += qSource_Cell_curr_[iCell];
    }
}
//=====================================================================================================================
void
porousLiIon_Cathode_dom1D::SetupThermoShop1(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
{
    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    updateElectrolyte(nv, solnElectrolyte_Curr);
    updateElectrode();
}
//=====================================================================================================================
void
porousLiIon_Cathode_dom1D::SetupThermoShop2(const NodalVars* const nvL, const doublereal* const solnElectrolyte_CurrL,
                                            const NodalVars* const nvR, const doublereal* const solnElectrolyte_CurrR,
                                            int type)
{
    double tempL = getPointTemperature(nvL, solnElectrolyte_CurrL);
    double tempR = getPointTemperature(nvR, solnElectrolyte_CurrR);
    temp_Curr_ = 0.5 * (tempL + tempR);
    /*
     * Get the pressure
     */
    pres_Curr_ = PressureReference_;

    size_t indexMFL = nvL->indexBulkDomainVar0(MoleFraction_Species);
    size_t indexMFR = nvR->indexBulkDomainVar0(MoleFraction_Species);

    mfElectrolyte_Soln_Curr_[0] = 0.5 * (solnElectrolyte_CurrL[indexMFL] +solnElectrolyte_CurrR[indexMFR]);
    mfElectrolyte_Soln_Curr_[1] = 0.5 * (solnElectrolyte_CurrL[indexMFL+1] +solnElectrolyte_CurrR[indexMFR+1]);
    mfElectrolyte_Soln_Curr_[2] = 0.5 * (solnElectrolyte_CurrL[indexMFL+2] +solnElectrolyte_CurrR[indexMFR+2]);
    double mf0 = std::max(mfElectrolyte_Soln_Curr_[0], 0.0);
    double mf1b = std::max(mfElectrolyte_Soln_Curr_[1], 0.0);
    double mf2b = std::max(mfElectrolyte_Soln_Curr_[2], 0.0);
    double mf1 = mf1b;
    double mf2 = mf2b;
    if (mf1b != mf2b) {
        mf1 = 0.5 * (mf1b + mf2b);
        mf2 = 0.5 * (mf1b + mf2b);
    }
    double tmp = mf0 + mf1 + mf2;

    mfElectrolyte_Thermo_Curr_[0] = mf0 / tmp;
    mfElectrolyte_Thermo_Curr_[1] = mf1 / tmp;
    mfElectrolyte_Thermo_Curr_[2] = mf2 / tmp;

    size_t indexVS = nvL->indexBulkDomainVar0(Voltage);
    double phiElectrolyteL = solnElectrolyte_CurrL[indexVS];
    indexVS = nvR->indexBulkDomainVar0(Voltage);
    double phiElectrolyteR = solnElectrolyte_CurrR[indexVS];
    phiElectrolyte_Curr_ = 0.5 * (phiElectrolyteL + phiElectrolyteR);

    if (type == 0) {
        porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ - 1] + porosity_Cell_[cIndex_cc_]);
    } else {
        porosity_Curr_ = 0.5 * (porosity_Cell_[cIndex_cc_ + 1] + porosity_Cell_[cIndex_cc_]);
    }
    /*
     *  Set the ThermoPhase states
     */
    ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);
    ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);
    //
    // Calculate the total concentration of the electrolyte kmol m-3 and store into concTot_Curr_
    //
    concTot_Curr_ = ionicLiquid_->molarDensity();
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
porousLiIon_Cathode_dom1D::updateElectrolyte(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
{ 
    /*
     * Get the temperature: Check to see if the temperature is in the solution vector.
     *   If it is not, then use the reference temperature
     */
    temp_Curr_ = getPointTemperature(nv, solnElectrolyte_Curr);
    /*
     * Get the pressure
     */
    pres_Curr_ = PressureReference_;

    getMFElectrolyte_soln(nv, solnElectrolyte_Curr);
    getVoltages(nv, solnElectrolyte_Curr);

    ionicLiquid_->setState_TPX(temp_Curr_, pres_Curr_, &mfElectrolyte_Thermo_Curr_[0]);
    ionicLiquid_->setElectricPotential(phiElectrolyte_Curr_);

    // Calculate the total concentration of the electrolyte kmol m-3.
    concTot_Curr_ = ionicLiquid_->molarDensity();
}
//=====================================================================================================================
void
porousLiIon_Cathode_dom1D::updateElectrode()
{
    Electrode* Electrode_ptr = Electrode_Cell_[cIndex_cc_];
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
porousLiIon_Cathode_dom1D::getVoltages(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    size_t indexVS = nv->indexBulkDomainVar0(Voltage);
    phiElectrolyte_Curr_ = solnElectrolyte_Curr[indexVS];
    phiElectrode_Curr_ = solnElectrolyte_Curr[indexVS + 1];
}
//=====================================================================================================================
void
porousLiIon_Cathode_dom1D::getMFElectrolyte_soln(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    size_t indexMF = nv->indexBulkDomainVar0(MoleFraction_Species);

    mfElectrolyte_Soln_Curr_[0] = solnElectrolyte_Curr[indexMF];
    mfElectrolyte_Soln_Curr_[1] = solnElectrolyte_Curr[indexMF + 1];
    mfElectrolyte_Soln_Curr_[2] = solnElectrolyte_Curr[indexMF + 2];
    double mf0  = std::max(mfElectrolyte_Soln_Curr_[0], 0.0);
    double mf1b = std::max(mfElectrolyte_Soln_Curr_[1], 0.0);
    double mf2b = std::max(mfElectrolyte_Soln_Curr_[2], 0.0);
    double mf1 = mf1b;
    double mf2 = mf2b;
    if (mf1b != mf2b) {
        mf1 = 0.5 * (mf1b + mf2b);
        mf2 = 0.5 * (mf1b + mf2b);
    }
    double tmp = mf0 + mf1 + mf2;

    mfElectrolyte_Thermo_Curr_[0] = mf0 / tmp;
    mfElectrolyte_Thermo_Curr_[1] = mf1 / tmp;
    mfElectrolyte_Thermo_Curr_[2] = mf2 / tmp;
}
//=====================================================================================================================
/*
 *  We assume that setupthermoshop1 has been called.
 *  This calculates the Heat capacity per cross-sectional area (Joules/K m2)
 *
 *   crossSectionalArea_
 */
double
porousLiIon_Cathode_dom1D::getCellHeatCapacity(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    Electrode* Electrode_ptr = Electrode_Cell_[cIndex_cc_];
    double cpMolar = ionicLiquid_->cp_mole();
    double lyteVol = porosity_Curr_ * xdelCell_Cell_[cIndex_cc_];
    double cpLyte =  lyteVol * concTot_Curr_ * cpMolar;
    double cpSolidTotal = Electrode_ptr->SolidHeatCapacityCV() ;
    double cpSolid =  cpSolidTotal / crossSectionalArea_;  
    return (cpSolid + cpLyte);
}
//=====================================================================================================================
void
porousLiIon_Cathode_dom1D::SetupTranShop(const double xdel, const int type)
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
    Tortuosity tort(bruggemannExp);
    for (int k = 0; k < nsp_; k++) {
        Vdiff_trCurr_[k] *= tort.tortuosityFactor(porosity_Curr_);
    }
    //Convert from diffusion velocity to diffusion flux
    for (int k = 0; k < nsp_; k++) {
        jFlux_trCurr_[k] = mfElectrolyte_Soln_Curr_[k] * concTot_Curr_ * Vdiff_trCurr_[k];
    }

    double volFSolid = (1.0 - porosity_Curr_);
    icurrElectrode_trCurr_ = -conductivityElectrode_ * pow(volFSolid, 1.5) * gradVElectrode_trCurr_;
}
//=====================================================================================================================
// Saving the solution on the domain in an xml node.
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
porousLiIon_Cathode_dom1D::saveDomain(Cantera::XML_Node& oNode,
                                      const Epetra_Vector* soln_GLALL_ptr,
                                      const Epetra_Vector* solnDot_GLALL_ptr,
                                      const double t,
                                      bool duplicateOnAllProcs)
{
    // get the NodeVars object pertaining to this global node
    GlobalIndices* gi = LI_ptr_->GI_ptr_;

    // Add a child for this domain
    Cantera::XML_Node& bdom = oNode.addChild("domain");

    // Number of equations per node
    int numEquationsPerNode = BDD_.NumEquationsPerNode;

    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType>& variableNameList = BDD_.VariableNameList;

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
        NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
        varContig[i] = nv->x0NodePos();
    }
    ctml::addNamedFloatArray(gv, "X0", varContig.size(), &(varContig[0]), "m", "length");

    for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
        VarType vt = variableNameList[iVar];
        i = 0;
        std::string nmm = vt.VariableName(200);
        for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
	    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
            int istart = nv->EqnStart_GbEqnIndex;
            varContig[i] = (*soln_GLALL_ptr)[istart + offset];
        }
        ctml::addNamedFloatArray(gv, nmm, varContig.size(), &(varContig[0]), "kmol/m3", "concentration");
    }

    if (PS_ptr->doHeatSourceTracking_) {
        std::string nmm = "qSource_Cell_curr_";
        ctml::addNamedFloatArray(gv, nmm, numNodes, &(qSource_Cell_curr_[0]), "Joule/s/m2", "");
    }

    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
        int iCell = iGbNode - firstGbNode;
        Electrode* ee = Electrode_Cell_[iCell];	
        ee->writeTimeStateFinal_toXML(bdom);
    }
}
//====================================================================================================================
//
//  This treatment assumes that the problem size stays constant. If this is not the case, the routine will
//  error exit. If we need to change the problem size, then we will need to reinitialize a lot more that just
//  the solution vector. This is best done after we have put in grid refinement.
//
//  Also, we assume that the number of variables stays the same. This may be fiddled with sooner than the number
//  of grid points. There are probably some interesting possibilities here with just initializing a subset of 
//  variables. However, right now, if the number of variables aren't equal and in the same order, the
//  routine will error exit.
//
//  Also we don't consider any interpolation in between time steps. We enter with an XML_Node specific
//  to a particular time step. And then populate the solution vector with the saved solution.
//
//  MP Implementation
//     We've punted on this for now. The basic strategy will be to identity which of the three situations
//     we are currently doing:
//     
//          1) global data into local-processor data structure
//          2) global data into global data structure
//          3) local-processor data into local-processor data structure
// 
//     We are currently set up for #1. However, that may change. Even #1 will fail
//     
void
porousLiIon_Cathode_dom1D::readDomain(const Cantera::XML_Node& SimulationNode,
				      Epetra_Vector * const soln_GLALL_ptr, Epetra_Vector * const solnDot_GLALL_ptr)
{
   // get the NodeVars object pertaining to this global node
    GlobalIndices *gi = LI_ptr_->GI_ptr_;

    string ids = id();
    Cantera::XML_Node *domainNode_ptr = SimulationNode.findNameID("domain", ids);

    // Number of equations per node
    int numEquationsPerNode = BDD_.NumEquationsPerNode;

    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType> &variableNameList = BDD_.VariableNameList;

    //! First global node of this bulk domain
    int firstGbNode = BDD_.FirstGbNode;

    //! Last Global node of this bulk domain
    int lastGbNode = BDD_.LastGbNode;
    int numNodes = lastGbNode - firstGbNode + 1;

    string iidd      = (*domainNode_ptr)["id"];
    string s_points  = (*domainNode_ptr)["points"];
    int points = atoi(s_points.c_str());
    if (points != numNodes) {
        printf("we have an unequal number of points\n");
        exit(-1);
    }
    string ttype    = (*domainNode_ptr)["type"];
    string snumVar  = (*domainNode_ptr)["numVariables"];
    int numVar = atoi(snumVar.c_str());
    if (numVar != numEquationsPerNode) {
       printf("we have an unequal number of equations\n");
       exit(-1);
    }
  //
    //  Go get the grid data XML node and read it in
    //
    const Cantera::XML_Node* gd_ptr = (*domainNode_ptr).findByName("grid_data");

    std::vector<double> varContig(numNodes);
    ctml::getFloatArray(*gd_ptr, varContig, true, "", "X0");
    int i = 0;
    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
      NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
      varContig[i] = nv->x0NodePos();
      nv->setupInitialNodePosition(varContig[i], 0.0);
    }

    for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
       VarType vt = variableNameList[iVar];
       i = 0;
       std::string nmm = vt.VariableName(200);
       ctml::getFloatArray(*gd_ptr, varContig, true, "", nmm);
       for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
          NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
          size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
          int istart = nv->EqnStart_GbEqnIndex;
          (*soln_GLALL_ptr)[istart + offset] =  varContig[i];
       }
    }

    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
        int iCell = iGbNode - firstGbNode;

        Electrode* ee = Electrode_Cell_[iCell];
        int cellNum = ee->electrodeCellNumber_;
        // ee->writeStateToXML();

        XML_Node *xCell = domainNode_ptr->findByAttr("cellNumber", int2str(cellNum), 1);
        if (!xCell) {
            throw m1d_Error("porousLiIon_Cathode_dom1D::readDomain() ERROR",
                            "Can't find cell number " + int2str(cellNum));
        }
        //  Read the state into the electrode object
        ee->loadTimeStateFinal(*xCell);
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
porousLiIon_Cathode_dom1D::writeSolutionTecplotHeader()
{
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid ; //only proc 0 should write

    if (doWrite) {

        //open tecplot file
        FILE* ofp;
        string sss = id();
        char filename[20];
        sprintf(filename,"%s%s",sss.c_str(),".dat");
        ofp = fopen(filename, "w");

        //write title and variable list
        fprintf(ofp, "TITLE = \"Solution on Domain %s\"\n",sss.c_str());

        // Number of equations per node
        int numVar = BDD_.NumEquationsPerNode;
        // Vector containing the variable names as they appear in the solution vector
        std::vector<VarType>& variableNameList = BDD_.VariableNameList;
        //! First global node of this bulk domain

        fprintf(ofp, "VARIABLES = ");
        fprintf(ofp, "\"x [m]\"  \n");

        for (int k = 0; k < numVar; k++) {
            VarType& vt = variableNameList[k];
            string name = vt.VariableName(15);
            fprintf(ofp, "\"%s\" \t", name.c_str());
        }
        fprintf(ofp, "\n");

        /////////////////////////////////////////
        //// end BulkDomain1D section
        //// start application specific output
        /////////////////////////////////////////

        // print potentials within each control volume
        fprintf(ofp, "\"Potential Drop Electrode-Electrolyte [V]\" \t");
        for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++) {
            fprintf(ofp, "\"Surf_%d Equilibrium Potential Drop Electrode-Electrolyte [V]\" \t", jSurf);
            fprintf(ofp, "\"Overpotential_%d [V]\" \t", jSurf);
        }
        fprintf(ofp, "\"Current Source [A/m^2]\" \t");
        fprintf(ofp, "\"Li+ flux [kmol/m^2/s]\" \t");
        fprintf(ofp, "\"Electrolyte Current [A/m^2]\" \t");
        fprintf(ofp, "\"Electrode Current [A/m^2]\" \t");
        fprintf(ofp, "\n");

        //  Print depth of discharge for each control volume
        fprintf(ofp, "\"Capacity discharged [A-s]\" \t");
        fprintf(ofp, "\"Depth of discharge [-]\" \t");
        fprintf(ofp, "\"Capacity remaining [A-s]\" \t");
        fprintf(ofp, "\"Initial capacity [A-s]\" \t");

        // print porosity, specific surface area, thickness for each control volume
        fprintf(ofp, "\"Porosity []\" \t");
        fprintf(ofp, "\"Surface Area [m2/m2] (area per area)\" \t");
        fprintf(ofp, "\"Specific current (per particle area) [A/m^2]\" \t");
        fprintf(ofp, "\"Control volume thickness [m]\" \t");
        fprintf(ofp, "\n");

        //set up Electrode objects so that we can plot moles of active materials
        Electrode* ee0 = Electrode_Cell_[0];
        std::vector<double> spmoles(nSpeciesElectrode_, 0.0);
        int solnPhase = ee0->solnPhaseIndex();
        int metalPhase =  ee0->metalPhaseIndex();
        int numVolPhasesE =  ee0->nVolPhases();

        // print mole numbers of active materials
        for (int vph = 0; vph < numVolPhasesE; vph++) {
            ThermoPhase* tp = &(ee0->volPhase(vph));
            int iph = ee0->getGlobalPhaseIndex(tp);
            if (iph == metalPhase || iph == solnPhase) {

            } else {
                int nspPhase = tp->nSpecies();
                int kStart =  ee0->getGlobalSpeciesIndex(iph, 0);
                for (int k = 0; k < nspPhase; k++) {
                    string sss = ee0->speciesName(kStart + k);
                    fprintf(ofp, "\"Moles %s [mol/m^2] (per electrode area)\" \t", sss.c_str());
                }
            }
        }
	//print thermal source terms
	// check dimensions!!
        fprintf(ofp, "\"qHeat_accum [J/m3]\" \t");
        fprintf(ofp, "\"qHeat_step [W/m3]\" \t");
        fprintf(ofp, "\"Joule_Lyte [W/m3]\" \t");
        fprintf(ofp, "\"Joule_Solid [W/m3]\" \t");
        fprintf(ofp, "\"electrodHeat [W/m3]\" \t");
        fprintf(ofp, "\"OverpotHeat [W/m3]\" \t");
        fprintf(ofp, "\"deltaSHeat [W/m3]\" \t");
        fprintf(ofp, "\n");

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
porousLiIon_Cathode_dom1D::writeSolutionTecplot(const Epetra_Vector* soln_GlAll_ptr,
                                                const Epetra_Vector* solnDot_GlAll_ptr,
                                                const double t)
{
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid ; //only proc 0 should write

    if (doWrite) {

        double deltaTime = t_final_ - t_init_;

        // get the NodeVars object pertaining to this global node
        GlobalIndices* gi = LI_ptr_->GI_ptr_;
        // Number of equations per node
        int numVar = BDD_.NumEquationsPerNode;
        //! First global node of this bulk domain
        int firstGbNode = BDD_.FirstGbNode;
        //! Last Global node of this bulk domain
        int lastGbNode = BDD_.LastGbNode;
        int numNodes = lastGbNode - firstGbNode + 1;
        std::vector<VarType>& variableNameList = BDD_.VariableNameList;


        //open tecplot file
        FILE* ofp;
        string sss = id();
        char filename[20];
        sprintf(filename,"%s%s",sss.c_str(),".dat");
        ofp = fopen(filename, "a");

        fprintf(ofp, "ZONE T = \"t = %g [s]\" I = %d SOLUTIONTIME = %19.13E\n", t, numNodes, t);

        for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++) {
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];

            //x-position
            fprintf(ofp, "%g \t", nv->xNodePos());

            for (int iVar = 0; iVar < numVar; iVar++) {
                VarType vt = variableNameList[iVar];
               	size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
                int istart = nv->EqnStart_GbEqnIndex;
                fprintf(ofp, "%g \t", (*soln_GlAll_ptr)[istart + offset]);
            }
            fprintf(ofp, "\n");

            /////////////////////////////////////////
            //// end BulkDomain1D section
            /////////////////////////////////////////

            int iCell = iGbNode - firstGbNode;
            // surface reaction data
            fprintf(ofp, "%g \t", deltaV_Cell_[iCell]);
            for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++) {
                fprintf(ofp, "%g \t", Ess_Surf_Cell_[nSurfsElectrode_ * iCell + jSurf]);
                fprintf(ofp, "%g \t", overpotential_Surf_Cell_[nSurfsElectrode_ * iCell + jSurf]);
            }
            fprintf(ofp, "%g \t", icurrRxn_Cell_[iCell]);
            fprintf(ofp, "%g \t", LiFlux_Cell_[iCell]);

            // current flux -- average of left and right fluxes
            fprintf(ofp, "%g \t", 0.5 * (icurrElectrolyte_CBL_[iCell] + icurrElectrolyte_CBR_[iCell]));
            fprintf(ofp, "%g \t", 0.5 * (icurrElectrode_CBL_[iCell] + icurrElectrode_CBR_[iCell]));

            fprintf(ofp, "\n");

            // capacity of individual control volumes
            fprintf(ofp, "%g \t", capacityDischargedPA_Cell_[iCell] * crossSectionalArea_);
            fprintf(ofp, "%g \t", depthOfDischargePA_Cell_[iCell]);
            fprintf(ofp, "%g \t", capacityLeftPA_Cell_[iCell] * crossSectionalArea_);
            fprintf(ofp, "%g \t", capacityPA_Cell_[iCell] * crossSectionalArea_);

            // print porosity, specific surface area, thickness for each control volume
            fprintf(ofp, "%g \t", porosity_Cell_[iCell]);
            fprintf(ofp, "%g \t", surfaceArea_Cell_[iCell]);
            fprintf(ofp, "%g \t", icurrInterfacePerSurfaceArea_Cell_[iCell]);
            fprintf(ofp, "%g \t", xdelCell_Cell_[iCell]);
            fprintf(ofp, "\n");

            // print the properties of the electrode particles
            Electrode* ee = Electrode_Cell_[iCell];
            std::vector<double> spmoles(nSpeciesElectrode_, 0.0);
            int solnPhase = ee->solnPhaseIndex();
            int metalPhase =  ee->metalPhaseIndex();
            int numVolPhasesE =  ee->nVolPhases();

            // print mole numbers of active materials
            ee->getMoleNumSpecies(DATA_PTR(spmoles));
            for (int vph = 0; vph < numVolPhasesE; vph++) {
                ThermoPhase* tp = &(ee->volPhase(vph));
                int iph = ee->getGlobalPhaseIndex(tp);
                if (iph == metalPhase || iph == solnPhase) {

                } else {
                    int nspPhase = tp->nSpecies();
                    int kStart =  ee->getGlobalSpeciesIndex(iph, 0);
                    for (int k = 0; k < nspPhase; k++) {
                        fprintf(ofp, "%g \t", spmoles[kStart + k] / crossSectionalArea_);
                    }
                }
            }
	    // print thermal source terms
	    fprintf(ofp, "%g \t", qSource_Cell_accumul_[iCell] / xdelCell_Cell_[iCell] );
	    if ( deltaTime > 1e-80 ) {
	      fprintf(ofp, "%g \t", qSource_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	      fprintf(ofp, "%g \t", jouleHeat_lyte_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	      fprintf(ofp, "%g \t", jouleHeat_solid_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	      fprintf(ofp, "%g \t", electrodeHeat_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	      fprintf(ofp, "%g \t", overPotentialHeat_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	      fprintf(ofp, "%g \t", deltaSHeat_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime );
	    } else {
	      fprintf(ofp, "0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t 0.0 \t " );
	    }
            fprintf(ofp, "\n");

        }
        fclose(ofp);

#undef DAKOTAOUT
#ifdef  DAKOTAOUT
        double firstOutputTime = 0.1;
        static double outputTime = 0.1;
        if (t > outputTime) {
            std::ofstream dfp;
            if (outputTime > firstOutputTime) {
                dfp.open("results.out", std::ios_base::app);
            } else {
                dfp.open("results.out", std::ios_base::out);
            }
            int indexVS = BDD_.VariableIndexStart_VarName[Voltage];
            double potentialE_left = (*soln_GlAll_ptr)[indexVS];
            double lithium_left = (*soln_GlAll_ptr)[1];

            NodalVars* nv = gi->NodalVars_GbNode[lastGbNode-1];
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
                                            * log(lithium_left / lithium_right) ;

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
porousLiIon_Cathode_dom1D::showSolution(const Epetra_Vector* soln_GlAll_ptr,
                                        const Epetra_Vector* solnDot_GlAll_ptr,
                                        const Epetra_Vector* soln_ptr,
                                        const Epetra_Vector* solnDot_ptr,
                                        const Epetra_Vector* solnOld_ptr,
                                        const Epetra_Vector_Owned* residInternal_ptr,
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
    std::vector<VarType>& variableNameList = BDD_.VariableNameList;
    int iBlock;
    int iGbNode;
    int n;
    int nPoints = BDD_.LastGbNode - BDD_.FirstGbNode + 1;

    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind1 = indent.c_str();
    char ind[120];
    strcpy(ind, ind1);
    ind[118] = '\0';
    doublereal v;
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
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
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
                doublereal x = nv->xNodePos();
                ss.print0("\n%s    %-10.4E ", ind, x);
                int istart = nv->EqnStart_GbEqnIndex;
                for (n = 0; n < 5; n++) {
                    int ivar = iBlock * 5 + n;
                    VarType vt = variableNameList[ivar];
		    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
		    if (offset == npos) {
                        throw m1d_Error("porousLiIon_Separator_dom1D::showSolution()", "cant find a variable");
                    }
		    v = (*soln_GlAll_ptr)[istart + offset];

                    ss.print0(" %-11.5E ", v);
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
            double oldVaxial = 0.0;
            double newVaxial = 0.0;

            for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
                doublereal x = nv->xNodePos();
                ss.print0("%s    %-10.4E ", ind, x);
                int istart = nv->EqnStart_GbEqnIndex;
                for (n = 0; n < nrem; n++) {
                    int ivar = iBlock * 5 + n;
                    VarType vt = variableNameList[ivar];
		    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
		    if (offset == npos) {
                        throw m1d_Error("porousLiIon_Separator_dom1D::showSolution()", "cant find a variable");
                    }
		    v = (*soln_GlAll_ptr)[istart + offset];
                    if (vt.VariableType == Velocity_Axial) {
                        newVaxial = v;
                        if (iGbNode ==  BDD_.LastGbNode) {
                            v = newVaxial;
                        } else {
                            v = 0.5 * (oldVaxial + newVaxial);
                        }
                        oldVaxial = newVaxial;
                    }
                    ss.print0(" %-11.5E ", v);
                }
                ss.print0("\n");
            }
        }
        drawline(indentSpaces, 80);
    }

    if (do0Write) {
        // -----------------------------------------------------------------------------------------------------------------
        // --             PRINT REACTION RATES within the cell --
        // -----------------------------------------------------------------------------------------------------------------
        ss.print0("\n");
        drawline0(indentSpaces, 80);
        ss.print0("%s        z       Delta_V        Ess    Overpotential icurrRxnCell", ind);
        ss.print0("\n");
        drawline0(indentSpaces, 80);
    }
    doublereal x;
    int iCell;
    for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
            iCell = iGbNode - BDD_.FirstGbNode;
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            x = nv->xNodePos();
            ss.print0("%s    %-10.4E ", ind, x);
            ss.print0("%11.4E ", deltaV_Cell_[iCell]);
            ss.print0("%11.4E ", Ess_Surf_Cell_[nSurfsElectrode_ * iCell]);
            ss.print0("%11.4E ", overpotential_Surf_Cell_[nSurfsElectrode_ * iCell]);
            ss.print0("%11.4E ", icurrRxn_Cell_[iCell]);
            ss.print0("\n");
        }
        print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    }


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
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
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
        ss.print0("%s        z      SurfaceArea  currPerSA   DeltaZ", ind);
        ss.print0("\n");
        drawline0(indentSpaces, 80);
    }
    for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
            iCell = iGbNode - BDD_.FirstGbNode;
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            x = nv->xNodePos();
            ss.print0("%s    %-10.4E ", ind, x);
            ss.print0("%11.4E ", surfaceArea_Cell_[iCell]);
            ss.print0("%11.4E ", surfaceArea_Cell_[iCell] / xdelCell_Cell_[iCell]);
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
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
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

    Electrode* ee0 = Electrode_Cell_[0];
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
            ThermoPhase* tp = &(ee0->volPhase(vph));
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
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            x = nv->xNodePos();
            ss.print0("%s    %-10.4E ", ind, x);
            Electrode* ee = Electrode_Cell_[iCell];
            ee->getMoleNumSpecies(DATA_PTR(spmoles));
            for (int vph = 0; vph < numVolPhasesE; vph++) {
                ThermoPhase* tp = &(ee->volPhase(vph));
                int iph = ee->getGlobalPhaseIndex(tp);
                if (iph == metalPhase || iph == solnPhase) {

                } else {
                    int nspPhase = tp->nSpecies();
                    int kStart =  ee->getGlobalSpeciesIndex(iph, 0);
                    for (int k = 0; k < nspPhase; k++) {
                        ss.print0("% -15.6E ", spmoles[kStart + k] / crossSectionalArea_);
                    }

                }
            }
            ss.print0("\n");
        }
        print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
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
    NodalVars* nvl;
    NodalVars* nvr;
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


    if (PS_ptr->doHeatSourceTracking_) {
	if (do0Write) {
	    // ----------------------------------------------------
	    // --         PRINT HEAT GENERATION TERMS for DOMAIN --
	    // ----------------------------------------------------

	    ss.print0("\n");
	    drawline0(indentSpaces, 100);
	    ss.print0("%s        z     qHeat_step    qHeat_accum Joule_Lyte Joule_Solid electrodHeat OverpotHeat deltaSHeat", ind);
	    ss.print0("\n");
	    drawline0(indentSpaces, 100);
	}
	//NodalVars* nvl;
	//NodalVars* nvr;
	for (iGbNode = BDD_.FirstGbNode; iGbNode <= BDD_.LastGbNode; iGbNode++) {
	    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
	    if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
		iCell = iGbNode - BDD_.FirstGbNode;
		NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
		x = nv->xNodePos();
		ss.print0("%s    %-10.4E ", ind, x);
		//Electrode* ee = Electrode_Cell_[iCell];

		ss.print0("% -11.4E ", qSource_Cell_curr_[iCell]);

		ss.print0("% -11.4E ",	qSource_Cell_accumul_[iCell]);
                ss.print0("% -11.4E ",  jouleHeat_lyte_Cell_curr_[iCell]);
		ss.print0("% -11.4E ",  jouleHeat_solid_Cell_curr_[iCell]);

 		ss.print0("% -11.4E ",  electrodeHeat_Cell_curr_[iCell]);
  		ss.print0("% -11.4E ",  overPotentialHeat_Cell_curr_[iCell]);
  		ss.print0("% -11.4E ",  deltaSHeat_Cell_curr_[iCell]);

		ss.print0("\n");
	    }
	    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
	}
    }
}
//=====================================================================================================================
// Set the underlying state of the system from the solution vector
/*
 *   Note this is an important routine for the speed of the solution.
 *   It would be great if we could supply just exactly what is changing here.
 *   This routine is always called at the beginning of the residual evaluation process.
 *
 *   This is a natural place to put any precalculations of nodal quantities that
 *   may be needed by the residual before its calculation.
 *
 *   Also, this routine is called with rdelta_t = 0. This implies that a step isn't being taken. However, the
 *   the initial conditions must be propagated.
 *
 * @param doTimeDependentResid
 * @param soln
 * @param solnDot
 * @param t
 * @param delta_t delta t. If zero then delta_t equals 0.
 * @param t_old
 */
void
porousLiIon_Cathode_dom1D::setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted *soln_ptr, 
						const Epetra_Vector_Ghosted *solnDot,
						const double t, const double delta_t, const double t_old)
{
    bool doAll = true;
    //bool doInit = false;
    //bool doFinal = true; // This is always true as you can only set the final Electrode object
    size_t  indexCent_EqnStart;
    if (doTimeDependentResid) {
	if (delta_t > 0.0) {
            doAll = false;
	    if (t <= t_old) {
		//doInit = true;
	    }
	}
    }
    const Epetra_Vector_Ghosted& soln = *soln_ptr;

    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;
	cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	/*
         *  ---------------- Get the index for the center node ---------------------------------
         *   Get the pointer to the NodalVars object for the center node
         *   Index of the first equation in the bulk domain of center node
         */
        NodalVars* nodeCent = cTmps.nvCent_;
        indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
       
	const double *solnCentStart = &(soln[indexCent_EqnStart]);
    
	/*
	 * Get the temperature: Check to see if the temperature is in the solution vector.
	 *   If it is not, then use the reference temperature
	 */
	temp_Curr_ = getPointTemperature(nodeCent, solnCentStart);
	/*
	 * Get the pressure
	 */
	pres_Curr_ = getPointPressure(nodeCent, solnCentStart);
	//
	//  fill in mfElectrolyte_Soln_Curr[]  mfElectrolyte_Thermo_Curr_[]
	//
	getMFElectrolyte_soln(nodeCent, solnCentStart);
	//
	//  Fill in phiElectroyte_Curr_ and phiElectrode_Curr_
	//
	getVoltages(nodeCent, solnCentStart);
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
//=====================================================================================================================
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
 *
 */
void
porousLiIon_Cathode_dom1D::initialConditions(const bool doTimeDependentResid,
                                             Epetra_Vector* soln_ptr,
                                             Epetra_Vector* solnDot,
                                             const double t,
                                             const double delta_t)
{
    Epetra_Vector& soln = *soln_ptr;

    int index_CentLcNode;
    NodalVars* nodeCent = 0;
    int indexCent_EqnStart;

    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        cIndex_cc_ = iCell;

        index_CentLcNode = Index_DiagLcNode_LCO[iCell];
        // pointer to the NodalVars object for the center node
        nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        // Index of the first equation in the bulk domain of center node
        indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];

        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
	size_t iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        size_t iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
        size_t iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
	size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
	size_t iVar_Pressure = nodeCent->indexBulkDomainVar0((size_t) Pressure_Axial);
        size_t iVar_Voltage_ED = iVar_Voltage + 1;
	//
	// Find the start of the solution at the current node
	//
	const double *solnCentStart = &(soln[indexCent_EqnStart]);

        soln[indexCent_EqnStart + iVAR_Vaxial] = 0.0;
	//
	// Set the temperature if it is part of the solution vector
	//	
	temp_Curr_ = PSinput.TemperatureReference_;
        if (iVar_Temperature != npos) {
            soln[indexCent_EqnStart + iVar_Temperature] = PSinput.TemperatureReference_;
        }
	//
	// Set the pressure if it is part of the solution vector
	//
	pres_Curr_ = PSinput.PressureReference_;
        if (iVar_Pressure != npos) {
            soln[indexCent_EqnStart + iVar_Pressure] = PSinput.PressureReference_;
        }

        /*
         * Get initial mole fractions from PSinput
         */
        int igECDMC  = PSinput.PhaseList_->globalSpeciesIndex("ECDMC");
        if (igECDMC < 0) {
            throw CanteraError("confused", "confused");
        }
        int igLip = PSinput.PhaseList_->globalSpeciesIndex("Li+");
        if (igLip < 0) {
            throw CanteraError("confused", "confused");
        }
        int igPF6m = PSinput.PhaseList_->globalSpeciesIndex("PF6-");
        if (igPF6m < 0) {
            throw CanteraError("confused", "confused");
        }

        soln[indexCent_EqnStart + iVar_Species + iECDMC_] = PSinput.electrolyteMoleFracs_[igECDMC];
        soln[indexCent_EqnStart + iVar_Species + iLip_  ] = PSinput.electrolyteMoleFracs_[igLip];
        soln[indexCent_EqnStart + iVar_Species + iPF6m_ ] = PSinput.electrolyteMoleFracs_[igPF6m];


        soln[indexCent_EqnStart + iVar_Voltage] = -0.07;

        //double icurr = PSinput.icurrDischargeSpecified_;
        double volt = PSinput.CathodeVoltageSpecified_;
        soln[indexCent_EqnStart + iVar_Voltage_ED] = volt;
	//
        //  Fill in phiElectroyte_Curr_ and phiElectrode_Curr_
        //
        getVoltages(nodeCent, solnCentStart);
	//
        //  fill in mfElectrolyte_Soln_Curr[]  mfElectrolyte_Thermo_Curr_[]
        //
        getMFElectrolyte_soln(nodeCent, solnCentStart);
	//
        // Electrode object will beupdated and we will compute the porosity
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
        // For cases where we aren't in time dependent mode set all states
        //
	ee->setInitStateFromFinal(true);
	ee->setFinalFinalStateFromFinal();
        //
        // Calculate the solid volume of the electrode
        //
	double solidVolCell = ee->SolidVol();
        double porosity = 1.0 - solidVolCell / (xdelCell_Cell_[iCell] * crossSectionalArea_);
        if (porosity <= 0.0) {
            throw m1d_Error("porousLiIon_Cathode_dom1D::initialConditions() ERROR",
                            "Calculated porosity, " + fp2str(porosity) + ", is less than zero\n"
                            "           There is an error in the setup of the anode");
        }    
        //
        // update porosity as computed from electrode input
        //
        porosity_Cell_[iCell] = porosity;
    }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void porousLiIon_Cathode_dom1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
                                              Epetra_Vector_Ghosted& atolVector,
                                              const Epetra_Vector_Ghosted* const atolV)
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
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
	int iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        int iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
        int iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        int iVar_Voltage_ED = iVar_Voltage + 1;

        /*
         * Set the atol value for the axial velocity
         *   arithmetically scaled -> so this is a characteristic value
         */
        double vax = soln[indexCent_EqnStart + iVAR_Vaxial];
        atolVector[indexCent_EqnStart + iVAR_Vaxial] = std::max(1.0E-4, 1.0E-1 * vax);

        /*
         * Set atol values for the species mole fractions
         */
        double val = atolDefault * 1.0E-2;
        if (val < 1.0E-12) {
            val = 1.0E-12;
        }
        for (int k = 0; k < nsp_; k++) {
            atolVector[indexCent_EqnStart + iVar_Species + k] = val;
        }
        /*
         * Set the atol value for the electrolyte voltage
         *      arithmetically scaled.-> so this is a characteristic value
         *         1 kcal gmol-1 = 0.05 volts
         */
        atolVector[indexCent_EqnStart + iVar_Voltage] = 0.05;
        atolVector[indexCent_EqnStart + iVar_Voltage_ED] = 0.05;
    }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void porousLiIon_Cathode_dom1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted& soln,
                                                      const Epetra_Vector_Ghosted& solnDot,
                                                      Epetra_Vector_Ghosted& atolVector,
                                                      const Epetra_Vector_Ghosted* const atolV)
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
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
	int iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        int iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
        int iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        int iVar_Voltage_ED = iVar_Voltage + 1;

        /*
         * Set the atol value for the axial velocity
         *   arithmetically scaled -> so this is a characteristic value
         */
        double vax = soln[indexCent_EqnStart + iVAR_Vaxial];
        atolVector[indexCent_EqnStart + iVAR_Vaxial] = std::max(1.0E-4, 1.0E-1 * vax);

        /*
         * Set atol values for the species mole fractions time derivatives
         */
        double val = atolDefault * 1.0E-1;
        if (val < 1.0E-6) {
            val = 1.0E-6;
        }
        for (int k = 0; k < nsp_; k++) {
            atolVector[indexCent_EqnStart + iVar_Species + k] = val;
        }

        /*
         * Set the atol value for the electrolyte voltage
         *      arithmetically scaled.-> so this is a characteristic value
         *         1 kcal gmol-1 = 0.05 volts
         */
        atolVector[indexCent_EqnStart + iVar_Voltage] = 0.05;
        atolVector[indexCent_EqnStart + iVar_Voltage_ED] = 0.05;
    }
}
//======================================================================================================================
void
porousLiIon_Cathode_dom1D::setAtolDeltaDamping(double atolDefault, double relcoeff,
                                               const Epetra_Vector_Ghosted& soln,
                                               Epetra_Vector_Ghosted& atolDeltaDamping,
                                               const Epetra_Vector_Ghosted* const atolV)
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
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];
        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */ 
	int iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        int iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
        int iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        int iVar_Voltage_ED = iVar_Voltage + 1;

        /*
         * Set the atol value for the axial velocity
         *   arithmetically scaled -> so this is a characteristic value
         */
        //double vax = soln[indexCent_EqnStart + iVAR_Vaxial];
        atolDeltaDamping[indexCent_EqnStart + iVAR_Vaxial] = 1.0E-4 * relcoeff;

        /*
         * Set atol values for the species mole fractions
         */
        double val = atolDefault * 1.0E3;
        if (val < 1.0E-4) {
            val = 1.0E-4;
        }
        for (int k = 0; k < nsp_; k++) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Species + k] = val * relcoeff;
        }

        /*
         * Set the atol value for the electrolyte voltage
         *      arithmetically scaled.-> so this is a characteristic value
         *         1 kcal gmol-1 = 0.05 volts
         */
        atolDeltaDamping[indexCent_EqnStart + iVar_Voltage] = 0.05 * relcoeff;
        atolDeltaDamping[indexCent_EqnStart + iVar_Voltage_ED] = 0.05 * relcoeff;
    }
}
//======================================================================================================================
void
porousLiIon_Cathode_dom1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff,
                                                       const Epetra_Vector_Ghosted& soln,
                                                       const Epetra_Vector_Ghosted& solnDot,
                                                       Epetra_Vector_Ghosted& atolDeltaDamping,
                                                       const Epetra_Vector_Ghosted* const atolV)
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
        NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[index_CentLcNode];

        /*
         *  Index of the first equation in the bulk domain of center node
         */
        int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
        /*
         * Offsets for the variable unknowns in the solution vector for the electrolyte domain
         */
	int iVAR_Vaxial  = nodeCent->indexBulkDomainVar0((size_t) Velocity_Axial);
        int iVar_Species = nodeCent->indexBulkDomainVar0((size_t) MoleFraction_Species);
        int iVar_Voltage = nodeCent->indexBulkDomainVar0((size_t) Voltage);
        int iVar_Voltage_ED = iVar_Voltage + 1;
        /*
         * Set the atol value for the axial velocity
         *   arithmetically scaled -> so this is a characteristic value
         */
        //double vax = soln[indexCent_EqnStart + iVAR_Vaxial];
        atolDeltaDamping[indexCent_EqnStart + iVAR_Vaxial] = 1.0E-4 * relcoeff;

        /*
         * Set atol values for the species mole fractions
         */
        double val = atolDefault * 1.0E5;
        if (val < 1.0E-1) {
            val = 1.0E-1;
        }
        for (int k = 0; k < nsp_; k++) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Species + k] = val * relcoeff;
        }

        /*
         * Set the atol value for the electrolyte voltage
         *      arithmetically scaled.-> so this is a characteristic value
         *         1 kcal gmol-1 = 0.05 volts
         */
        atolDeltaDamping[indexCent_EqnStart + iVar_Voltage] = 0.05 * relcoeff;
        atolDeltaDamping[indexCent_EqnStart + iVar_Voltage_ED] = 0.05 * relcoeff;
    }
}
//=====================================================================================================================
/**
 * Method to check for precipitation of the salts.
 * Returns index of offending cation or -1 if no precipitation
 */
int
porousLiIon_Cathode_dom1D::checkPrecipitation()
{
    //no precipitation
    return -1;
}
//=====================================================================================================================
void
porousLiIon_Cathode_dom1D::err(const char* msg)
{
    printf("porousLiIon_Cathode_dom1D: function not implemented: %s\n", msg);
    exit(-1);
}
//=====================================================================================================================
double porousLiIon_Cathode_dom1D::capacityPA(int platNum) const
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
double porousLiIon_Cathode_dom1D::capacityDischargedPA(int platNum) const
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
double porousLiIon_Cathode_dom1D::capacityLeftPA(int platNum, double voltsMax, double voltsMin) const
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
double porousLiIon_Cathode_dom1D::depthOfDischargePA(int platNum) const
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
double porousLiIon_Cathode_dom1D::depthOfDischargeStartingPA(int platNum) const
{
    double dodStarting = 0.0;
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; iCell++) {
	dodStarting += Electrode_Cell_[iCell]->depthOfDischargeStarting(platNum);
    }
    return dodStarting / crossSectionalArea_;
}
//=====================================================================================================================
// Reset the counters that keep track of the amount of discharge to date
void porousLiIon_Cathode_dom1D::resetCapacityDischargedToDate()
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
double porousLiIon_Cathode_dom1D::openCircuitPotentialQuick() const
{
    Electrode* ee = Electrode_Cell_[NumLcCells-1];
    double ocv = ee->openCircuitVoltage(nSurfsElectrode_ - 1);
    return ocv;
}
//=====================================================================================================================
//  WARNING -> will fail for multiprocessors, becauase the accumulation of information within eval_PostProc will fail.
double porousLiIon_Cathode_dom1D::effResistanceLayer(double &potAnodic, double &potCathodic,
						     double &voltOCV, double &current)
{
    static double resistance = 0.0;
    potAnodic = potentialAnodic_;
    potCathodic = potentialCathodic_;
    current = icurrElectrode_CBR_[NumLcCells-1];
    voltOCV = openCircuitPotentialQuick();
    //
    //  Calculate the effective resistance of the separator layer by taking the potential drop across the domain
    //  and dividing it by the current.
    //
    resistance = 0.0;
    if (fabs(current) > 1.0E-200) {
        resistance = (voltOCV - (potCathodic - potAnodic)) / current;
    }
    return resistance;
}

//=====================================================================================================================
} //namespace m1d
//=====================================================================================================================
