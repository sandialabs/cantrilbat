/**
 * @fie m1d_porousLiIon_Anode_dom1D.cpp
 */
/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_porousLiIon_Anode_dom1D.h"

#include "LiIon_PorousBat.h"
#include "m1d_BDT_porAnode_LiIon.h"
#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_DomainLayout_LiIon_PorousBat.h"
#include "m1d_ProblemStatementCell.h"
#include "m1d_SurfDomainDescription.h"
#include "m1d_SurDomain1D.h"
#include "m1d_globals.h"

#include "Electrode.h"
#include "Electrode_Factory.h"

#include "cantera/transport/Tortuosity.h"

using namespace std;
#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

//
//  Necessary expediency until we model dUdT correctly and fully, matching to experiment
//
// #define DELTASHEAT_ZERO 1

namespace m1d
{
//===================================================================================================================================
static void
drawline(int sp, int ll)
{
    for (int i = 0; i < sp; i++) {
        ZZCantera::writelog(" ");
    }
    for (int i = 0; i < ll; i++) {
        ZZCantera::writelog("-");
    }
    ZZCantera::writelog("\n");
}
//===================================================================================================================================
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

//===================================================================================================================================
porousLiIon_Anode_dom1D::porousLiIon_Anode_dom1D(BDT_porAnode_LiIon* bdd_anode_ptr) :
    porousElectrode_dom1D(bdd_anode_ptr),
    BDT_anode_ptr_(bdd_anode_ptr),
    nph_(0), nsp_(0),
    icurrInterfacePerSurfaceArea_Cell_(0),
    concTot_Cell_(0), concTot_Cell_old_(0),
    capacityDischargedPA_Cell_(0),
    depthOfDischargePA_Cell_(0),
    capacityLeftPA_Cell_(0),
    capacityPA_Cell_(0),
    Fleft_cc_(0.0), Fright_cc_(0.0), Vleft_cc_(0.0),
    Vcent_cc_(0.0), Vright_cc_(0.0), VElectrodeLeft_cc_(0.0), VElectrodeCent_cc_(0.0), VElectrodeRight_cc_(0.0),
    t_final_(0.0),
    t_init_(0.0),
    Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0), spCharge_(0),
    phiElectrode_Curr_(0.0),
    conductivityElectrode_(1.0E6), gradT_trCurr_(0.0),
    gradV_trCurr_(0.0), gradVElectrode_trCurr_(0.0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
    icurrElectrode_trCurr_(0.0),
    nSpeciesElectrode_(0), nSurfsElectrode_(0),
    electrodeSpeciesMoleDelta_Cell_(0),
    icurrInterface_Cell_(0), phaseMoleTransfer_(0),
    solnMoleFluxInterface_Cell_(0), icurrElectrode_CBL_(0), icurrElectrode_CBR_(0), icurrElectrolyte_CBL_(0),
    icurrElectrolyte_CBR_(0), deltaV_Cell_(0), 
    Ess_Surf_Cell_(0), 
    overpotential_Surf_Cell_(0), icurrRxn_Cell_(0),
    LiFlux_Cell_(0),
    iECDMC_(-1),
    iLip_(-1),
    iPF6m_(-1),
    solnTemp(0)
{

    //BDT_ptr_ = static_cast<BDT_porAnode_LiIon*>(&bdd);
    //BDT_porAnode_LiIon* fa = dynamic_cast<BDT_porAnode_LiIon*>(&bdd);
    if (!BDT_ptr_) {
        throw m1d_Error("confused", "confused");
    }
    /*
     * This is a shallow pointer copy. The BDT object owns the ionicLiquid_ object
     */
    ionicLiquid_ = BDT_ptr_->ionicLiquid_;
    /*
     *  This is a shallow pointer copy. The BDT object owns the transport object
     */
    trans_ = BDT_ptr_->trans_;
    /*
     *  This is a shallow pointer copy. The BDT object owns the Electrode object
     */
    nsp_ = BDT_ptr_->nSpeciesElectrolyte_;
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

    conductivityElectrode_ = PSinput.conductivityAnode_;

}
//==================================================================================================================================
porousLiIon_Anode_dom1D::porousLiIon_Anode_dom1D(const porousLiIon_Anode_dom1D& r) :
    porousElectrode_dom1D(r.BDT_anode_ptr_),
    BDT_anode_ptr_(r.BDT_anode_ptr_),
    nph_(0), nsp_(0),
    icurrInterfacePerSurfaceArea_Cell_(0), 
    concTot_Cell_(0), concTot_Cell_old_(0),
    capacityDischargedPA_Cell_(0),
    depthOfDischargePA_Cell_(0),
    capacityLeftPA_Cell_(0),
    capacityPA_Cell_(0),
    Fleft_cc_(0.0), Fright_cc_(0.0), Vleft_cc_(0.0),
    Vcent_cc_(0.0), Vright_cc_(0.0), VElectrodeLeft_cc_(0.0), VElectrodeCent_cc_(0.0), VElectrodeRight_cc_(0.0),
    t_final_(0.0),  t_init_(0.0),
    Xleft_cc_(0), Xcent_cc_(0), Xright_cc_(0), spCharge_(0),
    phiElectrode_Curr_(0.0), 
    conductivityElectrode_(1.0E6), gradT_trCurr_(0.0),
    gradV_trCurr_(0.0), gradVElectrode_trCurr_(0.0), gradX_trCurr_(0), Vdiff_trCurr_(0), jFlux_trCurr_(0),
    icurrElectrode_trCurr_(0.0),
    nSpeciesElectrode_(0), nSurfsElectrode_(0),
    electrodeSpeciesMoleDelta_Cell_(0),  icurrInterface_Cell_(0), phaseMoleTransfer_(0),
    solnMoleFluxInterface_Cell_(0), icurrElectrode_CBL_(0), icurrElectrode_CBR_(0), icurrElectrolyte_CBL_(0),
    icurrElectrolyte_CBR_(0), deltaV_Cell_(0), 
    Ess_Surf_Cell_(0), 
    overpotential_Surf_Cell_(0), 
    icurrRxn_Cell_(0),  LiFlux_Cell_(0),
    iECDMC_(-1), iLip_(-1),  iPF6m_(-1),
    solnTemp(0)
{
    porousLiIon_Anode_dom1D::operator=(r);
    conductivityElectrode_ = PSinput.conductivityAnode_;
}
//===================================================================================================================================
porousLiIon_Anode_dom1D::~porousLiIon_Anode_dom1D()
{
    for (int iCell = 0; iCell <NumLcCells; iCell++) {
        delete Electrode_Cell_[iCell];
        Electrode_Cell_[iCell] = 0;
    }
}
//===================================================================================================================================
porousLiIon_Anode_dom1D&
porousLiIon_Anode_dom1D::operator=(const porousLiIon_Anode_dom1D& r)
{
    if (this == &r) {
        return *this;
    }
    // Call the parent assignment operator
    porousElectrode_dom1D::operator=(r);

    BDT_anode_ptr_ = r.BDT_anode_ptr_;
    nph_ = r.nph_;
    nsp_ = r.nsp_;
    icurrInterfacePerSurfaceArea_Cell_ = r.icurrInterfacePerSurfaceArea_Cell_;
    concTot_Cell_ = r.concTot_Cell_;
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
    phiElectrode_Curr_ = r.phiElectrode_Curr_;
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
    iLip_  = r.iLip_;
    iPF6m_ = r.iPF6m_;
    solnTemp = r.solnTemp;

    throw CanteraError("", "not implemented");

    return *this;
}
//===================================================================================================================================
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
porousLiIon_Anode_dom1D::domain_prep(LocalNodeIndices* li_ptr)
{
    /*
     * First call the parent domain prep to get the node information.
     * Also all arrays defined by parent objects are sized appropriately.
     * Also, an initial attempt is made to calculate the porosity and the volumeFraction_Phases and molaNumber_Phases
     * is made using the reference temperature and pressure to calculate molar volumes.
     */
    porousElectrode_dom1D::domain_prep(li_ptr);

    //BDT_porAnode_LiIon* BDD_LiSi_Anode = dynamic_cast<BDT_porAnode_LiIon*>(&(BDD_));
    if (!BDT_anode_ptr_) {
        throw CanteraError(" porousLiIon_Anode_dom1D::domain_prep()", "bad dynamic cast ");
    }

    ZZCantera::Electrode* ee = BDT_anode_ptr_->Electrode_;
    nSpeciesElectrode_ = ee->nSpecies();
    nSurfsElectrode_ = ee->nSurfaces();

    //
    //  Placeholder for adding info about the solidSkeleton mole numbers in each cell
    //
    if (solidSkeleton_) {
    }

    /*
     * Resize vectors associated with this object
     */
    icurrInterfacePerSurfaceArea_Cell_.resize(NumLcCells, 0.0);
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
    spCharge_.resize(nsp_, 0.0);

    for (int k = 0; k < nsp_; k++) {
        spCharge_[k] = ionicLiquid_->charge(k);
    }

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

    mfElectrolyte_Soln_Cell_old_.resize(3, NumLcCells, 0.0);

    /*
     *  Set the velocity basis of the transport object. Initially, we are using
     *  mole-averaged velocities as the basis.
     */
    trans_->setVelocityBasis(ivb_);

    /*
     *  Set up the electrode object at each cell
     */
    instantiateElectrodeCells();

}
//==================================================================================================================================
//  An electrode object must be created and initialized for every cell in the domain
/*
 * Create electrode objects for every cell.
 * Correct the volume and number of moles of
 * active material within each of these electrode
 * objects to correspond to the discretized volume.
 *
 *  Everything that is used here comes from the anode.inp file.
 */
void
porousLiIon_Anode_dom1D::instantiateElectrodeCells()
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


    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        // Note that we need to create a Factory method to instantiate the desired electrode type.
        // int nSurPhases = PSinput.anode_input_->m_pl->nSurPhases();

        DomainLayout* DL = BDD_ptr_->DL_ptr_;

        DomainLayout_LiIon_PorousBat* dlc = dynamic_cast<DomainLayout_LiIon_PorousBat*>(DL);
        ProblemStatementCell* psc_ptr = dlc->pscInput_ptr_;
        ELECTRODE_KEY_INPUT* ai = psc_ptr->anode_input_;

        ZZCantera::Electrode* ee  = newElectrodeObject(ai->electrodeModelName);
        if (!ee) {
            throw  m1d_Error("porousLiIon_Anode_dom1D::instantiateElectrodeCells()",
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
        ee->setID(0, iCell);
        /*
         *  Set the initial subcycle deltaT policy
         */
        ee->choiceDeltaTsubcycle_init_ = 1;

        int retn = ee->electrode_model_create(PSinput.anode_input_);
        if (retn == -1) {
            printf("exiting with error\n");
            exit(-1);
        }

        retn = ee->setInitialConditions(PSinput.anode_input_);
        if (retn == -1) {
            printf("exiting with error\n");
            exit(-1);
        }
        /*
         *    Turn on the creation of xml state files. We need this for restart
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

	//
        // Compute total electrode volume
	//
        double totalElectrodeGrossVolume = -1.0;;
        double electrodeGrossArea = -1.0;
        double porosity = -1.0;
	double nonElectrodeVF = -1.0;
	//
        // Calculate the volume fraction of other phases than the electrode.
        // - Here we assume that the cells are all uniform. we can do this because we are in domain_prep()
	//
	double vfOther = volumeFractionOther(0);
        //
	// Calculate the gross electrode area
	//

        if (PSinput.anode_input_->electrodeGrossArea > 0.0) {
            electrodeGrossArea = PSinput.anode_input_->electrodeGrossArea;
        } else if (PSinput.anode_input_->electrodeGrossDiameter > 0.0) {
            electrodeGrossArea = Pi * 0.25 * PSinput.anode_input_->electrodeGrossDiameter *
                                 PSinput.anode_input_->electrodeGrossDiameter;
        }
	/*
         *  Save the cross-sectional area of the electrode to use throughout the code. It does not change within
         *  this calculation
         */
        if (electrodeGrossArea > 0.0) {
	    if (!doubleEqual(crossSectionalArea_, electrodeGrossArea, 1.0E-30, 9)) {
		throw m1d_Error("porousLiIon_Anode_dom1D::instantiateElectrodeCells())",
				"crossSectional Area, " + fp2str(crossSectionalArea_) + 
				", differs from cross sectional area input from anode file, " + fp2str(electrodeGrossArea));
	    }
	    electrodeGrossArea = crossSectionalArea_;
        } else {
	    electrodeGrossArea = crossSectionalArea_;
	}
        /*
         * If we have input the porosity directly, then we will honor that and calculate the electrode solid volume.
         * If we have input the electrode solid volume, then we will honor that calculation and calculate the resultant porosity.
         */
	if (PSinput.anode_input_->porosity > 0.0) {
	    porosity = PSinput.anode_input_->porosity;
	} else {
	    if (PSinput.anode_input_->electrodeGrossThickness > 0.0 && electrodeGrossArea > 0.0) {
		totalElectrodeGrossVolume = electrodeGrossArea * PSinput.anode_input_->electrodeGrossThickness;
		double totalSolidGrossVolume =  ee->SolidVol();
		porosity = 1.0 - totalSolidGrossVolume / totalElectrodeGrossVolume - vfOther;
		//printf("Anode porosity is %f with %g m^3 solid gross volume and %g m^3 electrode gross volume.\n",
		//       porosity, ee->SolidVol(), totalElectrodeGrossVolume);
		if (porosity <= 0.0) {
		    throw CanteraError("porousLiIon_Anode_dom1D::instantiateElectrodeCells()", "Computed porosity is not positive.");
		}
	    } else {
		throw m1d_Error("porousLiIon_Anode_dom1D::instantiateElectrodeCells()","Unknown Porosity");
            }
	}
  
	nonElectrodeVF = porosity + vfOther; 
	if (nonElectrodeVF < 0.0) {
	    throw m1d_Error("porousLiIon_Anode_dom1D::instantiateElectrodeCells()", "nonElectrodeVF < 0");
	}        if (nonElectrodeVF > 1.0) {
 	    throw m1d_Error("porousLiIon_Anode_dom1D::instantiateElectrodeCells()", "nonElectrodeVF > 1");
	}

	/*
	 * Reset the moles in the electrode using the computed porosity
	 * and the electrodeWidth FOR THIS node.
	 */
	ee->setElectrodeSizeParams(electrodeGrossArea, xdelCell_Cell_[iCell], nonElectrodeVF);

        // update porosity as computed from electrode input
        porosity_Cell_[iCell] = porosity;
	porosity_Cell_old_[iCell] = porosity;

	nVol_zeroStress_Electrode_Cell_[iCell] = ee->SolidVol();
	nVol_zeroStress_Electrode_Old_Cell_[iCell] = nVol_zeroStress_Electrode_Cell_[iCell];


        Electrode_Cell_[iCell] = ee;
    }
}
//==================================================================================================================================
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
porousLiIon_Anode_dom1D::advanceTimeBaseline(const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                                             const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr,
                                             const double t, const double t_old)
{
    const Epetra_Vector& soln = *soln_ptr;
    bool assumedAdvance = (t > t_old);
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
	//
	//  Advance the value of the old cell porosity
	//
        porosity_Cell_old_[iCell] = porosity_Curr_;
	//
	//  Advance the value of the old cell temperature
	//
	Temp_Cell_old_[iCell] = temp_Curr_;

        double* mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);
	for (size_t k = 0; k < (size_t) nsp_; ++k) {
            mfElectrolyte_Soln_old[k] = mfElectrolyte_Soln_Curr_[k];
        }
    
        /*
         * Tell the electrode object to accept the current step and prep for the next step.
         *
         * We might at this point do a final integration to make sure we nailed the conditions of the last step.
         * However, we will hold off at implementing this right now
         */
        Electrode* ee = Electrode_Cell_[iCell];
        ee->resetStartingCondition(t, assumedAdvance);
	/*
	 *  Store the Amount of Li element in each cell, for later use by conservation
	 *  checking routines
	 */
	size_t neSolid = ee->nElements();
        for (size_t elem = 0; elem < neSolid; ++elem) {
            elem_Solid_Old_Cell_(elem,iCell) = ee->elementSolidMoles("Li");
        }

        ee->updateState();

        ee->setInitStateFromFinal(true);

	if (energyEquationProbType_ == 3) { 
	    double volCellNew = xdelCell_Cell_[iCell];
	    // double volElectrodeCell = solidVolCell / crossSectionalArea_;
	    double solidEnthalpy = ee->SolidEnthalpy() / crossSectionalArea_;
	    double solidEnthalpyNew = solidEnthalpy;
	  
	    double lyteMolarEnthalpyNew = ionicLiquid_->enthalpy_mole();
	    double volLyteNew = porosity_Curr_ * volCellNew;
	    double lyteEnthalpyNew =  lyteMolarEnthalpyNew * concTot_Curr_ * volLyteNew;
	    double nEnthalpyInertNew = 0.0;
	    for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
		ExtraPhase* ep = ExtraPhaseList_[jPhase];
                ThermoPhase* tp = ep->tp_ptr;
		tp->setState_TP(temp_Curr_, pres_Curr_);
		double mEnth = tp->enthalpy_mole();
		double mn = moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + jPhase];
		nEnthalpyInertNew += mEnth * mn;
	    }
	    double nEnthalpy_New  = solidEnthalpyNew + lyteEnthalpyNew + nEnthalpyInertNew;
	    if (! checkDblAgree( nEnthalpy_New, nEnthalpy_New_Cell_[iCell] ) ) {
                    throw m1d_Error("porousLiIon_Anode_dom1D::advanceTimeBaseline",
				    "Disagreement on new enthalpy calc");
	    }
	    
	    nEnthalpy_Old_Cell_[iCell] = nEnthalpy_New_Cell_[iCell];
	    nEnthalpy_Electrode_Old_Cell_[iCell] = nEnthalpy_Electrode_New_Cell_[iCell];
	}

 	nVol_zeroStress_Electrode_Old_Cell_[iCell] = nVol_zeroStress_Electrode_Old_Cell_[iCell];

        if (numExtraCondensedPhases_ > 0) {
	for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; jPhase++) {
	    moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + jPhase] =  
	      moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase];
	    volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + jPhase] =  
	      volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase];  
	}
        }
    }
}

//==================================================================================================================================
// Revert the Residual object's conditions to the conditions at the start of the global time step
/*
 *  If there was a solution in t_final, this is wiped out and replaced with the solution at t_init_init.
 *
 *  Virtual from m1d_domain.h
 */
void
porousLiIon_Anode_dom1D::revertToInitialGlobalTime()
{
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        /*
         * Tell the electrode object to forget about any recent step and revert to t_init_init conditions
         */
        Electrode* ee = Electrode_Cell_[iCell];
        ee->revertToInitialTime(true);
    }
}


//===================================================================================================================================
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
porousLiIon_Anode_dom1D::residEval(Epetra_Vector& res,
                                   const bool doTimeDependentResid,
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
    residType_Curr_ = residType;
    //int index_RightLcNode;
    //int index_LeftLcNode;
    int index_CentLcNode;

    // -- example of how to get a particular residual id
    // if (counterResBaseCalcs_ == 65) {
    // }

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
    double newStuffTC;
    double oldStuffTC;

    //  Electrolyte mole fluxes - this is c V dot n at the boundaries of the cells
    double moleFluxRight = 0.0;
    double moleFluxLeft = 0.0;

    // Thermal Fluxes
    double fluxTright = 0.0;
    double fluxTleft = 0.0;
    double fluxL_JHPhi = 0.0;
    double fluxR_JHPhi = 0.0;
    double fluxL_JHelec = 0.0;
    double fluxR_JHelec = 0.0;
    double enthConvRight = 0.0;
    double enthConvLeft = 0.0;
    
    // Calculated current at the current collector that is conservative
    double icurrElectrode_LBcons = 0.0;

    // Flux of current in the electrode phase at the right and left cell boundaries
    double fluxVElectrodeRight = 0.0;
    double fluxVElectrodeLeft = 0.0;
    
    // scratch pad to calculate the expansion of the matrix
    std::vector<double> xratio(NumLcCells,0.0); 

    // scratch pad to calculate the new positions due to the strain of the elements
    std::vector<double> new_node_pos(NumLcCells,0.0);
  
    // average grad grad stress
    //    double avg_delta_matrix_pressure = 0;

  // for the anode material, we use the following values:
    // initial porosity  of the pyrolitic graphite  matrix_poro = 0.4
    // Volume frac of binder                        binder_vf   = 0.064
    
    // from http://www-ferp.ucsd.edu/LIB/PROPS/PANOS/c.html
    // v poisson's ratio             ~~0.5 
    // K  Bulk Modulus of solid material = 33 MPa
    // E Younds Modulus                  = 4800 MPa  = stress/strain
    // modification ala Simulation of elastic moduli of porous materials
    // CREWES Research Report — Volume 13 (2001) 83
    // Simulation of elastic moduli for porous materials
    // Charles P. Ursenbach 
    // @ 40% porosity of graphite, K = 0.226 * K0
    // @ 30% porosity,             K = 0.276 * K0
    // youngs mod == E =  3K(1-2v) 
    // shear mod G = 3*K(1-2v)/2(1+v)
    //
    // for the porous matrix, assume grain locking, and use the same values. 
    //    double Thermal_Expansion = 50e-6;
#ifdef MECH_MODEL
    double poisson = 0.45;
    double BulkMod = 33.0E6 * 0.226; // hard wired untill we get the porosity and Chi values. 
    double Eyoung = 3*BulkMod*(1.0 - 2*poisson);
    //    double G = 3*BulkMod*(1-2*poisson)/(2*(1+poisson));
    int pMECH_MODEL_extra = 1;
#endif

    //mole fraction fluxes
    std::vector<double> fluxXright(nsp_, 0.0);
    std::vector<double> fluxXleft(nsp_, 0.0);
    std::vector<double> resLocal_Species(nsp_, 0.0);
    double fluxL = 0.0;
    //double fluxR = 0.0;

    const Epetra_Vector& soln = *soln_ptr;
    double* mf_old;

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
        if (IOwnLeft && iCell == NumLcCells - 1) {
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
         * for this flux at the right cell edge of the previous cell. We don't need to redo it
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
            Vleft_cc_ = soln[indexLeft_EqnStart +  nodeTmpsLeft.Offset_Voltage ];
            VElectrodeLeft_cc_ = soln[indexLeft_EqnStart +  nodeTmpsLeft.Offset_Voltage  + 1];
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
                Xright_cc_[k] = soln[indexRight_EqnStart +  nodeTmpsRight.Offset_MoleFraction_Species + k];
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
	// We've turned this off as it is not necessary, if we are doing everything correctly
	// This is turned on for iCell = 0
        if (doLeftFluxCalc) {
            if (nodeLeft == 0) {

		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
                /*
                 *  We are here if we are at the left node boundary and we need a flux condition. The default now is to
                 *  set the flux to zero. This is good if we expect the negate/equalize the fluxes at a domain boundary,
		 *  since both sides are set to zero.
		 *
		 *  The other case is when we are specifying a flux boundary condition. In that case,
		 *  we set the flux here to zero. Then, in the surface domain we add the specification of
		 *  the flux to the residual of the node at the boundary.
                 */
		Fleft_cc_ = 0.0;
                moleFluxLeft = 0.0;
		fluxTleft = 0.0;
		fluxL_JHPhi = 0.0;
		fluxL_JHelec = 0.0;                 // Note this is non-zero, but we set the term elsewhere
		enthConvLeft = 0.0;
		// we will calculate this here
                icurrElectrolyte_CBL_[iCell] = 0.0;
                fluxVElectrodeLeft = 0.0;
                icurrElectrode_CBL_[iCell] = 0.0;   // Note this is non-zero, but we set the term elsewhere
                for (int k = 0; k < nsp_; k++) {
                    fluxXleft[k] = 0.0;
                }
		SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));
		//
		// Calculate the current at the boundary that is conservative.
		//   This is really  icurrElectrode_CBL_[iCell]. However, we leave that 0 because of the
		//   need to impose boundary conditions.
		//
		icurrElectrode_LBcons = icurrInterface_Cell_[iCell]  + icurrElectrode_CBR_[iCell];
		jFlux_EnthalpyPhi_metal_trCurr_ =  - icurrElectrode_LBcons / Faraday;
		jFlux_EnthalpyPhi_metal_trCurr_ *= EnthalpyPhiPM_metal_Curr_[0];
		fluxL_JHelec = jFlux_EnthalpyPhi_metal_trCurr_;
            } else {

                /*
                 *  Establish the environment at the left cell boundary
                 */
                SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent,  &(soln[indexCent_EqnStart]), 0);
                /*
                 *  Calculate the transport properties at the left cell boundary
                 */
                SetupTranShop(xdelL, 0);
                /*
                 * Calculate the flux at the left boundary for each equation
                 */
                moleFluxLeft = Fleft_cc_ * concTot_Curr_;
		fluxTleft = heatFlux_Curr_; 
		fluxL_JHPhi = jFlux_EnthalpyPhi_Curr_;
		fluxL_JHelec = jFlux_EnthalpyPhi_metal_trCurr_;
		enthConvLeft = moleFluxLeft * EnthalpyMolar_lyte_Curr_;
		/*
                 * Calculate the flux of species and the flux of charge
                 *   - the flux of charge must agree with the flux of species
                 */
                icurrElectrolyte_CBL_[iCell] = 0.0;
                fluxVElectrodeLeft = icurrElectrode_trCurr_;
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
                icurrElectrolyte_CBL_[iCell] *= (ZZCantera::Faraday);
            }
        } else {  // !doLeftFluxCalc
            /*
             * Copy the fluxes from the stored right side
             */
            moleFluxLeft = moleFluxRight;
	    fluxTleft = fluxTright;
	    fluxL_JHPhi = fluxR_JHPhi;
	    fluxL_JHelec = fluxR_JHelec;
	    enthConvLeft = enthConvRight;
            icurrElectrolyte_CBL_[iCell] = icurrElectrolyte_CBR_[iCell - 1];
            fluxVElectrodeLeft = fluxVElectrodeRight;
            icurrElectrode_CBL_[iCell] = icurrElectrode_CBR_[iCell - 1];
            for (int k = 0; k < nsp_; k++) {
                fluxXleft[k] = fluxXright[k];
            }
        }
        /*
         * ------------------- CALCULATE FLUXES AT THE RIGHT BOUNDARY -------------------------------
         *
         */
        if (nodeRight == 0) {
           
            AssertTrace(iCell == NumLcCells-1);
           
            SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
	    /*
             *  We are here if we are at the right node boundary and we need a flux
             *  condition. The default treatment is to set the flux to zero. The basic idea is that
	     *  if the right flux in this domain and the left flux in the next domain are zero,
	     *  then we have cancelled out the fluxes in the middle of a control volume that straddles
	     *  the two domains. We could put in a more sophisticated treatment later, but we would
	     *  have to ensure that treatment cancels fluxes in the middle.
             */
	    Fright_cc_ = 0.0;
            moleFluxRight = 0.0;
	    fluxTright = 0.0;
	    fluxR_JHPhi = 0.0;  
	    fluxR_JHelec = 0.0;
	    enthConvRight = 0.0;
            icurrElectrolyte_CBR_[iCell] = 0.0;
            fluxVElectrodeRight = 0.0;
            icurrElectrode_CBR_[iCell] = 0.0;   // Note this is nonzero. We set it elsewhere.
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
             * Calculate the flux at the right boundary for each equation
             * This is equal to
             *       Conc * Vaxial * phi
             */
            moleFluxRight = Fright_cc_ * concTot_Curr_;
	    /*
             * Calculate the heat flux - all of the types
             */
	    fluxTright = heatFlux_Curr_;
	    fluxR_JHPhi = jFlux_EnthalpyPhi_Curr_;
	    fluxR_JHelec = jFlux_EnthalpyPhi_metal_trCurr_;
	    enthConvRight = moleFluxRight * EnthalpyMolar_lyte_Curr_;

            /*
             * Calculate the flux of species and the flux of charge
             *   - the flux of charge must agree with the flux of species
             */
            icurrElectrolyte_CBR_[iCell] = 0.0;
            fluxVElectrodeRight = icurrElectrode_trCurr_;
            icurrElectrode_CBR_[iCell] = icurrElectrode_trCurr_;
            for (int k = 0; k < nsp_; k++) {
                fluxXright[k] = jFlux_trCurr_[k];

                icurrElectrolyte_CBR_[iCell] += fluxXright[k] * spCharge_[k];
                if (Fright_cc_ > 0.0) {
                    fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_;
                } else {
                    fluxXright[k] += Fright_cc_ * mfElectrolyte_Thermo_Curr_[k] * concTot_Curr_;
                }
            }
            icurrElectrolyte_CBR_[iCell] *= (ZZCantera::Faraday);

        }

#ifdef DEBUG_RESID
        if (doTimeDependentResid) {

            if (residType == Base_ResidEval) {
                printf(" Cell = %d, Totalflux_K+ = %10.3e,  Totalflux_Cl- = %10.3e \n", iCell, fluxXright[1], fluxXright[2]);
                printf("           Old porosity = %10.3e,  New porosity = %10.3e \n", porosity_Cell_old_[iCell], porosity_Cell_[iCell]);
                printf("           Vmolal = %10.3e, jd_Li+ = %10.3e  jd_K+ = %10.3e jd_Cl- = %10.3e\n", Fright_cc_,
                       jFlux_trCurr_[0], jFlux_trCurr_[1], jFlux_trCurr_[2]);
                printf("           Vmolal = %10.3e, vd_Li+ = %10.3e  vd_K+ = %10.3e vd_Cl- = %10.3e\n", Fright_cc_,
                       Vdiff_trCurr_[0], Vdiff_trCurr_[1], Vdiff_trCurr_[2]);
            }
        }
#endif

        /*
         * --------------- ADD FLUX TERMS INTO THE RESIDUALS --------------------------------------
         */
#ifdef DEBUG_RESID
        if (iCell == 9 && residType == Base_ResidEval) {
            if (iCell == NumLcCells - 1) {
                // printf("we are here fluxFLeft = %g\n", fluxFleft);
            }
        }
#endif

        /*
         *  Total continuity equation - fluxFright and fluxFleft represent the total mole
         *                              fluxes coming and going from the cell.
         *                    R =   d rho d t + del dot (rho V) = 0
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Electrolyte_Continuity] += (moleFluxRight - moleFluxLeft);

        /*
         * Momentum Equation - When there is a pressure unknown
         *            
         */
        if (PS_ptr->hasPressureEquation_) {
              // res[indexCent_EqnStart + nodeTmpsCenter.RO_Momentum_Axial] += permeability * (pressureRight - pressureLeft) / deltaX_right;
        }

        /*
         * Species continuity equation
         */
        for (int k = 0; k < nsp_; k++) {
            if (k != iECDMC_ && k != iPF6m_) {
                res[indexCent_EqnStart + nodeTmpsCenter.RO_Species_Eqn_Offset + k] += (fluxXright[k] - fluxXleft[k]);
            }
        }
        /*
         *   Current conservation equation in the electrolyte
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation] += 
	  icurrElectrolyte_CBR_[iCell] - icurrElectrolyte_CBL_[iCell];

        /*
         *   Current conservation equation - electrode
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation + 1] += (fluxVElectrodeRight - fluxVElectrodeLeft);

	/*
         *  Energy Equation
         */
        if (PS_ptr->energyEquationProbType_ == 3) {
	    if (nodeLeft == 0) {
		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
		icurrElectrode_LBcons = icurrInterface_Cell_[iCell]  + icurrElectrode_CBR_[iCell];
		jFlux_EnthalpyPhi_metal_trCurr_ =  - icurrElectrode_LBcons / Faraday;
		metalPhase_->setState_TP(temp_Curr_, pres_Curr_);
		metalPhase_->setElectricPotential(phiElectrode_Curr_);
		metalPhase_->getPartialMolarEnthalpies(&EnthalpyPhiPM_metal_Curr_[0]);
		EnthalpyPhiPM_metal_Curr_[0] -= Faraday * phiElectrode_Curr_;
		jFlux_EnthalpyPhi_metal_trCurr_ *= EnthalpyPhiPM_metal_Curr_[0];
		fluxL_JHelec = jFlux_EnthalpyPhi_metal_trCurr_;
	    }

	    AssertTrace(nodeTmpsCenter.RO_Enthalpy_Conservation != npos);
	    res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += (fluxTright - fluxTleft);
	    res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += (fluxR_JHPhi - fluxL_JHPhi);
	    res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += (fluxR_JHelec - fluxL_JHelec);
	    res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += (enthConvRight - enthConvLeft);

        }

        /*
         *   ------------------- ADD SOURCE TERMS TO THE CURRENT CELL CENTER --------------------------------------
         */

        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));

        /*
         *    Source terms for the species production rate of Li+.
         */
        for (int k = 0; k < nsp_; k++) {
            if (k != iECDMC_ && k != iPF6m_) {
                res[indexCent_EqnStart + nodeTmpsCenter.RO_Species_Eqn_Offset + k] -=
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
         *         units are coul /m2 /s
         */
        res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation] -= icurrInterface_Cell_[iCell] ;

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
                icurrElectrode_CBL_[iCell] = icurrInterface_Cell_[iCell] + icurrElectrode_CBR_[iCell];
            }
            if (residType == Base_ShowSolution || residType == Base_ResidEval) {
                DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Current_Conservation] = 
		  res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation];
                DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Current_Conservation + 1] =
		  res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation + 1];

		if (PS_ptr->energyEquationProbType_ == 3) {
		    DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation] = fluxL_JHelec;
		}
            }
      
        }
        /*
         * Special section if we own the right node of the domain. If we do
         * it will be cell 0. Store currents for later use.
         * These are the correct currents that work for the global balances
         */
        if (IOwnRight && iCell == (NumLcCells - 1)) {
            if (residType == Base_ResidEval) {
                DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Current_Conservation] = 
		  res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation];
                DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Current_Conservation + 1] =
		  res[indexCent_EqnStart + nodeTmpsCenter.RO_Current_Conservation + 1];
		if (PS_ptr->energyEquationProbType_ == 3) {
		    DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation] = 
		      res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation + 1];
		}
            }
            if (residType == Base_ShowSolution || residType == Base_ResidEval) {
                icurrElectrolyte_CBR_[iCell] =  icurrInterface_Cell_[iCell]  + icurrElectrolyte_CBL_[iCell];
            }
        }
        if (residType == Base_ShowSolution) {
            deltaV_Cell_[iCell] = Electrode_ptr->voltage();
            for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++) {
                Ess_Surf_Cell_[nSurfsElectrode_ * iCell + jSurf] = Electrode_ptr->openCircuitVoltage(jSurf);
                overpotential_Surf_Cell_[nSurfsElectrode_ * iCell + jSurf] = Electrode_ptr->overpotential(jSurf);
            }
            icurrRxn_Cell_[iCell] = icurrInterface_Cell_[iCell];

            LiFlux_Cell_[iCell] = jFlux_trCurr_[iLip_];
	}
      


        /*
         *   ------------------ ADD IN THE TIME DEPENDENT TERMS ----------------------------------------------------
         */
        if (doTimeDependentResid) {

#ifdef DEBUG_HKM_NOT
            if (residType == Base_ResidEval) {
                printf(" Cell = %d, Totalflux_Li+_r = %10.3e,  = %10.3e, Totalflux_Li+_l ", 
		       iCell, fluxXright[iLip_], fluxXleft[iLip_]);
            }
#endif
            newStuffTC = concTot_Curr_ * porosity_Curr_ * xdelCell;
            double newStuffSpecies0 = Xcent_cc_[iLip_] * newStuffTC;

            /*
             *   .................... Calculate quantities needed at the previous time step
             */

            mf_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);

            oldStuffTC = concTot_Cell_old_[iCell] * porosity_Cell_old_[iCell] * xdelCell;
            double oldStuffSpecies0 = mf_old[iLip_] * oldStuffTC;
            double tmp = (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;

#ifdef DEBUG_HKM_RESID
            if (residType == Base_ResidEval) {
                printf(" deltaT term = %10.3e BulkSum = %10.3e\n", tmp, tmp + (fluxXright[iLip_] - fluxXleft[iLip_]));
            }
#endif
            /*
             *   .................... Add these terms in the residual
             */
            /*
             *  Add in the time term for species 0
             */
            for (int k = 0; k < nsp_; k++) {
                if (k != iECDMC_ && k != iPF6m_) {
                    res[indexCent_EqnStart + nodeTmpsCenter.RO_Species_Eqn_Offset + k] += tmp;
                }
            }

            /*
             *   Add in the time term for the total continuity equation
             *         note: The current problem will have this term equally zero always.
             *               However, we put it in here for the next problem.
             */
            res[indexCent_EqnStart + nodeTmpsCenter.RO_Electrolyte_Continuity] += (newStuffTC - oldStuffTC) * rdelta_t;

        

	    if  (energyEquationProbType_ == 3) {
		//
		//  Need to count up the enthalpy over the electrode and the electrolyte on a per cross-sectional area basis
		//
		// Get the Solid enthalpy in Joules / m2
		//
		double solidEnthalpyNew = Electrode_ptr->SolidEnthalpy() / crossSectionalArea_;
		//
		// Calculate the enthalpy in the electrolyte Joules / kmol  (kmol/m3) (m) = Joules / m2
		// 
		double lyteMolarEnthalpyNew = ionicLiquid_->enthalpy_mole();
		double volLyteNew = porosity_Curr_ * xdelCell;
		double lyteEnthalpyNew =  lyteMolarEnthalpyNew * concTot_Curr_ * volLyteNew;
		// 
		// Calculate the enthalpy of other phases that may be present
		//
		double nEnthalpyInertNew = 0.0;
		for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
		    ExtraPhase* ep = ExtraPhaseList_[jPhase];
		    ThermoPhase* tp = ep->tp_ptr;
		    tp->setState_TP(temp_Curr_, pres_Curr_);
		    double mEnth = tp->enthalpy_mole();
		    double mn = moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + jPhase];
		    nEnthalpyInertNew += mEnth * mn;
		}
		//
		// Calculate and store the total enthalpy in the cell at the current conditions
		//
		nEnthalpy_New_Cell_[iCell] = solidEnthalpyNew + lyteEnthalpyNew +  nEnthalpyInertNew;
		/*
		 *   .................... Calculate quantities needed at the previous time step
		 *         These are already storred in nEnthalpy_Old_Cell_
		 *
		 *  Do the energy time derivative
		 */
		double tmp = (nEnthalpy_New_Cell_[iCell] - nEnthalpy_Old_Cell_[iCell]) * rdelta_t;
		res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += tmp;

	    } else if (energyEquationProbType_ == 4) {
		//
		//  Need to count up the heat capacity over the electrode and the electrolyte on a per cross-sectional area basis
		//
		double cp =  CpMolar_total_Cell_[iCell];
		double tempOld = Temp_Cell_old_[iCell];
		res[indexCent_EqnStart + nodeTmpsCenter.RO_Enthalpy_Conservation] += cp * (temp_Curr_ - tempOld) * rdelta_t;

	    }
	}
	//
	//  Section to find the residual contributions and the
	//       Axial velocity at the right domain boundary
	//
	if (residType == Base_ShowSolution) {
	    if (IOwnRight) {
		if (iCell == NumLcCells - 1) {
		    SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
		    SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));

		    double res_Cont = ((newStuffTC - oldStuffTC) * rdelta_t + ( - moleFluxLeft )
				       - solnMoleFluxInterface_Cell_[NumLcCells - 1]);
		    moleFluxRight = - res_Cont;
		    Fright_cc_ = moleFluxRight / concTot_Curr_;
		    DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity] = Fright_cc_;

		    for (size_t k = 0; k < (size_t) nsp_; k++) {
			resLocal_Species[k] = (fluxXright[k] - fluxXleft[k]);
			double newStuffSpecies0 =  Xcent_cc_[k] * newStuffTC;
			double oldStuffSpecies0 =  mf_old[k] * oldStuffTC;
			resLocal_Species[k] += (newStuffSpecies0 - oldStuffSpecies0) * rdelta_t;
			resLocal_Species[k] -= 
			  electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * iCell + indexStartEOelectrolyte + k] * rdelta_t /
			  crossSectionalArea_;
			//
			// Save the local domain residual to a vector for later validation studies
			//
			DomainResidVectorRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] = resLocal_Species[k];

			TotalFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] = - resLocal_Species[k];
			fluxXright[k] = - resLocal_Species[k];
			DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k] =  
			  TotalFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k]
			  -  Fright_cc_ * mfElectrolyte_Soln_Curr_[k] * concTot_Curr_;
		    }
		    if  (energyEquationProbType_ == 3) {
			double resLocal_Enth = (fluxTright - fluxTleft);
			resLocal_Enth += (fluxR_JHPhi - fluxL_JHPhi);
			resLocal_Enth += (fluxR_JHelec - fluxL_JHelec);
			resLocal_Enth += (enthConvRight - enthConvLeft);
			resLocal_Enth += (nEnthalpy_New_Cell_[iCell] - nEnthalpy_Old_Cell_[iCell]) * rdelta_t;
			DomainResidVectorRightBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation] = resLocal_Enth;
		    }
		}
	    }
	}

#ifdef MECH_MODEL
	if (solidMechanicsProbType_ > 0) {
 	    valCellTmps& valTmps = valCellTmpsVect_Cell_[iCell];
	    double thick_lc_now = -9e9;
	    double gross_vol_now = -9e9;
	    if(iCell == 0 ) {
	      // thick_lc_now =  (  nodeRight->x0NodePos()+soln[indexRight_EqnStart + nodeTmpsRight.Offset_Displacement_Axial ] - 
	      // 			 nodeCent->x0NodePos()+soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ] );
	      thick_lc_now =  (  nodeRight->x0NodePos() - nodeCent->x0NodePos() );
  	      gross_vol_now =  Electrode_Cell_[iCell]->SolidVol()/(1.0-calcPorosity(iCell)) + 
		0.5*  Electrode_Cell_[iCell+1]->SolidVol()/(1.0-calcPorosity(iCell+1));

	      //	      std::cout << " 00iCell "<<iCell<<" grossVol "<< gross_vol_now<<" thickness "<<thick_lc_now<<std::endl;
	    }
	    else  if (nodeLeft && nodeRight){
	      // thick_lc_now = 0.5* ( (nodeRight->x0NodePos()+soln[indexRight_EqnStart + nodeTmpsRight.Offset_Displacement_Axial ]) -
	      // 			    ( nodeLeft->x0NodePos()+soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Displacement_Axial ]));
	      thick_lc_now = 0.5* ( nodeRight->x0NodePos() - nodeLeft->x0NodePos());
	      gross_vol_now = 0.5*( Electrode_Cell_[iCell]->SolidVol()/(1.0-calcPorosity(iCell)) +
                		    Electrode_Cell_[iCell+1]->SolidVol()/(1.0-calcPorosity(iCell+1)));

	      if(iCell == NumLcCells-2) {
		gross_vol_now+= 0.5*Electrode_Cell_[iCell+1]->SolidVol()/(1.0-calcPorosity(iCell+1));
	      }
	      //	      std::cout << " biCell "<<iCell<<" grossVol "<< gross_vol_now<<" thickness "<<thick_lc_now<<std::endl;
	    }
	    else if (nodeLeft && (!nodeRight) ) {
	      // thick_lc_now =  ( (nodeCent->x0NodePos()+soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ]) -
	      // 			(nodeLeft->x0NodePos()+soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Displacement_Axial ]));
	      thick_lc_now =  nodeCent->x0NodePos() - nodeLeft->x0NodePos();
	      gross_vol_now =  Electrode_Cell_[iCell]->SolidVol()/(1.0-calcPorosity(iCell)) + 
	 	           0.5*Electrode_Cell_[iCell-1]->SolidVol()/(1.0-calcPorosity(iCell-1));

	      //	      std::cout << " nriCell "<<iCell<<" grossVol "<< gross_vol_now<<" thickness "<<thick_lc_now<<std::endl;
	    }
	    else {
	      std::cout.flush();
	      throw m1d_Error("porousLiIon_Anode_dom1D:: ERROR",
			      "thick_lc_now Contact Developer, unexpected nodeRight nodeLeft null / non-null error." );
	    }

	    xratio[iCell]=1.0;
	  
	    // All on or ChemEx expansion is turned on. 
	    if ( (Domain1D::ChemEx | Domain1D::All)  & solidMechanicsProbType_) {
	      xratio[iCell] *= gross_vol_now/(thick_lc_now * crossSectionalArea_);
	    }

	    // the divergence of the pressure == the trace of the STRESS tensor
	    // HOWEVER my understanding is that the pressure variable is the fluid pressure, not the solid matrix pressure. 
	    // SO the below result may need to be modified by the non-solid-volume ratio. 
	    // Also!!!
	    // With the assumption that the time step is much much greater than the sound speed across the whole battery layer;
	    // then this is quasi static: the """correct""" way to calculate the pressure induced strain would be to average the pressure
	    // across the whole battery, and adjust the strain in each cell so that average would be equal to the right boundary pressure. 
	    // While the implicit solution generated by calculating the residuals as below will give the same result, the convergance
	    // of the solution may be greatly lenghtened by not explicitly calculating the quasi-static solution. 

	    size_t iVar_Pressure = nodeCent->indexBulkDomainVar0((size_t) Pressure_Axial);
	    double pressure_STRESS = 0;
	    if (iVar_Pressure != npos) {
	      if(nodeRight && ! nodeLeft)
	        pressure_STRESS = -(1.0/BulkMod)*(soln[indexRight_EqnStart + iVar_Pressure]-soln[indexCent_EqnStart + iVar_Pressure]);
	      else if(nodeLeft && ! nodeRight)
		pressure_STRESS = -(1.0/BulkMod)*(soln[indexLeft_EqnStart + iVar_Pressure]-soln[indexCent_EqnStart + iVar_Pressure]);
	      else if (nodeLeft && nodeRight) {
		pressure_STRESS = -(1.0/BulkMod)*(
						  (soln[indexLeft_EqnStart + iVar_Pressure]-soln[indexCent_EqnStart + iVar_Pressure])+
						  (soln[indexRight_EqnStart + iVar_Pressure]-soln[indexCent_EqnStart + iVar_Pressure])
						  );
	      }	 
	    } // pressure exists
	    
	    double pressure_strain = pressure_STRESS/Eyoung;
	    if ( (Domain1D::FluidPr | Domain1D::All) & solidMechanicsProbType_) {
	      if(iCell ==1) xratio[iCell-1]*= (1.0+pressure_strain); 
	      xratio[iCell] *= (1.0+pressure_strain);
	    }
	    
	    // nodeTmpsCenter.Offset_Solid_Stress_Axial = nodeCent->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
	    // if(nodeLeft) 
	    //   nodeTmpsLeft.Offset_Solid_Stress_Axial   = nodeLeft->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
	    // else
	    //   nodeTmpsLeft.Offset_Solid_Stress_Axial   = nodeTmpsCenter.Offset_Solid_Stress_Axial;
	    // if(nodeRight)
	    //   nodeTmpsRight.Offset_Solid_Stress_Axial   = nodeRight->indexBulkDomainVar0((size_t) Solid_Stress_Axial);
	    // else
	    //   nodeTmpsRight.Offset_Solid_Stress_Axial   = nodeTmpsCenter.Offset_Solid_Stress_Axial;
	    
	    // double left_matrix_stress = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Solid_Stress_Axial] ;
	    // double center_matrix_stress = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Solid_Stress_Axial] ;
	    // double right_matrix_stress = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Solid_Stress_Axial] ;
	    // double matrix_pressure_left = - (left_matrix_stress- center_matrix_stress);
	    // double matrix_pressure_right= - (center_matrix_stress - right_matrix_stress);
	    // double matrix_LP_center = matrix_pressure_left - matrix_pressure_right;

	    // xratio[iCell-1] *= (1.0+ matrix_LP_center/Eyoung);
	    // if(iCell == NumLcCells-1)
	    //   xratio[iCell] *= (1.0+ (0.5*matrix_LP_center/Eyoung));
	    // avg_delta_matrix_pressure += matrix_pressure_left;
	  
	  // since we have the half control volumes at the right and left hand boundaries the divisor is NumLcCells-1
	  //	  avg_delta_matrix_pressure /= (NumLcCells-1);
	}

#endif
    }

#ifdef MECH_MODEL  
    if (solidMechanicsProbType_ > 0) {
	//  node[0] is pinned so it never moves, hence start at iCell=1
	/*
	 *  ------------------------------ LOOP OVER CELL -------------------------------------------------
	 *  Loop over the number of Cells in this domain on this processor
	 *  This loop is done from left to right.
	 */
	for (int iCell = 1; iCell < NumLcCells; iCell++) {
	    cIndex_cc_ = iCell;
	    cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	    NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
	    NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

	    NodalVars* nodeCent  = cTmps.nvCent_;
	    indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;

	    NodalVars* nodeLeft  = cTmps.nvLeft_;
	    indexLeft_EqnStart = nodeTmpsLeft.index_EqnStart;

	    nodeRight = cTmps.nvRight_;
	    indexRight_EqnStart = nodeTmpsRight.index_EqnStart;

	    nodeTmpsCenter.Offset_Displacement_Axial = nodeCent->indexBulkDomainVar0((size_t) Displacement_Axial);

	    if (nodeLeft != 0) {
		nodeTmpsLeft.Offset_Displacement_Axial   = nodeLeft->indexBulkDomainVar0((size_t) Displacement_Axial);
	    } else {
		nodeTmpsLeft.Offset_Displacement_Axial = -1;
	    }
	 
	    if(iCell ==1) {
	      new_node_pos[0] = nodeLeft->x0NodePos(); 
	      // Left most node is pinned at zero displacement
	    }
	    double delta_0 = (nodeCent->x0NodePos() + soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial]) 
	      - (nodeLeft->x0NodePos() + soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Displacement_Axial]);
	    double new_delta = delta_0 *  xratio[iCell-1]; 
	    new_node_pos[iCell] = new_node_pos[iCell-1] + new_delta;

	    // stress
	    // double left_matrix_stress = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Solid_Stress_Axial] ;
	    // double center_matrix_stress = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Solid_Stress_Axial] ;
	    // double lc_pressure = -(left_matrix_stress-center_matrix_stress);
	    // res[indexCent_EqnStart + nodeTmpsCenter.Offset_Solid_Stress_Axial] = left_matrix_stress + (avg_delta_matrix_pressure-lc_pressure); 
	    // if(soln[indexCent_EqnStart +  nodeTmpsCenter.Offset_Solid_Stress_Axial ]!=0.0) 
	    //   std::cout << " Anode::residEval iCell "<<iCell
	    //             <<" stress_axial "<<soln[indexCent_EqnStart +  nodeTmpsCenter.Offset_Solid_Stress_Axial ]
            //             <<std::endl;
	} // end of iCell loop

      for (int iCell = 0; iCell < NumLcCells; iCell++) {
	  cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	  NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	  NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;	
	  NodalVars* nodeCent = cTmps.nvCent_;
	  NodalVars* nodeLeft = cTmps.nvLeft_;
	  indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
	  nodeTmpsCenter.Offset_Displacement_Axial   = nodeCent->indexBulkDomainVar0((size_t) Displacement_Axial);

	  res[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ] = 
	    soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial]
	    -(new_node_pos[iCell]  - nodeCent->x0NodePos());

	  if (iCell == 0) { 
	    res[indexLeft_EqnStart + nodeTmpsLeft.Offset_Displacement_Axial]  = 
	      soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial] - 0.0;
	  }
	  if (pMECH_MODEL_extra) {
	    if(iCell==0)  std::cout << " anode::residEval sol,res,r-1.0,x-x0:: ";
	    std::cout << " ( "<<soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ]<<" "
		      << res[indexCent_EqnStart + nodeTmpsCenter.Offset_Displacement_Axial ] <<" "
		      << xratio[iCell]-1.0<<" "
		      <<(new_node_pos[iCell]  - nodeCent->x0NodePos())<<" ), ";
	    if(iCell == NumLcCells-1) std::cout<<endl;

	  }
      }
      //impose that -grad trace Solid_Stress == the average across this part of the battery.
      
    }
#endif
}
//==================================================================================================================================
// Utility function to calculate quantities before the main residual routine.
/*
 *  This is used for a loop over nodes. All calculated quantities must be internally stored.
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
porousLiIon_Anode_dom1D::residEval_PreCalc(const bool doTimeDependentResid,
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

    porousFlow_dom1D::residEval_PreCalc(doTimeDependentResid, soln_ptr, solnDot_ptr, solnOld_ptr,
                                        t, rdelta_t, residType, solveType); 
    residType_Curr_ = residType;
    const Epetra_Vector& soln = *soln_ptr;

    maxElectrodeSubIntegrationSteps_ = 0;

    t_final_ = t;
    if (rdelta_t < 1.0E-200) {
        t_init_ = t;
    } else {
        // We want an infinitly small time step 
        if (solveType == TimeDependentInitial) {
           t_init_ = t;
        } else {
           t_init_ = t - 1.0/rdelta_t;
        }
    }

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
        int indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
        /*
         * Calculate the cell width
         */
	xdelCell_Cell_[iCell] = cTmps.xCellBoundaryR_ - cTmps.xCellBoundaryL_;

        for (int k = 0; k < nsp_; k++) {
            Xcent_cc_[k] = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_MoleFraction_Species + k];
        }
        Vcent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage];
        VElectrodeCent_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Voltage + 1];
	//
        //  Alter the porosity here for some problem types. Other problem types involve calculating an equation 
        //      -> we do this here when the porosity is not a formal variable.
        //
        if (porosityEquationProbType_  & Porosity_EqnType_Status::CalculatedOutOfEqnSystem) {
            if (porosityEquationProbType_  & Porosity_EqnType_Status::PartOfMechanics) {
               // do something different
            } else {
                porosity_Cell_[iCell] = calcPorosity(iCell);
            }
        }
        /*
         * Setup the thermo
         */
        SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
        /*
         *  ----------- Integrate the electrode object during the current time step ------------------------------------------
         */
        int numSubcycles = calcElectrode();
        maxElectrodeSubIntegrationSteps_ = std::max(maxElectrodeSubIntegrationSteps_, numSubcycles);
        numElectrodeSubCycles_Cell_[iCell] = numSubcycles;
	//
	//  Get the total lyte, total solid and total cell heat capacity of the current cell
	//  Store them in vectors for later use.
	//
	if (energyEquationProbType_) {
	    CpMolar_total_Cell_[iCell] = getCellHeatCapacity(nodeCent, &(soln[indexCent_EqnStart]));

	    EnthalpyMolar_lyte_Cell_[iCell] = ionicLiquid_->enthalpy_mole();
	}
    } // loop over cell
}
//==================================================================================================================================
//
//  Calculate Various quantities after the solution has been obtained
//
void
porousLiIon_Anode_dom1D::eval_PostSoln(
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
            potentialAnodic_ = VElectrodeCent_cc_;
        }
        if (iCell == NumLcCells-1) {
            potentialCathodic_ = Vcent_cc_;
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
	    icurrElectrolyte_CBL_[iCell] *= (ZZCantera::Faraday);
	    icurrElectrode_CBL_[iCell] = icurrElectrode_trCurr_;
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
            icurrElectrolyte_CBR_[iCell] *= (ZZCantera::Faraday);
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
	// Add in the electrode contribution
#ifdef DELTASHEAT_ZERO
        overPotentialHeat_Cell_curr_[iCell] =
	  Electrode_ptr->getIntegratedThermalEnergySourceTerm_overpotential() / crossSectionalArea_;
        deltaSHeat_Cell_curr_[iCell] = 0.0;
        electrodeHeat_Cell_curr_[iCell] =  overPotentialHeat_Cell_curr_[iCell];
#else
        electrodeHeat_Cell_curr_[iCell] = Electrode_ptr->getIntegratedThermalEnergySourceTerm() / crossSectionalArea_;
        overPotentialHeat_Cell_curr_[iCell] = 
	  Electrode_ptr->getIntegratedThermalEnergySourceTerm_overpotential() / crossSectionalArea_;
        deltaSHeat_Cell_curr_[iCell] = 
	  Electrode_ptr->getIntegratedThermalEnergySourceTerm_reversibleEntropy() / crossSectionalArea_;
#endif
        qSource_Cell_curr_[iCell] += electrodeHeat_Cell_curr_[iCell];
        qSource_Cell_accumul_[iCell] += qSource_Cell_curr_[iCell];
    }
}
//==================================================================================================================================
//
//  Calculate Various quantities after the solution has been obtained
//
void
porousLiIon_Anode_dom1D::eval_HeatBalance(const int ifunc,
					    const double t,
					    const double deltaT,
					    const Epetra_Vector *soln_ptr,
					    const Epetra_Vector *solnDot_ptr,
					    const Epetra_Vector *solnOld_ptr,
                                            struct globalHeatBalVals& dVals)
{
    NodalVars* nodeCent = 0;
    NodalVars* nodeLeft = 0;
    NodalVars* nodeRight = 0;
    globalHeatBalValsBat* dValsB_ptr = dynamic_cast<globalHeatBalValsBat*>(&dVals);
    int indexCent_EqnStart, indexLeft_EqnStart, indexRight_EqnStart;
    const Epetra_Vector& soln = *soln_ptr;
    double xdelL; // Distance from the center node to the left node
    double xdelR; // Distance from the center node to the right node
#ifdef DEBUG_PRINT_CELL_TABLES
    int doPrint = 1;
#else
    int doPrint = 0;
#endif
    int doTimes = 2;
    int nColsTable = 173;
    double 	icurrElectrode_LBcons;
    //
    // Find the pointer for the Right separator
    //
    SurfDomainDescription *rightS = BDD_ptr_->RightSurf;
    BulkDomainDescription *rightDD = rightS->RightBulk;
    BulkDomain1D *rightD_sep = rightDD->BulkDomainPtr_;

    SurfDomainDescription *leftS = BDD_ptr_->LeftSurf;
    SurDomain1D *leftSD_collector = leftS->SurDomain1DPtr_;

    //double phiIcurrL = 0.0;
    //double phiIcurrR = 0.0;
    double moleFluxLeft = 0.0;
    double moleFluxRight = 0.0;
    for (int itimes = 0; itimes < doTimes; itimes++) {
	if (doPrint) {
	    drawline(1, nColsTable);
	    if (itimes == 0) {
		printf("                     ANODE Enthalpy Cell Balance \n\n");
		printf("Cell|   NewnEnth       OldnEnth        deltanEnth | ");
		printf("    fluxTLeft     FluxTRight | ");
		printf("    fluxL_JHPhi  fluxR_JHPhi | ");
		printf("   fluxL_JHelec fluxR_JHelec | ");
		printf("  enthConvLeft enthConvRight | ");
	
		printf("    Resid_Add     Residual  |");
		printf("\n");
	    }
	    if (itimes == 1) {
		printf("\n\n                Analysys of Heat Source Terms:   ANODE \n");
		printf("      ");
		printf("              JOULE HEATING    |  Conduction    |        REACTIONS     \n");
		printf("Cell| ");
		printf("  lyteJHSource   solidJHSource |  HeatSource    |  RxnHeatSrcTerm   OverPotHeat  deltaSHeat");
		printf("\n");
	    }
	}
	double fluxTright = 0.0, fluxTleft = 0.0, fluxR_JHPhi = 0.0, fluxL_JHPhi = 0.0, enthConvRight = 0.0,
	  enthConvLeft = 0.0, resid, residAdd = 0.0;
	double fluxL_JHelec = 0.0, fluxR_JHelec = 0.0;
	dVals.oldNEnthalpy = 0.0;
	dVals.newNEnthalpy = 0.0;
	dVals.totalHeatCapacity = 0.0;
	for (int iCell = 0; iCell < NumLcCells; iCell++) {
	    cIndex_cc_ = iCell;

	    cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	    NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
	    NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;

	    //Electrode* Electrode_ptr = Electrode_Cell_[iCell];
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
		potentialAnodic_ = VElectrodeCent_cc_;
	    }
	    if (iCell == NumLcCells-1) {
		potentialCathodic_ = Vcent_cc_;
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
		fluxTleft = heatFlux_Curr_;
		fluxL_JHPhi = jFlux_EnthalpyPhi_Curr_;
		fluxL_JHelec = jFlux_EnthalpyPhi_metal_trCurr_;

		Fleft_cc_ = soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Velocity_Axial];
		moleFluxLeft = Fleft_cc_ * concTot_Curr_;
		enthConvLeft = moleFluxLeft * EnthalpyMolar_lyte_Curr_;

		/*
		 * Calculate the flux of species and the flux of charge
		 *   - the flux of charge must agree with the flux of species
		 */
		icurrElectrolyte_CBL_[iCell] = 0.0;
		for (int k = 0; k < nsp_; k++) {
		    icurrElectrolyte_CBL_[iCell] += jFlux_trCurr_[k] * spCharge_[k];
		}
		icurrElectrolyte_CBL_[iCell] *= (ZZCantera::Faraday);
		icurrElectrode_CBL_[iCell] = icurrElectrode_trCurr_;
	    } else {
		//
		//  Setup shop at the left node
		//
		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
		SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));
		//phiIcurrL = icurrElectrolyte_CBL_[iCell] * phiElectrolyte_Curr_;
		fluxTleft = 0.0;
		fluxL_JHelec = 0.0;
		fluxL_JHPhi = 0.0;
		// Calculate the current at the boundary that is conservative.
		//   This is really  icurrElectrode_CBL_[iCell]. However, we leave that 0 because of the
		//   need to impose boundary conditions.
		//
		icurrElectrode_CBL_[iCell] = icurrInterface_Cell_[iCell]  + icurrElectrode_CBR_[iCell];
		icurrElectrode_LBcons = icurrInterface_Cell_[iCell]  + icurrElectrode_CBR_[iCell];
		jFlux_EnthalpyPhi_metal_trCurr_ =  - icurrElectrode_LBcons / Faraday;
		jFlux_EnthalpyPhi_metal_trCurr_ *= EnthalpyPhiPM_metal_Curr_[0];
		fluxL_JHelec = jFlux_EnthalpyPhi_metal_trCurr_;
		//phiIcurrL = icurrElectrolyte_CBL_[iCell] * phiElectrolyte_Curr_;
         
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
		fluxTright = heatFlux_Curr_;
		fluxR_JHPhi = jFlux_EnthalpyPhi_Curr_;
		fluxR_JHelec = jFlux_EnthalpyPhi_metal_trCurr_;
		Fright_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Velocity_Axial];
		moleFluxRight = Fright_cc_ * concTot_Curr_;
		enthConvRight = moleFluxRight * EnthalpyMolar_lyte_Curr_;
		/*
		 * Calculate the flux of species and the flux of charge
		 *   - the flux of charge must agree with the flux of species
		 */
		icurrElectrolyte_CBR_[iCell] = 0.0;
		for (int k = 0; k < nsp_; k++) {
		    icurrElectrolyte_CBR_[iCell] += jFlux_trCurr_[k]* spCharge_[k];
		}
		icurrElectrolyte_CBR_[iCell] *= (ZZCantera::Faraday);
	    } else {
		enthConvRight = 0.0;
		fluxR_JHPhi = 0.0;
		fluxR_JHelec = 0.0;
		fluxTright = 0.0;
	    } 

	    if (doPrint) {
		double deltanEnth = nEnthalpy_New_Cell_[iCell] - nEnthalpy_Old_Cell_[iCell];
		if (itimes == 0) {
		    resid = 0.0;
		    resid = deltanEnth + deltaT *( fluxTright - fluxTleft);
		    resid += deltaT * (fluxR_JHPhi - fluxL_JHPhi + enthConvRight - enthConvLeft);
		    resid += deltaT * (fluxR_JHelec - fluxL_JHelec);
		    
		    residAdd = 0.0;
		    if (iCell == 0) {
			residAdd = deltaT *
			  leftSD_collector->DomainResidVector_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation];
		    }
		    if (iCell == NumLcCells - 1) {
			residAdd = deltaT *
			  rightD_sep->DomainResidVectorLeftBound_LastResid_NE[nodeTmpsCenter.RO_Enthalpy_Conservation];
		    }
		    resid += residAdd;
		    
		    printf("%3d |  % 12.6E  % 12.6E  % 12.5E |", iCell, nEnthalpy_New_Cell_[iCell], 
			   nEnthalpy_Old_Cell_[iCell], deltanEnth);
		    printf("   % 12.5E  % 12.5E |",   deltaT *fluxTleft,    deltaT *fluxTright);
		    printf("   % 12.5E  % 12.5E |",   deltaT *fluxL_JHPhi,  deltaT *fluxR_JHPhi);
		    printf("   % 12.5E  % 12.5E |",   deltaT *fluxL_JHelec, deltaT *fluxR_JHelec);
		    printf("   % 12.5E  % 12.5E |",   deltaT *enthConvLeft, deltaT *enthConvRight);
		    printf("  % 12.5E  % 12.5E |", residAdd, resid);
		    printf("  \n");
		}
		if (itimes == 1) {

		    printf("%3d |  % 12.6E  % 12.6E  | ", iCell, jouleHeat_lyte_Cell_curr_[iCell],
			   jouleHeat_solid_Cell_curr_[iCell ]);
		    //double fluxL_JH = fluxL_JHPhi - phiIcurrL;
		    // double fluxR_JH = fluxR_JHPhi - phiIcurrR;
		    //double Einflux = deltaT * (enthConvLeft - enthConvRight + fluxL_JH -  fluxR_JH);
	
		    double HeatIn = deltaT * (fluxTleft - fluxTright);
		
		    double rxnEnt = electrodeHeat_Cell_curr_[iCell]; 
		    double nui    = overPotentialHeat_Cell_curr_[iCell];
		    double idels  = deltaSHeat_Cell_curr_[iCell];
		    printf(" % 12.6E | % 12.6E % 12.6E  % 12.6E ", 
			   HeatIn, rxnEnt, nui, idels);
		    printf("\n");

		}
	    }
	    dVals.totalHeatCapacity += CpMolar_total_Cell_[iCell];
	    //
	    //  Count up the total old and new cell enthalpies
	    //
	    dVals.oldNEnthalpy += nEnthalpy_Old_Cell_[iCell];
	    dVals.newNEnthalpy += nEnthalpy_New_Cell_[iCell];

	    if (iCell == 0) {
		dValsB_ptr->currentLeft = icurrElectrode_CBL_[iCell];
		dValsB_ptr->enthFluxOut = enthConvLeft;
		dValsB_ptr->phiSolid  = phiElectrode_Curr_;
		dValsB_ptr->enthalpyIVfluxRight = icurrElectrode_CBL_[iCell] *  phiElectrode_Curr_;
		dValsB_ptr->JHelecLeft = fluxL_JHelec;

	    }
	}
    }
}
//==================================================================================================================================
// Evaluate species and element balances
void
porousLiIon_Anode_dom1D::eval_SpeciesElemBalance(const int ifunc,
						 const double t,
						 const double deltaT,
						 const Epetra_Vector *soln_ptr,
						 const Epetra_Vector *solnDot_ptr,
						 const Epetra_Vector *solnOld_ptr,
						 class globalHeatBalVals& dVals)
{
    NodalVars* nodeCent = 0;
    NodalVars* nodeLeft = 0;
    NodalVars* nodeRight = 0;
    globalHeatBalValsBat* dValsB_ptr = dynamic_cast<globalHeatBalValsBat*>(&dVals);
    int indexCent_EqnStart, indexLeft_EqnStart, indexRight_EqnStart;
    const Epetra_Vector& soln = *soln_ptr;
    double xdelL; // Distance from the center node to the left node
    double xdelR; // Distance from the center node to the right node
    double xdelCell;
#ifdef DEBUG_PRINT_CELL_TABLES
    int doPrint = 1;
#else
    int doPrint = 0;
#endif

    // Initially, we'll limit it to Li transport, PF6- bal, and solvent balance
    int doTimes = 1;
    int nColsTable = 175;
  
    //
    // Find the pointer for the left Anode
    // SurfDomainDescription *leftS = BDD_ptr_->LeftSurf;
  
    //
    // Find the pointer for the left Cathode
    SurfDomainDescription *rightS = BDD_ptr_->RightSurf;
    BulkDomainDescription *rightDD = rightS->RightBulk;
    BulkDomain1D *rightD_sep = rightDD->BulkDomainPtr_;

    double moleFluxLeft = 0.0;
    double moleFluxRight = 0.0;
    double residAdd = 0.0;

    ZZCantera::Array2D elem_Lyte_New_Cell(nsp_, NumLcCells, 0.0);
    ZZCantera::Array2D elem_Lyte_Old_Cell(nsp_, NumLcCells, 0.0);

    std::vector<double>& elemLi_Lyte_New = dValsB_ptr->elem_Lyte_New;
    std::vector<double>& elemLi_Lyte_Old = dValsB_ptr->elem_Lyte_Old;
    std::vector<double>& elem_Solid_New = dValsB_ptr->elem_Solid_New;
    std::vector<double>& elem_Solid_Old = dValsB_ptr->elem_Solid_Old;

    std::vector<double> icurrCBL_Cell(NumLcCells, 0.0);
    std::vector<double> icurrCBR_Cell(NumLcCells, 0.0);

    std::vector<double> convRight(nsp_, 0.0);
    std::vector<double> convLeft(nsp_, 0.0);
    std::vector<double> jFluxRight(nsp_, 0.0);
    std::vector<double> jFluxLeft(nsp_, 0.0);
    std::vector<double> species_Lyte_New_Curr(nsp_, 0.0);
    std::vector<double> species_Lyte_Old_Curr(nsp_, 0.0);
    std::vector<double> species_Lyte_Src_Curr(nsp_, 0.0);

    std::vector<double>& species_Lyte_New_Total = dValsB_ptr->species_Lyte_New_Total;
    std::vector<double>& species_Lyte_Old_Total = dValsB_ptr->species_Lyte_Old_Total;   
    std::vector<double> species_convRight  = dValsB_ptr->species_convRight;
    std::vector<double> species_convLeft   = dValsB_ptr->species_convLeft;
    std::vector<double> species_jFluxRight = dValsB_ptr->species_jFluxRight;
    std::vector<double> species_jFluxLeft  = dValsB_ptr->species_jFluxLeft;
    std::vector<double>& species_Lyte_Src_Total = dValsB_ptr->species_Lyte_Src_Total;   
   
    ZZCantera::Array2D species_Lyte_New_Cell(nsp_, NumLcCells, 0.0);
    ZZCantera::Array2D species_Lyte_Old_Cell(nsp_, NumLcCells, 0.0);
    ZZCantera::Array2D res_Species(nsp_,NumLcCells, 0.0);

    /*
     *   Find the species index for the first species in the electrode object pertaining to the electrolyte
     */
    int sf = Electrode_Cell_[0]->solnPhaseIndex();
    int indexStartEOelectrolyte = Electrode_Cell_[0]->getGlobalSpeciesIndex(sf, 0);

    for (int itimes = 0; itimes < doTimes; itimes++) {
	if (doPrint) {
	    if (itimes == 0) {
		drawline(1, nColsTable);
		printf(" \n\n           ANODE ::          Species Cell balances \n");
		printf("Cell Spec|     NewnBal       OldnBal     deltanSpecies | ");
		printf("  - SrcTerm   | ");
		printf("    fluxLeft       FluxRight | ");
		printf("    convLeft    convRight | ");
		printf("  Residual-additions  Residual  |  Current_CBL  Current CBR ");
		printf("\n");
		drawline(1, nColsTable);
	    }
	  
	}
	double resid;

	for (int iCell = 0; iCell < NumLcCells; iCell++) {
	    cIndex_cc_ = iCell;

	    cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
	    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
	    NodeTmps& nodeTmpsLeft   = cTmps.NodeTmpsLeft_;
	    NodeTmps& nodeTmpsRight  = cTmps.NodeTmpsRight_;
	    double icurrCBL = 0.0;
	    double icurrCBR = 0.0;
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


	    xdelL = cTmps.xdelL_;
	    xdelR = cTmps.xdelR_;
	    xdelCell = cTmps.xdelCell_;


	    double* mfElectrolyte_Soln_old = mfElectrolyte_Soln_Cell_old_.ptrColumn(iCell);

	    /*
	     * ------------------------ Get the indexes for the right node ------------------------------------
	     */
	    nodeRight = cTmps.nvRight_;
	    indexRight_EqnStart = nodeTmpsRight.index_EqnStart;
	    for (size_t k = 0; k < (size_t) nsp_; k++) {
		Xright_cc_[k] = soln[indexRight_EqnStart + nodeTmpsRight.Offset_MoleFraction_Species + k];
	    }
	    Vright_cc_ = soln[indexRight_EqnStart + nodeTmpsRight.Offset_Voltage];

	    if (nodeLeft != 0) {
		/*
		 *  Establish the environment at the left cell boundary
		 */
		SetupThermoShop2(nodeLeft, &(soln[indexLeft_EqnStart]), nodeCent, &(soln[indexCent_EqnStart]), 0);	    
		SetupTranShop(xdelL, 0);
	
		Fleft_cc_ =soln[indexLeft_EqnStart + nodeTmpsLeft.Offset_Velocity_Axial];
		moleFluxLeft = Fleft_cc_ * concTot_Curr_;
		for (size_t k = 0; k < (size_t) nsp_; k++) {
		    convLeft[k] =  moleFluxLeft *  mfElectrolyte_Soln_Curr_[k];
		    jFluxLeft[k] = jFlux_trCurr_[k];
		    icurrCBL += jFluxLeft[k] *spCharge_[k] * Faraday;
		}
	    } else {
		//
		//  Setup shop at the left node
		//
		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
		SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));

		Fleft_cc_ = DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
		moleFluxLeft = Fleft_cc_ * concTot_Curr_;
		for (size_t k = 0; k <  (size_t) nsp_; k++) {
		    convLeft[k] = 0.0;
		}
		for (size_t k = 0; k < (size_t) nsp_; k++) {
		    jFlux_trCurr_[k] = - DiffFluxLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
		    jFluxLeft[k] = 0.0;
		    icurrCBL += jFlux_trCurr_[k] *spCharge_[k] * Faraday;
		}
	    }

	    if (nodeRight != 0) {
		/*
		 *  Establish the environment at the right cell boundary
		 */
		SetupThermoShop2(nodeCent, &(soln[indexCent_EqnStart]), nodeRight, &(soln[indexRight_EqnStart]), 1);
		SetupTranShop(xdelR, 1);

		Fright_cc_ = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Velocity_Axial];
		moleFluxRight = Fright_cc_ * concTot_Curr_;
		for (size_t k = 0; k <  (size_t) nsp_; k++) {
		    convRight[k] = moleFluxRight * mfElectrolyte_Soln_Curr_[k];
		    jFluxRight[k] = jFlux_trCurr_[k];
		    icurrCBR += jFluxRight[k] *spCharge_[k] * Faraday;
		}
	    } else {
		//
		//  Setup shop at the right node
		//
		SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));	
		SetupThermoShop1Extra(nodeCent, &(soln[indexCent_EqnStart]));

		Fright_cc_ = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
		moleFluxRight = Fright_cc_ * concTot_Curr_;
		for (size_t k = 0; k <  (size_t) nsp_; k++) {
		    convRight[k] = 0.0;
		}
		for (size_t k = 0; k <  (size_t) nsp_; k++) {
		    jFlux_trCurr_[k] = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
		    jFluxRight[k] = 0.0;
		    icurrCBR += jFlux_trCurr_[k] * spCharge_[k] * Faraday;
		}
	    }

	    //
	    // Return to the center node
	    //
	    SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
	    for (size_t k = 0; k < (size_t) nsp_; k++) {
		elem_Lyte_New_Cell(k, iCell) = porosity_Curr_ * concTot_Curr_ * xdelCell * mfElectrolyte_Soln_Curr_[k];
		elem_Lyte_Old_Cell(k, iCell) = (porosity_Cell_old_[iCell] * concTot_Cell_old_[iCell] *
						xdelCell  * mfElectrolyte_Soln_old[k]);
	    
		elemLi_Lyte_New[k] += elem_Lyte_New_Cell(k, iCell);
		elemLi_Lyte_Old[k] += elem_Lyte_Old_Cell(k, iCell);
	    }
	    //const std::vector<double>& species_ElectrodeMoles_new = Electrode_ptr->getMoleNumSpecies();
            double eLiMoles = Electrode_ptr->elementSolidMoles("Li");
            // printf(" 0 = %g\n", species_ElectrodeMoles_new[0]);
            //printf(" eLi = %g\n", eLiMoles );
            elem_Solid_New[iLip_] += eLiMoles;
            int eIndex = Electrode_ptr->elementIndex("Li");
            elem_Solid_Old[iLip_] += elem_Solid_Old_Cell_(eIndex, iCell);


	    for (size_t k = 0; k < (size_t) nsp_; k++) {
		species_Lyte_New_Curr[k] = porosity_Curr_ * concTot_Curr_* xdelCell  * mfElectrolyte_Soln_Curr_[k];
		species_Lyte_Old_Curr[k] = (porosity_Cell_old_[iCell] * concTot_Cell_old_[iCell]* xdelCell  * 
					    mfElectrolyte_Soln_old[k]);
		species_Lyte_New_Total[k] += species_Lyte_New_Curr[k];
		species_Lyte_Old_Total[k] += species_Lyte_Old_Curr[k];

		species_Lyte_Src_Curr[k] = 
                    electrodeSpeciesMoleDelta_Cell_[nSpeciesElectrode_ * iCell + indexStartEOelectrolyte + k] /
                    crossSectionalArea_;
		species_Lyte_Src_Total[k] += species_Lyte_Src_Curr[k];

		if (iCell == 0) {
		    species_convLeft[k] = convLeft[k];
		    species_jFluxLeft[k] =  jFluxLeft[k];	 
		}
		if (iCell == NumLcCells-1) {
		    species_convRight[k]  = convRight[k];
		    species_jFluxRight[k] = jFluxRight[k];
		}
	    }

	    if (doPrint) {
		if (itimes == 0) {
		    resid = 0.0;
		    for (size_t k = 0 ; k <  (size_t) nsp_; k++) {
			double deltanSp = species_Lyte_New_Curr[k] - species_Lyte_Old_Curr[k];
			double convSpRight = convRight[k];
			double convSpLeft = convLeft[k];
			double jFluxSpRight =  jFluxRight[k];
			double jFluxSpLeft =  jFluxLeft[k];
			resid = deltanSp;
			resid -= species_Lyte_Src_Curr[k];
			resid += deltaT * ( convSpRight - convSpLeft );
			resid += deltaT * ( jFluxSpRight - jFluxSpLeft );

			residAdd = 0.0;
			if (iCell == 0) {
			    residAdd = 0.0;
			    //residAdd = deltaT * 
			    //  leftD_anode->DomainResidVectorRightBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
			}
			if (iCell == NumLcCells - 1) {
			    residAdd = deltaT *
			      rightD_sep->DomainResidVectorLeftBound_LastResid_NE[nodeTmpsCenter.RO_Species_Eqn_Offset + k];
			}
			resid += residAdd;
			printf("%3d, %3d |  % 12.6E  % 12.6E  % 12.5E |", iCell, (int) k, 
			       species_Lyte_New_Total[k], species_Lyte_Old_Total[k], deltanSp);
			printf("  %12.5E |", - species_Lyte_Src_Curr[k]);
			printf("  % 12.5E  % 12.5E |",  - deltaT * jFluxSpLeft,  deltaT * jFluxSpRight);
			printf("  % 12.5E  % 12.5E |",  - deltaT * convSpLeft,  deltaT * convSpRight);
			printf("  % 12.5E  % 12.5E |", residAdd,  resid);
			printf("  % 12.5E  % 12.5E |", icurrCBL, icurrCBR);
			printf("  \n");
		    }
		}
	    }

	}
	if (doPrint) {
	    drawline(1, nColsTable);
	    printf("\n");
	}
    }
}
//==================================================================================================================================
//   Main routine to integrate the electrode for the current residual
/*
 *  This routine has been moved to the precalc category.
 */
int
porousLiIon_Anode_dom1D::calcElectrode()
{
    Electrode* Electrode_ptr = Electrode_Cell_[cIndex_cc_];
    int numSubcycles = 0;
    bool doInstanteous = false;
    double deltaT = t_final_ - t_init_;
    if (deltaT == 0.0) {
        doInstanteous = true;
    }
    Electrode_ptr->updateState();
    //
    // The above statement set the final state. Here we set the init and init_init state
    // to be equal to the final state
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
	if (residType_Curr_ == Base_ShowSolution) {
	    if (cIndex_cc_ == 6) {
		//printf("we are here\n");
	    }
	}
        numSubcycles = Electrode_ptr->integrate(deltaT);
#ifdef DEBUG_HKM
	if (ProblemResidEval::s_printFlagEnv > 0 && BDD_ptr_->Residual_printLvl_ > 8) {
	    if (numSubcycles > 15) {
		printf("      anode::calcElectrode(): problematic integration ( %d) : Cell = %d, counterNumberIntegrations_ = %d\n",
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
     * To get the surface area per m**2 we need to divide by the the cross sectional area of the electrode
     */
    surfaceArea_Cell_[cIndex_cc_] = saExternal / (crossSectionalArea_);

    //
    //  Get the extrinsic volume of the electrode at the current T, and isotropic P without any stress terms
    //
    if (! (porosityEquationProbType_ & Porosity_EqnType_Status::Constant)) { 
        nVol_zeroStress_Electrode_Cell_[cIndex_cc_] = Electrode_ptr->SolidVol();
    }

    /*
     *  Calculate the current per external surface area of electrode. To do this calculation
     *  we divide the total extrinsic current by the total extrinsic surface area
     */
    icurrInterfacePerSurfaceArea_Cell_[cIndex_cc_] = icurrInterface_Cell_[cIndex_cc_] * crossSectionalArea_ /
                                                     saExternal;

    if (energyEquationProbType_ == 3) {
	nEnthalpy_Electrode_New_Cell_[cIndex_cc_] = Electrode_ptr->SolidEnthalpy();
    }

    if (residType_Curr_ == Base_ShowSolution) {
        capacityDischargedPA_Cell_[cIndex_cc_] = Electrode_ptr->capacityDischarged() / crossSectionalArea_;
        depthOfDischargePA_Cell_[cIndex_cc_] = Electrode_ptr->depthOfDischarge() / crossSectionalArea_;
        capacityLeftPA_Cell_[cIndex_cc_] = Electrode_ptr->capacityLeft() / crossSectionalArea_;
        capacityPA_Cell_[cIndex_cc_]= Electrode_ptr->capacity() / crossSectionalArea_;
    }

    return numSubcycles;
}
//==================================================================================================================================
//  Setup shop at a particular point in the domain, calculating intermediate quantites
//  and updating Cantera's objects
/*
 *  All member data with the suffix, _Curr_, are updated by this function.
 *
 * @param solnElectrolyte_Curr  Current value of the solution vector
 * @param type                  Type of call
 *                              0 - at the current cell center
 */
void
porousLiIon_Anode_dom1D::SetupThermoShop1(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
{
    porosity_Curr_ = porosity_Cell_[cIndex_cc_];
    updateElectrolyte(nv, solnElectrolyte_Curr);
    updateElectrode();

    metalPhase_->setState_TP(temp_Curr_, pres_Curr_);
    metalPhase_->setElectricPotential(phiElectrode_Curr_);
    metalPhase_->getPartialMolarEnthalpies(&EnthalpyPhiPM_metal_Curr_[0]);
    EnthalpyPhiPM_metal_Curr_[0] -= Faraday * phiElectrode_Curr_;
}
//==================================================================================================================================
void
porousLiIon_Anode_dom1D::SetupThermoShop1Extra(const NodalVars* const nv, 
					       const doublereal* const solnElectrolyte_Curr)
{
    //
    // Calculate the EnthalpyPhi values at the CV interface and store these in  EnthalpyPhiPM_lyte_Curr_[]
    //
    ionicLiquid_->getPartialMolarEnthalpies(&EnthalpyPM_lyte_Curr_[0]);
    EnthalpyMolar_lyte_Curr_ = 0.0;
    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
        double z = ionicLiquid_->charge(k);
        EnthalpyPhiPM_lyte_Curr_[k] = EnthalpyPM_lyte_Curr_[k] + Faraday * z * phiElectrolyte_Curr_;
	EnthalpyMolar_lyte_Curr_ += mfElectrolyte_Soln_Curr_[k] * EnthalpyPM_lyte_Curr_[k];
    }
}
//==================================================================================================================================
//  Setup the thermo shop at a particular point in the domain, calculating intermediate quantites
//  and updating Cantera's objects
/*
 *  This routine will set up shop at a point intermediate to a left and a right point
 *  All member data with the suffix, _Curr_, are updated by this function.
 *
 * @param solnElectrolyte_CurrL  Current value of the solution vector at the left side
 * @param solnElectrolyte_CurrR  Current value of the solution vector at the right side
 * @param type                  Type of call
 *                              0 - at the left cell boundary
 *                              1 - at the right cell boundary
 */
void
porousLiIon_Anode_dom1D::SetupThermoShop2(const NodalVars* const nvL, const doublereal* const solnElectrolyte_CurrL,
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

    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
        mfElectrolyte_Soln_Curr_[k] = 0.5 * (solnElectrolyte_CurrL[indexMFL+k] +solnElectrolyte_CurrR[indexMFR+k]);
    }
    calcMFElectrolyte_Thermo(&mfElectrolyte_Soln_Curr_[0], &mfElectrolyte_Thermo_Curr_[0]); 

    size_t indexVS = nvL->indexBulkDomainVar0(Voltage);
    double phiElectrolyteL = solnElectrolyte_CurrL[indexVS];
    double phiElectrodeL = solnElectrolyte_CurrL[indexVS+1];
    indexVS = nvR->indexBulkDomainVar0(Voltage);
    double phiElectrolyteR = solnElectrolyte_CurrR[indexVS];
    double phiElectrodeR = solnElectrolyte_CurrR[indexVS+1];
    phiElectrolyte_Curr_ = 0.5 * (phiElectrolyteL + phiElectrolyteR);
    phiElectrode_Curr_ = 0.5 * (phiElectrodeL + phiElectrodeR);


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

    metalPhase_->setState_TP(temp_Curr_, pres_Curr_);
    metalPhase_->setElectricPotential(phiElectrode_Curr_);
    metalPhase_->getPartialMolarEnthalpies(&EnthalpyPhiPM_metal_Curr_[0]);
    EnthalpyPhiPM_metal_Curr_[0] -= Faraday * phiElectrode_Curr_;

    //
    // Calculate the total concentration of the electrolyte kmol m-3 and store into concTot_Curr_
    //
    concTot_Curr_ = ionicLiquid_->molarDensity();

    ionicLiquid_->getPartialMolarEnthalpies(&EnthalpyPM_lyte_Curr_[0]);
    EnthalpyMolar_lyte_Curr_ = 0.0;
    for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
        double z = ionicLiquid_->charge(k);
        EnthalpyPhiPM_lyte_Curr_[k] = EnthalpyPM_lyte_Curr_[k] + Faraday * z * phiElectrolyte_Curr_;
	EnthalpyMolar_lyte_Curr_ += mfElectrolyte_Soln_Curr_[k] * EnthalpyPM_lyte_Curr_[k];
    }

    //
    // Calculate the matrix thermal conductivity from a series resistance on the two sides
    //
    double tmp;
    if (type == 0) {
        tmp = 1.0 / thermalCond_Cell_[cIndex_cc_ - 1] + 1.0 / thermalCond_Cell_[cIndex_cc_];
    } else {
        tmp = 1.0 / thermalCond_Cell_[cIndex_cc_ + 1] + 1.0 / thermalCond_Cell_[cIndex_cc_];
    }
    thermalCond_Curr_ = 2.0 / tmp;

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
porousLiIon_Anode_dom1D::updateElectrolyte(const NodalVars* const nv, const doublereal* const solnElectrolyte_Curr)
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

    // Calculate the total concentration of the electrolyte kmol m-3.
    concTot_Curr_ = ionicLiquid_->molarDensity();

}
//=====================================================================================================================
//
//
void
porousLiIon_Anode_dom1D::updateElectrode()
{
    Electrode* Electrode_ptr = Electrode_Cell_[cIndex_cc_];
    /*
     * set the properties in the Electrode object
     *  -> temperature and pressure
     *  -> voltages of the phases
     */
    Electrode_ptr->setState_TP(temp_Curr_, pres_Curr_);
    Electrode_ptr->setVoltages(phiElectrode_Curr_, phiElectrolyte_Curr_);
    Electrode_ptr->setElectrolyteMoleNumbers(&(mfElectrolyte_Thermo_Curr_[0]), false);

    /*
     * Set the internal objects within the electrode for the final state
     */
    Electrode_ptr->updateState();
}
//=====================================================================================================================
// Retrieves the voltages from the solution vector and puts them into local storage
/*
 * @param solnElectrolyte start of the solution vector at the current node
 */
void
porousLiIon_Anode_dom1D::getVoltages(const NodalVars* const nv, const double* const solnElectrolyte)
{  
    size_t indexVS = nv->indexBulkDomainVar0(Voltage);
    phiElectrolyte_Curr_ = solnElectrolyte[indexVS];
    phiElectrode_Curr_ = solnElectrolyte[indexVS + 1];
}
//=====================================================================================================================
/*
 *  We assume that SetupThermoShop1 has been called.
 *  This calculates the Heat Capacity per cross-sectional area (Joules/K m2)
 *
 */
double
porousLiIon_Anode_dom1D::getCellHeatCapacity(const NodalVars* const nv, const double* const solnElectrolyte_Curr)
{
    Electrode* Electrode_ptr = Electrode_Cell_[cIndex_cc_];
    double cpMolar = ionicLiquid_->cp_mole();
    double lyteVol = porosity_Curr_ * xdelCell_Cell_[cIndex_cc_];
    CpMolar_lyte_Cell_[cIndex_cc_]  =  lyteVol * concTot_Curr_ * cpMolar;
    double cpSolidTotal = Electrode_ptr->SolidHeatCapacityCV();
    CpMolar_solid_Cell_[cIndex_cc_] =  cpSolidTotal / crossSectionalArea_;
    double cptotal = CpMolar_lyte_Cell_[cIndex_cc_] + CpMolar_solid_Cell_[cIndex_cc_];
    return cptotal;
}
//=====================================================================================================================
void
porousLiIon_Anode_dom1D::SetupTranShop(const double xdel, const int type)
{

    /*
     * Determine diffusion velocities
     */

    //set gradients
    gradT_trCurr_ = 0.0;
    valCellTmps& valTmps = valCellTmpsVect_Cell_[cIndex_cc_];

    if (type == 0) {
        // Left boundary
        gradV_trCurr_ = (Vcent_cc_ - Vleft_cc_) / xdel;
        gradVElectrode_trCurr_ = (VElectrodeCent_cc_ - VElectrodeLeft_cc_) / xdel;
	gradT_trCurr_ = (valTmps.Temperature.center - valTmps.Temperature.left) / xdel;

    } else {
        // Right boundary
        gradV_trCurr_ = (Vright_cc_ - Vcent_cc_) / xdel;
        gradVElectrode_trCurr_ = (VElectrodeRight_cc_ - VElectrodeCent_cc_) / xdel;  
	gradT_trCurr_ = (valTmps.Temperature.right - valTmps.Temperature.center) / xdel;
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
    // Convert from diffusion velocity to diffusion flux
    for (int k = 0; k < nsp_; k++) {
        jFlux_trCurr_[k] = mfElectrolyte_Soln_Curr_[k] * concTot_Curr_ * Vdiff_trCurr_[k];
    }

    jFlux_EnthalpyPhi_Curr_ = 0.0;
    for (size_t k = 0; k < (size_t) nsp_; ++k) {
        jFlux_EnthalpyPhi_Curr_ += jFlux_trCurr_[k] * EnthalpyPhiPM_lyte_Curr_[k];
    }

    //
    // Calculate the heat flux
    //
    heatFlux_Curr_ = - thermalCond_Curr_ * gradT_trCurr_;

    double volFSolid = (1.0 - porosity_Curr_);
    icurrElectrode_trCurr_ = - conductivityElectrode_ * pow(volFSolid, 1.5) * gradVElectrode_trCurr_;

    //
    //  Calculate =  j_electron * hphi_electron
    //
    jFlux_EnthalpyPhi_metal_trCurr_ =  - icurrElectrode_trCurr_ / Faraday;
    jFlux_EnthalpyPhi_metal_trCurr_ *= EnthalpyPhiPM_metal_Curr_[0];
   

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
porousLiIon_Anode_dom1D::saveDomain(ZZCantera::XML_Node& oNode,
                                    const Epetra_Vector* soln_GLALL_ptr,
                                    const Epetra_Vector* solnDot_GLALL_ptr,
                                    const double t,
                                    bool duplicateOnAllProcs)
{
    // get the NodeVars object pertaining to this global node
    GlobalIndices* gi = LI_ptr_->GI_ptr_;

    // Add an XML child for this domain. We will add the data to this child
    ZZCantera::XML_Node& bdom = oNode.addChild("domain");

    // Number of equations per node
    int numEquationsPerNode = BDD_ptr_->NumEquationsPerNode;

    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType>& variableNameList = BDD_ptr_->VariableNameList;

    //! First global node of this bulk domain
    int firstGbNode = BDD_ptr_->FirstGbNode;

    //! Last Global node of this bulk domain
    int lastGbNode = BDD_ptr_->LastGbNode;
    int numNodes = lastGbNode - firstGbNode + 1;

    bdom.addAttribute("id", id());
    bdom.addAttribute("points", numNodes);
    bdom.addAttribute("type", "bulk");
    bdom.addAttribute("numVariables", numEquationsPerNode);

    // Dump out the coordinates
    ZZCantera::XML_Node& gv = bdom.addChild("grid_data");

    std::vector<double> varContig(numNodes);

    int i = 0;
    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
        NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
        varContig[i] = nv->x0NodePos();
    }
    ZZctml::addNamedFloatArray(gv, "X0", varContig.size(), &(varContig[0]), "m", "length");

    for (int iVar = 0; iVar < numEquationsPerNode; iVar++) {
        VarType vt = variableNameList[iVar];
        i = 0;
        std::string nmm = vt.VariableName(200);
        for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            int istart = nv->EqnStart_GbEqnIndex;
	    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
	    if (offset == npos) {
		throw m1d_Error("porousLiIon_Separator_dom1D::saveDomain()", "cant find a variable");
	    }
            varContig[i] = (*soln_GLALL_ptr)[istart + offset];
        }
        ZZctml::addNamedFloatArray(gv, nmm, varContig.size(), &(varContig[0]), "kmol/m3", "concentration");
    }

    if (PS_ptr->doHeatSourceTracking_) {
        std::string nmm = "qSource_Cell_curr_";
        ZZctml::addNamedFloatArray(gv, nmm, numNodes, &(qSource_Cell_curr_[0]), "Joule/s/m2", "");
    }

    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++, i++) {
        int iCell = iGbNode - firstGbNode;
        Electrode* ee = Electrode_Cell_[iCell];
        ee->writeTimeStateFinal_toXML(bdom);
    }
}
//===================================================================================================================================
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
porousLiIon_Anode_dom1D::readDomain(const ZZCantera::XML_Node& SimulationNode,
				    Epetra_Vector * const soln_GLALL_ptr, Epetra_Vector * const solnDot_GLALL_ptr,
                                    double globalTimeRead)
{
    // get the NodeVars object pertaining to this global node
    GlobalIndices *gi = LI_ptr_->GI_ptr_;

    string ids = id();
    ZZCantera::XML_Node *domainNode_ptr = SimulationNode.findNameID("domain", ids);

    // Number of equations per node
    int numEquationsPerNode = BDD_ptr_->NumEquationsPerNode;

    // Vector containing the variable names as they appear in the solution vector
    std::vector<VarType> &variableNameList = BDD_ptr_->VariableNameList;

    //! First global node of this bulk domain
    int firstGbNode = BDD_ptr_->FirstGbNode;

    //! Last Global node of this bulk domain
    int lastGbNode = BDD_ptr_->LastGbNode;
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
    const ZZCantera::XML_Node* gd_ptr = (*domainNode_ptr).findByName("grid_data");

    std::vector<double> varContig(numNodes);
    ZZctml::getFloatArray(*gd_ptr, varContig, true, "", "X0");
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
       ZZctml::getFloatArray(*gd_ptr, varContig, true, "", nmm);
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
	    throw m1d_Error("porousLiIon_Anode_dom1D::readDomain() ERROR",
                            "Can't find cell number " + int2str(cellNum));
	}
	//  Read the state into the electrode object
	double timeE = ee->loadTimeStateFinal(*xCell);
	if (globalTimeRead >= 0.0) {
	    if (globalTimeRead != timeE) {
		ee->setTime(globalTimeRead);
            }
	}
    }
}




//===================================================================================================================================
// Method for writing the header for the bulk domain to a tecplot file.
void
porousLiIon_Anode_dom1D::writeSolutionTecplotHeader()
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
        int numVar = BDD_ptr_->NumEquationsPerNode;
        // Vector containing the variable names as they appear in the solution vector
        std::vector<VarType>& variableNameList = BDD_ptr_->VariableNameList;
        //! First global node of this bulk domain

        fprintf(ofp, "VARIABLES = ");
        fprintf(ofp, "\"x [m]\"  \n");

        for (int k = 0; k < numVar; k++) {
            VarType& vt = variableNameList[k];
            string name = vt.VariableName(15);
            fprintf(ofp, "\"%s\" \n", name.c_str());
        }
    
        /////////////////////////////////////////
        //// end BulkDomain1D section
        //// start application specific output
        /////////////////////////////////////////

        // print potentials within each control volume
        fprintf(ofp, "\"Potential Drop Electrode-Electrolyte [V]\" \n");
        for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++) {
            fprintf(ofp, "\"Surf_%d Equilibrium Potential Drop Electrode-Electrolyte [V]\" \n", jSurf);
            fprintf(ofp, "\"Overpotential_%d [V]\" \n", jSurf);
        }
        fprintf(ofp, "\"Current Source [A/m^2]\" \n");
        fprintf(ofp, "\"Li+ flux [kmol/m^2/s]\" \n");
        fprintf(ofp, "\"Electrolyte Current [A/m^2]\" \n");
        fprintf(ofp, "\"Electrode Current [A/m^2]\" \n");

        //  Print depth of discharge for each control volume
        fprintf(ofp, "\"CV Capacity discharged [A-s/m^2]\" \n");
        fprintf(ofp, "\"Depth of discharge in CV [A-s/m^2]\" \n");
        fprintf(ofp, "\"CV Capacity remaining [A-s/m^2]\" \n");
        fprintf(ofp, "\"Initial CV capacity [A-s/m^2]\" \n");

        // print porosity, specific surface area, thickness for each control volume
        fprintf(ofp, "\"Porosity []\" \n");
        fprintf(ofp, "\"Surface Area [m2/m2] (area per area)\" \n");
        fprintf(ofp, "\"Specific current (per particle area) [A/m^2]\" \n");
        fprintf(ofp, "\"Control volume thickness [m]\" \n");
	for (size_t i = 0; i < numExtraCondensedPhases_; i++) {
	    ExtraPhase* ep = ExtraPhaseList_[i];
	    fprintf(ofp, "\"volFrac %s []\" \t", ep->phaseName.c_str());
	}

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
                    fprintf(ofp, "\"MoleDens %s [kmol/m^3]\" \n", sss.c_str());
                }
            }
        }
	//print thermal source terms
	// check dimensions!!
        fprintf(ofp, "\"qHeat_accum [J/m3]\" \n");
        fprintf(ofp, "\"qHeat_step [W/m3]\" \n");
        fprintf(ofp, "\"Joule_Lyte [W/m3]\" \n");
        fprintf(ofp, "\"Joule_Solid [W/m3]\" \n");
        fprintf(ofp, "\"electrodHeat [W/m3]\" \n");
        fprintf(ofp, "\"OverpotHeat [W/m3]\" \n");
        fprintf(ofp, "\"deltaSHeat [W/m3]\" \n");


        fprintf(ofp, "\n");
        fclose(ofp);
    }

}
//======================================================================================================================
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
porousLiIon_Anode_dom1D::writeSolutionTecplot(const Epetra_Vector* soln_GlAll_ptr,
                                              const Epetra_Vector* solnDot_GlAll_ptr,
                                              const double t)
{
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid ; //only proc 0 should write

#define NEW_TECPLOTA
    if (doWrite) {

        double deltaTime = t_final_ - t_init_;

        // get the NodeVars object pertaining to this global node
        GlobalIndices* gi = LI_ptr_->GI_ptr_;
        // Number of equations per node
        int numVar = BDD_ptr_->NumEquationsPerNode;
        //! First global node of this bulk domain
        int firstGbNode = BDD_ptr_->FirstGbNode;
        //! Last Global node of this bulk domain
        int lastGbNode = BDD_ptr_->LastGbNode;
        int numNodes = lastGbNode - firstGbNode + 1;
#ifndef NEW_TECPLOTA
	std::vector<VarType>& variableNameList = BDD_ptr_->VariableNameList;
#endif
	std::vector<double> vars(numNodes);
        //
        //use SetupThermoShop methods to prepare for thermo model evaluation
        //

        //open tecplot file
        FILE* ofp;
        string sss = id();
        char filename[20];
        sprintf(filename,"%s%s",sss.c_str(),".dat");
        ofp = fopen(filename, "a");
#ifndef NEW_TECPLOTA
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
	    //
	    // volume fractions of control volume
	    //
	    for (size_t i = 0; i < numExtnumExtraCondensedPhases_; i++) {
		fprintf(ofp, "%g \t", volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + i];
	    }
            fprintf(ofp, "\n");

            // print the properties of the electrode particles
            Electrode* ee = Electrode_Cell_[iCell];
            std::vector<double> spmoles(nSpeciesElectrode_, 0.0);
            int solnPhase = ee->solnPhaseIndex();
            int metalPhase =  ee->metalPhaseIndex();
            int numVolPhasesE =  ee->nVolPhases();
	    //
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
                        fprintf(ofp, "%g \t", spmoles[kStart + k] / (crossSectionalArea_ * xdelCell_Cell_[iCell]));
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
#else
	fprintf(ofp, "ZONE T = \"ANOD, %10.4E\", I = %d, SOLUTIONTIME = %12.6E\n", t, numNodes, t);
	fprintf(ofp, "ZONETYPE = ORDERED\n");
	fprintf(ofp, "DATAPACKING = BLOCK\n");
	fprintf(ofp, "STRANDID = 1\n");
        // ----------------------------------------------------------------------------------------------------------------------
	//
	// Print the nodal position field
	//
	for (size_t iGbNode = (size_t) firstGbNode, i = 0; iGbNode <= (size_t) lastGbNode; ++iGbNode, ++i) {
	    NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
	    vars[i] = nv->xNodePos();
	}
	fwriteTecplotVector(ofp, vars, 13, 10);
	//
	// Print all of the solution variables
	//
	for (size_t iVar = 0; iVar < (size_t) numVar; iVar++) {
	    const VarType& vt = BDD_ptr_->VariableNameList[iVar];
	    for (size_t iGbNode = (size_t) firstGbNode, i = 0; iGbNode <= (size_t) lastGbNode; ++iGbNode, ++i) {
		NodalVars *nv = gi->NodalVars_GbNode[iGbNode];
		size_t istart = nv->EqnStart_GbEqnIndex;
		size_t ioffset = nv->indexBulkDomainVar(vt);
		if (ioffset == npos) {
		    throw m1d_Error("", "index prob");
		}
		vars[i] = (*soln_GlAll_ptr)[istart + ioffset];
	    }
	    fwriteTecplotVector(ofp, vars, 13, 10);
	}
	//
	//  Voltage of each CV = phi_Metal - phi_lyte (units = volts)
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = deltaV_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	for (int jSurf = 0; jSurf < nSurfsElectrode_ ; jSurf++) {
	    //
	    //  Open circuit potential
	    //
	    for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
		vars[iCell] = Ess_Surf_Cell_[nSurfsElectrode_ * iCell + jSurf];
	    }
	    fwriteTecplotVector(ofp, vars, 13, 10);
	    //
	    // Overpotential of each surface reaction
	    //
	    for (size_t iCell = 0; iCell < (size_t)  NumLcCells;  ++iCell) {
		vars[iCell] = overpotential_Surf_Cell_[nSurfsElectrode_ * iCell + jSurf];
	    }
	    fwriteTecplotVector(ofp, vars, 13, 10);
	}
	//
	//   Units are amps / m2 and the area is the cross sectional area of the electrode
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = icurrRxn_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Li+ flux [kmol/m^2/s]
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = LiFlux_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	//  Current electrolyte -  Units are amps / m2 and the area is the cross sectional area of the electrode
	//        (these are CB quantities, but are being output at nodes
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    if (iCell == 0) {
		vars[iCell] = icurrElectrolyte_CBL_[iCell];
	    } else if (iCell == (size_t) (NumLcCells - 1)) {
		vars[iCell] = icurrElectrolyte_CBR_[iCell];
	    } else {
		vars[iCell] = 0.5 * (icurrElectrolyte_CBL_[iCell] + icurrElectrolyte_CBR_[iCell]);
	    }
	}
	fwriteTecplotVector(ofp, vars, 13);

	//
	//  Current in electrode phase -  Units are amps / m2 and the area is the cross sectional area of the electrode
	//        (these are CB quantities, but are being output at nodes
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    if (iCell == 0) {
		vars[iCell] = icurrElectrode_CBL_[iCell];
	    } else if (iCell == (size_t) (NumLcCells - 1)) {
		vars[iCell] = icurrElectrode_CBR_[iCell];
	    } else {
		vars[iCell] = 0.5 * (icurrElectrode_CBL_[iCell] + icurrElectrode_CBR_[iCell]);
	    }
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	//  Capacity discharged in control volume per cross sectional area (amps-sec / m2)
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = capacityDischargedPA_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	//  Depth of discharge in cross sectional area (amps-sec / m2)
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = depthOfDischargePA_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	//  Capacity left in control volume per area (amps-sec / m2)
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = capacityLeftPA_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	//  Capacity in control volume per cross sectional area (amps-sec / m2)
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = capacityPA_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Porosity of control volume
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = porosity_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Surface area of control volume  units = m2 / m2. surface area per cross sectional area in each CV
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = surfaceArea_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Current per interfacial surface area in each CV units - amps / m2
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = icurrInterfacePerSurfaceArea_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Thickness of each control volume units = m
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = xdelCell_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// volume fractions of control volume
	//
	for (size_t i = 0; i < numExtraCondensedPhases_; i++) {
	    for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
		vars[iCell] = volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + i];
	    }
	    fwriteTecplotVector(ofp, vars, 13);
        }

	//
	// Gather the species moles in the electrode object
	//  Output moles of species -> 
	//
	Electrode* ee;
	std::vector<double> spmoles_Cell(nSpeciesElectrode_*NumLcCells, 0.0);
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    ee = Electrode_Cell_[iCell];
	    ee->getMoleNumSpecies(&spmoles_Cell[nSpeciesElectrode_*iCell]);
	}
	int solnPhase = ee->solnPhaseIndex();
	int metalPhase =  ee->metalPhaseIndex();
	size_t numVolPhasesE =  ee->nVolPhases();
	for (size_t vph = 0; vph < numVolPhasesE; vph++) {
	    ThermoPhase* tp = &(ee->volPhase(vph));
	    int iph = ee->getGlobalPhaseIndex(tp);
	    if (iph == metalPhase || iph == solnPhase) {
		
	    } else {
		int nspPhase = tp->nSpecies();
		int kStart = ee->getGlobalSpeciesIndex(iph, 0);
		for (size_t k = 0; k < (size_t) nspPhase; k++) {
		    for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
                        if (spmoles_Cell.size() -1 < nSpeciesElectrode_*iCell + kStart + k) {
                            exit(-1); 
                        }
                        if (kStart < 0) {
                            exit(-1);
                        }
			vars[iCell] = 
                          spmoles_Cell[nSpeciesElectrode_*iCell + kStart + k] / (crossSectionalArea_ * xdelCell_Cell_[iCell]);
		    }
		    fwriteTecplotVector(ofp, vars, 13);
		}
	    }
	}
	//
	// Total heat accumulated units = Joule/s/m3
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = qSource_Cell_accumul_[iCell] / xdelCell_Cell_[iCell];
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Total heat this step, units = Joule/s/m3
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = qSource_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Joule heating electrolyte, this step, units = Joule/s/m3
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = jouleHeat_lyte_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Joule heating solid, this step, units = Joule/s/m3
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = jouleHeat_solid_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Electrode Heat, this step, units = Joule/s/m3
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = electrodeHeat_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Irreversible electrode Heat, this step, units = Joule/s/m3
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = overPotentialHeat_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	}
	fwriteTecplotVector(ofp, vars, 13);
	//
	// Reversible electrode Heat, this step, units = Joule/s/m3
	//
	for (size_t iCell = 0; iCell < (size_t) NumLcCells;  ++iCell) {
	    vars[iCell] = deltaSHeat_Cell_curr_[iCell] / xdelCell_Cell_[iCell] / deltaTime;
	}
	fwriteTecplotVector(ofp, vars, 13);
#endif


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
            int indexVS = BDD_ptr_->VariableIndexStart_VarName[Voltage];
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
porousLiIon_Anode_dom1D::showSolution(const Epetra_Vector* soln_GlAll_ptr,
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

    std::vector<VarType>& variableNameList = BDD_ptr_->VariableNameList;
    int iBlock;
    int iGbNode;
    int n;
    int nPoints = BDD_ptr_->LastGbNode - BDD_ptr_->FirstGbNode + 1;

    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind = indent.c_str();
    doublereal v;
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
    // Number of points in each vector
    string sss = id();

    stream0 ss;
    print0_sync_start(0, ss, * (LI_ptr_->Comm_ptr_));
    if (do0Write) {
        drawline0(indentSpaces, 80);
        ss.print0("%s  Solution on Bulk Domain %12s : Number of variables = %d\n", ind, sss.c_str(), NumDomainEqns);
        ss.print0("%s                                         : Number of Nodes = %d\n", ind, nPoints);
        ss.print0("%s                                         : Beginning pos %g\n", ind, BDD_ptr_->Xpos_start);
        ss.print0("%s                                         : Ending    pos %g\n", ind, BDD_ptr_->Xpos_end);

#ifdef MECH_MODEL
	int firstGbNode = BDD_ptr_->FirstGbNode;
	int lastGbNode = BDD_ptr_->LastGbNode;
	if (solidMechanicsProbType_ > 0) {
	    for (int iGbNode = firstGbNode; iGbNode <= lastGbNode; iGbNode++) {
		NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
		double diff = nv->xNodePos() -nv->x0NodePos();
		ss.print0("%s                                         : Node[ %d ] Displacement %g\n", ind, iGbNode,diff);
	    }
	}
#endif 
    }


    // \todo remove this. 
#ifdef DEBUG_MECH_MODEL 
    if (do0Write) for (int iCell = 0; iCell < NumLcCells; iCell++) {
      Electrode* ee = Electrode_Cell_[iCell];
      ss.print0("[%d ]  SolidVolume %g \n",iCell,ee->SolidVol());
    }
#endif 
    

    if (do0Write) {
        for (iBlock = 0; iBlock < nn; iBlock++) {
            drawline0(indentSpaces, 80);
            ss.print0("%s iGbNode  z   ", ind);
            for (n = 0; n < 5; n++) {
                int ivar = iBlock * 5 + n;
                VarType vt = variableNameList[ivar];
                string name = vt.VariableName(15);
                ss.print0(" %15s", name.c_str());
            }
            ss.print0("\n");
            drawline(indentSpaces, 80);
            double oldVaxial = 0.0;
            double newVaxial = 0.0;

            for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {

                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
                doublereal x = nv->xNodePos();
                ss.print0("\n%s %4d  % -11.4E", ind, iGbNode, x);
                int istart = nv->EqnStart_GbEqnIndex;
                for (n = 0; n < 5; n++) {
                    int ivar = iBlock * 5 + n;
                    VarType vt = variableNameList[ivar];
		    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
                    if (offset == npos) {
                        throw m1d_Error("porousLiIon_Separator_dom1D::showSolution()", "cant find a variable");
                    }
		    v = (*soln_GlAll_ptr)[istart + offset];

                    /*
                     *  Special case the axial velocity because it's not located at the nodes.
                     */
                    if (vt.VariableType == Velocity_Axial) {
                        newVaxial = v;
                        if (iGbNode == BDD_ptr_->FirstGbNode) {
                            // we've applied a boundary condition here!
                            v = 0.0;
			} else if (iGbNode == BDD_ptr_->LastGbNode) {
			    cellTmps& cTmps = cellTmpsVect_Cell_[NumLcCells-1];
	                    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
                            // we've applied a boundary condition here!
                            v = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
			    newVaxial = v;
                        } else {
                            v = 0.5 * (oldVaxial + newVaxial);
                        }
                        oldVaxial = newVaxial;
                    }

                    ss.print0(" % -11.4E", v);
                }
            }
            ss.print0("\n");
        }

        int nrem = NumDomainEqns - 5 * nn;
        if (nrem > 0) {
            drawline(indentSpaces, 80);
            ss.print0("%s  iGbNode  z   ", ind);
            for (n = 0; n < nrem; n++) {
                int ivar = nn * 5 + n;
                VarType vt = variableNameList[ivar];
                string name = vt.VariableName(15);
                ss.print0(" %15s", name.c_str());
            }
            ss.print0("\n");
            drawline0(indentSpaces, 80);
            double oldVaxial = 0.0;
            double newVaxial = 0.0;

            for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
                doublereal x = nv->xNodePos();
                ss.print0("%s %4d   % -11.4E ", ind, iGbNode, x);
                int istart = nv->EqnStart_GbEqnIndex;
                for (n = 0; n < nrem; n++) {
                    int ivar = iBlock * 5 + n;
                    VarType vt = variableNameList[ivar];
		    size_t offset = nv->indexBulkDomainVar(vt.VariableType, vt.VariableSubType);
                    if (offset == npos) {
                        throw m1d_Error("porousLiIon_Separator_dom1D::showSolution()", "cant find a variable");
                    }
		    v = (*soln_GlAll_ptr)[istart + offset];
		    /*
                     *  Special case the axial velocity because it's not located at the nodes.
                     */
		    if (vt.VariableType == Velocity_Axial) {
                        newVaxial = v;
                        if (iGbNode == BDD_ptr_->FirstGbNode) {
                            // we've applied a boundary condition here!
                            v = 0.0;
			} else if (iGbNode == BDD_ptr_->LastGbNode) {
			    cellTmps& cTmps = cellTmpsVect_Cell_[NumLcCells-1];
	                    NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
                            // we've applied a boundary condition here!
                            v = DiffFluxRightBound_LastResid_NE[nodeTmpsCenter.RO_Electrolyte_Continuity];
			    newVaxial = v;
                        } else {
                            v = 0.5 * (oldVaxial + newVaxial);
                        }
                        oldVaxial = newVaxial;
                    }

                    ss.print0(" % -11.4E ", v);
                }
                ss.print0("\n");
            }
        }
        drawline(indentSpaces, 80);
    }
    print0_sync_end(0, ss, * (LI_ptr_->Comm_ptr_));
    print0_sync_start(0, ss, * (LI_ptr_->Comm_ptr_));
    if (do0Write) {
        // ----------------------------------------------------
        // --             print surface potentials within the cell --
        // ----------------------------------------------------
        ss.print0("\n");
        drawline0(indentSpaces, 80);
        ss.print0("%s  iGbNode  z       Delta_V        Ess    Overpotential icurrRxnCell", ind);
        ss.print0("\n");
        drawline0(indentSpaces, 80);
    }
    print0_sync_end(0, ss, * (LI_ptr_->Comm_ptr_));
    doublereal x;
    int iCell;
    for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
            iCell = iGbNode - BDD_ptr_->FirstGbNode;
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            x = nv->xNodePos();
            ss.print0("%s %4d  % -11.4E ", ind, iGbNode, x);
            ss.print0("% 11.4E ", deltaV_Cell_[iCell]);
            ss.print0("% 11.4E ", Ess_Surf_Cell_[nSurfsElectrode_ * iCell]);
            ss.print0("% 11.4E ", overpotential_Surf_Cell_[nSurfsElectrode_ * iCell]);
            ss.print0("% 11.4E ", icurrRxn_Cell_[iCell]);
            ss.print0("\n");
        }
        print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    }

    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (do0Write) {
        // -----------------------------------------------------------------------------------------------------------------
        // --             print porosity within the cell --
        // -----------------------------------------------------------------------------------------------------------------
        ss.print0("\n");
        drawline0(indentSpaces, 80);
        ss.print0("%s  iGbNode  z      Porosity    VF_Electrode", ind);
	for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
	    ExtraPhase* ep = ExtraPhaseList_[jPhase];
	    ss.print0("%-11.11s ", ep->phaseName.c_str());
	}
        ss.print0("\n");
        drawline0(indentSpaces, 80);
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
            iCell = iGbNode - BDD_ptr_->FirstGbNode;
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
	    double volCell =  crossSectionalArea_ * xdelCell_Cell_[iCell];
            x = nv->xNodePos();
            ss.print0("%s %4d % -11.4E ", ind, iGbNode, x);
            ss.print0("%11.4E ", porosity_Cell_[iCell]);
	    double vole = nVol_zeroStress_Electrode_Cell_[iCell];
	    double vfE = vole / volCell;
	    ss.print0("% -11.4E ", vfE);
	    for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
		ss.print0("% -11.4E ", volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + jPhase]);
	    }
            ss.print0("\n");
        }
        print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    }

    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (do0Write) {
        // -----------------------------------------------------------------------------------------------------------------
        // --             PRINT SURFACE REACTION DETAILS ABOUT EACH CELL --
        // -----------------------------------------------------------------------------------------------------------------
        ss.print0("\n");
        drawline0(indentSpaces, 80);
        ss.print0("%s   iGbNode  z     SurfaceArea  SurfAreaDens  currPerSA  DeltaZ", ind);
        ss.print0("\n");
        drawline0(indentSpaces, 80);
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
            iCell = iGbNode - BDD_ptr_->FirstGbNode;
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            x = nv->xNodePos();
            ss.print0("%s %4d  % -10.4E ", ind, iGbNode, x);
            ss.print0("% 11.4E ", surfaceArea_Cell_[iCell]);
            ss.print0("% 11.4E ", surfaceArea_Cell_[iCell] / xdelCell_Cell_[iCell]);
            ss.print0("% 11.4E ", icurrInterfacePerSurfaceArea_Cell_[iCell]);
            ss.print0("% 11.4E ", xdelCell_Cell_[iCell]);
            ss.print0("\n");
        }
        print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    }

    print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
    if (do0Write) {
        // -----------------------------------------------------------------------------------------------------------------
        // --             PRINT DEPTH OF DISCHARGE VALUES FOR EACH CELL --
        // -----------------------------------------------------------------------------------------------------------------
        ss.print0("\n");
        drawline0(indentSpaces, 80);
        ss.print0("%s        z    capDischarged DepthDischged capLeft    capZeroDOD  numSubSteps", ind);
        ss.print0("\n");
        ss.print0("%s       (m)      (coul/m2)   (coul/m2)   (coul/m2)   (coul/m2)", ind);
        ss.print0("\n");
        drawline0(indentSpaces, 80);
    }
    print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
    for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
            iCell = iGbNode - BDD_ptr_->FirstGbNode;
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            x = nv->xNodePos();

            ss.print0("%s    % -11.4E ", ind, x);
            ss.print0("% -11.4E ", capacityDischargedPA_Cell_[iCell]);
            ss.print0("% -11.4E ", depthOfDischargePA_Cell_[iCell]);
            ss.print0("% -11.4E ", capacityLeftPA_Cell_[iCell]);
            ss.print0("% -11.4E ", capacityPA_Cell_[iCell]);
	    ss.print0("  %6d  ", numElectrodeSubCycles_Cell_[iCell]);

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
        ss.print0("%s   iGbNode  z       ", ind);

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
    for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
            iCell = iGbNode - BDD_ptr_->FirstGbNode;
            NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
            x = nv->xNodePos();
            ss.print0("%s %4d  %-10.4E ", ind, iGbNode, x);
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

#ifdef MECH_MODEL
    //    if (do0Write) {
    //   drawline0(indentSpaces, 80);
    //   ss.print0("%s    CellBound    SolidStressAxial ", ind);
    //   ss.print0("\n");
    //   drawline(indentSpaces, 80);

    //   const Epetra_Vector& soln = *soln_ptr;
     
    //   for (int iCell = 1; iCell < NumLcCells; iCell++) {
    //   	int index_CentLcNode = Index_DiagLcNode_LCO[iCell];
    //   	int indexCent_EqnStart = LI_ptr_->IndexLcEqns_LcNode[index_CentLcNode];
    //   	cellTmps& cTmps          = cellTmpsVect_Cell_[iCell];
    //   	NodeTmps& nodeTmpsCenter = cTmps.NodeTmpsCenter_;
    //   	indexCent_EqnStart = nodeTmpsCenter.index_EqnStart;
    //   	double Solid_Stress_Axial = soln[indexCent_EqnStart + nodeTmpsCenter.Offset_Solid_Stress_Axial];
    //   	ss.print0("%s %d-%d   % -11.4E ",ind,iCell-1,iCell, Solid_Stress_Axial);
    //   	ss.print0("\n");
    //   }
      
    // }
#endif 

    if (do0Write) {
        // ----------------------------------------------------
        // --             PRINT FLUXES AT THE CELL BOUNDARIES --
        // ----------------------------------------------------
        ss.print0("\n");
        drawline0(indentSpaces, 80);
        ss.print0("%s    CellBound    z     IcurrElectrolyte IcurrElectrode ", ind);
        ss.print0("\n");
        drawline(indentSpaces, 80);
    }



    NodalVars* nvl;
    NodalVars* nvr;

    for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
        print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
        if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {

            ss.print0("%s    ", ind);
            if (iGbNode == BDD_ptr_->FirstGbNode) {

                iCell = 0;
                nvr = gi->NodalVars_GbNode[BDD_ptr_->FirstGbNode];
                x = nvr->xNodePos();
                ss.print0("Lft-0     % -11.4E ", x);
                ss.print0("% -11.4E ", icurrElectrolyte_CBL_[0]);
                ss.print0("% -11.4E \n", icurrElectrode_CBL_[0]);
                ss.print0("%s    ", ind);
                iCell = iGbNode - BDD_ptr_->FirstGbNode;
                nvl = gi->NodalVars_GbNode[iGbNode];
                nvr = gi->NodalVars_GbNode[iGbNode + 1];
                x = 0.5 * (nvl->xNodePos() + nvr->xNodePos());
                ss.print0("%3d-%-3d   % -11.4E ", iCell, iCell + 1, x);
                ss.print0("% -11.4E ", icurrElectrolyte_CBR_[iCell]);
                ss.print0("% -11.4E ", icurrElectrode_CBR_[iCell]);
            } else if (iGbNode == BDD_ptr_->LastGbNode) {
                iCell = BDD_ptr_->LastGbNode - BDD_ptr_->FirstGbNode;
                nvr = gi->NodalVars_GbNode[BDD_ptr_->LastGbNode];
                x = nvr->xNodePos();
                ss.print0("%3d-Rgt   % -11.4E ", iCell, x);
                ss.print0("% -11.4E ", icurrElectrolyte_CBR_[iCell]);
                ss.print0("% -11.4E ", icurrElectrode_CBR_[iCell]);
            } else {
                iCell = iGbNode - BDD_ptr_->FirstGbNode;
                nvl = gi->NodalVars_GbNode[iGbNode];
                nvr = gi->NodalVars_GbNode[iGbNode + 1];
                x = 0.5 * (nvl->xNodePos() + nvr->xNodePos());
                ss.print0("%3d-%-3d   % -11.4E ", iCell, iCell + 1, x);
                ss.print0("% -11.4E ", icurrElectrolyte_CBR_[iCell]);
                ss.print0("% -11.4E ", icurrElectrode_CBR_[iCell]);
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
        for (iGbNode = BDD_ptr_->FirstGbNode; iGbNode <= BDD_ptr_->LastGbNode; iGbNode++) {
            print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
            if (iGbNode >= FirstOwnedGbNode && iGbNode <= LastOwnedGbNode) {
                iCell = iGbNode - BDD_ptr_->FirstGbNode;
                NodalVars* nv = gi->NodalVars_GbNode[iGbNode];
                x = nv->xNodePos();
                ss.print0("%s   % -11.4E ", ind, x);

                ss.print0("% -11.4E ",  qSource_Cell_curr_[iCell]);

                ss.print0("% -11.4E ",  qSource_Cell_accumul_[iCell]);
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
porousLiIon_Anode_dom1D::setStateFromSolution(const bool doTimeDependentResid, const Epetra_Vector_Ghosted *soln_ptr, 
					      const Epetra_Vector_Ghosted *solnDot,
					      const double t, const double delta_t, const double t_old)
{
    bool doAll = true;
    //bool doInit = false;
    //bool doFinal = true; // This is always true as you can only set the final Electrode object
    size_t  indexCent_EqnStart;
    if (doTimeDependentResid) {
	doAll = false;
	if (delta_t > 0.0) {
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
	//
	// Find the start of the solution at the current node
	//
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
	// For cases where we aren't in time dependent mode set all states
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
 */
void
porousLiIon_Anode_dom1D::initialConditions(const bool doTimeDependentResid,  Epetra_Vector* soln_ptr,
                                           Epetra_Vector* solnDot,  const double t,
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
	NodalVars* nodeLeft = nodeCent;
	NodalVars* nodeRight = nodeCent;
	if (iCell !=0) {
	    int index_LeftLcNode = Index_LeftLcNode_LCO[iCell];
	    nodeLeft =  LI_ptr_->NodalVars_LcNode[ index_LeftLcNode];
	}
	if ((int) iCell !=NumLcCells - 1) {
	    int index_RightLcNode = Index_RightLcNode_LCO[iCell];
	    nodeRight = LI_ptr_->NodalVars_LcNode[index_RightLcNode];
	}
	double cellThickness = 0.5*(nodeRight->xNodePos() - nodeLeft->xNodePos());

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
	// Set the temperature if it is part of the solution vector: We get the temperature
	//     from the Reference Temperature card in the input deck.
	//	
	temp_Curr_ = PSinput.TemperatureReference_;
        if (iVar_Temperature != npos) {
            soln[indexCent_EqnStart + iVar_Temperature] = PSinput.TemperatureReference_;
        }
	Temp_Cell_old_[iCell] = temp_Curr_;

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
	for (size_t k = 0; k < BDT_ptr_->nSpeciesElectrolyte_; ++k) {
	    soln[indexCent_EqnStart + iVar_Species + k] = PSinput.electrolyteMoleFracs_[k];
	}

	//
	// Set some arbitrary values for the two voltages
	//
        soln[indexCent_EqnStart + iVar_Voltage] = -0.07;
        soln[indexCent_EqnStart + iVar_Voltage_ED] = 0.0;

	SetupThermoShop1(nodeCent, &(soln[indexCent_EqnStart]));
	//
        //  Fill in phiElectroyte_Curr_ and phiElectrode_Curr_
        //
        getVoltages(nodeCent, solnCentStart);
	//
        //  fill in mfElectrolyte_Soln_Curr[]  mfElectrolyte_Thermo_Curr_[]
        //
        getMFElectrolyte_soln(nodeCent, solnCentStart);
        //
        // Need to update the Electrode objects with the state of the solution
        //
        Electrode* ee = Electrode_Cell_[iCell];
	//
        //  Set the temperature and pressure and voltages in the final_ state
	//      (NOTE, we are not setting the internal state here! This may have to change)
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
        // Porosity Setup - We have essentially already done this. However, just to make sure that
	//                  we don't have unknown errors, we'll calculate the porosity given the
	//                  current conditions and the store it in the porosity_Cell_[] vector
        //
	double solidVolCell = ee->SolidVol();
        int offS = 0;
        double vfe = solidVolCell / (xdelCell_Cell_[iCell] * crossSectionalArea_);
        porosity_Curr_ = 1.0 - vfe;
	//double domainThickness = BDT_ptr_->Xpos_end - BDT_ptr_->Xpos_start;

        for (size_t k = 0; k < ExtraPhaseList_.size(); ++k) {
	    ExtraPhase* ep = ExtraPhaseList_[k];
	    ThermoPhase* tp = ep->tp_ptr;
	    tp->setState_TP(temp_Curr_, pres_Curr_);
	    double mvp = tp->molarVolume();
	    volumeFraction_Phases_Cell_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction;
	    volumeFraction_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction;
	    moleNumber_Phases_Cell_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction* cellThickness * mvp;
	    moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + offS + k] = ep->volFraction* cellThickness * mvp;
	    porosity_Curr_ -= ep->volFraction;
	}

        if (porosity_Curr_ <= 0.0) {
            throw m1d_Error("porousLiIon_Anode_dom1D::initialConditions() ERROR",
                            "Calculated porosity, " + fp2str(porosity_Curr_) + ", is less than zero\n"
                            "           There is an error in the setup of the anode");
        }    
	//
        // update porosity as computed from electrode input
	//
        porosity_Cell_[iCell] = porosity_Curr_;
        porosity_Cell_old_[iCell] = porosity_Curr_;
          
	if (energyEquationProbType_ == 3) { 
	    double volCellNew = xdelCell_Cell_[iCell];
	    // double volElectrodeCell = solidVolCell / crossSectionalArea_;
	    double solidEnthalpy = ee->SolidEnthalpy() / crossSectionalArea_;
	    double solidEnthalpyNew = solidEnthalpy;
	  
	    double lyteMolarEnthalpyNew = ionicLiquid_->enthalpy_mole();
	    double volLyteNew = porosity_Curr_ * volCellNew;
	    double lyteEnthalpyNew =  lyteMolarEnthalpyNew * concTot_Curr_ * volLyteNew;
	    double nEnthalpyInertNew = 0.0;
	    for (size_t jPhase = 0; jPhase < numExtraCondensedPhases_; ++jPhase) {
                ExtraPhase* ep = ExtraPhaseList_[jPhase];
                ThermoPhase* tp = ep->tp_ptr;
                tp->setState_TP(temp_Curr_, pres_Curr_);
                double mEnth = tp->enthalpy_mole();
                double mn = moleNumber_Phases_Cell_old_[numExtraCondensedPhases_ * iCell + jPhase];
                nEnthalpyInertNew += mEnth * mn;
            }


	    double nEnthalpy_New  = solidEnthalpyNew + lyteEnthalpyNew +   nEnthalpyInertNew;

	    nEnthalpy_Old_Cell_[iCell] = nEnthalpy_New; 
	    nEnthalpy_New_Cell_[iCell] = nEnthalpy_New;

	    nEnthalpy_Electrode_New_Cell_[iCell] = solidEnthalpy;
	    nEnthalpy_Electrode_Old_Cell_[iCell] = solidEnthalpy;

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
void porousLiIon_Anode_dom1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
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
        atolVector[indexCent_EqnStart + iVar_Voltage] = 1.0E-6;
        atolVector[indexCent_EqnStart + iVar_Voltage_ED] = 1.0E-6;
	/*
         * Set the tolerance on the temperature
         */
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        if (iVar_Temperature != npos) {
            atolVector[indexCent_EqnStart + iVar_Temperature] = 1.0E-7;
        }

    }
}
//===================================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void porousLiIon_Anode_dom1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted& soln,
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
        atolVector[indexCent_EqnStart + iVar_Voltage] = 1.0E-6;
        atolVector[indexCent_EqnStart + iVar_Voltage_ED] = 1.0E-6;
	/*
         * Set the tolerance on the temperature
         */
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        if (iVar_Temperature != npos) {
            atolVector[indexCent_EqnStart + iVar_Temperature] = 1.0E-7;
        }
    }
}
//======================================================================================================================
void
porousLiIon_Anode_dom1D::setAtolDeltaDamping(double atolDefault, double relcoeff,
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
	/*
         * Set the tolerance on the temperature
         */
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        if (iVar_Temperature != npos) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Temperature] = 5.0;
        }

	/*
         * Set the tolerance on the displacment
	 *   -> turning it off for the moment.
         */
        size_t iVar_disp = nodeCent->indexBulkDomainVar0((size_t) Displacement_Axial);
        if (iVar_disp != npos) {
            atolDeltaDamping[indexCent_EqnStart + iVar_disp] = 5.0 * relcoeff;
        }

	/*
         * Set the tolerance on the Velocity_Radial
         */
        size_t iVar_Vrad = nodeCent->indexBulkDomainVar0((size_t) Velocity_Radial);
        if (iVar_Vrad != npos) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Vrad] = 5.0 * relcoeff;
        }


    }
}
//======================================================================================================================
void
porousLiIon_Anode_dom1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff,
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
	/*
         * Set the tolerance on the temperature
         */
        size_t iVar_Temperature = nodeCent->indexBulkDomainVar0((size_t) Temperature);
        if (iVar_Temperature != npos) {
            atolDeltaDamping[indexCent_EqnStart + iVar_Temperature] = 5.0;
        }
    }
}
//=====================================================================================================================
/**
 * Method to check for precipitation of the salts.
 * Returns index of offending cation or -1 if no precipitation
 */
int
porousLiIon_Anode_dom1D::checkPrecipitation()
{
    return -1;
}
//=====================================================================================================================
void
porousLiIon_Anode_dom1D::err(const char* msg)
{
    printf("porousLiIon_Anode_dom1D: function not implemented: %s\n", msg);
    exit(-1);
}
//=====================================================================================================================
double porousLiIon_Anode_dom1D::capacityPA(int platNum) const
{
    double totalCapacity = 0.0;
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
	totalCapacity += ee->capacity(platNum);
        capacityPA_Cell_[iCell] = ee->capacity() / crossSectionalArea_;
    }
    totalCapacity /= crossSectionalArea_;
    return totalCapacity;
}
//=====================================================================================================================
double porousLiIon_Anode_dom1D::capacityDischargedPA(int platNum) const
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
double porousLiIon_Anode_dom1D::capacityLeftPA(int platNum, double voltsMax, double voltsMin) const
{
    double totalCapacity = 0.0;
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
	totalCapacity += ee->capacityLeft(platNum, voltsMax, voltsMin);
        capacityLeftPA_Cell_[iCell] = ee->capacityLeft() / crossSectionalArea_;
    }
    totalCapacity /= crossSectionalArea_;
    return totalCapacity;
}
//=====================================================================================================================
double porousLiIon_Anode_dom1D::depthOfDischargePA(int platNum) const
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
// Initial starting depth of discharge in coulombs per cross sectional area
/*
 *   When there is capacity lost, this number may be modified.
 *
 *   @param platNum  Plateau number. Default is -1 which treats all plateaus as a single entity.
 *                   If positive or zero, each plateau is treated as a separate entity.
 */
double porousLiIon_Anode_dom1D::depthOfDischargeStartingPA(int platNum) const
{
    double dodStarting = 0.0;
    for (int iCell = 0; iCell < NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
	dodStarting += ee->depthOfDischargeStarting(platNum);
    }
    dodStarting /= crossSectionalArea_;
    return dodStarting;
}
//=====================================================================================================================
void porousLiIon_Anode_dom1D::resetCapacityDischargedToDate()
{
    for (size_t iCell = 0; iCell < (size_t) NumLcCells; iCell++) {
        Electrode* ee = Electrode_Cell_[iCell];
        ee->resetCapacityDischargedToDate();
    }
}
//=====================================================================================================================
// Return a value for the open circuit potential without doing a formally correct calculation
/*
 *  Currently this is defined as the open circuit potential on the outside electrode.
 *
 *   @return return the open circuit potential 
 */
double porousLiIon_Anode_dom1D::openCircuitPotentialQuick() const
{
    Electrode* ee = Electrode_Cell_[0];
    double ocv = ee->openCircuitVoltage(nSurfsElectrode_ - 1);
    return ocv;
}
//=====================================================================================================================
//  WARNING -> will fail for multiprocessors, becauase the accumulation of information within eval_PostProc will fail.
double porousLiIon_Anode_dom1D::effResistanceLayer(double &potAnodic, double &potCathodic,
			                           double &voltOCV, double &current)
{
    static double resistance = 0.0;
    potAnodic = potentialAnodic_;
    potCathodic = potentialCathodic_;
    current = icurrElectrode_CBL_[0];
    voltOCV = openCircuitPotentialQuick();
    //
    //  Calculate the effective resistance of the separator layer by taking the potential drop across the domain
    //  and dividing it by the current.
    //
    resistance = 0.0;
    if (fabs(current) > 1.0E-200) {
        resistance = ((potAnodic - potCathodic) - voltOCV) / current;
    }
    return resistance;
}
//=====================================================================================================================
} //namespace m1d
//=====================================================================================================================


