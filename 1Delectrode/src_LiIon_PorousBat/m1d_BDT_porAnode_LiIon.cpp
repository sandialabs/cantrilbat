/**
 * @file m1d_BDT_porAnode_LiKCl.cpp
 */

/*
 *  $Id: m1d_BDT_porAnode_LiIon.cpp 552 2013-03-01 21:25:03Z hkmoffa $
 */

#include "m1d_BDT_porAnode_LiIon.h"
#include "m1d_porousLiIon_Anode_dom1D.h"
#include "Electrode_Factory.h"
#include "m1d_ProblemStatementCell.h"
extern m1d::ProblemStatementCell PSinput;

using namespace std;
using namespace Cantera;

namespace m1d
{

//====================================================================================================================
BDT_porAnode_LiIon::BDT_porAnode_LiIon(DomainLayout* dl_ptr) :
    BulkDomainDescription(dl_ptr),
    ionicLiquid_(0),
    trans_(0),
    m_position(0),
    Electrode_(0),
    Capacity_initial_(0.0),
    CapacityDischarged_initial_(0.0),
    CapacityLeft_initial_(0.0)
{
    int eqnIndex = 0;
    IsAlgebraic_NE.resize(7,0);
    IsArithmeticScaled_NE.resize(7,0);
    /*
     * Store a copy of the electrolyte ThermoPhase object
     */
    int iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.electrolytePhase_);
    if (iph < 0) {
        throw CanteraError("BDT_porAnode_LiIon::BDT_porAnode_LiIon()",
                           "Can't find the phase in the phase list: " + PSinput.electrolytePhase_);
    }
    ThermoPhase* tmpPhase = & (PSinput.PhaseList_)->thermo(iph);
    ionicLiquid_ = tmpPhase->duplMyselfAsThermoPhase();

    /*
     *  Create and Store a pointer to the Transport Manager
     */
    trans_ = Cantera::newTransportMgr("Simple", ionicLiquid_, 1);

    /*
     *  Find the hook into the input for the cathode electrode model
     */
    ProblemStatementCell* psc_ptr = &PSinput;
    ELECTRODE_KEY_INPUT* ai = psc_ptr->anode_input_;

    /*
     *  Use the ElectrodeModelName value as input to the electrode factory to create the electrode
     */
    Electrode_  = newElectrodeObject(ai->electrodeModelName);
    if (!Electrode_) {
        throw  m1d_Error("BDT_porAnode_LiIon::instantiateElectrodeCells()",
                         "Electrode factory method failed");
    }
    ELECTRODE_KEY_INPUT* ai_new = newElectrodeKeyInputObject(ai->electrodeModelName);
    string commandFile = ai->commandFile_;
    BEInput::BlockEntry* cfA = new BEInput::BlockEntry("command_file");
    /*
     *  Parse the complete child input file
     */
    int retn = ai_new->electrode_input_child(commandFile, cfA);
    if (retn == -1) {
        throw  m1d_Error("BDT_porAnode_LiIon::BDT_porAnode_LiKCl()",
                         "Electrode input child method failed");
    }
    /*
     * Switch the pointers around so that the child input file is returned.
     * Delete the original pointer.
     */
    delete ai;
    PSinput.anode_input_ = ai_new;

    /*
     *   Initialize the electrode model using the input from the ELECTRODE_KEY_INPUT object
     */
    retn = Electrode_->electrode_model_create(PSinput.anode_input_);
    if (retn == -1) {
        throw CanteraError("BDT_porAnode_LiIon::BDT_porAnode_LiIon()",
                           "Error initializing the anode electrode object");
    }
    retn = Electrode_->setInitialConditions(PSinput.anode_input_);
    if (retn == -1) {
        throw CanteraError("BDT_porAnode_LiIon::BDT_porAnode_LiIon()",
                           "Electrode::setInitialConditions() failed");
    }

    delete cfA;


    Capacity_initial_ = Electrode_->capacity();

    CapacityDischarged_initial_ = Electrode_->capacityDischarged();

    CapacityLeft_initial_ = Electrode_->capacityLeft();

    /*
     *  Create a vector of Equation Names
     *  This is the main place to specify the ordering of the equations within the code
     */
    EquationNameList.clear();
    VariableNameList.clear();

    // Continuity is used to solve for bulk velocity
    // Note that this is a single phase continuity so phase change will result in a source term
    //         Equation 0: = Continuity         variable 0 = Axial Velocity

    EquationNameList.push_back(EqnType(Continuity, 0, "Continuity: Bulk Velocity"));
    VariableNameList.push_back(VarType(Velocity_Axial, 0, 0));
    IsAlgebraic_NE[eqnIndex] = 1;
    IsArithmeticScaled_NE[eqnIndex] = 1;
    eqnIndex++;

    // List of species in the electrolyte
    const std::vector<std::string>& namesSp = ionicLiquid_->speciesNames();
    int nsp = ionicLiquid_->nSpecies();

    /*
     *  Loop over the species in the electrolyte phase. Each gets its own equation.
     *  Here, we hard code the mole fraction sum equation to the solvent, ECDMC, and we
     *  hardcode the charge conservation equation to PF6m. All other species are assigned
     *  the species conservation equation.
     */
    int iMFS = -1;
    int iCN = -1;
    for (int k = 0; k < nsp; k++) {
        if (namesSp[k] == "ECDMC") {
            iMFS = k;
            if (k != 0) {
              throw m1d_Error("  ", "ECDMC hardwired to species 0 not true: " + int2str(k));
            } 
            VariableNameList.push_back(VarType(MoleFraction_Species, k, (namesSp[k]).c_str()));
            EquationNameList.push_back(EqnType(MoleFraction_Summation, 0));
            IsAlgebraic_NE[1 + k] = 2;
        } else if (namesSp[k] == "PF6-") {
            iCN = k;
            if (k != 2) {
              throw m1d_Error("  ", "PF6- hardwired to species 2 not true: " + int2str(k));
            }
            VariableNameList.push_back(VarType(MoleFraction_Species, k, (namesSp[k]).c_str()));
            EquationNameList.push_back(EqnType(ChargeNeutrality_Summation, 0));
            IsAlgebraic_NE[1 + k] = 2;
        } else {
            VariableNameList.push_back(VarType(MoleFraction_Species, k, (namesSp[k]).c_str()));
            EquationNameList.push_back(EqnType(Species_Conservation, k, (namesSp[k]).c_str()));
            IsAlgebraic_NE[1 + k] = 0;
        }
        eqnIndex++;
    }
    if (iMFS < 0) {
        throw CanteraError("sep", "no ECDMC");
    }
    if (iCN < 0) {
        throw CanteraError("sep", "no PF6-");
    }



    //   Current conservation is used to solve for electrostatic potential
    //           Equation 4: Current Conservation - Electrolyte   Variable 4: Volts_Electrolyte
    EquationNameList.push_back(EqnType(Current_Conservation, 0, "Current Conservation"));
    VariableNameList.push_back(VarType(Voltage, 0, "Electrolyte"));
    IsAlgebraic_NE[eqnIndex] = 1;
    IsArithmeticScaled_NE[eqnIndex] = 1;
    eqnIndex++;

    // Current conservation is used to solve for electrostatic potential
    //  Equation 5: Current Conservation - Cathode   Variable 5: Volts_cathode

    EquationNameList.push_back(EqnType(Current_Conservation, 1, "Anode Current Conservation"));
    VariableNameList.push_back(VarType(Voltage, 1, "AnodeVoltage"));
    IsAlgebraic_NE[eqnIndex] = 1;
    IsArithmeticScaled_NE[eqnIndex] = 1;
    eqnIndex++;

    //  Temperature field is a given function across the domain
    //          
    if (PSinput.Energy_equation_prob_type_ == 2) {
	EqnType td = EqnType(Dirichlet_Specification, 0, "Temperature_Dirichlet");
	EquationNameList.push_back(td);
	VariableNameList.push_back(VarType(Temperature, 0, "Temperature"));
	IsAlgebraic_NE[eqnIndex] = 1;
	IsArithmeticScaled_NE[eqnIndex] = 0;
	eqnIndex++;
    }

    //  Enthalpy conservation is used to solve for the temperature
    //           EquationNameList.push_back(EqnType(Enthalpy_conservation, 0, "Enthalpy Conservation"));
    //           This is not an algebraic equation. It has a time derivative.
    if (PSinput.Energy_equation_prob_type_ == 3) {
	EqnType ht = EqnType(Enthalpy_Conservation, 0, "Enthalpy Conservation");
	EquationNameList.push_back(ht);
	VariableNameList.push_back(VarType(Temperature, 0, "Temperature"));
	IsAlgebraic_NE[eqnIndex] = 0;
	IsArithmeticScaled_NE[eqnIndex] = 0;
	eqnIndex++;
    }

    //  Thermal conservation is used to solve for the temperature
    //           
    if (PSinput.Energy_equation_prob_type_ == 4) {
	EqnType ht = EqnType(Thermal_Conservation, 0, "Thermal_Conservation");
	EquationNameList.push_back(ht);
	VariableNameList.push_back(VarType(Temperature, 0, "Temperature"));
	IsAlgebraic_NE[eqnIndex] = 0;
	IsArithmeticScaled_NE[eqnIndex] = 0;
	eqnIndex++;
    }

}
//=====================================================================================================================
BDT_porAnode_LiIon::BDT_porAnode_LiIon(const BDT_porAnode_LiIon& r) :
    BulkDomainDescription(r.DL_ptr_), ionicLiquid_(0), trans_(0), m_position(0), Electrode_(0)
{
    *this = r;
}
//=====================================================================================================================
BDT_porAnode_LiIon::~BDT_porAnode_LiIon()
{
    /*
     * Delete objects that we own
     */
    safeDelete(ionicLiquid_);
    safeDelete(trans_);
    safeDelete(Electrode_);
}
//=====================================================================================================================
BDT_porAnode_LiIon&
BDT_porAnode_LiIon::operator=(const BDT_porAnode_LiIon& r)
{
    if (this == &r) {
        return *this;
    }
    BulkDomainDescription::operator=(r);

    delete ionicLiquid_;
    ionicLiquid_ = (r.ionicLiquid_)->duplMyselfAsThermoPhase();

    delete trans_;
    trans_ = Cantera::newTransportMgr("Simple", ionicLiquid_, 1);

    m_position = r.m_position;
    delete Electrode_;
    // ok this is wrong and needs to be changed
    Electrode_ = r.Electrode_->duplMyselfAsElectrode();

    exit(-1);

    return *this;
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 *
 * @return  Returns a pointer to the object that will calculate the residual efficiently
 */
BulkDomain1D*
BDT_porAnode_LiIon::mallocDomain1D()
{
    BulkDomainPtr_ = new porousLiIon_Anode_dom1D(*this);
    return BulkDomainPtr_;
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
