/**
 * @file m1d_SDD_FlatAnode.cpp
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_SDD_FlatAnode.h"
#include "m1d_SurDomain_FlatLiSiAnode.h"
#include "m1d_defs.h"

#include "m1d_ProblemStatementCell.h"

#include "Electrode.h"
#include "Electrode_input.h"
#include "Electrode_Factory.h"
#include "zuzax/transport.h"      // transport properties
#include "zuzax/thermo.h"      // transport properties
#include "zuzax/thermo/IonsFromNeutralVPSSTP.h"  // ion properties
#include "zuzax/base/ctexceptions.h"

extern m1d::ProblemStatementCell PSinput;

#include  <string>

using namespace std;
using namespace Zuzax;

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
SDD_FlatAnode::SDD_FlatAnode(DomainLayout* dl_ptr, int pos) :
    SDD_Mixed(dl_ptr),
    m_position(pos), 
    ElectrodeA_(nullptr),
    ionicLiquidIFN_(nullptr)
{
/*
    ElectrodeA_ = new Zuzax::Electrode_SuccessiveSubstitution();

    ELECTRODE_KEY_INPUT* electrodeA_input = new ELECTRODE_KEY_INPUT();

    std::string commandFileA = "anode.inp";
    // Initialize a block input structure for the command file
    BEInput::BlockEntry* cfA = new BEInput::BlockEntry("command_file");

     // Go get the problem description from the input file
    electrodeA_input->printLvl_ = 5;
    int retn = electrodeA_input->electrode_input(commandFileA, cfA);

    if (retn == -1) {
        printf("exiting with error\n");
        exit(-1);
    }

    retn = ElectrodeA_->electrode_model_create(electrodeA_input);
    if (retn == -1) {
        printf("exiting with error\n");
        exit(-1);
    }

    retn = ElectrodeA_->setInitialConditions(electrodeA_input);
    if (retn == -1) {
        printf("exiting with error\n");
        exit(-1);
    }
*/
    size_t iph = (PSinput.PhaseList_)->globalPhaseIndex(PSinput.electrolytePhase_);
    if (iph == npos) {
        //throw ZuzaxError("SDD_FlatCathode::SDD_FlatCathode", "Can't find the phase in the phase list: " + PSinput.electrolytePhase_);
        throw Zuzax::ZuzaxError("SDD_FlatAnode::SDD_Anode", "Can't find the phase in the phase list: " );
    }
    thermo_t_double* ionicLiquid_ = & (PSinput.PhaseList_)->thermo(iph);
    ionicLiquidIFN_ = dynamic_cast<Zuzax::IonsFromNeutralVPSSTP *>( ionicLiquid_->duplMyselfAsThermoPhase() );
    ionicLiquid_ = (ThermoPhase*) ionicLiquidIFN_;

    ELECTRODE_KEY_INPUT *ai = PSinput.anode_input_;
    ElectrodeA_ = newElectrodeObject(ai->electrodeModelName);
    if (!ElectrodeA_) {
        throw  m1d_Error("SDD_FlatAnode::SDD_Anode", "Electrode factory method failed");
    }
    ELECTRODE_KEY_INPUT *ai_new = newElectrodeKeyInputObject(ai->electrodeModelName);
    std::string commandFile = ai->commandFile_;
    BEInput::BlockEntry *cfA = new BEInput::BlockEntry("command_file");


    //  Parse the complete child input file
    int retn = ai_new->electrode_input_child(commandFile, cfA);
    if (retn == -1) {
        throw  m1d_Error("SDD_FlatAnode::SDD_FlatAnode", "Electrode input child method failed");
    }
    /*
     * Switch the pointers around so that the child input file is returned.
     * Delete the original pointer.
     */
    delete ai;
    PSinput.anode_input_ = ai_new;

    retn = ElectrodeA_->electrode_model_create(PSinput.anode_input_);
    if (retn == -1) {
        throw  m1d_Error("SDD_FlatCathode::SDD_FlatAnode", "Electrode model create method failed");
    }
    retn = ElectrodeA_->setInitialConditions(PSinput.anode_input_);
    if (retn == -1) {
        throw  m1d_Error("SDD_FlatCathode::SDD_FlatAnode", "setInitialConditions method failed");
    }

    /*
     *  Add an equation for this surface domain
     *    For the anode we will install a boundary condition on it of a constant voltage
     *    This is the voltage datum for the system.
     */
    EquationNameList.push_back(EqnType(Current_Conservation, 1, "Anode Current Conservation"));
    VariableNameList.push_back(VarType(Voltage, 1, "AnodeVoltage"));

    safeDelete(cfA);
}
//=====================================================================================================================
SDD_FlatAnode::SDD_FlatAnode(const SDD_FlatAnode& r) :
    SDD_Mixed(r.DL_ptr_),
    m_position(0),
    ElectrodeA_(nullptr),
    ionicLiquidIFN_(nullptr)
{
    *this = r;
}
//=====================================================================================================================
SDD_FlatAnode::~SDD_FlatAnode()
{
    safeDelete(ElectrodeA_);
}
//=====================================================================================================================
SDD_FlatAnode& SDD_FlatAnode::operator=(const SDD_FlatAnode& r)
{
    if (this == &r) {
        return *this;
    }

    SDD_Mixed::operator=(r);
    m_position = r.m_position;


    delete ElectrodeA_;
    ElectrodeA_ =  (r.ElectrodeA_)->duplMyselfAsElectrode();

    delete  ionicLiquidIFN_;
    ionicLiquidIFN_ = new IonsFromNeutralVPSSTP(*(r.ionicLiquidIFN_));

    return *this;
}
//=====================================================================================================================
void
SDD_FlatAnode::SetEquationsVariablesList()
{
    EquationNameList.clear();
    VariableNameList.clear();
    /*
     *  Add an equation for this surface domain
     *    For the anode we will install a boundary condition on it of a constant voltage
     *    This is the voltage datum for the system.
     */
    EquationNameList.push_back(EqnType(Current_Conservation, 1, "Anode Current Conservation"));
    VariableNameList.push_back(VarType(Voltage, 1, "AnodeVoltage"));
}
//=====================================================================================================================
// Set the equation description
/*
 *  This routine is responsible for setting the variables:
 *    - NumEquationsPerNode
 *    - VariableNameList
 *    - EquationNameList
 *    - EquationIndexStart_EqName
 */
void
SDD_FlatAnode::SetEquationDescription()
{
    /*
     * Set the policy for connecting bulk domains
     * This really isn't set yet.
     */
    setRLMapping(1);
    /*
     * Fill in the rest of the information
     */
    SurfDomainDescription::SetEquationDescription();

    /*
     *  Add the dirichlet condition
     */
    EqnType e1 = EquationNameList[0];
    VarType v1 = VariableNameList[0];
    addDirichletCondition(e1, v1, 0.0);
}
//=====================================================================================================================
// Malloc and Return the object that will calculate the residual efficiently
/*
 * @return  Returns a pointer to the object that will calculate the residual
 *          efficiently
 */
SurDomain1D*
SDD_FlatAnode::mallocDomain1D()
{
    SurDomain_FlatLiSiAnode* s1d = new SurDomain_FlatLiSiAnode(*this, 1);
    return s1d;
}

//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================

