/**
 * @file m1d_SD_CathodeCollector.cpp
 *    Definitions for the Cathode Collector Surface Domain region
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_SDD_CathodeCollector.h"
#include "m1d_SurDomain_CathodeCollector.h"
#include "m1d_exception.h"

#include "m1d_ProblemStatementCell.h"
#include "m1d_BC_Battery.h"
#include "m1d_CanteraElectrodeGlobals.h"

#include "m1d_BDD_porousElectrode.h"


//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
SDD_CathodeCollector::SDD_CathodeCollector(DomainLayout* dl_ptr, int pos, const char* domainName) :
    SDD_Mixed(dl_ptr,domainName),
    m_position(pos),
    voltageVarBCType_(0),
    icurrCathodeSpecified_(0.0),
    voltageCathodeSpecified_(1.9),
    cathodeCCThickness_(0.0),
    extraResistanceCathode_(0.0),
    ResistanceLoad_(0.0),
    VoltageLoad_(0.0)
{

    voltageVarBCType_ = PSCinput_ptr->cathodeBCType_;
    icurrCathodeSpecified_ = PSCinput_ptr->icurrDischargeSpecified_;

    cathodeCCThickness_ = PSCinput_ptr->cathodeCCThickness_;
    extraResistanceCathode_ = PSCinput_ptr->extraCathodeResistance_;
    cathodeTempBCType_ = PSCinput_ptr->cathodeTempBCType_;
    ResistanceLoad_ = PSCinput_ptr->ResistanceLoad_;
    VoltageLoad_ = PSCinput_ptr->VoltageLoad_;

    /*
     *  Add an equation for this surface domain
     *    For the cathode we will install a boundary condition on it of either a constant voltage or a constant current.
     */
}
//==================================================================================================================================
SDD_CathodeCollector::SDD_CathodeCollector(const SDD_CathodeCollector& r) :
    SDD_Mixed(r.DL_ptr_),
    m_position(0),
    voltageVarBCType_(0),
    icurrCathodeSpecified_(0.0),
    voltageCathodeSpecified_(1.9),
    cathodeCCThickness_(0.0),
    extraResistanceCathode_(0.0) ,
    ResistanceLoad_(0.0) ,
    VoltageLoad_(0.0)
{
    *this = r;
}
//==================================================================================================================================
SDD_CathodeCollector::~SDD_CathodeCollector()
{
}
//==================================================================================================================================
SDD_CathodeCollector&
SDD_CathodeCollector::operator=(const SDD_CathodeCollector& r)
{
    if (this == &r) {
        return *this;
    }
    SDD_Mixed::operator=(r);

    m_position = r.m_position;
    voltageVarBCType_ = r.voltageVarBCType_;
    icurrCathodeSpecified_ = r.icurrCathodeSpecified_;
    voltageCathodeSpecified_ = r.voltageCathodeSpecified_;
    cathodeCCThickness_ = r.cathodeCCThickness_;
    extraResistanceCathode_ = r.extraResistanceCathode_;
    ResistanceLoad_ = r.ResistanceLoad_;
    VoltageLoad_ = r.VoltageLoad_;

    return *this;
}
//==================================================================================================================================
/*
 *  This routine is responsible for setting the variables:
 *    - NumEquationsPerNode
 *    - VariableNameList
 *    - EquationNameList
 *    - EquationIndexStart_EqName
 */
void
SDD_CathodeCollector::SetEquationDescription()
{
    /*
     * Set the policy for connecting bulk domains
     * This really isn't set yet.
     */
    setRLMapping(0);
    /*
     * Fill in the rest of the information
     */
    SurfDomainDescription::SetEquationDescription();

    // BDD_porousElectrode* bddPE = dynamic_cast< BDD_porousElectrode*>(LeftBulk);
    //ThermoPhase* ionicLiquid = bddPE->ionicLiquid_;
    //size_t nSpeciesElectrolyte = ionicLiquid->nSpecies();
    //size_t iMFS_index_ = bddPE->iMFS_index_;
    //size_t iCN_index_ = bddPE->iCN_index_;


    /*
     *  If we are just fixing the voltage at the cathode, we can set the plain Dirichlet condition here.
     *  If we are setting the current, we will add in a residual equation by hand in residEval().
     */
    voltageCathodeSpecified_ = PSCinput_ptr->CathodeVoltageSpecified_;
    double (*timeDepFunction)(double) = PSCinput_ptr->TimeDepFunction_;
    BoundaryCondition* BC_timeDep = PSCinput_ptr->BC_TimeDep_;
    EqnType e1(Current_Conservation, 2, "Cathode Current Conservation");
    VarType v1(Voltage, 2, "CathodeVoltage");
    double area = PSCinput_ptr->cathode_input_->electrodeGrossArea;
    if (voltageVarBCType_ == 10) {
        delete BC_timeDep;
        BC_timeDep = new BC_cathodeCC(cathodeCCThickness_, extraResistanceCathode_, area, voltageCathodeSpecified_);
    }
    if (voltageVarBCType_ == 11) {
        delete BC_timeDep;
        BC_timeDep = new BC_cathodeCCLoad(cathodeCCThickness_, extraResistanceCathode_, area,
                                          voltageCathodeSpecified_, ResistanceLoad_, VoltageLoad_);
    }

    switch (voltageVarBCType_) {
    case 0:
        addDirichletCondition(e1, v1, voltageCathodeSpecified_);
        break;
    case 1:
        addFluxCondition(e1, v1, icurrCathodeSpecified_);
        break;
    case 2:
        addDirichletCondition(e1, v1, voltageCathodeSpecified_, timeDepFunction);
        break;
    case 3:
        addFluxCondition(e1, v1, icurrCathodeSpecified_, timeDepFunction);
        break;
    case 4:
    case 6:
    case 8:
        addDirichletCondition(e1, v1, voltageVarBCType_, BC_timeDep);
        break;
    case 5:
    case 7:
    case 9:
        addFluxCondition(e1, v1, voltageVarBCType_, BC_timeDep);
        break;
    case 10:
        addRobinCondition(e1, v1, BC_timeDep, 10);
        break;
    case 11:
        addRobinCondition(e1, v1, BC_timeDep, 11);
        break;

    default:
        throw m1d_Error("SDD_CathodeCollector::SetEquationDescription",
                        "voltageVarBCType not 0, 1, 2 for Dirichlet, Neumann, or Diriclet with sin oscillation");
    }

    /*
     *  All of the other boundary conditions default to zero flux at the interface
     *      This includes:
     *
     *           flux of Li+
     *           flux of K+
     *           flux of Cl-
     *           flux of current
     *
     *  Because of the staggered grid, a dirichlet condition must be set on the axial velocity. Or else the
     *  last axial velocity unknown is not represented in the solution vector, leading to a singular matrix.
     *
     */
#ifdef DEBUG_OLD_CC_FLOW_CONDITION
    EqnType e2(Continuity, 0, "Continuity: Bulk Velocity");
    VarType v2(Velocity_Axial, 0, "Axial_Velocity");
    addDirichletCondition(e2, v2, 0.0);
#endif

    /*
     * Temperature boundary condition
     */
    EqnType et(Enthalpy_Conservation, 0, "");
    VarType vt(Temperature, 0, "");
    if (cathodeTempBCType_ != -1) {
        if (cathodeTempBCType_ == 0) {
            double tref =  PSCinput_ptr->cathodeTempRef_;
            addDirichletCondition(et, vt, tref);
        }
        if (cathodeTempBCType_ == 10) {
            double ht = PSCinput_ptr->cathodeHeatTranCoeff_;
            double tref =  PSCinput_ptr->cathodeTempRef_;
            double area = PSCinput_ptr->crossSectionalArea_;
            BoundaryCondition* BC_timeDep = new BC_heatTransfer(ht, tref, area);
            addRobinCondition(et, vt, BC_timeDep, cathodeTempBCType_);
        }
    }
}
//==================================================================================================================================
SurDomain1D* SDD_CathodeCollector::mallocDomain1D()
{
    SurDomain_CathodeCollector* s1d = new SurDomain_CathodeCollector(*this, 1);
    SurDomain1DPtr_ = s1d;
    return s1d;
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

