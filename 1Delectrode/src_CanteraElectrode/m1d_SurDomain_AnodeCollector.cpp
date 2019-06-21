/**
 * @file m1d_SurDomain_Cu2S.cpp
 *  object to calculate the  surface domains in the Cu2S problem
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_SurDomain_AnodeCollector.h"

#include "m1d_NodalVars.h"
#include "m1d_GlobalIndices.h"
#include "m1d_BulkDomain1D.h"
#include "m1d_DomainLayout.h" 
#include "m1d_Comm.h"
#include "m1d_SDD_AnodeCollector.h"
#include "m1d_BatteryResidEval.h"

using namespace std;

//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
static void drawline0(stream0& ss, int sp, int ll)
{
    for (int i = 0; i < sp; i++) {
        ss.print0(" ");
    }
    for (int i = 0; i < ll; i++) {
        ss.print0("-");
    }
    ss.print0("\n");
}
//==================================================================================================================================
SurDomain_AnodeCollector::SurDomain_AnodeCollector(SurfDomainDescription& sdd, int problemType) :
    SurBC_Dirichlet(sdd),
    bedd_(0),
    phiElectrolyte_(0.0),
    phiAnode_(0.0),
    phiAnodeCC_(0.0),
    icurrCollector_(0.0),
    CCThickness_(0.0),
    extraAnodeResistance_(0.0)
{
    //! Determine bedd_
    //       use SurfDomainDescription &SDD_;
    //                    LeftBulk or RightBulk
    bedd_ = SDD_.RightBulk;
    if (!bedd_) {
        bedd_ = SDD_.LeftBulk;
    }
    if (!bedd_) {
        throw m1d_Error("SurDomain_FlatLiSiAnode::SurDomain_FlatLiSiAnode", "Can't find adjoining bulk electrolyte domain");
    }
    const SDD_AnodeCollector* sdta_ptr = dynamic_cast<const SDD_AnodeCollector*>(&SDD_);
    CCThickness_ = sdta_ptr->anodeCCThickness_;
}
//==================================================================================================================================
SurDomain_AnodeCollector::SurDomain_AnodeCollector(const SurDomain_AnodeCollector& r) :
    SurBC_Dirichlet(r.SDD_), 
    bedd_(0),
    phiElectrolyte_(0.0),
    phiAnode_(0.0),
    phiAnodeCC_(0.0),
    icurrCollector_(0.0),
    CCThickness_(0.0),
    extraAnodeResistance_(0.0)
{
    operator=(r);
}
//==================================================================================================================================
SurDomain_AnodeCollector::~SurDomain_AnodeCollector()
{
}
//==================================================================================================================================
SurDomain_AnodeCollector&
SurDomain_AnodeCollector::operator=(const SurDomain_AnodeCollector& r)
{
    if (this == &r) {
        return *this;
    }
    SurBC_Dirichlet::operator=(r);

    SpecFlag_NE = r.SpecFlag_NE;
    Value_NE = r.Value_NE;
    bedd_ = r.bedd_;
    phiElectrolyte_ = r.phiElectrolyte_;
    phiAnode_ = r.phiAnode_;
    phiAnodeCC_ = r.phiAnodeCC_;
    icurrCollector_ = r.icurrCollector_;
    CCThickness_ =r.CCThickness_;
    extraAnodeResistance_ = r.extraAnodeResistance_;

    return *this;
}
//==================================================================================================================================
void
SurDomain_AnodeCollector::domain_prep(LocalNodeIndices* li_ptr)
{
    /*
     * First call the parent domain prep to get the node information
     *   -  Index_LcNode
     *   -  NumEqns
     *   - NodalVarPtr
     */
    SurBC_Dirichlet::domain_prep(li_ptr);

    for (size_t j = 0; j < NumNodeEqns; j++) {
        if (BC_Type_NE[j] == 10) {
            //SpecFlag_NE[j] = 0;
        }
    }
}
//==================================================================================================================================
void
SurDomain_AnodeCollector::residEval(Epetra_Vector& res,
                                    const bool doTimeDependentResid,
                                    const Epetra_Vector* const soln_ptr,
                                    const Epetra_Vector* const solnDot_ptr,
                                    const Epetra_Vector* const solnOld_ptr,
                                    const double t,
                                    const double rdelta_t,
                                    const Zuzax::ResidEval_Type residType,
                                    const Zuzax::Solve_Type solveType)
{
    residType_Curr_ = residType;
    int ieqn;
    const Epetra_Vector& soln = *soln_ptr;
    incrementCounters(residType);
    /*
     *  Quick return if we don't own the node that the boundary condition
     *  is to be applied on.
     */
    if (NumOwnedNodes == 0) {
        return;
    }
    /*
     *  Figure out the equation start for this node
     *   We start at the start of the equations for this node
     *   because we will be applying dirichlet conditions on the bulk
     *   equations.
     */
    int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
    //
    //  Store the entrance residual for later processing in balancing applications
    //
    if (residType == Zuzax::ResidEval_Type::Base_ResidEval || residType == Zuzax::ResidEval_Type::Base_ShowSolution) {
        size_t ieqnTemp = NodalVarPtr->Offset_VarType[Temperature];
        if (ieqnTemp != npos) {
            TempCollector_ = soln[index_EqnStart + ieqnTemp];
        } else {
            TempCollector_ = TemperatureReference_;
        }
    }

    /*
     * get the offsets for the BulkDomain and the surface domain.
     */
    int offsetBD = NodalVarPtr->OffsetIndex_BulkDomainEqnStart_BDN[0];

    int EQ_Current_offset_BD = offsetBD + bedd_->EquationIndexStart_EqName[Current_Conservation];
    int EQ_Current_offset_ED = EQ_Current_offset_BD + 1;
    /*
     *  Loop over the equations that the boundary conditions are going to be applied to
     *    -> This takes care of the current surface domain equation
     */
    double res_contrib = 0.0;
    for (size_t i = 0; i < NumNodeEqns; i++) {
        if (SpecFlag_NE[i]) {
	    res_contrib = 0.0;
            ieqn = index_EqnStart + i;
            double solnVal = soln[ieqn];
            double val = Value_NE[i];
            switch (BC_Type_NE[i]) {
            case 0:
                /*
                 *  For Dirichlet equations, replace the equation
                 */
                res_contrib = val - solnVal;
                res[ieqn] = res_contrib;
                break;
            case 1:
                /*
                 *  For flux boundary conditions, subtract from equation indicating a flux out
                 */
		res_contrib = -val;
                res[ieqn] += res_contrib;
                // CAL: WARNING m1d_SurDomain_CathodeCollector.cpp has += val
                break;
            case 2:
                /*
                 *  For Dirichlet boundary conditions with oscillation, replace the equation
                 */
		res_contrib =  val * TimeDep_NE[i](t) - solnVal;
                res[ieqn] = res_contrib;
                break;
            case 3:
                /*
                 *  For flux boundary conditions with oscillation, replace the equation
                 */
		res_contrib = val * TimeDep_NE[i](t);
                res[ieqn] += res_contrib;
                break;
            case 4: // voltage BCconstant
            case 6: // voltage BCsteptable
            case 8: // voltage BClineartable
                /*
                 *  For time dependent Dirichlet boundary condition using BoundaryCondition class
                 */
                res_contrib = BC_TimeDep_NE[i]->value(t) - solnVal;
                res[ieqn] = res_contrib;
                break;
            case 5: // current BCconstant
            case 7: // current BCsteptable
            case 9: // current BClineartable
                /*
                 *  For time dependent flux boundary condition using BoundaryCondition class
                 */
		res_contrib = BC_TimeDep_NE[i]->value(t);
                res[ieqn] += res_contrib;
                break;

            case 10: // Anode collector plate at constant voltage
                /*
                *     current_dens dot n = (V_cathode - V_cathCC) / Resistance_CC
                *
                *     A note about signs
                *     This is the residual equation for enthalpy. The temperature time derivative is always positive
                *     in all residuals. We only include the heat conduction term with an expression for flux in the equation
                *            flux = - lambda del T.
                *
                *      (1)    Normal temp residual = dH/dtDelV  + flux_Right - flux_left = 0.0
                *
                *     The boundary condition for the cathode side, which is on the +x right side is
                *
                *       (2)  flux_Right = - h (T_ref - T_N) = h ( T_N - T_ref ) = RobinBoundaryCondition
                *
                *     Therefore, we substitute (2) into (1) to get the signs correct, which is a + for cathode.
                *
                *     On the Anode side
                *       (3) - flux_Left = - h (T_ref - T_0) = h ( T_0 - T_ref ) = RobinBoundaryCondition
                *
                *     Therefore, we substitute (3) into (1) to get the signs correct, which is a + for anode
                *
                */
		res_contrib = BC_TimeDep_NE[i]->valueAtTime(t, solnVal);
                res[ieqn] += res_contrib;
                break;

            default:
                throw m1d_Error("SurDomain_AnodeCollector::residEval",
                                "BC_Type_NE[i] 0-9 for Dirichlet, Neumann, and Time Dependence");
            }
	    if (residType == Zuzax::ResidEval_Type::Base_ShowSolution) { 
		DomainResidVector_LastResid_NE[i] = res_contrib;
	    }
        }
    }
#ifdef DEBUG_HKM
    if (doTimeDependentResid) {
        if (residType == Zuzax::ResidEval_Type::Base_ResidEval) {
        }
    }
#endif

    getVoltages(&(soln[index_EqnStart]));
    /*
     * get the consistent currents
     */
    BulkDomain1D* bd = bedd_->BulkDomainPtr_;
    icurrCollector_ = bd->DiffFluxLeftBound_LastResid_NE[EQ_Current_offset_ED];
}
//==================================================================================================================================
void
SurDomain_AnodeCollector::getVoltages(const double* const solnElectrolyte)
{
    int indexVS = bedd_->VariableIndexStart_VarName[Voltage];
    phiElectrolyte_ = solnElectrolyte[indexVS];
    phiAnode_ = solnElectrolyte[indexVS + 1];
}
//==================================================================================================================================
void SurDomain_AnodeCollector::showSolution(const Epetra_Vector* soln_GlAll_ptr, const Epetra_Vector* solnDot_GlAll_ptr,
                                            const Epetra_Vector* soln_ptr, const Epetra_Vector* solnDot_ptr,
                                            const Epetra_Vector* solnOld_ptr, const Epetra_Vector_Owned* residInternal_ptr, 
					    const double t,
                                            const double rdelta_t, int indentSpaces, bool duplicateOnAllProcs)
{
    int locGbNode = SDD_.LocGbNode;
    // int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind = indent.c_str();
    // get the NodeVars object pertaining to this global node
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
    NodalVars* nv = gi->NodalVars_GbNode[locGbNode];
    int eqnStart = nv->EqnStart_GbEqnIndex;
    //std::vector<VarType> &variableNameList = SDD_.VariableNameList;
    std::vector<VarType>& variableNameListNode = nv->VariableNameList_EqnNum;
    int numVar = nv->NumEquations;
    std::string sss = id();
    BulkDomain1D* bd = bedd_->BulkDomainPtr_;
    int offsetBD = NodalVarPtr->OffsetIndex_BulkDomainEqnStart_BDN[0];
    int EQ_Current_offset_BD = offsetBD + bedd_->EquationIndexStart_EqName[Current_Conservation];
    int EQ_Current_offset_ED = EQ_Current_offset_BD + 1;
    icurrCollector_ = bd->DiffFluxLeftBound_LastResid_NE[EQ_Current_offset_ED];

    stream0 ss;
    print0_sync_start(0, ss, * (LI_ptr_->Comm_ptr_));
    if (doWrite) {
        drawline0(ss, indentSpaces, 80);
        ss.print0("%s  Solution on Surface Domain %10s : Number of variables = %d\n", ind, sss.c_str(), numVar);
        ss.print0("%s                                           : Number of boundary conditions = %d\n", ind, NumBCs);
        double x0 = nv->x0NodePos();
        ss.print0("%s                                           : Node %d at pos %g\n", ind, locGbNode, x0);
        drawline0(ss, indentSpaces, 80);
        ss.print0("%s     VariableName         Value        DirichletCondition\n", ind);
        drawline0(ss, indentSpaces + 2, 60);
        int jDir = 0;
        for (int k = 0; k < numVar; k++) {
            VarType& vt = variableNameListNode[k];
            std::string name = vt.VariableName(20);
            double sval = (*soln_GlAll_ptr)[eqnStart + k];
            ss.print0("%s   %-20s   %-10.4E ", ind, name.c_str(), sval);
	    if (SpecFlag_NE[k] != 0) {
		if (BC_Type_NE[k] == 0) {
		    ss.print0(" (Dir %d val = %-10.4E)", jDir, Value_NE[jDir]);
		} else {
		    ss.print0(" (BC %d Type %d val = %-10.4E)", jDir, BC_Type_NE[k], Value_NE[jDir]);
		}
                jDir++;
            }
            ss.print0("\n");
            if (vt.VariableType == Voltage && vt.VariableSubType == 1) {
                if (BC_Type_NE[k] == 10) {
                    ss.print0("%s   Volts(AnodeCurrentCollector) %-10.4E (thickness = %-10.4E)\n", ind, phiAnodeCC_, CCThickness_);
                }
            }
        }
	drawline0(ss, indentSpaces + 2, 60);
	ss.print0("%s   %-20s   %-10.4E \n", ind, "Current", icurrCollector_);

        drawline0(ss, indentSpaces + 2, 60);
        drawline0(ss, indentSpaces, 80);
    }
    print0_sync_end(0, ss, * (LI_ptr_->Comm_ptr_));
}
//==================================================================================================================================
void
SurDomain_AnodeCollector::eval_HeatBalance(const int ifunc,
                                             const double t,
                                             const double deltaT,
                                             const Epetra_Vector *soln_ptr,
                                             const Epetra_Vector *solnDot_ptr,
                                             const Epetra_Vector *solnOld_ptr,
                                             struct globalHeatBalVals& dVals)
 {

     globalHeatBalValsBat* dValsBat_ptr = dynamic_cast< globalHeatBalValsBat *>(& dVals);
     // 
     // We may add heat capacity later. However, not yet implemented.
     //
     dVals.totalHeatCapacity = 0.0;
     double resTempContrib = 0.0;
     const Epetra_Vector& soln = *soln_ptr;
   
     size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
     DomainLayout* dl = SDD_.DL_ptr_;
     ProblemResidEval* pb = dl->problemResid_;
     int timeRegion = pb->m_currentTimeRegion;
     //
     // get the offset of the temperature
     //
     size_t ieqnTemp = NodalVarPtr->Offset_VarType[Temperature];

     /*
      * get the offsets for the BulkDomain and the surface domain.
      */
     size_t offsetBD = NodalVarPtr->OffsetIndex_BulkDomainEqnStart_BDN[0];
     size_t EQ_Current_offset_BD = offsetBD + bedd_->EquationIndexStart_EqName[Current_Conservation];
     size_t EQ_Current_offset_ED = EQ_Current_offset_BD + 1;


     getVoltages(& (soln[index_EqnStart]));

     if (! SpecFlag_NE[ieqnTemp]) {
	 dVals.HeatFluxLeft = 0.0;
     } else {
   
   
	 size_t ieqn = index_EqnStart + ieqnTemp;
	 double solnValTemp = soln[ieqn];
	 double val = Value_NE[ieqnTemp];
	     
	 switch (BC_Type_NE[ieqnTemp]) {
	 case 0:
	     resTempContrib = val - solnValTemp;
	     break;	 
	 case 1:
	     resTempContrib= val;
	     dVals.HeatFluxLeft = resTempContrib;
	     break;
	 case 2:

	     resTempContrib = val * TimeDep_NE[ieqnTemp](t) - solnValTemp;
	     dVals.HeatFluxLeft = resTempContrib;
	     break;
	 case 3:
 
	     resTempContrib = val * TimeDep_NE[ieqnTemp](t);
	     dVals.HeatFluxLeft = resTempContrib;
	     break;

	 case 4: // voltage BCconstant
	 case 6: // voltage BCsteptable
	 case 8: // voltage BClineartable
	     /*
	      *  For time dependent Dirichlet boundary condition using BoundaryCondition class
	      *
	      *  Replace the equation for current at the boundary with a Dirichlet condition
	      * for the electrode
	      *  voltage involving a function
	      */
	     resTempContrib = BC_TimeDep_NE[ieqnTemp]->value(t, timeRegion) - solnValTemp;
	     break;

	 case 5: // current BCconstant
	 case 7: // current BCsteptable
	 case 9: // current BClineartable
	     /*
	      *  For time dependent flux boundary condition using BoundaryCondition class
	      *
	      *  For flux boundary conditions we must supply the value of  = icurr dot n_out
	      *  icurr is the current density (coul / (sec m2)) dotted into the outward facing normal
	      */
	     resTempContrib = BC_TimeDep_NE[ieqnTemp]->value(t, timeRegion);
	     dVals.HeatFluxLeft = resTempContrib;
	     break;

	 case 10:  // Collector constant current with collector plat resistance
	     /*
	      *     current_dens dot n = (V_cathode - V_cathCC) / Resistance_CC
	      *
	      *     A note about signs
	      *     This is the residual equation for enthalpy. The temperature time derivative is always positive
	      *     in all residuals. We only include the heat conduction term with an expression for flux in the equation
	      *            flux = - lambda del T.
	      *
	      *      (1)    Normal temp residual = dH/dtDelV  + flux_Right - flux_left = 0.0
	      *
	      *     The boundary condition for the cathode side, which is on the +x right side is
	      *
	      *       (2)  flux_Right = - h (T_ref - T_N) = h ( T_N - T_ref ) = RobinBoundaryCondition
	      *
	      *     Therefore, we substitute (2) into (1) to get the signs correct, which is a + for cathode.
	      *
	      *     On the Anode side
	      *       (3) - flux_Left = - h (T_ref - T_0) = h ( T_0 - T_ref ) = RobinBoundaryCondition
	      *
	      *     Therefore, we substitute (3) into (1) to get the signs correct, which is a + for anode
	      *
	      */
	     resTempContrib = BC_TimeDep_NE[ieqnTemp]->valueAtTime(t, solnValTemp, timeRegion);
	     dVals.HeatFluxLeft = resTempContrib;
	     break;

	 default:
	     throw m1d_Error("SurDomain_CathodeCollector::residEval",
			     "BC_Type_NE[i] 0-9 for Dirichlet, Neumann, and Time Dependence");
	 }
     }
     /*
     * Get the consistent currents
     */
     //SDT_CathodeCollector* SDD_cathode_ptr = dynamic_cast<SDT_CathodeCollector*>(&SDD_);
 
    //
    //  Retrieve the value of icurrElectrode_CBR_[iCell] from the cathode domain calculation
    //  We will use this as the official current coming out of the cathode.
    // 
    BulkDomain1D* bd = bedd_->BulkDomainPtr_;
    icurrCollector_ = bd->DiffFluxLeftBound_LastResid_NE[EQ_Current_offset_ED];
    double jfluxR = bd->DiffFluxLeftBound_LastResid_NE[ ieqnTemp ];
    dValsBat_ptr->currentLeft = icurrCollector_;
    dValsBat_ptr->JHelecLeft = jfluxR;
    dValsBat_ptr->phiSolid = phiAnode_;

 }
//==================================================================================================================================
void
SurDomain_AnodeCollector::initialConditions(const bool doTimeDependentResid, Epetra_Vector* const soln_ptr,
                                            Epetra_Vector* const solnDot, const double t, const double delta_t)
{
    size_t ieqn;
    Epetra_Vector& soln = *soln_ptr;
    /*
     *  Quick return if we don't own the node that the boundary condition
     *  is to be applied on.
     */
    if (NumOwnedNodes == 0) {
        return;
    }
    /*
     *  Figure out the equation start for this node
     *   We start at the start of the equations for this node
     *   because we will be applying Dirichlet conditions on the bulk
     *   equations.
     */
    size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
    /*
     *  Loop over the equations that the boundary conditions are going to be applied to
     */
    for (size_t i = 0; i < NumNodeEqns; i++) {
        if (SpecFlag_NE[i]) {
            /*
             *  For Dirichlet equations, replace the solution
             */
            ieqn = index_EqnStart + i;
            double val = Value_NE[i];
            if (BC_Type_NE[i] == 0) {
                soln[ieqn] = val;
            } else if (BC_Type_NE[i] == 2) {
                soln[ieqn] = val * TimeDep_NE[i](t);
            } else if (BC_Type_NE[i] == 4 || BC_Type_NE[i] == 6 || BC_Type_NE[i] == 8) {
                soln[ieqn] = BC_TimeDep_NE[i]->value(t);
            }
        }
    }
}
//==================================================================================================================================
void
SurDomain_AnodeCollector::saveDomain(Zuzax::XML_Node& oNode, const Epetra_Vector* const soln_GLALL_ptr,
                                     const Epetra_Vector* const solnDot_GLALL_ptr, const double t, bool duplicateOnAllProcs)
{
    // Find the global node number of the node where this domain resides
    size_t locGbNode = SDD_.LocGbNode;

    // get the NodeVars object pertaining to this global node
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
    NodalVars* nv = gi->NodalVars_GbNode[locGbNode];
    size_t eqnStart = nv->EqnStart_GbEqnIndex;
    //XML_Node& inlt = o.addChild("inlet");
    Zuzax::XML_Node& inlt = oNode.addChild("domain");
    size_t numVar = nv->NumEquations;
    inlt.addAttribute("id", id());
    inlt.addAttribute("points", 1);
    inlt.addAttribute("type", "surface");
    inlt.addAttribute("numVariables", numVar);
    double x0pos = nv->x0NodePos();
    double xpos = nv->xNodePos();
    double xfrac = nv->xFracNodePos();
    ztml::addFloat(inlt, "X0", x0pos, "", "", Zuzax::Undef, Zuzax::Undef);
    ztml::addFloat(inlt, "X", xpos, "", "", Zuzax::Undef, Zuzax::Undef);
    ztml::addFloat(inlt, "Xfraction", xfrac, "", "", Zuzax::Undef, Zuzax::Undef);

    for (size_t k = 0; k < numVar; k++) {
        double sval = (*soln_GLALL_ptr)[eqnStart + k];
        std::string nm = nv->VariableName(k);
        VarType vv = nv->VariableNameList_EqnNum[k];
        std::string type = VarType::VarMainName(vv.VariableType);
        ztml::addFloat(inlt, nm, sval, "", "", Zuzax::Undef, Zuzax::Undef);
    }
    std::string nm = "Volts(AnodeCCVoltage)";
    ztml::addFloat(inlt, nm, phiAnodeCC_, "", "", Zuzax::Undef, Zuzax::Undef);
}
//==================================================================================================================================
} 
//----------------------------------------------------------------------------------------------------------------------------------
