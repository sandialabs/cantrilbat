/**
 * @file m1d_SurDomain_Cu2S.cpp
 *  object to calculate the  surface domains in the Cu2S problem
 */


#include "m1d_SurDomain_CathodeCollector.h"

#include "m1d_NodalVars.h"
#include "m1d_SDD_Mixed.h"

#include "m1d_GlobalIndices.h"
#include "m1d_BulkDomainDescription.h"
#include "m1d_BulkDomain1D.h"

#include "m1d_SDD_CathodeCollector.h"

#include "m1d_DomainLayout.h" 
#include "m1d_Comm.h"
#include "m1d_materials.h"

#include "m1d_ProblemStatementCell.h"
#include "m1d_CanteraElectrodeGlobals.h"
#include "m1d_BatteryResidEval.h"

using namespace std;

//==================================================================================================================================
namespace m1d
{
//==================================================================================================================================
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
SurDomain_CathodeCollector::SurDomain_CathodeCollector(SurfDomainDescription& sdd, int problemType) :
    SurBC_Dirichlet(sdd),
    bedd_(0),
    phiElectrolyte_(0.0),
    phiCathode_(0.0),
    phiCathodeCC_(0.0),
    phiLoad_(0.0),
    icurrCollector_(0.0),
    CCThickness_(0.0),
    extraCathodeResistance_(0.0)
{
    //! Determine bedd_
    //       use SurfDomainDescription &SDD_;
    //                    LeftBulk or RightBulk
    bedd_ = SDD_.RightBulk;
    if (!bedd_) {
        bedd_ = SDD_.LeftBulk;
    }
    if (!bedd_) {
        throw m1d_Error("SurDomain_FlatLiSiCathode::SurDomain_FlatLiSiCathode",
                        "Can't find adjoining bulk electrolyte domain");
    }
}
//=================================================================================================================================
SurDomain_CathodeCollector::SurDomain_CathodeCollector(const SurDomain_CathodeCollector& r) :
    SurBC_Dirichlet(r.SDD_),
    bedd_(0),
    phiElectrolyte_(0.0),
    phiCathode_(0.0),
    phiCathodeCC_(0.0),
    phiLoad_(0.0),
    icurrCollector_(0.0),
    CCThickness_(0.0),
    extraCathodeResistance_(0.0)
{
    operator=(r);
}
//=====================================================================================================================
// Destructor
SurDomain_CathodeCollector::~SurDomain_CathodeCollector()
{
}
//=====================================================================================================================
// Assignment operator
/*
 * @param r      Object to be copied into the current object
 * @return       Returns a changeable reference to the current object
 */
SurDomain_CathodeCollector&
SurDomain_CathodeCollector::operator=(const SurDomain_CathodeCollector& r)
{
    if (this == &r) {
        return *this;
    }
    SurBC_Dirichlet::operator=(r);

    SpecFlag_NE                    = r.SpecFlag_NE;
    Value_NE                       = r.Value_NE;
    bedd_                          = r.bedd_;
    phiElectrolyte_                = r.phiElectrolyte_;
    phiCathode_                    = r.phiCathode_;
    phiCathodeCC_                  = r.phiCathodeCC_;
    phiLoad_                       = r.phiLoad_;
    icurrCollector_                = r.icurrCollector_;
    CCThickness_                   = r.CCThickness_;
    extraCathodeResistance_        = r.extraCathodeResistance_;

    return *this;
}
//===================================================================================================================================
// Prepare all of the indices for fast calculation of the residual
/*
 *  Here we collect all of the information necessary to
 *  speedily implement SpecFlag_NE and Value_NE within the
 *  residual calculation.
 *  We transfer the information from SDD_Dirichlet structure to
 * this structure for quick processing.
 */
void SurDomain_CathodeCollector::domain_prep(LocalNodeIndices* li_ptr)
{
    /*
     * First call the parent domain prep to get the node information
     *   -  Index_LcNode
     *   -  NumEqns
     *   - NodalVarPtr
     *   and figure out:
     *    - SpecFlag_NE[i]
     *    - Value_NE[i];
     */
    SurBC_Dirichlet::domain_prep(li_ptr);

    phiCathodeCC_   =  PSCinput_ptr->CathodeVoltageSpecified_;
    phiCathode_ =  phiCathodeCC_;
    phiLoad_ =  phiCathodeCC_;

    // TODO: In SurBC_Dirichlet::domain_prep() is is presumed that the
    // voltage/current BC is specified by a Dirichlet condition (specified voltage)
    // Here check to see if a current is specified (PSinput.cathodeBCType_
    // being an odd integer indicates current specification) and
    // if so set SpecFlag_NE(voltageEqn) = 0 so that we know it is a Neumann condition
}
//==================================================================================================================================
// Basic function to calculate the residual for the domain.
/*
 *  We calculate the additions and/or replacement of the
 *  residual here for the equations that this Dirichlet condition
 *  is responsible for.
 *
 * @param res           Output vector containing the residual
 * @param doTimeDependentResid  boolean indicating whether the time
 *                         dependent residual is requested
 * @param soln         Solution vector at which the residual should be
 *                     evaluated
 */
void SurDomain_CathodeCollector::residEval(Epetra_Vector& res, const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                                           const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr, const double t,
                                           const double rdelta_t, const Zuzax::ResidEval_Type residType,
                                           const Zuzax::Solve_Type solveType)
{
    residType_Curr_ = residType;
    size_t ieqn;
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
     * Find the current time region
     */
    DomainLayout* dl = SDD_.DL_ptr_;
    ProblemResidEval* pb = dl->problemResid_;
    int timeRegion = pb->m_currentTimeRegion;

    /*
     *  Figure out the equation start for this node
     *   We start at the start of the equations for this node because we will be applying dirichlet conditions on the bulk equations.
     */
    size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
    //
    //  Store the entrance residual for later processing in balancing applications
    //
    if (residType == Zuzax::ResidEval_Type::Base_ResidEval || residType == Zuzax::ResidEval_Type::Base_ShowSolution) {
	for (size_t i = 0; i < (size_t) NumNodeEqns; i++) {
	    Resid_BeforeSurDomain_NE[i] = res[index_EqnStart + i];
	}
	size_t ieqnTemp = NodalVarPtr->Offset_VarType[Temperature];
        if (ieqnTemp != npos) {
            TempCollector = soln[index_EqnStart + ieqnTemp];
        } else {
            TempCollector = TemperatureReference_; 
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
		res_contrib = val;
                res[ieqn] += res_contrib;
                break;

            case 2:
                /*
                 *  For Dirichlet boundary conditions with oscillation, replace the equation
                 */
		res_contrib  = val * TimeDep_NE[i](t) - solnVal;
                res[ieqn] = res_contrib;
                break;

            case 3:
                /*
                 *  For flux boundary conditions with oscillation, replace the equation
                 */
		res_contrib =  val * TimeDep_NE[i](t);
                res[ieqn] += res_contrib;
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
		res_contrib = BC_TimeDep_NE[i]->value(t, timeRegion) - solnVal;
                res[ieqn] = res_contrib;
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
		res_contrib += BC_TimeDep_NE[i]->value(t, timeRegion);
                res[ieqn] += res_contrib;
                break;

            case 10:
            case 11:  // Collector constant current with collector plat resistance
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
		res_contrib = BC_TimeDep_NE[i]->valueAtTime(t, solnVal, timeRegion);
                res[ieqn] += res_contrib;
                break;

            default:
                throw m1d_Error("SurDomain_CathodeCollector::residEval",
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
    //
    //  Retrieve the values of phiElectrolyte_ and phiCathode_ from the current solution vector
    //
    getVoltages(& (soln[index_EqnStart]));
    /*
     * Get the consistent currents
     */
    SDD_CathodeCollector* SDD_cathode_ptr = dynamic_cast<SDD_CathodeCollector*>(&SDD_);
    BulkDomain1D* bd = bedd_->BulkDomainPtr_;
    //
    //  Retrieve the value of icurrElectrode_CBR_[iCell] from the cathode domain calculation
    //  We will use this as the official current coming out of the cathode.
    //
    icurrCollector_ = bd->DiffFluxRightBound_LastResid_NE[EQ_Current_offset_ED];

    double resistivity = resistivity_aluminum(298.);
    double denom = resistivity * SDD_cathode_ptr->cathodeCCThickness_ +
                   SDD_cathode_ptr->extraResistanceCathode_ * crossSectionalArea_;
    double phiCathodeCCcalc = phiCathode_ - icurrCollector_ * denom;

    if ((SDD_cathode_ptr->voltageVarBCType_ != 10)) {
        phiCathodeCC_  = phiCathodeCCcalc;
    }

   

    if (energyEquationProbType_ == 3) {
	size_t offsetTemp = NodalVarPtr->indexBulkDomainVar0((size_t) Temperature);
        if (offsetTemp != npos) {
	    res_contrib =  icurrCollector_ * denom;
	    res[index_EqnStart + offsetTemp] -= res_contrib;
	    DomainResidVector_LastResid_NE[offsetTemp] -= res_contrib;
        }
    }
    if (residType_Curr_ == Zuzax::ResidEval_Type::Base_ShowSolution) {
	//if (SDD_cathode_ptr->voltageVarBCType_ == 11) {
	if (SDD_cathode_ptr->voltageVarBCType_ != 10) {
	    phiCathodeCC_ =  phiCathodeCCcalc;
	}
	phiLoad_ = phiCathodeCCcalc -  icurrCollector_ * (SDD_cathode_ptr->ResistanceLoad_ * crossSectionalArea_);
	//}
    }

#undef DAKOTAOUT
#ifdef  DAKOTAOUT
    double firstOutputTime = 0.061;
    static double outputTime = 0.06;
    if (t > outputTime) {
        std::ofstream dfp;
        if (outputTime > firstOutputTime) {
            dfp.open("results.out", std::ios_base::app);
        } else {
            dfp.open("results.out", std::ios_base::out);
        }
        std::cout << phiCathode_ << "  t" << outputTime << std::endl;
        dfp << phiCathode_ << "  t" << outputTime << std::endl;
        outputTime *= 10.0;
    }
#endif
#undef DAKOTAOUT
}
//==================================================================================================================================
void SurDomain_CathodeCollector::getVoltages(const double* const solnElectrolyte)
{
    int indexVS = bedd_->VariableIndexStart_VarName[Voltage];
    phiElectrolyte_ = solnElectrolyte[indexVS];
    phiCathode_ = solnElectrolyte[indexVS + 1];
}
//==================================================================================================================================
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
void SurDomain_CathodeCollector::showSolution(const Epetra_Vector* soln_GlAll_ptr,
                                              const Epetra_Vector* solnDot_GlAll_ptr,
                                              const Epetra_Vector* soln_ptr, const Epetra_Vector* solnDot_ptr,
                                              const Epetra_Vector* solnOld_ptr, const Epetra_Vector_Owned* residInternal_ptr, 
					      const double t,
                                              const double rdelta_t, int indentSpaces, bool duplicateOnAllProcs)
{
    int locGbNode = SDD_.LocGbNode;
    // int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
    string indent = "";
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
    string sss = id();
    //
    //  Retrieve the value of icurrElectrode_CBR_[iCell] from the cathode domain calculation
    //  We will use this as the official current coming out of the cathode.
    //
    BulkDomain1D* bd = bedd_->BulkDomainPtr_;
    int offsetBD = NodalVarPtr->OffsetIndex_BulkDomainEqnStart_BDN[0];
    int EQ_Current_offset_BD = offsetBD + bedd_->EquationIndexStart_EqName[Current_Conservation];
    int EQ_Current_offset_ED = EQ_Current_offset_BD + 1;
    icurrCollector_ = bd->DiffFluxRightBound_LastResid_NE[EQ_Current_offset_ED];

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
            string name = vt.VariableName(20);
            double sval = (*soln_GlAll_ptr)[eqnStart + k];
            ss.print0("%s   %-20s   % -10.4E ", ind, name.c_str(), sval);
            if (SpecFlag_NE[k] != 0) {
		if (BC_Type_NE[k] == 0) {
		    ss.print0(" (Dir %d val = %-10.4E)", jDir, Value_NE[jDir]);
		} else {
		    ss.print0(" (BC %d Type %d val = %-10.4E)", jDir, BC_Type_NE[k], Value_NE[jDir]);
		}
                jDir++;
            }
            ss.print0("\n");
	    SDD_CathodeCollector* SDD_cathode_ptr = dynamic_cast<SDD_CathodeCollector*>(&SDD_);
            if (vt.VariableType == Voltage && vt.VariableSubType == 2) {
                if (BC_Type_NE[k] == 10 || BC_Type_NE[k] == 11 || 
		    (CCThickness_ != 0.0 ) || (SDD_cathode_ptr->extraResistanceCathode_ != 0.0) || 
		    (SDD_cathode_ptr->ResistanceLoad_ != 0.0)) {
			ss.print0("%s   Volts(CathodeCC)        %-10.4E  (thickness = %-10.4E resistCC = % -10.4E ohm m2)\n",
				  ind, phiCathodeCC_, CCThickness_, SDD_cathode_ptr->extraResistanceCathode_ * crossSectionalArea_);
		}
		if ((BC_Type_NE[k] == 10) || (BC_Type_NE[k] == 11) || (SDD_cathode_ptr->ResistanceLoad_ != 0.0)) {
		    ss.print0("%s   Volts(Load)            % -10.4E  (resisLoad =  % -10.4E ohm m2, voltLoad =  % -10.4E volts)\n",
			      ind, phiLoad_, SDD_cathode_ptr->ResistanceLoad_ * crossSectionalArea_,
			      SDD_cathode_ptr->VoltageLoad_);
		}
            }
        }
	drawline0(ss, indentSpaces + 2, 60);
	ss.print0("%s   %-20s    %-10.4E \n", ind, "Current", icurrCollector_);


        drawline0(ss, indentSpaces + 2, 60);
        drawline0(ss, indentSpaces, 80);
    }
    print0_sync_end(0, ss, * (LI_ptr_->Comm_ptr_));
}
//==================================================================================================================================
void
SurDomain_CathodeCollector::eval_HeatBalance(const int ifunc,
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
	 dVals.HeatFluxRight = 0.0;
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
	     dVals.HeatFluxRight = resTempContrib;
	     break;
	 case 2:

	     resTempContrib = val * TimeDep_NE[ieqnTemp](t) - solnValTemp;
	     dVals.HeatFluxRight = resTempContrib;
	     break;
	 case 3:
 
	     resTempContrib = val * TimeDep_NE[ieqnTemp](t);
	     dVals.HeatFluxRight = resTempContrib;
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
	     dVals.HeatFluxRight = resTempContrib;
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
	     dVals.HeatFluxRight = resTempContrib;
	     break;

	 default:
	     throw m1d_Error("SurDomain_CathodeCollector::residEval",
			     "BC_Type_NE[i] 0-9 for Dirichlet, Neumann, and Time Dependence");
	 }
     }
     /*
      * Get the consistent currents
      */
     //SDD_CathodeCollector* SDD_cathode_ptr = dynamic_cast<SDD_CathodeCollector*>(&SDD_);
 
     //
     //  Retrieve the value of icurrElectrode_CBR_[iCell] from the cathode domain calculation
     //  We will use this as the official current coming out of the cathode.
     // 
     BulkDomain1D* bd = bedd_->BulkDomainPtr_;
     icurrCollector_ = bd->DiffFluxRightBound_LastResid_NE[EQ_Current_offset_ED];
     double jfluxR = bd->DiffFluxRightBound_LastResid_NE[ ieqnTemp ];
     dValsBat_ptr->currentRight = icurrCollector_;
     dValsBat_ptr->JHelecRight =   jfluxR;
     dValsBat_ptr->phiSolid = phiCathode_;

     SDD_CathodeCollector* SDD_cathode_ptr = dynamic_cast<SDD_CathodeCollector*>(&SDD_);
     double resistivity = resistivity_aluminum(298.);
     double denom = resistivity * SDD_cathode_ptr->cathodeCCThickness_ + SDD_cathode_ptr->extraResistanceCathode_ * crossSectionalArea_;

     if (energyEquationProbType_ == 3) {
	 size_t offsetTemp = NodalVarPtr->indexBulkDomainVar0((size_t) Temperature);
	 if (offsetTemp != npos) {
	     double res_contrib = icurrCollector_ * denom;
	     dValsBat_ptr->sourceTermExtra = res_contrib;
	 }
     }

 }
//==================================================================================================================================
// Generate the initial conditions
/*
 *   For surface Dirichlet conditions, we impose the t = 0- condition.
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
void SurDomain_CathodeCollector::initialConditions(const bool doTimeDependentResid, Epetra_Vector* const soln_ptr,
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
     * Find the current time region
     */
    DomainLayout* dl = SDD_.DL_ptr_;
    ProblemResidEval* pb = dl->problemResid_;
    int timeRegion = pb->m_currentTimeRegion;

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
                soln[ieqn] = BC_TimeDep_NE[i]->value(t, timeRegion);
            }
        }
    } //end for loop over equations

    // Store either the cathode voltage or current as specified by the BC
    //
    //  The cathode voltage variable index is 1 greater than the electrolyte voltage index
    size_t indexCathVoltage = bedd_->VariableIndexStart_VarName[Voltage] + 1;
    ieqn = index_EqnStart + indexCathVoltage;

    //  if we are setting the voltage, then use the value computed above for phiCathode_
    if (BC_Type_NE[indexCathVoltage] == 0 || BC_Type_NE[indexCathVoltage] == 2 || BC_Type_NE[indexCathVoltage] == 4
            || BC_Type_NE[indexCathVoltage] == 6 || BC_Type_NE[indexCathVoltage] == 8) {
        phiCathode_ = soln[ieqn];
    }
    phiCathodeCC_ = phiCathode_;

    //  For Neumann BC it is convenient to store the current fixed current specified by "Discharge Current" keyword
    //  within this object in icurrCollector_. It's also equal to  Value_NE[indexCathVoltage], which is fixed.
    //
    if (BC_Type_NE[indexCathVoltage] == 1) {
        icurrCollector_ = Value_NE[indexCathVoltage];
    }
    //  Time dependent function to specify the current.
    //  Store the initial value of icurrCollector_ within the variable icurrCollector_.
    //
    else if (BC_Type_NE[indexCathVoltage] == 3) {
        icurrCollector_ = Value_NE[indexCathVoltage] * TimeDep_NE[indexCathVoltage](t);
    }
    //
    //  Boundary condition class (constant, linear or step changes) : specification of current directly
    //  Calculate the initial current and store it in icurrCollector_
    //
    else if (BC_Type_NE[indexCathVoltage] == 5 || BC_Type_NE[indexCathVoltage] == 7 || BC_Type_NE[indexCathVoltage] == 9) {
        icurrCollector_ = BC_TimeDep_NE[indexCathVoltage]->value(t, timeRegion);
    }
}
//=================================================================================================================
void SurDomain_CathodeCollector::saveDomain(ZZCantera::XML_Node& oNode, const Epetra_Vector* soln_GLALL_ptr,
                                            const Epetra_Vector* solnDot_GLALL_ptr, const double t, bool duplicateOnAllProcs)
{
    // const double* s = soln_GLALL_ptr + loc();
    // Find the number of global equations on this domain, whether it's local or not
    //int numEquationsGb = SDD_.NumEquationsPerNode;
    // Find the global node number of the node where this domain resides
    int locGbNode = SDD_.LocGbNode;

    // get the NodeVars object pertaining to this global node
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
    NodalVars* nv = gi->NodalVars_GbNode[locGbNode];
    int eqnStart = nv->EqnStart_GbEqnIndex;
    //XML_Node& inlt = o.addChild("inlet");
    ZZCantera::XML_Node& inlt = oNode.addChild("domain");
    int numVar = nv->NumEquations;
    inlt.addAttribute("id", id());
    inlt.addAttribute("points", 1);
    inlt.addAttribute("type", "surface");
    inlt.addAttribute("numVariables", numVar);
    double x0pos = nv->x0NodePos();
    double xpos = nv->xNodePos();
    double xfrac = nv->xFracNodePos();
    ZZctml::addFloat(inlt, "X0", x0pos, "", "", ZZCantera::Undef, ZZCantera::Undef);
    ZZctml::addFloat(inlt, "X", xpos, "", "", ZZCantera::Undef, ZZCantera::Undef);
    ZZctml::addFloat(inlt, "Xfraction", xfrac, "", "", ZZCantera::Undef, ZZCantera::Undef);

    for (int k = 0; k < numVar; k++) {
        double sval = (*soln_GLALL_ptr)[eqnStart + k];
        string nm = nv->VariableName(k);
        VarType vv = nv->VariableNameList_EqnNum[k];
        string type = VarType::VarMainName(vv.VariableType);
        ZZctml::addFloat(inlt, nm, sval, "", "", ZZCantera::Undef, ZZCantera::Undef);
    }
    string nm = "Volts(CurrentCCVoltage)";
    ZZctml::addFloat(inlt, nm, phiCathodeCC_, "", "", ZZCantera::Undef, ZZCantera::Undef);

}
//===================================================================================================================================
} /* End of Namespace */
//===================================================================================================================================
