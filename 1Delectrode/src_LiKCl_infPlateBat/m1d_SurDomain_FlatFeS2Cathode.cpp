/**
 * @file m1d_SurDomain_Cu2S.cpp
 *  object to calculate the  surface domains in the Cu2S problem
 */


#include "m1d_SurDomain_FlatFeS2Cathode.h"

#include "m1d_NodalVars.h"
#include "m1d_SDD_Mixed.h"

#include "m1d_exception.h"
#include "m1d_GlobalIndices.h"
#include "m1d_BulkDomainDescription.h"

#include "m1d_SDD_FlatCathode.h"

#include "Epetra_Comm.h"
#include "Epetra_Vector.h"

#include "m1d_Comm.h"

using namespace std;

//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
static void
drawline(int sp, int ll)
{
  for (int i = 0; i < sp; i++) {
    Zuzax::writelog(" ");
  }
  for (int i = 0; i < ll; i++) {
    Zuzax::writelog("-");
  }
  Zuzax::writelog("\n");
}
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================

SurDomain_FlatFeS2Cathode::SurDomain_FlatFeS2Cathode(SurfDomainDescription &sdd, int problemType) :
  ReactingSurDomain(sdd), ElectrodeC_(0), electrolyteThermo_(0), bedd_(0), mfElectrolyte_Soln(0),
      mfElectrolyte_Thermo(0), phiElectrolyte_(0.0), phiCathode_(0.0), NumNodeEqns(0), SpecFlag_NE(0), Value_NE(0),
      electrodeSpeciesProdRates_(0), phaseMoleFlux_(0), surfaceArea_(1.0), concTot_Curr_(0.0), voltageVarBCType_(0),
  icurrCathodeSpecified_(0.0),
  icurrCollector_(0.0)
{

  SDD_FlatCathode *fc = dynamic_cast<SDD_FlatCathode *> (&sdd);
  if (!fc) {
    throw m1d_Error("confused", "confused");
  }
  ElectrodeC_ = fc->ElectrodeC_;
  voltageVarBCType_ = fc->voltageVarBCType_;
  icurrCathodeSpecified_ = fc->icurrCathodeSpecified_;

  int iph = ElectrodeC_->solnPhaseIndex();
  electrolyteThermo_ = &(ElectrodeC_->thermo(iph));

  //! Determine bedd_
  //       use SurfDomainDescription &SDD_;
  //                    LeftBulk or RightBulk
  bedd_ = SDD_.LeftBulk;
  if (!bedd_) {
    bedd_ = SDD_.RightBulk;
  }
  if (!bedd_) {
    throw m1d_Error("SurDomain_FlatLiSiAnode::SurDomain_FlatLiSiAnode", "Can't find adjoining bulk electrolyte domain");
  }
}
//=====================================================================================================================
SurDomain_FlatFeS2Cathode::SurDomain_FlatFeS2Cathode(const SurDomain_FlatFeS2Cathode &r) :
  ReactingSurDomain(r.SDD_), ElectrodeC_(0), electrolyteThermo_(0), bedd_(0), mfElectrolyte_Soln(0),
      mfElectrolyte_Thermo(0), phiElectrolyte_(0.0), phiCathode_(0.0), NumNodeEqns(0), SpecFlag_NE(0), Value_NE(0),
      electrodeSpeciesProdRates_(0), phaseMoleFlux_(0), surfaceArea_(1.0), concTot_Curr_(0.0), voltageVarBCType_(0),
  icurrCathodeSpecified_(0.0),
  icurrCollector_(0.0)
{
  operator=(r);
}
//=====================================================================================================================
// Destructor
SurDomain_FlatFeS2Cathode::~SurDomain_FlatFeS2Cathode()
{
}
//=====================================================================================================================
// Assignment operator
/*
 * @param r      Object to be copied into the current object
 * @return       Returns a changeable reference to the current object
 */
SurDomain_FlatFeS2Cathode &
SurDomain_FlatFeS2Cathode::operator=(const SurDomain_FlatFeS2Cathode &r)
{
  if (this == &r) {
    return *this;
  }
  ReactingSurDomain::operator=(r);

  ElectrodeC_ = r.ElectrodeC_;
  electrolyteThermo_ = r.electrolyteThermo_;
  bedd_ = r.bedd_;
  mfElectrolyte_Soln = r.mfElectrolyte_Soln;
  mfElectrolyte_Thermo = r.mfElectrolyte_Thermo;
  phiElectrolyte_ = r.phiElectrolyte_;
  phiCathode_ = r.phiCathode_;
  NumNodeEqns = r.NumNodeEqns;
  SpecFlag_NE = r.SpecFlag_NE;
  Value_NE = r.Value_NE;
  electrodeSpeciesProdRates_ = r.electrodeSpeciesProdRates_;
  phaseMoleFlux_ = r.phaseMoleFlux_;
  surfaceArea_ = r.surfaceArea_;
  concTot_Curr_ = r.concTot_Curr_;
  voltageVarBCType_ = r.voltageVarBCType_;
  icurrCathodeSpecified_ = r.icurrCathodeSpecified_;
  icurrCollector_ = r.icurrCollector_;

  return *this;
}
//==================================================================================================================================
// Prepare all of the indices for fast calculation of the residual
/*
 *  Here we collect all of the information necessary to
 *  speedily implement SpecFlag_NE and Value_NE within the
 *  residual calculation.
 *  We transfer the information from SDT_Dirichlet structure to
 * this structure for quick processing.
 */
void
SurDomain_FlatFeS2Cathode::domain_prep(LocalNodeIndices *li_ptr)
{
  /*
   * First call the parent domain prep to get the node information
   *   -  Index_LcNode
   *   -  NumEqns
   *   - NodalVarPtr
   */
  SurDomain1D::domain_prep(li_ptr);
  /*
   *  Now figure out:
   *    - SpecFlag_NE[i]
   *    - Value_NE[i];
   */

  /*
   *  Make sure that the SurfDomainType associated with this domain
   *  is a straight surface Dirichlet condition
   */
  SDD_Mixed *sdt = dynamic_cast<SDD_Mixed *> (&SDD_);
  AssertThrow(sdt, "bad cast");
  /*
   * Zero out the current domain's setup section.
   */
  NumNodeEqns = NodalVarPtr->NumEquations;

  SpecFlag_NE.resize(NumNodeEqns);
  Value_NE.resize(NumNodeEqns);
  for (int j = 0; j < NumNodeEqns; j++) {
    SpecFlag_NE[j] = 0;
    Value_NE[j] = 0.0;
  }
  NumBCs = 0;

  /*
   * Loop over the Dirichlet conditions defined in the SurfDomainType
   * Structure. Transfer the information to this structure for quick
   * processing.
   */
  for (int i = 0; i < sdt->NumConditions; i++) {
    VarType vtDir = sdt->VariableID[i];
    VAR_TYPE vtmDir = vtDir.VariableType;
    VAR_TYPE_SUBNUM stDir = vtDir.VariableSubType;
    // If the subtype is -1, then we apply to all equations of the
    // variable main type
    /*
     *
     */
    if (stDir == -1) {
      for (int j = 0; j < NumNodeEqns; j++) {
        VarType eqt = NodalVarPtr->VariableNameList_EqnNum[j];
        if ((vtmDir == Variable_Type_Any) || (eqt.VariableType == vtmDir)) {
          NumBCs++;
          SpecFlag_NE[j] = 1;
          Value_NE[j] = sdt->Value[i];
        } else {
          SpecFlag_NE[j] = 0;
        }
      }
    } else {
      for (int j = 0; j < NumNodeEqns; j++) {
        VarType vtNode = NodalVarPtr->VariableNameList_EqnNum[j];
        if (vtDir == vtNode) {
          NumBCs++;
          SpecFlag_NE[j] = 1;
          Value_NE[j] = sdt->Value[i];
        } else {
          SpecFlag_NE[j] = 0;
        }
      }
    }
  }

  mfElectrolyte_Soln.resize(3, 0.0);
  mfElectrolyte_Thermo.resize(3, 0.0);
  electrodeSpeciesProdRates_.resize(30, 0.0);
  phaseMoleFlux_.resize(30, 0.0);
}
//=====================================================================================================================
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
void
SurDomain_FlatFeS2Cathode::residEval(Epetra_Vector &res,
                                     const bool doTimeDependentResid,
                                     const Epetra_Vector *soln_ptr,
                                     const Epetra_Vector *solnDot_ptr,
                                     const Epetra_Vector *solnOld_ptr,
                                     const double t,
                                     const double rdelta_t,
                                     const Zuzax::ResidEval_Type residType,
				     const Zuzax::Solve_Type solveType)
{
  residType_Curr_ = residType;
  int ieqn;
  const Epetra_Vector &soln = *soln_ptr;
  incrementCounters(residType);
  /*
   *  Quick return if we don't own the node that the boundary condition
   *  is to be applied on.
   */
  if (NumOwnedNodes == 0) {
    return;
  }

  /*
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr_->NodalVars_LcNode[Index_LcNode];

  /*
   *  Figure out the equation start for this node
   *   We start at the start of the equations for this node
   *   because we will be applying dirichlet conditions on the bulk
   *   equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  /*
   *  Figure out the surface domain on which these boundary conditions
   *  will be applied. There is one controlling one
   */
  // NodalVars *nv = NodalVarPtr;
  /*
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr->NodalVars_LcNode[Index_LcNode];

  /*
   *
   * Update Zuzax with the solution vector
   */
  updateDependencies(soln_ptr, t, residType);

  // Zuzax::ReactingSurDomain *rSurDomain = ElectrodeC_->m_rSurDomain;

  /*
   * get the offsets for the BulkDomain and the surface domain.
   */
  int offsetBD = NodalVarPtr->OffsetIndex_BulkDomainEqnStart_BDN[0];
  int offsetSD = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[0];

  int EQ_Current_offset_BD = offsetBD + bedd_->EquationIndexStart_EqName[Current_Conservation];
  //int EQ_Current_offset_SD = offsetSD + SDD_.EquationIndexStart_EqName[Current_Conservation];

  int EQ_TCont_offset_BD = offsetBD + bedd_->EquationIndexStart_EqName[Continuity];
  int EQ_Species_offset_BD = offsetBD + bedd_->EquationIndexStart_EqName[Species_Conservation];

  int indexCent_EqnStart_BD = index_EqnStart + offsetBD;
  int iVAR_Vaxial_BD = bedd_->VariableIndexStart_VarName[Velocity_Axial];
  /*  Get the equation offsets at the node
   */

  /*
   *   Calculate the state of the system on that controlling surface domain
   */

  /*
   * Calculate the rates of production of all species in the Electrode
   */
  double icurr = ElectrodeC_->getNetSurfaceProductionRatesCurrent(0, &electrodeSpeciesProdRates_[0]);

  /*
   * Get the phase mole flux
   */
  ElectrodeC_->getPhaseMoleFlux(0, &phaseMoleFlux_[0]);
  int sf = ElectrodeC_->solnPhaseIndex();
  double solnMoleFlux = phaseMoleFlux_[sf];
  int speciesIndex0 = ElectrodeC_->globalSpeciesIndex(sf, 0);

  /*
   *  Loop over the equations that the boundary conditions are going to be applied to
   *    -> This takes care of the current surface domain equation
   */
  for (int i = 0; i < NumNodeEqns; i++) {
    if (SpecFlag_NE[i]) {
      /*
       *  For Dirichlet equations, replace the equation
       */
      ieqn = index_EqnStart + i;
      double solnVal = soln[ieqn];
      double val = Value_NE[i];
      res[ieqn] = val - solnVal;
    }
  }
  /*
   * Boundary conditions on top
   */
  double area_cvb = 1.0;
  if (coordinateSystemType_ == Cylindrical_Coordinates) {
    //  area_cvb = Pi * m_cellBound[m_nodes - 1] * m_xbot0;
    //  whole program needs to be fixed for Cylindrical and spherical coordinates
    //area_cvb = CellArea_;
    exit(-1);
  } else if (coordinateSystemType_ == Cylindrical_Coordinates) {
    //area_cvb = 2 * Pi * m_cellBound[m_nodes - 1] * m_cellBound[m_nodes - 1];
    //area_cvb = CellArea_;
    exit(-1);
  }
  surfaceArea_ = 0.1;

  // ------------------ BULK DOMAIN -----------------------------------
  /*
   *  Total Continuity Equation
   *     Mole flux of soln into the soln phase creates a source term
   *     Note, this term is left out of the solution on purpose.
   */
  double Fright_cc_ = soln[indexCent_EqnStart_BD + iVAR_Vaxial_BD];

  double porosity_Curr_ = 0.5;

  res[index_EqnStart + EQ_TCont_offset_BD] = Fright_cc_ * concTot_Curr_ * porosity_Curr_ + solnMoleFlux * surfaceArea_
      * area_cvb;

  /*
   * Species 0 Conservation equation
   */
  res[index_EqnStart + EQ_Species_offset_BD] -= electrodeSpeciesProdRates_[speciesIndex0] * surfaceArea_ * area_cvb;
  /*
   * Sum MF equation - no surface contribution
   */
  /*
   *  Electroneutrality equation - no surface contribution
   */
  /*
   *  Current equation within the electrolyte
   */
  res[index_EqnStart + EQ_Current_offset_BD] -= icurr * surfaceArea_ * area_cvb;

  // ------------------ SURFACE DOMAIN -----------------------------------
  /*
   *   Surface domain's voltage is handled two ways. The first way is a Dirichlet condition. This
   *   is not handled here.
   *   The second way is handled here. The current is set to a specified value. This creates an
   *   equation below.
   */
  if (voltageVarBCType_ == 1) {
    int EQ_Current_offset_SD = offsetSD + SDD_.EquationIndexStart_EqName[Current_Specification];

    res[index_EqnStart + EQ_Current_offset_SD] = icurr * surfaceArea_ - icurrCathodeSpecified_;
  }
  // }
  /*
   *  For the steady state equations case, we are counting on the fact that the mesh movement residual
   *  within the volume domain handles the case.
   */

  if (doTimeDependentResid) {

#ifdef DEBUG_HKM
    if (residType == Zuzax::ResidEval_Type::Base_ResidEval) {
      double tmp = -electrodeSpeciesProdRates_[speciesIndex0] * surfaceArea_ * area_cvb;

      //printf(" Cell = %d, CathdSource = %10.3e, ResTotal = %10.3e\n", 9, tmp,
      //       res[index_EqnStart + EQ_Species_offset_BD]);
      tmp += 1.0;
    }
#endif
  }

  //
  // Store the current for use in printouts
  //
  icurrCollector_ = icurr * surfaceArea_;
}

//=====================================================================================================================
// utility routine to update the objects used to calculate quantities at the surface
/*
 *
 * @param soln_ptr     solution vector at which the residual should be
 *                     evaluated
 * @param t            time
 * @param residType    Residual evaluation type
 */
void
SurDomain_FlatFeS2Cathode::updateDependencies(const Epetra_Vector *soln_ptr,
                                              const double t,
                                              const Zuzax::ResidEval_Type residType)
{
  /*
   *  Figure out the equation start for this node
   *  We start at the start of the equations for this node
   *  because we will be applying dirichlet conditions on the bulk
   *  equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  // Find the offset for the bulk domain in the solution vector.
  int offsetElectrolyteBulk = NodalVarPtr->OffsetIndex_BulkDomainEqnStart_BDN[0];
  AssertTrace(offsetElectrolyteBulk == 0);

  const double *solnElectrolyte = &((*soln_ptr)[index_EqnStart + offsetElectrolyteBulk]);

  int offsetSolid = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[0];
  const double *solnSolid = &((*soln_ptr)[index_EqnStart + offsetSolid]);

  /*
   * Get the temperature: Check to see if the temperature is in the solution vector.
   *   If it is not, then use the reference temperature
   */
  double Temp = TemperatureReference_;
  int iTemp = SDD_.VariableIndexStart_VarName[Temperature];
  if (iTemp >= 0) {
    Temp = solnElectrolyte[iTemp];
  }

  /*
   * Get the pressure
   */
  double pres = PressureReference_;

  /*
   * Get the electrolyte mole fractions and store them in local storage
   */
  getMFElectrolyte_soln(solnElectrolyte);

  /*
   *  Get the voltages for the metal and solution and store them in
   *  local storage
   */
  getVoltages(solnElectrolyte, solnSolid);

  /*
   * Get the species vector for all bulk and species phases
   */
  //NodalVars *nv = NodalVarPtr;

  /*
   * set the properties in the Electrode object
   *  -> temperature and pressure
   *  -> voltages of the phases
   */
  ElectrodeC_->setState_TP(Temp, pres);
  ElectrodeC_->setVoltages(phiCathode_, phiElectrolyte_);
  ElectrodeC_->setElectrolyteMoleNumbers(&(mfElectrolyte_Thermo[0]), false);
  /*
   * Set anything else here for the electrode object. However, for this initial
   * implementation there is nothing more
   */

  /*
   * Set the internal objects within the electrode
   */
  //ElectrodeC_->downloadMP();
  ElectrodeC_->updateState();
  // Calculate the total concentration of the electrolyte kmol m-3.

  concTot_Curr_ = electrolyteThermo_->molarDensity();
}
//=====================================================================================================================
void
SurDomain_FlatFeS2Cathode::getVoltages(const double * const solnElectrolyte, const double * const solnSolid)
{
  int indexVS = bedd_->VariableIndexStart_VarName[Voltage];
  phiElectrolyte_ = solnElectrolyte[indexVS];

  int indexVE = SDD_.VariableIndexStart_VarName[Voltage];
  phiCathode_ = solnSolid[indexVE];
}
//=====================================================================================================================
void
SurDomain_FlatFeS2Cathode::getMFElectrolyte_soln(const double * const solnBulk)
{
  int indexMF = bedd_->VariableIndexStart_VarName[MoleFraction_Species];
  mfElectrolyte_Soln[0] = solnBulk[indexMF];
  mfElectrolyte_Soln[1] = solnBulk[indexMF + 1];
  mfElectrolyte_Soln[2] = solnBulk[indexMF + 2];
  double mf0 = std::max(mfElectrolyte_Soln[0], 0.0);
  double mf1 = std::max(mfElectrolyte_Soln[1], 0.0);
  double tmp = mf0 + mf1;

  mfElectrolyte_Thermo[0] = (mf0) * 0.5 / tmp;
  mfElectrolyte_Thermo[1] = (mf1) * 0.5 / tmp;
  mfElectrolyte_Thermo[2] = 0.5;
}
//=====================================================================================================================
// Base class for saving the solution on the domain in an xml node.
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
SurDomain_FlatFeS2Cathode::saveDomain(Zuzax::XML_Node& oNode,
                                      const Epetra_Vector *soln_GLALL_ptr,
                                      const Epetra_Vector *solnDot_GLALL_ptr,
                                      const double t,
                                      bool duplicateOnAllProcs)
{
  // const double* s = soln_GLALL_ptr + loc();
  // Find the number of global equations on this domain, whether it's local or not
  //int numEquationsGb = SDD_.NumEquationsPerNode;
  // Find the global node number of the node where this domain resides
  int locGbNode = SDD_.LocGbNode;

  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  NodalVars *nv = gi->NodalVars_GbNode[locGbNode];
  int eqnStart = nv->EqnStart_GbEqnIndex;
  //XML_Node& inlt = o.addChild("inlet");
  Zuzax::XML_Node& inlt = oNode.addChild("domain");
  int numVar = nv->NumEquations;
  inlt.addAttribute("id", id());
  inlt.addAttribute("points", 1);
  inlt.addAttribute("type", "surface");
  inlt.addAttribute("numVariables", numVar);
  double x0pos = nv->x0NodePos();
  double xpos = nv->xNodePos();
  ztml::addFloat(inlt, "X0", x0pos, "", "", Zuzax::Undef, Zuzax::Undef);
  ztml::addFloat(inlt, "X", xpos, "", "", Zuzax::Undef, Zuzax::Undef);

  for (int k = 0; k < numVar; k++) {
    double sval = (*soln_GLALL_ptr)[eqnStart + k];
    string nm = nv->VariableName(k);
    VarType vv = nv->VariableNameList_EqnNum[k];
    string type = VarType::VarMainName(vv.VariableType);
    ztml::addFloat(inlt, nm, sval, "", "", Zuzax::Undef, Zuzax::Undef);
  }
}
//=====================================================================================================================
// Class for writing the solution on the domain to a logfile.
/*
 *
 * @param soln_GlALL_ptr       Pointer to the Global-All solution vector
 * @param solnDot_GlALL_ptr    Pointer to the Global-All solution dot vector
 * @param soln_ptr             Pointer to the solution vector
 * @param solnDot_ptr          Pointer to the time-derivative of the solution vector
 * @param solndOld_ptr         Pointer to the solution vector at the old time step
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
SurDomain_FlatFeS2Cathode::showSolution(const Epetra_Vector *soln_GlAll_ptr,
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
  char buf[132];
  int locGbNode = SDD_.LocGbNode;
  //int mypid = LI_ptr_->Comm_ptr_->MyPID();
  // Changed to whether I own the node
  bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
  string indent = "";
  for (int i = 0; i < indentSpaces; i++) {
    indent += " ";
  }
  stream0 ss;
  print0_sync_start(0, ss, *(LI_ptr_->Comm_ptr_));
  if (doWrite) {
    const char *ind = indent.c_str();
    // get the NodeVars object pertaining to this global node
    GlobalIndices *gi = LI_ptr_->GI_ptr_;
    NodalVars *nv = gi->NodalVars_GbNode[locGbNode];
    int eqnStart = nv->EqnStart_GbEqnIndex;
    std::vector<VarType> &variableNameListNode = nv->VariableNameList_EqnNum;
    int numVar = nv->NumEquations;
    int sf = ElectrodeC_->solnPhaseIndex();
    int speciesIndex0 = ElectrodeC_->globalSpeciesIndex(sf, 0);

    updateDependencies(soln_ptr, t);

    /*
     * Calculate the rates of production of all species in the Electrode
     */
    double icurr = ElectrodeC_->getNetSurfaceProductionRatesCurrent(0, &electrodeSpeciesProdRates_[0]);

    /*
     * Get the phase mass flux
     */
    ElectrodeC_->getPhaseMoleFlux(0, &phaseMoleFlux_[0]);

    string sss = id();

    drawline(indentSpaces, 80);
    sprintf(buf, "%s  Solution on Surface Domain %10s : Number of variables = %d\n", ind, sss.c_str(), numVar);
    Zuzax::writelog(buf);
    sprintf(buf, "%s                                           : Number of boundary conditions = %d\n", ind, NumBCs);
    Zuzax::writelog(buf);
    double x0 = nv->x0NodePos();
    sprintf(buf, "%s                                           : Node %d at pos %g\n", ind, locGbNode, x0);
    Zuzax::writelog(buf);
    drawline(indentSpaces, 80);
    sprintf(buf, "%s     VariableName         Value        DirichletCondition\n", ind);
    Zuzax::writelog(buf);
    drawline(indentSpaces + 2, 60);
    int jDir = 0;
    for (int k = 0; k < numVar; k++) {
      VarType &vt = variableNameListNode[k];
      string name = vt.VariableName(15);
      double sval = (*soln_GlAll_ptr)[eqnStart + k];
      sprintf(buf, "%s   %-15s   %-10.4E ", ind, name.c_str(), sval);
      Zuzax::writelog(buf);
      if (SpecFlag_NE[k] != 0) {
        sprintf(buf, " (Dir %d val = %-10.4E)", jDir, Value_NE[jDir]);
        Zuzax::writelog(buf);
        jDir++;
      }
      sprintf(buf, "\n");
      Zuzax::writelog(buf);
    }
    drawline(indentSpaces + 2, 60);

    double deltaV = ElectrodeC_->voltage();
    sprintf(buf, "%s   Delta Voltage = %g volts\n", ind, deltaV);
    Zuzax::writelog(buf);
    double Ess = ElectrodeC_->openCircuitVoltage(0);
    sprintf(buf, "%s   Ess Voltage = %g volts\n", ind, Ess);
    Zuzax::writelog(buf);
    double op = ElectrodeC_->overpotential(0);
    sprintf(buf, "%s   Overpotential = %g volts\n", ind, op);
    Zuzax::writelog(buf);
    sprintf(buf, "%s   Current       = %g amps/m2\n", ind, icurr * surfaceArea_);
    Zuzax::writelog(buf);
    sprintf(buf, "%s   Li+ Production Rate = %g kmol/m2/s\n", ind, electrodeSpeciesProdRates_[speciesIndex0]
        * surfaceArea_);
    Zuzax::writelog(buf);

    drawline(indentSpaces, 80);
  }
  print0_sync_end(0, ss, *(LI_ptr_->Comm_ptr_));
}

//=====================================================================================================================
// Method for writing the header for the surface domain to a tecplot file.
void
SurDomain_FlatFeS2Cathode::writeSolutionTecplotHeader()
{
  int locGbNode = SDD_.LocGbNode;
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid ; //only proc 0 should write

  if (doWrite) {

    // get the NodeVars object pertaining to this global node
    GlobalIndices *gi = LI_ptr_->GI_ptr_;
    NodalVars *nv = gi->NodalVars_GbNode[locGbNode];
    std::vector<VarType> &variableNameListNode = nv->VariableNameList_EqnNum;
    int numVar = nv->NumEquations;
    
    //open tecplot file
    FILE* ofp;
    ofp = fopen( "FlatFeS2Cathode.dat", "w");

    //write title and variable list
    string sss = id();
    fprintf( ofp, "TITLE = \"Solution on Domain %s\"\n",sss.c_str() );
    fprintf( ofp, "VARIABLES = ");

    fprintf( ofp, "\"t [s]\"  \n" );
    fprintf( ofp, "\"x [m]\"  \n" );
    fprintf( ofp, "\"<greek>D</greek> Voltage [V]\" \n" );
    fprintf( ofp, "\"Open circuit voltage [V]\" \n" );
    fprintf( ofp, "\"Overpotential [V]\" \n" );
    fprintf( ofp, "\"Current [amps/m^2]\" \n" );
    fprintf( ofp, "\"Li+ Prod Rate [kmol/m2/s]\"\n" );

    for (int k = 0; k < numVar; k++) {
      VarType &vt = variableNameListNode[k];
      string name = vt.VariableName(15);
      fprintf( ofp, "\"%s\" \t", name.c_str() );
    }

    fprintf(ofp, "\n" );
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
SurDomain_FlatFeS2Cathode::writeSolutionTecplot(const Epetra_Vector *soln_GlAll_ptr,
						const Epetra_Vector *solnDot_GlAll_ptr,
						const double t )
{
  int locGbNode = SDD_.LocGbNode;
  int mypid = LI_ptr_->Comm_ptr_->MyPID();
  bool doWrite = !mypid ; //only proc 0 should write

  // get the NodeVars object pertaining to this global node
  GlobalIndices *gi = LI_ptr_->GI_ptr_;
  NodalVars *nv = gi->NodalVars_GbNode[locGbNode];
  int eqnStart = nv->EqnStart_GbEqnIndex;

  int numVar = nv->NumEquations;
  int sf = ElectrodeC_->solnPhaseIndex();
  size_t speciesIndex0 = ElectrodeC_->globalSpeciesIndex(sf, 0);

  updateDependencies(soln_GlAll_ptr, t);


  if (doWrite) {

    FILE* ofp;
    ofp = fopen( "FlatFeS2Cathode.dat", "a");

    /*
     * Calculate the rates of production of all species in the Electrode
     */
    double icurr = ElectrodeC_->getNetSurfaceProductionRatesCurrent(0, &electrodeSpeciesProdRates_[0]);

    //time
    fprintf( ofp, "%g \t", t );

    //x-position
    double x0 = nv->x0NodePos();
    fprintf( ofp, "%g \t", x0 );

    //Delta Voltage
    double deltaV = ElectrodeC_->voltage();
    fprintf( ofp, "%g \t", deltaV );

    //Open circuit voltage
    double Ess = ElectrodeC_->openCircuitVoltage(0);
    fprintf( ofp, "%g \t", Ess );

    //Overpotential
    double op = ElectrodeC_->overpotential(0);
    fprintf( ofp, "%g \t", op );

    //Current 
    fprintf( ofp, "%g \t", icurr );

    //Li+ Production Rate 
    fprintf( ofp, "%g \t", electrodeSpeciesProdRates_[speciesIndex0] );

// Write general variables
    for (int k = 0; k < numVar; k++) {
      double sval = (*soln_GlAll_ptr)[eqnStart + k];
      fprintf(ofp, "%g \t", sval );

    }

    fprintf(ofp, "\n" );
    fclose(ofp);
  }

}
//=====================================================================================================================
// Generate the initial conditions
/*
 *   For surface dirichlet conditions, we impose the t = 0- condition.
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
SurDomain_FlatFeS2Cathode::initialConditions(const bool doTimeDependentResid,
                                             Epetra_Vector *soln_ptr,
                                             Epetra_Vector *solnDot,
                                             const double t,
                                             const double delta_t)
{
  int ieqn;
  Epetra_Vector &soln = *soln_ptr;
  /*
   *  Quick return if we don't own the node that the boundary condition
   *  is to be applied on.
   */
  if (NumOwnedNodes == 0) {
    return;
  }

  /*
   *   Figure out the node to apply the equation on
   *   It is one node -> This can be done using base class level
   */
  // NodalVars* nodeCent = LI_ptr->NodalVars_LcNode[Index_LcNode];

  /*
   *  Figure out the equation start for this node
   *   We start at the start of the equations for this node
   *   because we will be applying dirichlet conditions on the bulk
   *   equations.
   */
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

  /*
   *  Figure out the surface domain on which these boundary conditions
   *  will be applied. There is one controlling one
   */

  /*
   *   Calculate the state of the system on that controlling surface domain
   */

  /*
   *  Loop over the equations that the boundary conditions are going to be applied to
   */
  for (int i = 0; i < NumNodeEqns; i++) {
    if (SpecFlag_NE[i]) {
      /*
       *  For Dirichlet equations, replace the solution
       */
      ieqn = index_EqnStart + i;
      double val = Value_NE[i];
      soln[ieqn] = val;
    }
  }
  if (voltageVarBCType_ == 1) {

    int offsetSD = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[0];
    int EQ_Current_offset_SD = offsetSD + SDD_.EquationIndexStart_EqName[Current_Specification];
    soln[index_EqnStart + EQ_Current_offset_SD] = 1.9;
  }
}
//=====================================================================================================================
//  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void
SurDomain_FlatFeS2Cathode::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
					 Epetra_Vector_Ghosted & atolVector,
					 const Epetra_Vector_Ghosted * const atolV)
{
  /*
   *  Quick return if we don't own the node that the boundary condition
   *  is to be applied on.
   */
  if (NumOwnedNodes == 0) {
    return;
  }
  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
  int offsetSD = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[0];
  int EQ_Current_offset_SD = SDD_.VariableIndexStart_VarName[Voltage];

  /*
   * Set the atol value for the electrolyte voltage
   *      arithmetically scaled.-> so this is a characteristic value
   *         1 kcal gmol-1 = 0.05 volts
   */
  atolVector[index_EqnStart + offsetSD + EQ_Current_offset_SD] = 1.0E-6;
}
//=====================================================================================================================
//!  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*!
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void 
SurDomain_FlatFeS2Cathode::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
					       const Epetra_Vector_Ghosted & solnDot,
					       Epetra_Vector_Ghosted & atolVector,
					       const Epetra_Vector_Ghosted * const atolV)
{
  /*
   *  Quick return if we don't own the node that the boundary condition
   *  is to be applied on.
   */
  if (NumOwnedNodes == 0) {
    return;
  }

  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
  int offsetSD = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[0];
  int EQ_Current_offset_SD = SDD_.VariableIndexStart_VarName[Voltage];

  /*
   * Set the atol value for the electrolyte voltage
   *      arithmetically scaled.-> so this is a characteristic value
   *         1 kcal gmol-1 = 0.05 volts
   */
  atolVector[index_EqnStart + offsetSD + EQ_Current_offset_SD] = 1.0E-6;  
}
//=====================================================================================================================
//!  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*!
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void 
SurDomain_FlatFeS2Cathode::setAtolDeltaDamping(double atolDefault,  double relcoeff, 
					     const Epetra_Vector_Ghosted & soln, 
					     Epetra_Vector_Ghosted & atolDeltaDamping,
					     const Epetra_Vector_Ghosted * const atolV)
{
  /*
   *  Quick return if we don't own the node that the boundary condition
   *  is to be applied on.
   */
  if (NumOwnedNodes == 0) {
    return;
  }

  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
  int offsetSD = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[0];
  int EQ_Current_offset_SD = SDD_.VariableIndexStart_VarName[Voltage];

  /*
   * Set the atol value for the electrolyte voltage
   *      arithmetically scaled.-> so this is a characteristic value
   *         1 kcal gmol-1 = 0.05 volts
   */
  atolDeltaDamping[index_EqnStart + offsetSD + EQ_Current_offset_SD] = 0.05;  
}
//=====================================================================================================================
//!  Fill the vector atolVector with the values from the DomainDescription for abs tol
/*!
 *  @param atoldDefault default value
 *  @param soln soln vector
 *  @param atolVector Reference for the atol vector to fill up
 */
void 
SurDomain_FlatFeS2Cathode::setAtolDeltaDamping_DAEInit(double atolDefault,  double relcoeff, 
						     const Epetra_Vector_Ghosted & soln, 
						     const Epetra_Vector_Ghosted & solnDot,
						     Epetra_Vector_Ghosted & atolDeltaDamping,
						     const Epetra_Vector_Ghosted * const atolV)
{
  /*
   *  Quick return if we don't own the node that the boundary condition
   *  is to be applied on.
   */
  if (NumOwnedNodes == 0) {
    return;
  }

  int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
  int offsetSD = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[0];
  int EQ_Current_offset_SD = SDD_.VariableIndexStart_VarName[Voltage];

  /*
   * Set the atol value for the electrolyte voltage
   *      arithmetically scaled.-> so this is a characteristic value
   *         1 kcal gmol-1 = 0.05 volts
   */
  atolDeltaDamping[index_EqnStart + offsetSD + EQ_Current_offset_SD] = 0.05;  
}

//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
