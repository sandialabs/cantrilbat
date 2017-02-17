/**
 * @file m1d_SurDomain1d.cpp
 *  Basic object to calculate the surface residuals for surface domains.
 */

#include "m1d_SurDomain1D.h"

#include "m1d_NodalVars.h"
#include "m1d_SurfDomainTypes.h"
#include "m1d_BulkDomain1D.h"

#include "m1d_GlobalIndices.h"
#include "m1d_SurfDomainDescription.h"
#include "m1d_DomainLayout.h"

#include "m1d_Comm.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
SurDomain1D::SurDomain1D(SurfDomainDescription& bdd) :
    Domain1D(),
    SDD_(bdd),
    NumOwnedNodes(0),
    FirstGbNode(-1),
    LastOwnedGbNode(-1),
    NumBCs(-1),
    IOwnLeft(false),
    IOwnRight(false),
    Index_LcNode(-1),
    NodalVarPtr(0),
    Index_NodalSD_(-1),
    LI_ptr_(nullptr),
    NumNodeEqns(0),
    IsAlgebraic_Node(0),
    IsArithmeticScaled_Node(0),
    NumDomainEqnsLeft_(0),
    NumDomainEqnsRight_(0),
    DiffFluxLeftBulkDomain_LastResid_NE(0),
    DiffFluxRightBulkDomain_LastResid_NE(0),
    TotalFluxLeftBulkDomain_LastResid_NE(0),
    TotalFluxRightBulkDomain_LastResid_NE(0),
    VarVectorLeftBulkDomain_LastResid_NE(0),
    VarVectorRightBulkDomain_LastResid_NE(0)
{
    NumDomainEqns = SDD_.NumEquationsPerNode;
}
//==================================================================================================================================
SurDomain1D::SurDomain1D(const SurDomain1D& r) :
    Domain1D(),
    SDD_(r.SDD_),
    NumOwnedNodes(0),
    FirstGbNode(-1),
    LastOwnedGbNode(-1),
    NumBCs(-1),
    IOwnLeft(false),
    IOwnRight(false),
    Index_LcNode(-1),
    NodalVarPtr(0),
    Index_NodalSD_(-1),
    LI_ptr_(nullptr),
    NumNodeEqns(0),
    IsAlgebraic_Node(0),
    IsArithmeticScaled_Node(0),
    NumDomainEqnsLeft_(0),
    NumDomainEqnsRight_(0),
    DiffFluxLeftBulkDomain_LastResid_NE(0),
    DiffFluxRightBulkDomain_LastResid_NE(0),
    TotalFluxLeftBulkDomain_LastResid_NE(0),
    TotalFluxRightBulkDomain_LastResid_NE(0),
    VarVectorLeftBulkDomain_LastResid_NE(0),
    VarVectorRightBulkDomain_LastResid_NE(0)
{
    SurDomain1D::operator=(r);
}
//==================================================================================================================================
SurDomain1D::~SurDomain1D()
{
}
//==================================================================================================================================
SurDomain1D&
SurDomain1D::operator=(const SurDomain1D& r)
{
    if (this == &r) {
        return *this;
    }
    Domain1D::operator=(r);

    SDD_ = r.SDD_;
    NumOwnedNodes = r.NumOwnedNodes;
    FirstGbNode = r.FirstGbNode;
    LastOwnedGbNode = r.LastOwnedGbNode;
    NumBCs = r.NumBCs;
    IOwnLeft = r.IOwnLeft;
    IOwnRight = r.IOwnRight;
    Index_LcNode = r.Index_LcNode;
    NodalVarPtr = r.NodalVarPtr;
    Index_NodalSD_ = r.Index_NodalSD_;
    LI_ptr_ = r.LI_ptr_;
    NumNodeEqns = r.NumNodeEqns;
    IsAlgebraic_Node = r.IsAlgebraic_Node;
    IsArithmeticScaled_Node = r.IsArithmeticScaled_Node;
    NumDomainEqnsLeft_  = r.NumDomainEqnsLeft_;
    NumDomainEqnsRight_ = r.NumDomainEqnsRight_;

    DiffFluxLeftBulkDomain_LastResid_NE = r.DiffFluxLeftBulkDomain_LastResid_NE;
    DiffFluxRightBulkDomain_LastResid_NE = r.DiffFluxRightBulkDomain_LastResid_NE;
    TotalFluxLeftBulkDomain_LastResid_NE = r.TotalFluxLeftBulkDomain_LastResid_NE;
    TotalFluxRightBulkDomain_LastResid_NE = r.TotalFluxRightBulkDomain_LastResid_NE;
    VarVectorLeftBulkDomain_LastResid_NE = r.VarVectorLeftBulkDomain_LastResid_NE;
    VarVectorRightBulkDomain_LastResid_NE = r.VarVectorRightBulkDomain_LastResid_NE;
    Resid_BeforeSurDomain_NE = r.Resid_BeforeSurDomain_NE;

    DomainResidVector_LastResid_NE = r.DomainResidVector_LastResid_NE;


    return *this;
}
//==================================================================================================================================
std::string SurDomain1D::id() const
{
    int id = SDD_.ID();
    if (m_id != "") {
        return m_id;
    } else {
        return std::string("SurDomain1D_") + ZZCantera::int2str(id);
    }
}
//==================================================================================================================================
void SurDomain1D::domain_prep(LocalNodeIndices* const li_ptr)
{
    LI_ptr_ = li_ptr;
    /*
     *  Find the surface domain on which to apply the bc.
     */
    /*
     * Find the global node where the boundary condition is applied
     */
    int gbnode = SDD_.LocGbNode;
    /*
     * Translate from global node number to local node number.
     * The answer may be -1
     */
    Index_LcNode = LI_ptr_->GbNodeToLcNode(gbnode);

    GlobalIndices* gi_ptr = LI_ptr_->GI_ptr_;

    /*
     * If the node doesn't exist on this processor we fill in some
     * empty information and return
     */
    if (Index_LcNode < 0) {
        NumOwnedNodes = 0;
        FirstGbNode = gbnode;
        LastOwnedGbNode = gbnode;
        IOwnLeft = false;
        IOwnRight = false;
        NumBCs = 0;
        NumDomainEqns = 0;
        NodalVarPtr = gi_ptr->NodalVars_GbNode[gbnode];
        return;
    } else {
        NumOwnedNodes = 1;
        FirstGbNode = gbnode;
        LastOwnedGbNode = gbnode;
        IOwnLeft = true;
        IOwnRight = true;
        /*
         * Get the pointer to the NodalVar structure
         */
        NodalVarPtr = LI_ptr_->NodalVars_LcNode[Index_LcNode];

        // Get the total number of equations defined at this node
        NumNodeEqns = NodalVarPtr->NumEquations;

        AssertTrace((gi_ptr->NodalVars_GbNode[gbnode]) == (LI_ptr_->NodalVars_LcNode[Index_LcNode]));

        /*
         * Download the number of equations that are defined on the domain
         */
        NumDomainEqns = SDD_.NumEquationsPerNode;

        /*
         * Determine Index_NodalSD_
         */
        for (int i = 0; i < NodalVarPtr->NumSurfDomains; i++) {
            int id = NodalVarPtr->SurfDomainIndex_SDN[i];
            if (id == SDD_.ID()) {
                Index_NodalSD_ = id;
                break;
            }
        }
        if (Index_NodalSD_ == -1) {
            throw m1d_Error("SurDomain1D::initialConditions", "logic error");
        }
        /*
         * Set the IsAlgebraic flag for the node list to "I don't know"
         */
        IsAlgebraic_Node.resize(NumNodeEqns, -1);

        /*
         * Set the IsArithmeticScaled flag for the node list to "I don't know"
         */
        IsArithmeticScaled_Node.resize(NumNodeEqns, -1);
        //
        //  Resize based on number of equations
        //
        Resid_BeforeSurDomain_NE.resize(NumNodeEqns, 0.0);
        DomainResidVector_LastResid_NE.resize(NumNodeEqns, 0.0);
    }

    // Get Pointer to the left domain's domain description
    DomainDescription* ldd = SDD_.LeftDomain;
    NumDomainEqnsLeft_ = 0;
    if (ldd) {
        BulkDomainDescription* lbdd = dynamic_cast<BulkDomainDescription*>(ldd);
        // Download the number of equations in this domain
        NumDomainEqnsLeft_ = lbdd->NumEquationsPerNode;

        DiffFluxLeftBulkDomain_LastResid_NE.resize(NumDomainEqnsLeft_, 0.0);
        TotalFluxLeftBulkDomain_LastResid_NE.resize(NumDomainEqnsLeft_, 0.0);
        VarVectorLeftBulkDomain_LastResid_NE.resize(NumDomainEqnsLeft_, 0.0);
    }

    // Get Pointer to the right domain's domain description
    DomainDescription* rdd = SDD_.RightDomain;
    if (rdd) {
        BulkDomainDescription* rbdd = dynamic_cast<BulkDomainDescription*>(rdd);
        // Download the number of equations in this domain
        NumDomainEqnsRight_ = rbdd->NumEquationsPerNode;

        DiffFluxRightBulkDomain_LastResid_NE.resize(NumDomainEqnsRight_, 0.0);
        TotalFluxRightBulkDomain_LastResid_NE.resize(NumDomainEqnsRight_, 0.0);
        VarVectorRightBulkDomain_LastResid_NE.resize(NumDomainEqnsRight_, 0.0);

        if (!ldd) {
            DiffFluxLeftBulkDomain_LastResid_NE.resize(NumDomainEqnsRight_, 0.0);
            TotalFluxLeftBulkDomain_LastResid_NE.resize(NumDomainEqnsRight_, 0.0);
            VarVectorLeftBulkDomain_LastResid_NE.resize(NumDomainEqnsRight_, 0.0);
        }
    } else {
        if (ldd) {
            DiffFluxRightBulkDomain_LastResid_NE.resize(NumDomainEqnsLeft_, 0.0);
            TotalFluxRightBulkDomain_LastResid_NE.resize(NumDomainEqnsLeft_, 0.0);
            VarVectorRightBulkDomain_LastResid_NE.resize(NumDomainEqnsLeft_, 0.0);
        }
    }
}
//==================================================================================================================================
void SurDomain1D::residEval(Epetra_Vector& res, const bool doTimeDependentResid, const Epetra_Vector* const soln_ptr,
                            const Epetra_Vector* const solnDot_ptr, const Epetra_Vector* const solnOld_ptr, const double t,
                            const double rdelta_t, ResidEval_Type_Enum residType, const Solve_Type_Enum solveType)
{
    /*
     *   We call an error routine because the base class residual calculation should
     *   never be called
     */
    residType_Curr_ = residType;

    DomainDescription* ldd = SDD_.LeftBulk;
    BulkDomain1D* lbd1d = 0;
    BulkDomainDescription* lbdd = 0;
    if (ldd) {
        lbdd = dynamic_cast<BulkDomainDescription*>(ldd);
        lbd1d = lbdd->BulkDomainPtr_;
        if (lbd1d) {
            for (int i = 0; i < NumDomainEqnsLeft_; i++) {
                DiffFluxLeftBulkDomain_LastResid_NE[i] = lbd1d->DiffFluxRightBound_LastResid_NE[i];
                TotalFluxLeftBulkDomain_LastResid_NE[i] = lbd1d->TotalFluxRightBound_LastResid_NE[i];
                VarVectorLeftBulkDomain_LastResid_NE[i] = lbd1d->VarVectorRightBound_LastResid_NE[i];
            }
        }
    }

    DomainDescription* rdd = SDD_.RightBulk;
    BulkDomain1D* rbd1d = 0;
    if (rdd) {
        BulkDomainDescription* rbdd = dynamic_cast<BulkDomainDescription*>(rdd);
        rbd1d = rbdd->BulkDomainPtr_;
        if (rbd1d) {
            for (int i = 0; i < NumDomainEqnsRight_; i++) {
                DiffFluxRightBulkDomain_LastResid_NE[i] = rbd1d->DiffFluxLeftBound_LastResid_NE[i];
                TotalFluxRightBulkDomain_LastResid_NE[i] = rbd1d->TotalFluxLeftBound_LastResid_NE[i];
                VarVectorRightBulkDomain_LastResid_NE[i] = rbd1d->VarVectorLeftBound_LastResid_NE[i];
            }
        }
    } else {
        if (lbd1d) {
            for (int i = 0; i < NumDomainEqnsRight_; i++) {
                DiffFluxRightBulkDomain_LastResid_NE[i] = lbd1d->DiffFluxRightBound_LastResid_NE[i];
                TotalFluxRightBulkDomain_LastResid_NE[i] = lbd1d->TotalFluxRightBound_LastResid_NE[i];
                VarVectorRightBulkDomain_LastResid_NE[i] = lbd1d->VarVectorRightBound_LastResid_NE[i];
            }
        }
    }
    if (!ldd) {
        if (rbd1d) {
            for (int i = 0; i < NumDomainEqnsLeft_; i++) {
                DiffFluxLeftBulkDomain_LastResid_NE[i] = rbd1d->DiffFluxLeftBound_LastResid_NE[i];
                TotalFluxLeftBulkDomain_LastResid_NE[i] = rbd1d->TotalFluxLeftBound_LastResid_NE[i];
                VarVectorLeftBulkDomain_LastResid_NE[i] = rbd1d->VarVectorLeftBound_LastResid_NE[i];
            }
        }
    }
    err("residEval()");
}
//==================================================================================================================================
// Transfer the bulk flux vectors to the surface flux vectors
/*
 *  This routine update the following right flux vectors from the neighboring right bulk domain
 *      DiffFluxRightBulkDomain_LastResid_NE[i]
 *	  TotalFluxRightBulkDomain_LastResid_NE[i]
 *	  VarVectorRightBulkDomain_LastResid_NE[i]
 *
 *  It then updates the corresponding left flux vectors from the left bulk domain
 *
 *  This routine is typically called from the residEval routine. However, it's modular
 *  enough to be carved out as its own routine.
 */
void SurDomain1D::updateBulkFluxVectors()
{
    DomainDescription* ldd = SDD_.LeftBulk;
    BulkDomain1D* lbd1d = 0;
    BulkDomainDescription* lbdd = 0;
    if (ldd) {
        lbdd = dynamic_cast<BulkDomainDescription*>(ldd);
        lbd1d = lbdd->BulkDomainPtr_;
        if (lbd1d) {
            for (int i = 0; i < NumDomainEqnsLeft_; i++) {
                DiffFluxLeftBulkDomain_LastResid_NE[i] = lbd1d->DiffFluxRightBound_LastResid_NE[i];
                TotalFluxLeftBulkDomain_LastResid_NE[i] = lbd1d->TotalFluxRightBound_LastResid_NE[i];
                VarVectorLeftBulkDomain_LastResid_NE[i] = lbd1d->VarVectorRightBound_LastResid_NE[i];
            }
        }
    }

    DomainDescription* rdd = SDD_.RightBulk;
    BulkDomain1D* rbd1d = 0;
    if (rdd) {
        BulkDomainDescription* rbdd = dynamic_cast<BulkDomainDescription*>(rdd);
        rbd1d = rbdd->BulkDomainPtr_;
        if (rbd1d) {
            for (int i = 0; i < NumDomainEqnsRight_; i++) {
                DiffFluxRightBulkDomain_LastResid_NE[i] = rbd1d->DiffFluxLeftBound_LastResid_NE[i];
                TotalFluxRightBulkDomain_LastResid_NE[i] = rbd1d->TotalFluxLeftBound_LastResid_NE[i];
                VarVectorRightBulkDomain_LastResid_NE[i] = rbd1d->VarVectorLeftBound_LastResid_NE[i];
            }
        }
    } else {
        if (lbd1d) {
            for (int i = 0; i < NumDomainEqnsRight_; i++) {
                DiffFluxRightBulkDomain_LastResid_NE[i] = lbd1d->DiffFluxRightBound_LastResid_NE[i];
                TotalFluxRightBulkDomain_LastResid_NE[i] = lbd1d->TotalFluxRightBound_LastResid_NE[i];
                VarVectorRightBulkDomain_LastResid_NE[i] = lbd1d->VarVectorRightBound_LastResid_NE[i];
            }
        }
    }
    if (!ldd) {
        if (rbd1d) {
            for (int i = 0; i < NumDomainEqnsLeft_; i++) {
                DiffFluxLeftBulkDomain_LastResid_NE[i] = rbd1d->DiffFluxLeftBound_LastResid_NE[i];
                TotalFluxLeftBulkDomain_LastResid_NE[i] = rbd1d->TotalFluxLeftBound_LastResid_NE[i];
                VarVectorLeftBulkDomain_LastResid_NE[i] = rbd1d->VarVectorLeftBound_LastResid_NE[i];
            }
        }
    }
}
//==================================================================================================================================
void SurDomain1D::initialConditions(const bool doTimeDependentResid, Epetra_Vector* const soln_ptr,
                                    Epetra_Vector* const solnDot_ptr, const double t, const double delta_t)
{
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
     *   because we will be applying dirichlet conditions on the bulk
     *   equations.
     */
    size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
    /*
     *   Find the offset of this domain on the node
     */
    size_t offsetSD = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[Index_NodalSD_];
    /*
     *   Find the number of equations actually owned by this surface domain
     */
    size_t numEqSD = SDD_.NumEquationsPerNode;
    /*
     *  Zero the solution components corresponding to this solution domain
     */
    size_t istart = index_EqnStart + offsetSD;
    if (solnDot_ptr) {
        Epetra_Vector& solnDot = *solnDot_ptr;
        for (size_t ieqn = 0; ieqn < numEqSD; ieqn++) {
            soln[istart + ieqn] = 0.0;
            solnDot[istart + ieqn] = 0.0;
        }
    } else {
        for (size_t ieqn = 0; ieqn < numEqSD; ieqn++) {
            soln[istart + ieqn] = 0.0;
        }
    }
}
//==================================================================================================================================
void SurDomain1D::saveDomain(ZZCantera::XML_Node& oNode, const Epetra_Vector* const soln_GLALL_ptr,
                             const Epetra_Vector* const solnDot_GLALL_ptr, const double t, bool duplicateOnAllProcs)
{
    // const double* s = soln_GLALL_ptr + loc();
    // Find the number of global equations on this domain, whether it's local or not
    //int numEquationsGb = SDD_.NumEquationsPerNode;
    // Find the global node number of the node where this domain resides
    //int locGbNode = SDD_.LocGbNode;

    size_t eqnStart = NodalVarPtr->EqnStart_GbEqnIndex;
    //XML_Node& inlt = o.addChild("inlet");
    ZZCantera::XML_Node& inlt = oNode.addChild("domain");
    size_t numVar = NodalVarPtr->NumEquations;
    inlt.addAttribute("id", id());
    inlt.addAttribute("points", 1);
    inlt.addAttribute("type", "surface");
    inlt.addAttribute("numVariables", numVar);
    double x0pos = NodalVarPtr->x0NodePos();
    double xpos =  NodalVarPtr->xNodePos();
    double xfrac = NodalVarPtr->xFracNodePos();
    ZZctml::addFloat(inlt, "X0", x0pos, "", "", ZZCantera::Undef, ZZCantera::Undef);
    ZZctml::addFloat(inlt, "X", xpos, "", "", ZZCantera::Undef, ZZCantera::Undef);
    ZZctml::addFloat(inlt, "Xfraction", xfrac, "", "", ZZCantera::Undef, ZZCantera::Undef);

    for (size_t k = 0; k < numVar; k++) {
        double sval = (*soln_GLALL_ptr)[eqnStart + k];
        std::string nm = NodalVarPtr->VariableName(k);
        VarType vv = NodalVarPtr->VariableNameList_EqnNum[k];
        std::string type = VarType::VarMainName(vv.VariableType);
        ZZctml::addFloat(inlt, nm, sval, "", "", ZZCantera::Undef, ZZCantera::Undef);
    }
}
//==================================================================================================================================
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
void
SurDomain1D::readDomain(const ZZCantera::XML_Node& simulationNode, Epetra_Vector* const soln_GLALL_ptr,
                        Epetra_Vector* const solnDot_GLALL_ptr, double globalTimeRead)
{
    /*
     *   Find the global node number of the node where this surface domain resides
     */
    //int locGbNode = SDD_.LocGbNode;

    // Number of equations per node
    //int numEquationsPerNode = SDD_.NumEquationsPerNode;
    std::string ids = id();
    ZZCantera::XML_Node* domainNode_ptr = simulationNode.findNameID("domain", ids);
    /*
     * get the NodeVars object pertaining to this global node
     */
    size_t numVar = NodalVarPtr->NumEquations;

    /*
     *  Get the global equation number index of the unknowns at this surface domain
     */
    size_t eqnStart = NodalVarPtr->EqnStart_GbEqnIndex;

    std::string iidd = (*domainNode_ptr)["id"];

    std::string s_points  = (*domainNode_ptr)["points"];
    size_t points = atoi(s_points.c_str());
    if (points != 1) {
        printf("we have an unequal number of points\n");
        exit(-1);
    }
    std::string ttype  = (*domainNode_ptr)["type"];
    std::string snumVar  = (*domainNode_ptr)["numVariables"];
    size_t numVarRstart = atoi(snumVar.c_str());
    if (numVarRstart != numVar) {
        printf("we have an unequal number of variables at the node\n");
        exit(-1);
    }

    double x0pos = ZZctml::getFloat(*domainNode_ptr, "X0", "toSI");
    double xpos = ZZctml::getFloat(*domainNode_ptr, "X", "toSI");
    double xfrac = ZZctml::getFloat(*domainNode_ptr, "Xfraction", "toSI");
    NodalVarPtr->setupInitialNodePosition(x0pos, xfrac);
    NodalVarPtr->changeNodePosition(xpos);

    for (size_t k = 0; k < numVar; k++) {
        double sval = (*soln_GLALL_ptr)[eqnStart + k];
        std::string nm = NodalVarPtr->VariableName(k);
        VarType vv = NodalVarPtr->VariableNameList_EqnNum[k];
        std::string type = VarType::VarMainName(vv.VariableType);
        sval = ZZctml::getFloat(*domainNode_ptr, nm, "toSI");
        (*soln_GLALL_ptr)[eqnStart + k] = sval;
    }
}
//==================================================================================================================================
void SurDomain1D::fillIsAlgebraic(Epetra_IntVector& isAlgebraic)
{
    // Find the global node number of the node where this domain resides
    size_t locGbNode = SDD_.LocGbNode;

    // Get the NodeVars object pertaining to this global node
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
    NodalVars* nv = gi->NodalVars_GbNode[locGbNode];
    AssertTrace(NodalVarPtr == nv);

    int mySDD_ID = SDD_.ID();

    // This node may not exist on this processor. If it doesn't, there is no entry in isAlgebraic[] to fill
    if (Index_LcNode < 0) {
        return;
    }

    size_t bmatch = NodalVarPtr->SurfDomainIndex_fromID[mySDD_ID];
    size_t indexSurfDomainOffset = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[bmatch];
    /*
     *   Figure out the equation start for this node
     *   We start at the start of the equations for this node
     *   because we will be applying dirichlet conditions on the bulk
     *   equations.
     */
    size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

    size_t numVar = SDD_.NumEquationsPerNode;
    for (size_t k = 0; k < numVar; k++) {
        int isA = SDD_.IsAlgebraic_NE[k];
        isAlgebraic[index_EqnStart + indexSurfDomainOffset + k] = isA;
    }
    // Override the bulk domains too if we have information that overrides the equations for the bulk domains
    for (size_t k = 0; k < static_cast<size_t>(NumNodeEqns); k++) {
        int isA = IsAlgebraic_Node[k];
        if (isA >= 0) {
            isAlgebraic[index_EqnStart + k] = isA;
        }
    }
}
//==================================================================================================================================
void 
SurDomain1D::fillIsArithmeticScaled(Epetra_IntVector& isArithmeticScaled)
{
    // Find the global node number of the node where this domain resides
    size_t locGbNode = SDD_.LocGbNode;

    // Get the NodeVars object pertaining to this global node
    GlobalIndices* gi = LI_ptr_->GI_ptr_;
    NodalVars* nv = gi->NodalVars_GbNode[locGbNode];
    AssertTrace(NodalVarPtr == nv);

    int mySDD_ID = SDD_.ID();

    // This node may not exist on this processor. If it doesn't, there is no entry in isAlgebraic[] to fill
    if (Index_LcNode < 0) {
        return;
    }

    size_t bmatch = NodalVarPtr->SurfDomainIndex_fromID[mySDD_ID];
    size_t indexSurfDomainOffset = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[bmatch];
    /*
     *   Figure out the equation start for this node
     *   We start at the start of the equations for this node
     *   because we will be applying dirichlet conditions on the bulk
     *   equations.
     */
    size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];

    size_t numVar = SDD_.NumEquationsPerNode;
    for (size_t k = 0; k < numVar; k++) {
        int isA = SDD_.IsArithmeticScaled_NE[k];
        isArithmeticScaled[index_EqnStart + indexSurfDomainOffset + k] = isA;
    }
    // Override the bulk domains too if we have information that overrides the equations for the bulk domains
    for (size_t k = 0; k < static_cast<size_t>(NumNodeEqns); k++) {
        int isA = IsArithmeticScaled_Node[k];
        if (isA >= 0) {
            isArithmeticScaled[index_EqnStart + k] = isA;
        }
    }
}
//===================================================================================================================================
void
SurDomain1D::calcDeltaSolnVariables(const double t, const Epetra_Vector& soln, const Epetra_Vector* const solnDot_ptr,
                                    Epetra_Vector& deltaSoln, const Epetra_Vector* const atolVector_ptr,
                                    const Solve_Type_Enum solveType, const Epetra_Vector* const solnWeights)
{

    int mySDD_ID = SDD_.ID();
    // This node may not exist on this processor. If it doesn't, there is no entry in isAlgebraic[] to fill
    if (Index_LcNode < 0) {
        return;
    }
    int bmatch = NodalVarPtr->SurfDomainIndex_fromID[mySDD_ID];
    if (bmatch != -1) {

        /*
         *   Figure out the equation start for this node
         *   We start at the start of the equations for this node
         *   because we will be applying dirichlet conditions on the bulk
         *   equations.
         */
        size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
        size_t indexSurfDomainOffset = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[bmatch];
        //std::vector<VarType>& vnl = SDD_.VariableNameList;
        size_t numVar = SDD_.NumEquationsPerNode;
        size_t index;
        double base;
        for (size_t k = 0; k < numVar; k++) {
            //const VarType& vt = vnl[k];
            index = index_EqnStart + indexSurfDomainOffset + k;
            base = soln[index];
            deltaSoln[index] = 1.0E-5 * fabs(base) + 1.0E-9;
        }
    } else {
        printf("we shouldn't be here\n");
        exit(-1);
    }
}
//===================================================================================================================================
void SurDomain1D::setAtolVector(double atolDefault, const Epetra_Vector_Ghosted& soln,
                                Epetra_Vector_Ghosted& atolVector, const Epetra_Vector_Ghosted* const atolV)
{
    int mySDD_ID = SDD_.ID();

    // This node may not exist on this processor. If it doesn't, there is no entry in isAlgebraic[] to fill
    if (Index_LcNode < 0) {
        return;
    }
    int bmatch = NodalVarPtr->SurfDomainIndex_fromID[mySDD_ID];
    if (bmatch != -1) {
        int indexSurfDomainOffset = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[bmatch];
        /*
         *   Figure out the equation start for this node
         *   We start at the start of the equations for this node
         *   because we will be applying dirichlet conditions on the bulk
         *   equations.
         */
        size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
        size_t numVar = SDD_.NumEquationsPerNode;
        for (size_t k = 0; k < numVar; k++) {
            // int isA = SDD_.IsArithmeticScaled_NE[k];
            atolVector[index_EqnStart + indexSurfDomainOffset + k] = atolDefault;
        }
    } else {
        printf("we shouldn't be here\n");
        exit(-1);
    }
}
//===================================================================================================================================
void 
SurDomain1D::setAtolDeltaDamping(double atolDefault, double relcoeff, const Epetra_Vector_Ghosted& soln,
                                 Epetra_Vector_Ghosted& atolDeltaDamping, const Epetra_Vector_Ghosted* const atolV)
{
    int mySDD_ID = SDD_.ID();

    // This node may not exist on this processor. If it doesn't, there is no entry in isAlgebraic[] to fill
    if (Index_LcNode < 0) {
        return;
    }
    int bmatch = NodalVarPtr->SurfDomainIndex_fromID[mySDD_ID];
    size_t indexSurfDomainOffset = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[bmatch];
    /*
     *   Figure out the equation start for this node
     *   We start at the start of the equations for this node
     *   because we will be applying dirichlet conditions on the bulk
     *   equations.
     */
    size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
    size_t numVar = SDD_.NumEquationsPerNode;
    for (size_t k = 0; k < numVar; k++) {
        atolDeltaDamping[index_EqnStart + indexSurfDomainOffset + k] = relcoeff * atolDefault;
    }
}
//===================================================================================================================================
void 
SurDomain1D::setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted& soln,
                                   const Epetra_Vector_Ghosted& solnDot, Epetra_Vector_Ghosted& atolVector_DAEInit,
                                   const Epetra_Vector_Ghosted* const atolV)
{
    int mySDD_ID = SDD_.ID();

    // This node may not exist on this processor. If it doesn't, there is no entry in isAlgebraic[] to fill
    if (Index_LcNode < 0) {
        return;
    }
    int bmatch = NodalVarPtr->SurfDomainIndex_fromID[mySDD_ID];
    size_t indexSurfDomainOffset = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[bmatch];
    /*
     *   Figure out the equation start for this node
     *   We start at the start of the equations for this node
     *   because we will be applying dirichlet conditions on the bulk
     *   equations.
     */
    size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
    size_t numVar = SDD_.NumEquationsPerNode;
    for (size_t k = 0; k < numVar; k++) {
        // int isA = SDD_.IsArithmeticScaled_NE[k];
        atolVector_DAEInit[index_EqnStart + indexSurfDomainOffset + k] = atolDefault;
    }
}
//==================================================================================================================================
void 
SurDomain1D::setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, const Epetra_Vector_Ghosted& soln,
                                         const Epetra_Vector_Ghosted& solnDot, Epetra_Vector_Ghosted& atolDeltaDamping,
                                         const Epetra_Vector_Ghosted* const atolV)
{
    int mySDD_ID = SDD_.ID();

    // This node may not exist on this processor. If it doesn't, there is no entry in isAlgebraic[] to fill
    if (Index_LcNode < 0) {
        return;
    }
    int bmatch = NodalVarPtr->SurfDomainIndex_fromID[mySDD_ID];
    size_t indexSurfDomainOffset = NodalVarPtr->OffsetIndex_SurfDomainEqnStart_SDN[bmatch];
    /*
     *   Figure out the equation start for this node
     *   We start at the start of the equations for this node
     *   because we will be applying dirichlet conditions on the bulk
     *   equations.
     */
    size_t index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
    size_t numVar = SDD_.NumEquationsPerNode;
    for (size_t k = 0; k < numVar; k++) {
        atolDeltaDamping[index_EqnStart + indexSurfDomainOffset + k] = relcoeff * atolDefault;
    }
}
//===================================================================================================================================
void SurDomain1D::writeSolutionTecplotHeader()
{
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid; //only proc 0 should write

    if (doWrite) {
        std::vector<VarType>& variableNameListNode = NodalVarPtr->VariableNameList_EqnNum;
        size_t numVar = NodalVarPtr->NumEquations;

        //open tecplot file
        FILE* ofp;
        const std::string sss = id();
        char filename[20];
        sprintf(filename, "%s%s", sss.c_str(), ".dat");
        ofp = fopen(filename, "w");

        //write title and variable list
        fprintf(ofp, "TITLE = \"Solution on Domain %s\"\n", sss.c_str());
        fprintf(ofp, "VARIABLES = ");

        fprintf(ofp, "\"t [s]\"  \n");
        fprintf(ofp, "\"x [m]\"  \n");

        for (size_t k = 0; k < numVar; k++) {
            VarType& vt = variableNameListNode[k];
            std::string name = vt.VariableName(15);
            fprintf(ofp, "\"%s\" \t", name.c_str());
        }
        fprintf(ofp, "\n");
        fclose(ofp);
    }
}
//===================================================================================================================================
void SurDomain1D::writeSolutionTecplot(const Epetra_Vector_GlAll* const soln_GlAll_ptr,
                                       const Epetra_Vector_GlAll* const solnDot_GlAll_ptr, const double t)
{
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid; //only proc 0 should write

    size_t eqnStart = NodalVarPtr->EqnStart_GbEqnIndex;
    size_t numVar = NodalVarPtr->NumEquations;

    if (doWrite) {
        //open tecplot file
        FILE* ofp;
        std::string sss = id();
        char filename[20];
        sprintf(filename, "%s%s", sss.c_str(), ".dat");
        ofp = fopen(filename, "a");
        //time
        fprintf(ofp, "%g \t", t);
        //x-position
        double x0 = NodalVarPtr->x0NodePos();
        fprintf(ofp, "%g \t", x0);

        // Write general variables
        for (size_t k = 0; k < numVar; k++) {
            double sval = (*soln_GlAll_ptr)[eqnStart + k];
            fprintf(ofp, "%g \t", sval);
        }
        fprintf(ofp, "\n");
        fclose(ofp);
    }
}
//===================================================================================================================================
// Extract the double value out of a solution vector for a particular
// variable defined at a node corresponding to the surface domain
/*
 *  @param  soln_ptr   Pointer to the ghosted solution vector
 *  @param  v1         VarType of the variable to be returned
 *
 *  @return           Returns the value of the variable. If the variable doesn't exist
 *                    on the processor this routine returns the value of
 *                    -1.0E300.
 */
double SurDomain1D::extractSolnValue(Epetra_Vector_Ghosted* const soln_ptr, VarType v1)
{
    double val = -1.0E300;
    /*
     *  Figure out the equation start for this node. Note, if this node isn't located on the processor
     *  return an error indicator
     */
    if (Index_LcNode < 0) {
        return val;
    }
    int index_EqnStart = LI_ptr_->IndexLcEqns_LcNode[Index_LcNode];
    if (index_EqnStart < 0) {
        return val;
    }
    int ioffset = -1;
    size_t numNodeEqns = NodalVarPtr->NumEquations;
    for (size_t j = 0; j < numNodeEqns; j++) {
        VarType vtNode = NodalVarPtr->VariableNameList_EqnNum[j];
        if (v1 == vtNode) {
            ioffset = j;
            break;
        }
    }
    if (ioffset < 0) {
        return val;
    }
    val = (*soln_ptr)[index_EqnStart + ioffset];
    return val;
}
//===================================================================================================================================
static void drawline(int sp, int ll)
{
    for (int i = 0; i < sp; i++) {
        ZZCantera::writelog(" ");
    }
    for (int i = 0; i < ll; i++) {
        ZZCantera::writelog("-");
    }
    ZZCantera::writelog("\n");
}
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
//=====================================================================================================================
void 
SurDomain1D::showSolution(const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
                          const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr, 
                          const Epetra_Vector* const solnOld_ptr, const Epetra_Vector_Owned* const residInternal_ptr,
                          const double t, const double rdelta_t, int indentSpaces, bool duplicateOnAllProcs)
{
#ifdef DO_OLD_WAY
    showSolution0All(soln_GlAll_ptr, solnDot_GlAll_ptr, soln_ptr, solnDot_ptr, solnOld_ptr,
                     residInternal_ptr,t, rdelta_t, indentSpaces, duplicateOnAllProcs);
    return;
#else
    int locGbNode = SDD_.LocGbNode;
    //int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind = indent.c_str();

    size_t eqnStart = NodalVarPtr->EqnStart_GbEqnIndex;
    std::vector<VarType>& variableNameListNode = NodalVarPtr->VariableNameList_EqnNum;
    size_t numVar = NodalVarPtr->NumEquations;
    std::string sss = id();
    stream0 ss;

    BulkDomain1D* bd = 0;
    BulkDomainDescription* bedd = SDD_.LeftBulk;
    double* diffFlux = 0;
    if (!bedd) {
        bedd = SDD_.RightBulk;
        if (bedd) {
            bd = bedd->BulkDomainPtr_;
            diffFlux = &(bd->DiffFluxLeftBound_LastResid_NE[0]);
        }
    } else {
        bd = bedd->BulkDomainPtr_;
        diffFlux = &(bd->DiffFluxRightBound_LastResid_NE[0]);
    }


    print0_sync_start(0, ss, * (LI_ptr_->Comm_ptr_));
    if (doWrite) {
        drawline0(ss, indentSpaces, 80);
        ss.print0("%s  Solution on Surface Domain %10s : Number of variables = %d\n", ind, sss.c_str(), static_cast<int>(numVar));
        ss.print0("%s                                           : Number of boundary conditions = %d\n", ind, NumBCs);
        double x0 = NodalVarPtr->x0NodePos();
        ss.print0("%s                                           : Node %d at pos %g\n", ind, locGbNode, x0);
        drawline0(ss, indentSpaces, 80);
        ss.print0("%s     VariableName         Value        DirichletCondition\n", ind);
        drawline0(ss, indentSpaces + 2, 60);
        for (size_t k = 0; k < numVar; k++) {
            VarType& vt = variableNameListNode[k];
            std::string name = vt.VariableName(15);
            double sval = (*soln_GlAll_ptr)[eqnStart + k];
            if (vt.VariableType == Velocity_Axial) {
                size_t offset = NodalVarPtr->Offset_VarTypeVector[(size_t) Velocity_Axial];
                sval = diffFlux[offset];
                ss.print0("%s   %-15s * % -10.4E ", ind, name.c_str(), sval);
            } else {
                ss.print0("%s   %-15s   % -10.4E ", ind, name.c_str(), sval);
            }
            ss.print0("\n");
        }
        drawline0(ss, indentSpaces + 2, 60);
        drawline0(ss, indentSpaces, 80);
    }
    print0_sync_end(0, ss, * (LI_ptr_->Comm_ptr_));
#endif
}
//==================================================================================================================================
void 
SurDomain1D::showSolutionVector(std::string& solnVecName, const Epetra_Vector* const solnVector_GlAll_ptr,
                                const Epetra_Vector* const solnVector_ptr, const double t, const double rdelta_t, int indentSpaces,
                                bool duplicateOnAllProcs, FILE* const of)
{

    int locGbNode = SDD_.LocGbNode;
    //int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind = indent.c_str();

    size_t eqnStart = NodalVarPtr->EqnStart_GbEqnIndex;
    std::vector<VarType>& variableNameListNode = NodalVarPtr->VariableNameList_EqnNum;
    size_t numVar = NodalVarPtr->NumEquations;
    std::string sss = id();
    stream0 ss(of);

    print0_sync_start(0, ss, * (LI_ptr_->Comm_ptr_));
    if (doWrite) {
        drawline0(ss, indentSpaces, 100);
        ss.print0("%s  %s Vector on Surface Domain %10s : Number of variables = %d\n", ind, solnVecName.c_str(), sss.c_str(),
                  numVar);
        ss.print0("%s                                           : Number of boundary conditions = %d\n", ind, NumBCs);
        double x0 = NodalVarPtr->x0NodePos();
        ss.print0("%s                                           : Node %d at pos %g\n", ind, locGbNode, x0);
        drawline0(ss, indentSpaces, 100);
        ss.print0("%s     VariableName   GblEqnInd      Value  \n", ind);
        drawline0(ss, indentSpaces + 2, 60);
        for (size_t k = 0; k < numVar; k++) {
            VarType& vt = variableNameListNode[k];
            std::string name = vt.VariableName(15);
            double sval = (*solnVector_GlAll_ptr)[eqnStart + k];
            ss.print0("%s   %-15.15s  %7d  % -11.4E ", ind, name.c_str(), eqnStart + k, sval);
            ss.print0("\n");
        }
        drawline0(ss, indentSpaces + 2, 60);
        drawline0(ss, indentSpaces, 100);
    }
    print0_sync_end(0, ss, * (LI_ptr_->Comm_ptr_));
}
//==================================================================================================================================
void 
SurDomain1D::showSolutionIntVector(std::string& solnVecName, const Epetra_IntVector* const solnIntVector_GlAll_ptr,
                                   const Epetra_IntVector* const solnIntVector_ptr, const double t, const double rdelta_t,
                                   int indentSpaces, bool duplicateOnAllProcs, FILE* of)
{

    int locGbNode = SDD_.LocGbNode;
    //int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = (NumOwnedNodes > 0) || duplicateOnAllProcs;
    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind = indent.c_str();
    // get the NodeVars object pertaining to this global node

    size_t eqnStart = NodalVarPtr->EqnStart_GbEqnIndex;
    std::vector<VarType>& variableNameListNode = NodalVarPtr->VariableNameList_EqnNum;
    size_t numVar = NodalVarPtr->NumEquations;
    std::string sss = id();
    stream0 ss(of);

    print0_sync_start(0, ss, * (LI_ptr_->Comm_ptr_));
    if (doWrite) {
        drawline0(ss, indentSpaces, 100);
        ss.print0("%s  %s Vector on Surface Domain %10s : Number of variables = %d\n", ind, solnVecName.c_str(), sss.c_str(),
                  numVar);
        ss.print0("%s                                           : Number of boundary conditions = %d\n", ind, NumBCs);
        double x0 = NodalVarPtr->x0NodePos();
        ss.print0("%s                                           : Node %d at pos %g\n", ind, locGbNode, x0);
        drawline0(ss, indentSpaces, 100);
        ss.print0("%s     VariableName   GblEqnInd      Value  \n", ind);
        drawline0(ss, indentSpaces + 2, 60);
        for (size_t k = 0; k < numVar; k++) {
            VarType& vt = variableNameListNode[k];
            std::string name = vt.VariableName(15);
            int sval = (*solnIntVector_GlAll_ptr)[eqnStart + k];
            ss.print0("%s   %-15.15s  %7d  % -11d ", ind, name.c_str(), eqnStart + k, sval);
            ss.print0("\n");
        }
        drawline0(ss, indentSpaces + 2, 60);
        drawline0(ss, indentSpaces, 100);
    }
    print0_sync_end(0, ss, * (LI_ptr_->Comm_ptr_));
}

//==================================================================================================================================
void 
SurDomain1D::showSolution0All(const Epetra_Vector* const soln_GlAll_ptr, const Epetra_Vector* const solnDot_GlAll_ptr,
                              const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
                              const Epetra_Vector* const solnOld_ptr, const Epetra_Vector_Owned* residInternal_ptr, const double t,
                              const double rdelta_t, int indentSpaces, bool duplicateOnAllProcs)
{
    char buf[132];
    int locGbNode = SDD_.LocGbNode;
    int mypid = LI_ptr_->Comm_ptr_->MyPID();
    bool doWrite = !mypid || duplicateOnAllProcs;
    std::string indent = "";
    for (int i = 0; i < indentSpaces; i++) {
        indent += " ";
    }
    const char* ind = indent.c_str();
    size_t eqnStart = NodalVarPtr->EqnStart_GbEqnIndex;
    std::vector<VarType>& variableNameListNode = NodalVarPtr->VariableNameList_EqnNum;
    size_t numVar = NodalVarPtr->NumEquations;
    std::string sss = id();
    if (doWrite) {
        drawline(indentSpaces, 80);
        sprintf(buf, "%s  Solution on Surface Domain %10s : Number of variables = %d\n", ind, sss.c_str(), static_cast<int>(numVar));
        ZZCantera::writelog(buf);
        sprintf(buf, "%s                                           : Number of boundary conditions = %d\n", ind, NumBCs);
        ZZCantera::writelog(buf);
        double x0 = NodalVarPtr->x0NodePos();
        sprintf(buf, "%s                                           : Node %d at pos %g\n", ind, locGbNode, x0);
        ZZCantera::writelog(buf);
        drawline(indentSpaces, 80);
        sprintf(buf, "%s     VariableName         Value        DirichletCondition\n", ind);
        ZZCantera::writelog(buf);
        drawline(indentSpaces + 2, 60);
        for (size_t k = 0; k < numVar; k++) {
            VarType& vt = variableNameListNode[k];
            std::string name = vt.VariableName(15);
            double sval = (*soln_GlAll_ptr)[eqnStart + k];
            sprintf(buf, "%s   %-15s   %-10.4E ", ind, name.c_str(), sval);
            ZZCantera::writelog(buf);
            sprintf(buf, "\n");
            ZZCantera::writelog(buf);
        }
        drawline(indentSpaces + 2, 60);
        drawline(indentSpaces, 80);
    }
}
//==================================================================================================================================
void SurDomain1D::err(const char* msg)
{
    printf("Domain1D: function not implemented: %s\n", msg);
    exit(-1);
}
//==================================================================================================================================
//==================================================================================================================================
//==================================================================================================================================

SurBC_Dirichlet::SurBC_Dirichlet(SurfDomainDescription& sdd) :
    SurDomain1D(sdd),
    SpecFlag_NE(0),
    Value_NE(0),
    TimeDep_NE(0),
    BC_TimeDep_NE(0),
    BC_Type_NE(0)
{

}
//=====================================================================================================================
SurBC_Dirichlet::SurBC_Dirichlet(const SurBC_Dirichlet& r) :
    SurDomain1D(r.SDD_),
    SpecFlag_NE(0),
    Value_NE(0),
    TimeDep_NE(0),
    BC_TimeDep_NE(0),
    BC_Type_NE(0)
{
    operator=(r);
}
//=====================================================================================================================
// Destructor
SurBC_Dirichlet::~SurBC_Dirichlet()
{
}
//=====================================================================================================================
// Assignment operator
/*
 * @param r      Object to be copied into the current object
 * @return       Returns a changeable reference to the current object
 */
SurBC_Dirichlet&
SurBC_Dirichlet::operator=(const SurBC_Dirichlet& r)
{
    if (this == &r) {
        return *this;
    }
    SurDomain1D::operator=(r);
    SpecFlag_NE = r.SpecFlag_NE;
    Value_NE = r.Value_NE;
    TimeDep_NE = r.TimeDep_NE;
    BC_Type_NE = r.BC_Type_NE;
    BC_TimeDep_NE = r.BC_TimeDep_NE;

    return *this;
}
//=====================================================================================================================
// Prepare all of the indices for fast calculation of the residual
/*
 *  Here we collect all of the information necessary to
 *  speedily implement SpecFlag_NE and Value_NE within the
 *  residual calculation.
 *  We transfer the information from SDT_Dirichlet structure to
 * this structure for quick processing.
 */
void SurBC_Dirichlet::domain_prep(LocalNodeIndices* li_ptr)
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
    SDT_Dirichlet* sdt = dynamic_cast<SDT_Dirichlet*>(&SDD_);
    AssertThrow(sdt, "bad cast");
    /*
     * Zero out the current domain's setup section.
     */
    NumNodeEqns = NodalVarPtr->NumEquations;

    SpecFlag_NE.resize(NumNodeEqns);
    Value_NE.resize(NumNodeEqns);
    TimeDep_NE.resize(NumNodeEqns);
    BC_TimeDep_NE.resize(NumNodeEqns);
    BC_Type_NE.resize(NumNodeEqns);
    IsAlgebraic_Node.resize(NumNodeEqns, -1);
    for (int j = 0; j < NumNodeEqns; j++) {
        SpecFlag_NE[j] = 0;
        Value_NE[j] = 0.0;
        TimeDep_NE[j] = 0;
        BC_TimeDep_NE[j] = 0;
        BC_Type_NE[j] = 0;
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
        bool ifound = false;
        // If the subtype is -1, then we apply to all equations of the
        // variable main type
        /*
         *
         */
        if (stDir == -1) {
            for (int j = 0; j < NumNodeEqns; j++) {
                VarType eqt = NodalVarPtr->VariableNameList_EqnNum[j];
                if ((vtmDir == Variable_Type_Any) || (eqt.VariableType == vtmDir)) {
                    if (SpecFlag_NE[j] == 1) {
                        throw m1d_Error("SurBC_Dirichlet::domain_prep",
                                        "Error: multiple boundary conditions applied to the same variable " + ZZCantera::int2str(j));
                    }
                    SpecFlag_NE[j] = 1;
                    Value_NE[j] = sdt->Value[i];
                    TimeDep_NE[j] = sdt->TimeDep[i];
                    BC_TimeDep_NE[j] = sdt->BC_TimeDep_[i];
                    BC_Type_NE[j] = sdt->BC_Type_[i];
                    IsAlgebraic_Node[j] = 1;
                    NumBCs++;
                    ifound = true;
                }
            }
        } else {
            for (int j = 0; j < NumNodeEqns; j++) {
                VarType vtNode = NodalVarPtr->VariableNameList_EqnNum[j];
                if (vtDir == vtNode) {
                    if (SpecFlag_NE[j] == 1) {
                        throw m1d_Error("SurBC_Dirichlet::domain_prep",
                                        "Error: multiple boundary conditions applied to the same variable " + ZZCantera::int2str(j));
                    }
                    SpecFlag_NE[j] = 1;
                    Value_NE[j] = sdt->Value[i];
                    TimeDep_NE[j] = sdt->TimeDep[i];
                    BC_TimeDep_NE[j] = sdt->BC_TimeDep_[i];
                    BC_Type_NE[j] = sdt->BC_Type_[i];
                    IsAlgebraic_Node[j] = 1;
                    NumBCs++ ;
                    ifound = true;
                }
            }
        }
        if (!ifound) {
            std::string ss = vtDir.VariableName(24);
            printf("SurBC_Dirichlet::domain_prep: didn't find equation for variable %s\n", ss.c_str());

        }
    }

}
//=====================================================================================================================
int SurBC_Dirichlet::changeDirichletConditionValue(VarType vtDir, double newVal)
{
    int num = 0;
    for (int j = 0; j < NumNodeEqns; j++) {
        VarType vtNode = NodalVarPtr->VariableNameList_EqnNum[j];
        if (vtDir == vtNode) {
            if (SpecFlag_NE[j] != 1) {
                throw m1d_Error("SurBC_Dirichlet::domain_prep",
                                "Error: multiple boundary conditions applied to the same variable " + ZZCantera::int2str(j));
            }
            Value_NE[j] = newVal;
            num++;
        }
    }
    return num;
}
//=====================================================================================================================

// Change the boundary condition applied to a variable
/*
 *   @param vtDir    Variable type class. Note, general matches are allowed with this parameter
 *   @param BC_Type  Type of the boundary condition
 *   @param value    Value of the dirichlet condition or flux - default 0.0
 *   @param BC_TimeDep BoundaryCondition Pointers for time dependent BC for BC_Tppe = 3,4,5
 *                   (default 0)
 *   @param TimeDep  Function pointer to a function that returns a double given a single parameter (the time).
 *                   Defaults to a NULL pointer.
 *
 *   @return  Returns the number of boundary conditions matched.
 *            A negative number means that an error has been encountered
 */
int SurBC_Dirichlet::changeBoundaryCondition(VarType vtDir, int BC_Type, double value, BoundaryCondition* BC_TimeDep,
                                             TimeDepFunctionPtr TimeDep)
{
    int num = 0;
    for (int j = 0; j < NumNodeEqns; j++) {
        VarType vtNode = NodalVarPtr->VariableNameList_EqnNum[j];
        if (vtDir == vtNode) {
            num++;
            if (BC_Type < 0) {
                if (SpecFlag_NE[j] != 0) {
                    SpecFlag_NE[j] = 0;
                    IsAlgebraic_Node[j] = -1;
                    NumBCs--;
                }
            } else {
                if (SpecFlag_NE[j] != 1) {
                    SpecFlag_NE[j] = 1;
                    IsAlgebraic_Node[j] = 1;
                    NumBCs++;
                }
                BC_Type_NE[j] = BC_Type;
                Value_NE[j] = value;
                TimeDep_NE[j] = TimeDep;
                BC_TimeDep_NE[j] = BC_TimeDep;
            }
        }
    }
    return num;
}
//=====================================================================================================================
//   Report on the boundary condition applied on the first match to VarType
/*
 *   Inputs:
 *   @param time     Current time for evaluating time dependent BC
 *   @param vtDir    Variable type class. Note, general matches are allowed with this parameter.
 *   Outputs
 *   @param BC_Type  Type of the boundary condition
 *   @param value    Value of the dirichlet condition or flux - default 0.0
 *   @param BC_TimeDep BoundaryCondition Pointers for time dependent BC for BC_Tppe = 3,4,5
 *                   (default 0)
 *   @param TimeDep  Function pointer to a function that returns a double given a single parameter (the time).
 *                   Defaults to a NULL pointer.
 *
 *   @return  Returns the number of boundary conditions matched.
 *            A negative number means that an error has been encountered
 */
int SurBC_Dirichlet::reportBoundaryCondition(double time, const VarType vtDir, int& BC_Type, double& value,
                                             BoundaryCondition*& BC_TimeDep, TimeDepFunctionPtr& TimeDep) const
{
    int num = 0;
    /*
     * Find the current time region
     */
    DomainLayout* dl = SDD_.DL_ptr_;
    ProblemResidEval* pb = dl->problemResid_;
    int timeRegion = pb->m_currentTimeRegion;

    for (int j = 0; j < NumNodeEqns; j++) {
        VarType vtNode = NodalVarPtr->VariableNameList_EqnNum[j];
        if (vtDir == vtNode) {
            BC_Type = BC_Type_NE[j];
            value = Value_NE[j];
            if (BC_Type_NE[j] == 2 || BC_Type_NE[j] == 3) {
                value *= TimeDep_NE[j](time);
            } else if (BC_Type_NE[j] >= 4 && BC_Type_NE[j] <= 9) {
                value = BC_TimeDep_NE[j]->value(time, timeRegion);
            }

            TimeDep = TimeDep_NE[j];
            BC_TimeDep = BC_TimeDep_NE[j];
            num++;
            return num;
        }
    }
    return num;
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
void SurBC_Dirichlet::residEval(Epetra_Vector& res, const bool doTimeDependentResid, const Epetra_Vector* soln_ptr,
                                const Epetra_Vector* solnDot_ptr, const Epetra_Vector* solnOld_ptr, const double t,
                                const double rdelta_t, const ResidEval_Type_Enum residType, const Solve_Type_Enum solveType)
{
    int ieqn;
    const Epetra_Vector& soln = *soln_ptr;
    /*
     *  Quick return if we don't own the node that the boundary condition
     *  is to be applied on.
     */
    if (NumOwnedNodes == 0) {
        return;
    }

    /*
     * Pull in the flux vectors into this object
     */
    updateBulkFluxVectors();
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
     * Find the current time region
     */
    DomainLayout* dl = SDD_.DL_ptr_;
    ProblemResidEval* pb = dl->problemResid_;
    int timeRegion = pb->m_currentTimeRegion;

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

            ieqn = index_EqnStart + i;
            double solnVal = soln[ieqn];
            double val = Value_NE[i];
            switch (BC_Type_NE[i]) {
            case 0:
                /*
                 *  For Dirichlet equations, replace the equation
                 */
                res[ieqn] = val - solnVal;
                break;
            case 1:
                /*
                 *  For flux boundary conditions, subtract from equation indicating a flux out
                 */
                res[ieqn] -= val;
                // CAL: WARNING m1d_SurDomain_CathodeCollector.cpp has += val
                break;
            case 2:
                /*
                 *  For Dirichlet boundary conditions with oscillation, replace the equation
                 */
                res[ieqn] = val * TimeDep_NE[i](t) - solnVal;
                break;
            case 3:
                /*
                 *  For flux boundary conditions with oscillation, replace the equation
                 */
                // CAL: WARNING m1d_SurDomain_CathodeCollector.cpp has += val
                res[ieqn] -= val * TimeDep_NE[i](t);
                break;
            case 4: // voltage BCconstant
            case 6: // voltage BCsteptable
            case 8: // voltage BClineartable
                /*
                 *  For time dependent Dirichlet boundary condition using BoundaryCondition class
                 */
                res[ieqn] = BC_TimeDep_NE[i]->value(t, timeRegion) - solnVal;
                break;
            case 5: // current BCconstant
            case 7: // current BCsteptable
            case 9: // current BClineartable
                /*
                 *  For time dependent flux boundary condition using BoundaryCondition class
                 */
                // CAL: WARNING m1d_SurDomain_CathodeCollector.cpp has += val
                res[ieqn] -= BC_TimeDep_NE[i]->value(t, timeRegion);
                break;
            case 10: // cathode collector
                val = Value_NE[i];
                res[ieqn] -= BC_TimeDep_NE[i]->valueAtTime(t, val, timeRegion);
                break;
            default:
                throw m1d_Error("SDT_CathodeCollector::SetEquationDescription",
                                "voltageVarBCType not 0-9 for Dirichlet, Neumann, and Time Dependence");
            }

        }
    }
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
void SurBC_Dirichlet::saveDomain(ZZCantera::XML_Node& oNode, const Epetra_Vector* soln_GLALL_ptr,
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
        std::string nm = nv->VariableName(k);
        VarType vv = nv->VariableNameList_EqnNum[k];
        std::string type = VarType::VarMainName(vv.VariableType);
        ZZctml::addFloat(inlt, nm, sval, "", "", ZZCantera::Undef, ZZCantera::Undef);
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
void SurBC_Dirichlet::showSolution(const Epetra_Vector* soln_GlAll_ptr, const Epetra_Vector* solnDot_GlAll_ptr,
                                   const Epetra_Vector* soln_ptr, const Epetra_Vector* solnDot_ptr,
                                   const Epetra_Vector* solnOld_ptr, const Epetra_Vector_Owned* residInternal_ptr, const double t,
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

    BulkDomain1D* bd = 0;
    BulkDomainDescription* bedd = SDD_.LeftBulk;
    double* diffFlux = 0;
    if (!bedd) {
        bedd = SDD_.RightBulk;
        if (bedd) {
            bd = bedd->BulkDomainPtr_;
            diffFlux = &(bd->DiffFluxLeftBound_LastResid_NE[0]);
        }
    } else {
        bd = bedd->BulkDomainPtr_;
        diffFlux = &(bd->DiffFluxRightBound_LastResid_NE[0]);
    }

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
            if (vt.VariableType == Velocity_Axial) {
                size_t offset = NodalVarPtr->Offset_VarTypeVector[(size_t) Velocity_Axial];
                sval = diffFlux[offset];
                ss.print0("%s   %-20s * % -10.4E ", ind, name.c_str(), sval);
            } else {
                ss.print0("%s   %-20s   % -10.4E ", ind, name.c_str(), sval);
            }
            if (SpecFlag_NE[k] != 0) {
                ss.print0(" (Dir %d val = %-10.4E)", jDir, Value_NE[jDir]);
                jDir++;
            }
            ss.print0("\n");
        }
        drawline0(ss, indentSpaces + 2, 60);
        drawline0(ss, indentSpaces, 80);
    }
    print0_sync_end(0, ss, * (LI_ptr_->Comm_ptr_));
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
void SurBC_Dirichlet::initialConditions(const bool doTimeDependentResid, Epetra_Vector* soln_ptr,
                                        Epetra_Vector* solnDot,
                                        const double t, const double delta_t)
{
    int ieqn;
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
            if (BC_Type_NE[i] == 0) {
                soln[ieqn] = val;
            } else if (BC_Type_NE[i] == 2) {
                soln[ieqn] = val * TimeDep_NE[i](t);
            } else if (BC_Type_NE[i] == 4 || BC_Type_NE[i] == 6 || BC_Type_NE[i] == 8) {
                soln[ieqn] = BC_TimeDep_NE[i]->value(t, timeRegion);
            }
        }
    }
}

//=====================================================================================================================
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================
