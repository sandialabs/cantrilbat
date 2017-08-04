/**
 *  @file m1d_EpetraJac.cpp
 *      Implementation file for class EpetraJac
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#include "m1d_EpetraJac.h"

#include "m1d_GlobalIndices.h"
#include "m1d_LocalNodeIndices.h"
#include "m1d_VBRIndices.h"
#include "m1d_Comm.h"

#include "m1d_solvers.h"

#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"

#include "BlockEntryGlobal.h"

#include "cantera/base/utilities.h"

#include "Epetra_CrsMatrix.h"

using namespace std;
using namespace BEInput;
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
EpetraJac::EpetraJac(ProblemResidEval& r) :
    m_factored(false),
    m_columnScaled(false),
    m_rowScaled(false),
    m_n(0),
    m_kl(0),
    m_ku(0),
    A_(0),
    Arow_(0),
    Comm_ptr_(0),
    m_BmAX(0),
    m_AX(0),
    deltaSoln_(0),
    m_resid(&r),
    m_rtol(1.0E-3),
    m_atol(1.0E-9),
    m_elapsed(0.0),
    m_ssdiag(0),
    Acrs_(0),
    m_isAlgebraic(0),
    m_nevals(0),
    m_age(100000),
    m_NumGbEqns(0),
    solverType_(Direct),
    directSolverName_("Umfpack"),
    m_GbBlockRows(0),
    GI_ptr_(0),
    LI_ptr_(0),
    LRN_VBR_ptr_(0),
    solveType_(Zuzax::Solve_Type::SteadyState_Solve)
{
    m_NumGbEqns = r.NumGbEqns();
    m_GbBlockRows = r.NumGbNodes();


    GI_ptr_ = r.GI_ptr_;
    LI_ptr_ = r.LI_ptr_;
    Comm_ptr_ = LI_ptr_->Comm_ptr_;

    bool CopyMode_ = false;
    LRN_VBR_ptr_ = new m1d::LocalRowNodeVBRIndices(Comm_ptr_, CopyMode_, GI_ptr_, LI_ptr_);

    LRN_VBR_ptr_->informProblem(&r);

}
//==================================================================================================================================
EpetraJac::~EpetraJac()
{
    safeDelete(A_);
    safeDelete(Arow_);
    safeDelete(LRN_VBR_ptr_);
    safeDelete(m_BmAX);
    safeDelete(m_AX);
    safeDelete(deltaSoln_);
    safeDelete(m_isAlgebraic);
}
//==================================================================================================================================
EpetraJac::EpetraJac(const EpetraJac& r) :
    EpetraJac(*(r.m_resid))
{
    operator=(r);
}
//==================================================================================================================================
EpetraJac& EpetraJac::operator=(const EpetraJac& r)
{
    if (this == &r) {
        return *this;
    }

    m_factored     = r.m_factored;
    m_columnScaled = r.m_columnScaled;
    m_rowScaled    = r.m_rowScaled;
    m_n            = r.m_n;
    m_kl           = r.m_kl;
    m_ku           = r.m_ku;

    // fix
    safeDelete(A_);
    if (r.A_) {
        A_   =  new Epetra_VbrMatrix( *(r.A_) );
    }
    // fix
    safeDelete(Arow_);
    if (r.Arow_) {
         Arow_  = new Epetra_VbrRowMatrix( *(r.Arow_));
    }

    // Assuming this is a shallow pointer
    Comm_ptr_ = r.Comm_ptr_;

    safeDelete(m_BmAX);
    if (r.m_BmAX) {
        m_BmAX = new Epetra_Vector(*(r.m_BmAX));
    }
    safeDelete(m_AX);
    if (r.m_AX) {
        m_AX = new Epetra_Vector(*(r.m_AX));
    }
    safeDelete(deltaSoln_);
    if (r.deltaSoln_) {
        deltaSoln_ = new Epetra_Vector(*(r.deltaSoln_));
    }
    m_resid = r.m_resid;
    m_rtol = r.m_rtol;
    m_atol = r.m_atol;
    m_elapsed = r.m_elapsed;

    safeDelete(m_ssdiag);
    if (r.m_ssdiag) {
        m_ssdiag = new Epetra_Vector(*(r.m_ssdiag));
    }
    safeDelete(Acrs_);
    if (r.Acrs_) {
        Acrs_ = new Epetra_CrsMatrix(*(r.Acrs_));
    }
    safeDelete(m_isAlgebraic)
    if (r.m_isAlgebraic) {
        m_isAlgebraic = new Epetra_IntVector(*(r.m_isAlgebraic));
    }
    m_nevals = r.m_nevals;
    m_age = r.m_age;
    m_NumGbEqns = r.m_NumGbEqns;
    solverType_ = r.solverType_;
    directSolverName_ = r.directSolverName_;
    m_GbBlockRows = r.m_GbBlockRows;
    // fix
    GI_ptr_ = r.GI_ptr_;
    // fix
    LI_ptr_ = r.LI_ptr_;

    LRN_VBR_ptr_ = r.LRN_VBR_ptr_;
    varDiffString_LcNode = r.varDiffString_LcNode;
    solveType_ = r.solveType_;
    // error exit until tested
    exit(-1);

    return *this;
}
//==================================================================================================================================
void EpetraJac::allocateMatrix()
{
    A_ = LRN_VBR_ptr_->alloc_VbrMatrix();

    /*
     *  Add a row matrix representation. This is actually used to solve the matrix
     */
    Arow_ = LRN_VBR_ptr_->alloc_VbrRowMatrix(A_);

    m_BmAX = new Epetra_Vector(A_->RangeMap());
    m_AX = new Epetra_Vector(A_->RangeMap());

    /*
     * Output to understand what's going on. The RowMatrixRowMap, RowMatrixColMap,
     * OperatorRangeMap, OperatorDomainMap, and Map()
     * maps are all point maps. What point maps mean is the block size is equal to one. Therefore, they
     * are similar to what the Epetra_CsrMatrix would produce.
     *
     * All of the maps are the same except for RowMatrixColMap. They correspond to the concept of
     * Epetra_Vector_Owned. RowMatrixColMap map corresponds to the concept of EpetraVector_Ghosted.
     */
    /*
    Epetra_Map * m_ArowMRow = new Epetra_Map(Arow_->RowMatrixRowMap());
    print0_epBlockMap(*m_ArowMRow);

    Epetra_Map * m_ArowMCol = new Epetra_Map(Arow_->RowMatrixColMap());
    print0_epBlockMap(*m_ArowMCol);

    Epetra_Map * m_ArowORange = new Epetra_Map(Arow_->OperatorRangeMap());
    print0_epBlockMap(*m_ArowORange);

    Epetra_Map * m_ArowODomain = new Epetra_Map(Arow_->OperatorDomainMap());
    print0_epBlockMap(*m_ArowODomain);

    Epetra_BlockMap * m_ArowMap = new Epetra_BlockMap(Arow_->Map());
    print0_epBlockMap(*m_ArowMap);

    const Epetra_BlockMap &rMap = A_->RangeMap();
    print0_epBlockMap(rMap);

    */

    int numBlockCols = A_->NumMyBlockCols();
    varDiffString_LcNode.resize(numBlockCols, "");
    int num = LI_ptr_->NumLcNodes;

    if (num != numBlockCols) {
        printf("don't know what's going on%d %d \n", num, numBlockCols);
        exit(-1);
    }

    //deltaSoln_ = new Epetra_Vector(A_->RangeMap());
    //const Epetra_BlockMap &rMap = A_->RangeMap();
    //print0_epBlockMap(rMap);
    //const Epetra_BlockMap &mMap = A_->RowMatrixColMap();
    //print0_epBlockMap(mMap);

    //deltaSoln_ = new Epetra_Vector(A_->RangeMap());
    deltaSoln_ = new Epetra_Vector(*(LI_ptr_->GbBlockNodeEqnstoLcBlockNodeEqnsColMap));


    m_isAlgebraic = new Epetra_IntVector(A_->RangeMap());
    m_resid->fillIsAlgebraic(*m_isAlgebraic);
}
//==================================================================================================================================
void EpetraJac::updateTransient(double rdt, int* mask)
{
    int n;
    for (n = 0; n < m_NumGbEqns; n++) {
        //  value(n, n) = m_ssdiag[n] - mask[n] * rdt;
    }
}
//==================================================================================================================================
void EpetraJac::incrementDiagonal(int j, double d)
{
    // m_ssdiag[j] += d;
    // value(j, j) = m_ssdiag[j];
}
//==================================================================================================================================
RecordTree_base*
EpetraJac::setupMDinput_pass1(BEInput::BlockEntry* Parent, RecordTree_base* dbb)
{
    BEinput_EpetraJac* db = 0;
    if (dbb) {
        db = dynamic_cast<BEinput_EpetraJac*>(dbb);
        if (!db) {
            exit(-1);
        }
    }

    // Find the block relevant to this member
    BEInput::BlockEntry* lsb = new  BEInput::BlockEntry("Linear Solver", 0);
    Parent->addSubBlock(lsb);

    if (!db) {
        db = new BEinput_EpetraJac();
    }
    if (lsb) {
        const char* cppkkm[2] = {"Direct", "Iterative"};
        LE_PickList* lepkm = new LE_PickList("Solver Type", (int*) &(db->I_solverType), cppkkm, 2, 1, "SolverType");
        lepkm->set_default("Direct");
        lsb->addLineEntry(lepkm);

        LE_OneStr* ledsn = new LE_OneStr("Direct Solver Name", &(db->S_directSolverName), 2, 1, 0, "DirectSolverName");
        ledsn->set_default("Umfpack");
        lsb->addLineEntry(ledsn);
    }
    return static_cast<m1d::RecordTree_base*>(db);
    //return db;
}
//==================================================================================================================================
Epetra_IntVector& EpetraJac::transientMask()
{
    return *m_isAlgebraic;
}
//==================================================================================================================================
void EpetraJac::matrixEval(const bool doTimeDependentResid, const Epetra_Vector* const solnBase_ptr,
                           const Epetra_Vector* const solnDotBase_ptr, const double t, const double rdelta_t,
                           const Zuzax::Solve_Type solveType)
{
    // Epetra_Vector *resBase = new Epetra_Vector(*(GI_ptr_->GbEqnstoOwnedLcEqnsMap));
    Epetra_Vector* resBase = new Epetra_Vector(*(LI_ptr_->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap));
    matrixEval1(doTimeDependentResid, solnBase_ptr, solnDotBase_ptr, resBase, t, rdelta_t, solveType);
    fillVbr();
    safeDelete(resBase);
}
//===================================================================================================================================
void EpetraJac::matrixResEval(const bool doTimeDependentResid, const Epetra_Vector* const solnBase_ptr,
                              const Epetra_Vector* const solnDotBase_ptr, Epetra_Vector* const resBase,
                              const double t, const double rdelta_t, Zuzax::Solve_Type solveType)
{
    matrixEval1(doTimeDependentResid, solnBase_ptr, solnDotBase_ptr, resBase, t, rdelta_t, solveType);
    fillVbr();
}
//==================================================================================================================================
/*
 *   solnBase contains ghost nodes
 */
void EpetraJac::matrixEval1(const bool doTimeDependentResid, const Epetra_Vector* const solnBase_ptr,
                            const Epetra_Vector* const solnDotBase_ptr, Epetra_Vector* const resBase,
                            const double t, const double rdelta_t, Zuzax::Solve_Type solveType)
{
    solveType_ = solveType;
    m_resid->residEval(resBase, doTimeDependentResid, solnBase_ptr, solnDotBase_ptr, t, rdelta_t, 
                       Zuzax::ResidEval_Type::Base_ResidEval, solveType);

    eval(doTimeDependentResid, solnBase_ptr, solnDotBase_ptr, *resBase, t, rdelta_t);
}
//===================================================================================================================================
void EpetraJac::eval(const bool doTimeDependentResid, const Epetra_Vector* const solnBase_ptr,
                     const Epetra_Vector* const solnDotBase_ptr, const Epetra_Vector& resBase,
                     double t, double rdelta_t)
{
    m_nevals++;
    double timeDim = 0.0;
    if (solveType_ == Zuzax::Solve_Type::DAESystemInitial_Solve) {
        timeDim = 1.0E-3;
    }
    //
    // Calculate the vector of deltas used in the numerical jacobian
    //
    m_resid->calcDeltaSolnVariables(t, *solnBase_ptr, solnDotBase_ptr, *deltaSoln_, solveType_);

    // To print out the solution deltas:
#ifdef DEBUG_JACOBIAN_CALCULATION
    string ss = "DeltaSolutionVars";
    m_resid->showSolutionVector(ss, t, rdelta_t, *deltaSoln_);
#endif


    // Set the age to zero
    m_age = 0;

    int numLcRowNodes = LRN_VBR_ptr_->NumLcRowNodes;
    const Epetra_Vector& solnBase = *solnBase_ptr;
    //Epetra_Vector *res = new Epetra_Vector(*(GI_ptr_->GbEqnstoOwnedLcEqnsMap));

    Epetra_Vector* res = new Epetra_Vector(*(LI_ptr_->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap));


    int numEqnColors = LI_ptr_->NumEqnColors;
    Epetra_Vector soln(*solnBase_ptr);
    Epetra_Vector eChange(*solnBase_ptr);
    //Epetra_Vector delta(*solnBase_ptr);

    Epetra_Vector* solnDot_ptr = 0;
    if (solnDotBase_ptr) {
        solnDot_ptr = new Epetra_Vector(*solnDotBase_ptr);
    }

    Epetra_Vector jacColn(*solnBase_ptr);
    jacColn.PutScalar(0.0);

    zeroMatrix();


    for (int iColor = 0; iColor < numEqnColors; iColor++) {

        for (int iLcEqn = 0; iLcEqn < LI_ptr_->NumLcEqns; iLcEqn++) {
            //  int jColor = (*(LI_ptr_->EqnColorMap))[iLcEqn];
            int jColor = (*(LI_ptr_->EqnColors))[iLcEqn];
            eChange[iLcEqn] = 0.0;
            //delta[iLcEqn] = 0.0;
            if (iColor == jColor) {
                //soln[iLcEqn] = m_resid->deltaSolnCompJac(*solnBase_ptr, iLcEqn);
                //delta[iLcEqn] = soln[iLcEqn] - solnBase[iLcEqn];
                soln[iLcEqn] += (*deltaSoln_)[iLcEqn];
                if (solveType_ != Zuzax::Solve_Type::DAESystemInitial_Solve || (*m_isAlgebraic)[iLcEqn] == 1) {
                    if (solnDotBase_ptr) {
                        //(*solnDot_ptr)[iLcEqn] += delta[iLcEqn] * rdelta_t;
                        (*solnDot_ptr)[iLcEqn] += (*deltaSoln_)[iLcEqn] * rdelta_t;
                    }
                } else {
                    soln[iLcEqn] = solnBase[iLcEqn];
                    //delta[iLcEqn] /=  timeDim;
                    (*deltaSoln_)[iLcEqn] /= timeDim;
                    // (*solnDot_ptr)[iLcEqn] += delta[iLcEqn];
                    (*solnDot_ptr)[iLcEqn] += (*deltaSoln_)[iLcEqn];
                }
                eChange[iLcEqn] = 1.0;
            }
        }

        m_resid->residEval(res, doTimeDependentResid, &soln, solnDot_ptr, t, rdelta_t, 
                           Zuzax::ResidEval_Type::JacDelta_ResidEval, solveType_);

        // loop over the block rows
        for (int iBlockRow = 0; iBlockRow < numLcRowNodes; iBlockRow++) {
            // Find out the global node number
            //int gbNode = m1d::LI_ptr->IndexGbNode_LcNode[iBlockRow];

            // Find out the number of block columns corresponding to this block row
            int numBlocks = LRN_VBR_ptr_->NumColBlocks_LcRNodes[iBlockRow];

            // Find out the block column indices for this block row
            /*
             *  bIndexes[j] is the local node number of the block column
             */
            int* bIndexes = LRN_VBR_ptr_->IndexLcNodeColBlock_LcRNodes[iBlockRow];
            //int *colSize_row = m1d::LRN_VBR_ptr->ColSizeColBlock_LcRNodes[iBlockRow];

            // Pull down the vector of dense matrices corresponding to the block row
            // Jacobian entries
            Epetra_SerialDenseMatrix** rowBlock = LRN_VBR_ptr_->BlockMatrices[iBlockRow];

            // Get the column size of each block column in the current block row.
            int* ColSizeColBlock = LRN_VBR_ptr_->ColSizeColBlock_LcRNodes[iBlockRow];

            // Loop over the block columns. At this point at most one of the block
            // column variables has changed.
            int ifound = -1;
            for (int jColIndex = 0; jColIndex < numBlocks; jColIndex++) {
                // Find the block row number
                int blockColLcNode = bIndexes[jColIndex];
                // Find the Jacobian block for this interaction
                Epetra_SerialDenseMatrix* rowColBlock = rowBlock[jColIndex];

                // Get the column size. This should be equal to the number of equations
                // located at the column node - Check it
                int colSize = ColSizeColBlock[jColIndex];
                AssertTrace(LI_ptr_->NumEqns_LcNode[blockColLcNode] == colSize);

                // Loop over the block column equations
                for (int je = 0; je < colSize; je++) {
                    // Get the local equation number value
                    int jLcEqn = je + LI_ptr_->IndexLcEqns_LcNode[blockColLcNode];
                    // discover which one changed
                    if (eChange[jLcEqn]) {
                        // Check to make sure that no other equation changed
                        if (ifound != -1) {
                            AssertThrow("Logic Error","");
                        }
                        ifound = jLcEqn;

                        //Write the block row to the block column.
                        // Get the delta
                        //double dd = 1.0 / delta[jLcEqn];
                        double dd = 1.0 / (*deltaSoln_)[jLcEqn];
                        // get the number of rows in the block
                        int numRowEqns = LRN_VBR_ptr_->RowSize_LcRNodes[iBlockRow];
                        int istart = LI_ptr_->IndexLcEqns_LcNode[iBlockRow];
                        for (int ie = 0; ie < numRowEqns; ie++) {
                            int iLcEqn = istart + ie;
                            double value = ((*res)[iLcEqn] - resBase[iLcEqn]) * dd;
#ifdef DEBUG_MODE
                            ZZCantera::checkFinite(value);
#endif
                            (*rowColBlock)(ie, je) = value;
                        }
                    }
                }
            }
        }
        for (int iLcEqn = 0; iLcEqn < LI_ptr_->NumLcEqns; iLcEqn++) {
            soln[iLcEqn] = solnBase[iLcEqn];
        }
        if (solnDotBase_ptr) {
            for (int iLcEqn = 0; iLcEqn < LI_ptr_->NumLcEqns; iLcEqn++) {
                (*solnDot_ptr)[iLcEqn] = (*solnDotBase_ptr)[iLcEqn];
            }
        }

    }
    m_factored = false;
    m_columnScaled = false;
    m_rowScaled = false;
    safeDelete(res);
    safeDelete(solnDot_ptr);
}
//==================================================================================================================================
void EpetraJac::fillVbr()
{
    bool copyMode = LRN_VBR_ptr_->CopyMode;
    if (!copyMode) {
        return;
    }
    int numLcRowNodes = LRN_VBR_ptr_->NumLcRowNodes;

    for (int iBlockRow = 0; iBlockRow < numLcRowNodes; iBlockRow++) {
        //int GbNode = m1d::LI_ptr->IndexGbNode_LcNode[iBlockRow];

        int numBlocks = LRN_VBR_ptr_->NumColBlocks_LcRNodes[iBlockRow];

        int* bIndexes = LRN_VBR_ptr_->IndexLcNodeColBlock_LcRNodes[iBlockRow];
        //int *colSize_row = m1d::LRN_VBR_ptr->ColSizeColBlock_LcRNodes[iBlockRow];
        Epetra_SerialDenseMatrix** rowBlock = LRN_VBR_ptr_->BlockMatrices[iBlockRow];

        if (copyMode) {
            A_->BeginReplaceMyValues(iBlockRow, numBlocks, bIndexes);
        }

        for (int jColIndex = 0; jColIndex < numBlocks; jColIndex++) {
            // int blockColIndex = bIndexes[jColIndex];

            Epetra_SerialDenseMatrix* rowColBlock = rowBlock[jColIndex];
            A_->SubmitBlockEntry(*rowColBlock);

        }

        if (copyMode) {
            A_->EndSubmitEntries();
        }
    }
}
//==================================================================================================================================
void EpetraJac::getRowScales(Epetra_Vector* const rowScales) const
{
    LRN_VBR_ptr_->calcInvRowSumScales(*rowScales);
}
//==================================================================================================================================
void EpetraJac::columnScale(const Epetra_Vector* const colScales)
{
    if (m_factored) {
        throw m1d_Error("EpetraJac::columnScale", "matrix is factored");
    }
    if (m_rowScaled) {
        throw m1d_Error("EpetraJac::columnScale", "matrix is rowScaled");
    }
    A_->RightScale(*colScales);
    m_columnScaled = true;
}
//==================================================================================================================================
void EpetraJac::rowScale(const Epetra_Vector* const rowScales)
{
    if (m_factored) {
        throw m1d_Error("EpetraJac::rowScale", "matrix is factored");
    }
    // multiply each row of the matrix by a single scale.
    A_->LeftScale(*rowScales);
    m_rowScaled = true;
}
//==================================================================================================================================
void EpetraJac::queryMatrixStructure(std::ostream& oo)
{
    LRN_VBR_ptr_->queryMatrixStructure(A_, oo);
}
//===================================================================================================================================
const Epetra_Vector& EpetraJac::deltaSolnJac() const
{
    return *deltaSoln_;
}
//===================================================================================================================================
void EpetraJac::zeroMatrix()
{
    // loop over the block rows
    int numLcRowNodes = LRN_VBR_ptr_->NumLcRowNodes;
    for (int iBlockRow = 0; iBlockRow < numLcRowNodes; iBlockRow++) {

        // Find out the global node number
        //int gbNode = m1d::LI_ptr->IndexGbNode_LcNode[iBlockRow];

        // Find out the number of block columns corresponding to this block row
        int numBlocks = LRN_VBR_ptr_->NumColBlocks_LcRNodes[iBlockRow];

        // Find out the block column indices for this block row
        /*
         *  bIndexes[j] is the local node number of the block column
         */
        // int *bIndexes = m1d::LRN_VBR_ptr->IndexLcNodeColBlock_LcRNodes[iBlockRow];
        //int *colSize_row = m1d::LRN_VBR_ptr->ColSizeColBlock_LcRNodes[iBlockRow];

        // Pull down the vector of dense matrices corresponding to the block row
        // Jacobian entries
        Epetra_SerialDenseMatrix** rowBlock = LRN_VBR_ptr_->BlockMatrices[iBlockRow];
        // Get the column size of each block column in the current block row.
        int* ColSizeColBlock = LRN_VBR_ptr_->ColSizeColBlock_LcRNodes[iBlockRow];
        int numRowEqns = LRN_VBR_ptr_->RowSize_LcRNodes[iBlockRow];
        // Loop over the block columns. At this point at most one of the block
        // column variables has changed.
        for (int jColIndex = 0; jColIndex < numBlocks; jColIndex++) {
            // Get the column size. This should be equal to the number of equations
            // located at the column node - Check it
            int colSize = ColSizeColBlock[jColIndex];
            // Find the Jacobian block for this interaction
            Epetra_SerialDenseMatrix* rowColBlock = rowBlock[jColIndex];
            // AssertTrace(LI_ptr->NumEqns_LcNode[blockColLcNode] == colSize);

            // Loop over the block column equations
            for (int je = 0; je < colSize; je++) {
                for (int ie = 0; ie < numRowEqns; ie++) {
                    (*rowColBlock)(ie, je) = 0.0;
                }
            }
        }
    }
    m_factored = false;
    m_columnScaled = false;
    m_rowScaled = false;
}
//==================================================================================================================================
/*
 *   Return a changeable pointer into the matrix given Global Block Row indices
 *   I don't think this will ever be used in practice, because of the speed. But, here it is.
 *   We will only return a valid pointer if the block row and local equation number is owned by this processor. If it is
 *   not owned or is an external row on this processor will will return the null pointer.
 *
 *  @param gbRow          Global row number
 *  @param lcRowIndex     local row index of the equation on that row
 *  @param gbCol          Global col number
 *  @param lcColIndex     local col index of the variable on that col
 *  @return    Returns a pointer to the position in the block matrix corresponding
 *            to that row and column. Will return 0 if this number is not on the current
 *            processor. Will return 0 if this number isn't in the sparce matrix stencil.
 */
double* EpetraJac::GbBlkValue(int gbRow, int lcRowIndex, int gbCol, int lcColIndex) const
{
    double* vptr = nullptr;

    // Find the local block row
    int iBlockRow = LI_ptr_->GbNodeToLcNode(gbRow);
    // First check to see that the row is defined on this machine. If it is not, we will return nullptr.
    if (iBlockRow >= 0) {
        // Next check to see that the block row is an external row for this processor. If it is, we will return a nullptr
        if (! LI_ptr_->IsExternal_LcNode[iBlockRow]) {

            // find the local column node index on the processor given the global column node, which corresponds to the global node number
            int lccol = LI_ptr_->GbNodeToLcNode(gbCol);
            // Note this should always work. So should add an error
            if (lccol >= 0) {

                // Get the number of blocks columns in the block row
                int numBlocks = LRN_VBR_ptr_->NumColBlocks_LcRNodes[iBlockRow];

                // Get a pointer for the block column indices for this block row
                /*
                 *  bIndexes[j] is the starting local node number of the block column
                 */
                int* bIndexes = LRN_VBR_ptr_->IndexLcNodeColBlock_LcRNodes[iBlockRow];
                //int *colSize_row = m1d::LRN_VBR_ptr->ColSizeColBlock_LcRNodes[iBlockRow];

                // Pull down the vector of dense matrices corresponding to the block row Jacobian entries
                Epetra_SerialDenseMatrix** rowBlock = LRN_VBR_ptr_->BlockMatrices[iBlockRow];

                // Loop over the block columns. At this point at most one of the block column variables has changed.
                int ifound = -1;
                for (int jColIndex = 0; jColIndex < numBlocks; jColIndex++) {
                    // Find the block row number
                    int blockColLcNode = bIndexes[jColIndex];
                    if (blockColLcNode == lccol) {
                        ifound = lccol;
                        break;
                    }
                }
                if (ifound == -1) {
                    throw m1d_Error("EpetraJac::GbBlkValue", "indexing error");
                }

#ifdef DEBUG_MODE
                // Get the column size of each block column in the current block row.
                int* ColSizeColBlock = LRN_VBR_ptr_->ColSizeColBlock_LcRNodes[iBlockRow];
                if (lcColIndex >= ColSizeColBlock[ifound] || lcColIndex < 0) {
                    throw m1d_Error("EpetraJac::GbBlkValue", "lcol indexing error");
                }
                if (lcColIndex >= ColSizeColBlock[ifound] || lcColIndex < 0) {
                    throw m1d_Error("EpetraJac::GbBlkValue", "lcol indexing error");
                }
                if (lcRowIndex >= LI_ptr_->NumEqns_LcNode[iBlockRow] || lcRowIndex < 0) {
                    throw m1d_Error("EpetraJac::GbBlkValue", "lrow indexing error");
                }
#endif
                // Find the Jacobian block for this interaction
                Epetra_SerialDenseMatrix* rowColBlock = rowBlock[ifound];
                vptr = &(*rowColBlock)(lcRowIndex, lcColIndex);
            }
#ifdef DEBUG_MODE
            else {
                throw m1d_Error("EpetraJac::GbBlkValue", "lcol indexing error");
            }
#endif
        }

    }
    return vptr;
}
//==================================================================================================================================
double& EpetraJac::value(const int iGlobalEqn, const int jGlobalEqn)
{
    int lcRowIndex, lcColIndex;
    int igbRow = GI_ptr_->GbEqnToGbNode(iGlobalEqn, lcRowIndex);
    int igbCol = GI_ptr_->GbEqnToGbNode(iGlobalEqn, lcColIndex);
#ifdef DEBUG_MODE
    double* pos = GbBlkValue(igbRow, lcRowIndex, igbCol, lcColIndex);
    if (pos == nullptr) {
        throw m1d_Error("EpetraJac::value()", "bad indecises");
    }
    return *pos;
#else
    return GbBlkValue(igbRow, lcRowIndex, igbCol, lcColIndex);
#endif
}
//==================================================================================================================================
double& EpetraJac::operator()(const int iGlobalEqn, const int jGlobalEqn)
{
    return value(iGlobalEqn, jGlobalEqn);
}
//==================================================================================================================================
double EpetraJac::value(const int iGlobalEqn, const int jGlobalEqn) const
{
    int lcRowIndex, lcColIndex;
    int igbRow = GI_ptr_->GbEqnToGbNode(iGlobalEqn, lcRowIndex);
    int igbCol = GI_ptr_->GbEqnToGbNode(iGlobalEqn, lcColIndex);
    const double* pos = GbBlkValue(igbRow, lcRowIndex, igbCol, lcColIndex);
    if (pos == nullptr) {
        return 0.0;
    }
    return *pos;
}
//==================================================================================================================================
int EpetraJac::nRows() const
{
    return GI_ptr_->NumGbEqns;
}
//==================================================================================================================================
/*
 * Multiply A*b and write result to \c prod.
 */
void EpetraJac::mult(const Epetra_Vector& b, Epetra_Vector& prod) const
{
    int err = A_->Multiply1(false, b, prod);
    if (err) {
        throw m1d_Error("EpetraJac::mult", "error code returned");
    }
}
//==================================================================================================================================
/*
 * Multiply b*A and write result to \c prod.
 */
void EpetraJac::leftMult(const Epetra_Vector& b, Epetra_Vector& prod) const
{
    int err = A_->Multiply1(true, b, prod);
    if (err) {
        throw m1d_Error("EpetraJac::mult", "error code returned");
    }
}
//==================================================================================================================================
int EpetraJac::factor()
{
    if (solverType_ == Direct) {

    } else if (solverType_ == Iterative) {

    }
    m_factored = true;
    return 0;
}
//==================================================================================================================================
int EpetraJac::solve(Epetra_Vector* const b, Epetra_Vector* const x, int& its, double& norm, bool doRes)
{
    int retn = 0;
    double residualL2;
    norm = 0.0;
    /*
     * Option Here to print out A_. It may be scaled
     */

    /*
     *  Option to convert to Csr matrix
     */

    /*
     * Set optional parameters on each of these solves
     */
    if (solverType_ == Direct) {
        retn = solve_amesos(Arow_, x, b, this);
        its = 1;
    } else if (solverType_ == Iterative) {
        retn = solve_aztecoo(Arow_, x, b, this);
    }
    /*
     * Find the norm of ||Ax-b|| if asked to do so
     */
    if (doRes) {
        A_->Multiply(false, *x, *m_AX);
        (*m_AX).Update(1.0, *b, -1.0);
        (*m_AX).Norm2(&residualL2);
        norm = residualL2;
    }
    m_factored = true;
    return retn;
}
//==================================================================================================================================
double EpetraJac::elapsedTime() const
{
    return m_elapsed;
}
//==================================================================================================================================
int EpetraJac::nEvals() const
{
    return m_nevals;
}
//==================================================================================================================================
int EpetraJac::age() const
{
    return m_age;
}
//==================================================================================================================================
void EpetraJac::setAge(int age)
{
    m_age = age;
}
//==================================================================================================================================
void EpetraJac::incrementAge()
{
    m_age++;
}
//==================================================================================================================================
int EpetraJac::nColumns() const
{
    return m_n;
}
//==================================================================================================================================
int EpetraJac::nSubDiagonals() const
{
    return m_kl;
}
//==================================================================================================================================
int EpetraJac::nSuperDiagonals() const
{
    return m_ku;
}
//==================================================================================================================================
int EpetraJac::ldim() const
{
    return 2 * m_kl + m_ku + 1;
}
//==================================================================================================================================
void EpetraJac::process_BEinput(const RecordTree_base* const dbb)
{
    const BEinput_EpetraJac* const db = dynamic_cast<const BEinput_EpetraJac*>(dbb);
    if (db) {
        solverType_ = db->I_solverType;
        directSolverName_ = db->S_directSolverName;
    }
}
//==================================================================================================================================
void EpetraJac::process_input(BlockEntry* const cf)
{
    // Find the block
    BlockEntry* sb = cf->match_block("Linear Solver");
    if (sb) {
        // Process the line entries and fill in the record
        LineEntry* le = sb->searchLineEntry("Solver Type");
        LE_PickList* lep = dynamic_cast<LE_PickList*>(le);
        solverType_ = (m1d::SolverType) lep->currentTypedValue();

        le = sb->searchLineEntry("Direct Solver Name");
        LE_OneStr* ledsn = dynamic_cast<LE_OneStr*>(le);
        directSolverName_ = ledsn->currentTypedValue();
    }
}
//==================================================================================================================================
ostream& operator<<(ostream& os, const EpetraJac& m)
{
    Epetra_VbrMatrix& mat = *m.A_;
    const Epetra_BlockMap& rMap = mat.RowMap();
    int MyPID = rMap.Comm().MyPID();

    if (MyPID == 0) {
        os << "\nNumber of Global Block Rows  = ";
        os << mat.NumGlobalBlockRows();
        os << endl;
        os << "Number of Global Block Cols  = ";
        os << mat.NumGlobalBlockCols();
        os << endl;
        os << "Number of Global Block Diags = ";
        os << mat.NumGlobalBlockDiagonals();
        os << endl;
        os << "Number of Global Blk Entries = ";
        os << mat.NumGlobalBlockEntries();
        os << endl;
        os << "Global Max Num Block Entries = ";
        os << mat.GlobalMaxNumBlockEntries();
        os << endl;
        os << "\nNumber of Global Rows        = ";
        os << mat.NumGlobalRows();
        os << endl;
        os << "Number of Global Cols        = ";
        os << mat.NumGlobalCols();
        os << endl;
        os << "Number of Global Diagonals   = ";
        os << mat.NumGlobalDiagonals();
        os << endl;
        os << "Number of Global Nonzeros    = ";
        os << mat.NumGlobalNonzeros();
        os << endl;
        os << "Global Maximum Num Entries   = ";
        os << mat.GlobalMaxNumNonzeros();
        os << endl;
        if (mat.LowerTriangular()) {
            os << " ** Matrix is Lower Triangular **";
        }
        os << endl;
        if (mat.UpperTriangular()) {
            os << " ** Matrix is Upper Triangular **";
        }
        os << endl;
        if (mat.NoDiagonal()) {
            os << " ** Matrix has no diagonal     **";
        }
        os << endl;
        os << endl;
    }

    os << "\nNumber of My Block Rows  = ";
    os << mat.NumMyBlockRows();
    os << endl;
    os << "Number of My Block Cols  = ";
    os << mat.NumMyBlockCols();
    os << endl;
    os << "Number of My Block Diags = ";
    os << mat.NumMyBlockDiagonals();
    os << endl;
    os << "Number of My Blk Entries = ";
    os << mat.NumMyBlockEntries();
    os << endl;
    os << "My Max Num Block Entries = ";
    os << mat.MaxNumBlockEntries();
    os << endl;
    os << "\nNumber of My Rows        = ";
    os << mat.NumMyRows();
    os << endl;
    os << "Number of My Cols        = ";
    os << mat.NumMyCols();
    os << endl;
    os << "Number of My Diagonals   = ";
    os << mat.NumMyDiagonals();
    os << endl;
    os << "Number of My Nonzeros    = ";
    os << mat.NumMyNonzeros();
    os << endl;
    os << "My Maximum Num Entries   = ";
    os << mat.MaxNumBlockEntries();
    os << endl;
    os << endl;

    os << flush;
    int NumBlockRows1 = mat.NumMyBlockRows();
    int MaxNumBlockEntries1 = mat.MaxNumBlockEntries();
    int* BlockIndices1 = new int[MaxNumBlockEntries1];
    Epetra_SerialDenseMatrix** Entries1;
    int RowDim1, NumBlockEntries1;
    int i, j;

    os.width(8);
    os << "   Processor ";
    os.width(10);
    os << "   Block Row Index ";
    os.width(10);
    os << "   Block Col Index";
    os.width(10);
    os << "   Values     ";
    os << endl;

    for (i = 0; i < NumBlockRows1; i++) {
        int BlockRow1 = mat.GRID(i); // Get global row number
        mat.ExtractGlobalBlockRowPointers(BlockRow1, MaxNumBlockEntries1, RowDim1, NumBlockEntries1, BlockIndices1, Entries1);

        for (j = 0; j < NumBlockEntries1; j++) {
            os.width(8);
            os << MyPID;
            os << "    ";
            os.width(10);
            os << BlockRow1;
            os << "    ";
            os.width(10);
            os << BlockIndices1[j];
            os << "    ";
            os.width(20);

            if (Entries1[j] == 0) {
                os << "Block Entry == NULL" << endl;
                continue;
            }

            Epetra_SerialDenseMatrix entry(View, Entries1[j]->A(), Entries1[j]->LDA(), RowDim1, Entries1[j]->N());
            int nrow = RowDim1;
            int ncol = entry.N();
            for (int i = 0; i < nrow; i++) {
                if (i != 0) {
                    os << "                                                  ";
                }
                for (int j = 0; j < ncol; j++) {
                    os.width(12);
                    os << entry[j][i] << "  ";
                }
                os << endl;
            }

        }
        os << endl;
    }
    /*
     * Cleanup section
     */
    delete[] BlockIndices1;

    return os;
}
//==================================================================================================================================
EpetraJac::BEinput_EpetraJac::BEinput_EpetraJac() :
    RecordTree_base(),
    I_solverType(Direct),
    S_directSolverName("Umfpack")
{
}
//==================================================================================================================================
EpetraJac::BEinput_EpetraJac::~BEinput_EpetraJac()
{
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
