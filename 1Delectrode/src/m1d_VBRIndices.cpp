/**
 * @file m1d_VBRIndices.cpp
 *
 */

/*
 *  $Id: m1d_VBRIndices.cpp 504 2013-01-07 22:32:48Z hkmoffa $
 */
#include "m1d_VBRIndices.h"
#include "m1d_GlobalIndices.h"
#include "m1d_ProblemResidEval.h"
#include "m1d_Comm.h"
#include "m1d_LocalNodeIndices.h"

using namespace std;

namespace m1d
{

//LocalRowNodeVBRIndices *LRN_VBR_ptr = 0;

//==========================================================================================================================
LocalRowNodeVBRIndices::LocalRowNodeVBRIndices(Epetra_Comm* comm_ptr,
                                               bool copyMode,
                                               GlobalIndices* gi_ptr,
                                               LocalNodeIndices* li_ptr) :
    Comm_ptr_(comm_ptr), CopyMode(copyMode), NumLcRowNodes(-1), GI_ptr_(gi_ptr), LI_ptr_(li_ptr),
    Func_ptr_(0)
{
}
//=======================================================================================================================
LocalRowNodeVBRIndices::~LocalRowNodeVBRIndices()
{
    freeBlockMatrices();
    for (int iBlockRow = 0; iBlockRow < NumLcRowNodes; iBlockRow++) {
        if (IndexLcNodeColBlock_LcRNodes[iBlockRow] != 0) {
            delete[] ColSizeColBlock_LcRNodes[iBlockRow];
        }
        if (IndexLcNodeColBlock_LcRNodes[iBlockRow] != 0) {
            delete[] IndexLcNodeColBlock_LcRNodes[iBlockRow];
        }
        if (Copy_IndexLcNodeColBlock_LcRNodes[iBlockRow] != 0) {
            delete[] Copy_IndexLcNodeColBlock_LcRNodes[iBlockRow];
        }
    }
}
//===========================================================================
LocalRowNodeVBRIndices::LocalRowNodeVBRIndices(const LocalRowNodeVBRIndices& r) :
    Comm_ptr_(r.Comm_ptr_), CopyMode(r.CopyMode), NumLcRowNodes(0), GI_ptr_(r.GI_ptr_), LI_ptr_(r.LI_ptr_)
{
    *this = r;
}
//===========================================================================
LocalRowNodeVBRIndices&
LocalRowNodeVBRIndices::operator=(const LocalRowNodeVBRIndices& r)
{
    if (this == &r) {
        return *this;
    }

    Comm_ptr_ = r.Comm_ptr_;
    CopyMode = r.CopyMode;
    /*
     * First delete the malloced memory. Note, this is duplicate code
     * wrt the destructor, so we may consolidate it later.
     */
    for (int iBlockRow = 0; iBlockRow < NumLcRowNodes; iBlockRow++) {
        if (IndexLcNodeColBlock_LcRNodes[iBlockRow] != 0) {
            delete[] ColSizeColBlock_LcRNodes[iBlockRow];
        }
        if (IndexLcNodeColBlock_LcRNodes[iBlockRow] != 0) {
            delete[] IndexLcNodeColBlock_LcRNodes[iBlockRow];
        }
        if (Copy_IndexLcNodeColBlock_LcRNodes[iBlockRow] != 0) {
            delete[] Copy_IndexLcNodeColBlock_LcRNodes[iBlockRow];
        }
    }

    // Ok this is wrong, will fix it later!
    NumLcRowNodes = r.NumLcRowNodes;
    NumColBlocks_LcRNodes = r.NumColBlocks_LcRNodes;
    RowSize_LcRNodes = r.RowSize_LcRNodes;

    IndexLcNodeColBlock_LcRNodes.resize(NumLcRowNodes);
    Copy_IndexLcNodeColBlock_LcRNodes.resize(NumLcRowNodes);
    ColSizeColBlock_LcRNodes.resize(NumLcRowNodes);
    for (int i = 0; i < NumLcRowNodes; i++) {
        IndexLcNodeColBlock_LcRNodes[i] = new int[NumColBlocks_LcRNodes[i]];
        Copy_IndexLcNodeColBlock_LcRNodes[i] = new int[NumColBlocks_LcRNodes[i]];
        ColSizeColBlock_LcRNodes[i] = new int[NumColBlocks_LcRNodes[i]];
        for (int j = 0; j < NumColBlocks_LcRNodes[i]; j++) {
            IndexLcNodeColBlock_LcRNodes[i][j] = r.IndexLcNodeColBlock_LcRNodes[i][j];
            Copy_IndexLcNodeColBlock_LcRNodes[i][j] = r.Copy_IndexLcNodeColBlock_LcRNodes[i][j];
            ColSizeColBlock_LcRNodes[i][j] = r.ColSizeColBlock_LcRNodes[i][j];
        }
    }

    /*
     *   This is a complicated enough operation that we put it into a
     *   subroutine
     */
    copyBlockMatrices(r.BlockMatrices);

    GI_ptr_ = r.GI_ptr_;
    LI_ptr_ = r.LI_ptr_;
    Func_ptr_  = r.Func_ptr_;

    return *this;
}
//===========================================================================

void
LocalRowNodeVBRIndices::freeBlockMatrices()
{
    for (int i = 0; i < NumLcRowNodes; i++) {
        if (BlockMatrices[i]) {
            int numBlocks = NumColBlocks_LcRNodes[i];
            Epetra_SerialDenseMatrix** rowBlocks = BlockMatrices[i];
            for (int j = 0; j < numBlocks; j++) {
                safeDelete(rowBlocks[j]);
            }
            delete[] BlockMatrices[i];
            BlockMatrices[i] = 0;
        }
    }
}
//===========================================================================
void
LocalRowNodeVBRIndices::initSize(int numLcRowNodes)
{
    NumLcRowNodes = numLcRowNodes;
    NumColBlocks_LcRNodes.resize(NumLcRowNodes);
    RowSize_LcRNodes.resize(NumLcRowNodes, 0);
    ColSizeColBlock_LcRNodes.resize(NumLcRowNodes, 0);
    IndexLcNodeColBlock_LcRNodes.resize(NumLcRowNodes, 0);
    Copy_IndexLcNodeColBlock_LcRNodes.resize(NumLcRowNodes, 0);
}
//===========================================================================
// determine the number of block columns in each block row
/*
 *   This determines the matrix stencil
 */
void
LocalRowNodeVBRIndices::determineBlockCol()
{
    int index;
    /*
     *
     * Here we restrict the stencil calculation to a very simple form
     * until we need and can develop a more extensive implementation
     */
    for (int i = 0; i < GI_ptr_->NumOwnedLcNodes; i++) {
        int GbNode = LI_ptr_->IndexGbNode_LcNode[i];
        RowSize_LcRNodes[i] = LI_ptr_->NumEqns_LcNode[i];
        NumColBlocks_LcRNodes[i] = 3;
        if (GbNode == 0) {
            NumColBlocks_LcRNodes[i] = 2;
        }
        if (GbNode == GI_ptr_->NumGbNodes - 1) {
            NumColBlocks_LcRNodes[i] = 2;
        }
    }

    /*
     *   Now we will create the basic block-diagonal matrix stencil
     */
    for (int i = 0; i < NumLcRowNodes; i++) {
        /*
         * Make sure that the vector is still filled with zeroes, so that
         * we don't have a memory leak
         */
        AssertTrace(IndexLcNodeColBlock_LcRNodes[i] == 0);
        /*
         *   Now malloc the int vector
         */
        IndexLcNodeColBlock_LcRNodes[i] = new int[NumColBlocks_LcRNodes[i]];
        Copy_IndexLcNodeColBlock_LcRNodes[i] = new int[NumColBlocks_LcRNodes[i]];
        ColSizeColBlock_LcRNodes[i] = new int[NumColBlocks_LcRNodes[i]];
        int* hh = ColSizeColBlock_LcRNodes[i];
        int* ih = IndexLcNodeColBlock_LcRNodes[i];

        /*
         * Do the left block
         */
        index = 0;

        int left = LI_ptr_->IDLeftLcNode_LcNode[i];
        if (left >= 0) {
            ih[index] = left;
            hh[index] = LI_ptr_->NumEqns_LcNode[left];
            index++;
        }

        /*
         * Do the diagonal block
         */
        ih[index] = i;
        hh[index] = LI_ptr_->NumEqns_LcNode[i];
        index++;
        /*
         * Do the right block
         */
        int right = LI_ptr_->IDRightLcNode_LcNode[i];
        if (right >= 0) {
            ih[index] = right;
            hh[index] = LI_ptr_->NumEqns_LcNode[right];
            index++;
        }
        /*
         * Do a sanity check on the process
         */
        AssertTrace(index == NumColBlocks_LcRNodes[i]);
    }
}
//===========================================================================
void
LocalRowNodeVBRIndices::mallocBlockMatrices()
{
    NumOwnedEqns = LI_ptr_->NumLcOwnedEqns;
    /*
     *   Now we will create the basic block-diagonal matrix stencil
     */
    BlockMatrices.resize(LI_ptr_->NumOwnedLcNodes, 0);
    for (int i = 0; i < GI_ptr_->NumOwnedLcNodes; i++) {
        int rowSize = RowSize_LcRNodes[i];
        int numBlocks = NumColBlocks_LcRNodes[i];
        BlockMatrices[i] = new Epetra_SerialDenseMatrix *[numBlocks];
        Epetra_SerialDenseMatrix** rowBlocks = BlockMatrices[i];
        for (int j = 0; j < numBlocks; j++) {
            int colSize = ColSizeColBlock_LcRNodes[i][j];
            rowBlocks[j] = new Epetra_SerialDenseMatrix(rowSize, colSize);
        }
    }
}
//===========================================================================
void
LocalRowNodeVBRIndices::copyBlockMatrices(const std::vector<Epetra_SerialDenseMatrix**>& c_BlockMatrices)
{

    int numOwnedLcNodes = c_BlockMatrices.size();
    AssertTrace(numOwnedLcNodes != 0);
    AssertTrace(numOwnedLcNodes == NumLcRowNodes);
    freeBlockMatrices();
    const Epetra_SerialDenseMatrix* c_block;
    BlockMatrices.resize(numOwnedLcNodes, 0);
    for (int i = 0; i < numOwnedLcNodes; i++) {
        int rowSize = RowSize_LcRNodes[i];
        int numBlocks = NumColBlocks_LcRNodes[i];
        BlockMatrices[i] = new Epetra_SerialDenseMatrix *[numBlocks];
        /*
         * Find the row size of the copied block and see that it is the same
         */

        Epetra_SerialDenseMatrix** rowBlocks = BlockMatrices[i];
        Epetra_SerialDenseMatrix** c_rowBlocks = c_BlockMatrices[i];
        for (int j = 0; j < numBlocks; j++) {
            int colSize = ColSizeColBlock_LcRNodes[i][j];
            c_block = c_rowBlocks[j];
            AssertTrace(c_block != 0);
            int c_colSize = c_block->N();
            int c_rowSize = c_block->M();
            AssertTrace(rowSize == c_rowSize);
            AssertTrace(colSize == c_colSize);
            /*
             *   Use the copy constructor to copy blocks at the lowest level
             */
            rowBlocks[j] = new Epetra_SerialDenseMatrix(*c_block);
        }
    }

}
//===========================================================================
Epetra_VbrMatrix*
LocalRowNodeVBRIndices::alloc_VbrMatrix()
{
    m1d::stream0 w0;

    m1d::GlobalIndices& GI = *GI_ptr_;

    /*
     *  Get a vector between the local nodes and the global node indices.
     */

    initSize(LI_ptr_->NumLcRowNodes);

    determineBlockCol();

    mallocBlockMatrices();

    /*
     *  Now create the matrix and vectors, which will have dimensions specified by our Map object.
     *
     *  Changing this from Copy to View caused the answer to change.
     */
    Epetra_VbrMatrix* A;
    //CopyMode = false;
    if (CopyMode) {
        A = new Epetra_VbrMatrix(Copy, *GI.GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap,
                                 *(LI_ptr_->GbBlockNodeEqnstoLcBlockNodeEqnsColMap), DATA_PTR(NumColBlocks_LcRNodes));
    } else {
        A = new Epetra_VbrMatrix(View, *GI.GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap,
                                 *(LI_ptr_->GbBlockNodeEqnstoLcBlockNodeEqnsColMap), DATA_PTR(NumColBlocks_LcRNodes));

    }

    /* NOW WE NEED TO FILL IN THE MATRIX STENCIL AND THEN FREEZE THE MATRIX */

    for (int iBlockRow = 0; iBlockRow < NumLcRowNodes; iBlockRow++) {

        int rowSize = RowSize_LcRNodes[iBlockRow];
        int colSize;
        int* bIndexes = IndexLcNodeColBlock_LcRNodes[iBlockRow];
        int numBlocks = NumColBlocks_LcRNodes[iBlockRow];
        int* colSize_row = ColSizeColBlock_LcRNodes[iBlockRow];
        Epetra_SerialDenseMatrix** rowBlock = BlockMatrices[iBlockRow];
        int* cp_bIndexes = Copy_IndexLcNodeColBlock_LcRNodes[iBlockRow];
        for (int j = 0; j < numBlocks; j++) {
            cp_bIndexes[j] = bIndexes[j];
        }
        A->BeginInsertMyValues(iBlockRow, numBlocks, cp_bIndexes);
        for (int j = 0; j < numBlocks; j++) {
            colSize = colSize_row[j];
            Epetra_SerialDenseMatrix& rcblock = *(rowBlock[j]);
            /*
             * Put some extra checks in here about making sure the sizes are consistent
             */
            AssertTrace(rowSize == rcblock.M());
            AssertTrace(colSize == rcblock.N());
            for (int jj = 0; jj < colSize; jj++) {
                double* dcol = rcblock[jj];
                for (int ii = 0; ii < rowSize; ii++) {
                    dcol[ii] = 0.0;
                }
            }
            A->SubmitBlockEntry(rcblock);
        }
        A->EndSubmitEntries();

    }

    A->FillComplete();

    A->PutScalar(0.0);

    return A;
}
//=====================================================================================================================
Epetra_VbrRowMatrix*
LocalRowNodeVBRIndices::alloc_VbrRowMatrix(Epetra_VbrMatrix* A)
{
    Epetra_VbrRowMatrix* B = new Epetra_VbrRowMatrix(A);
    return B;
}
//=====================================================================================================================
// Calculate the row sum scales
/*
 *  This calculates the sum of the abs value of all entries
 *  on a row and returns its inverse in a vector.
 *
 * @param rowScales  returns the row sum scales
 */
void
LocalRowNodeVBRIndices::calcInvRowSumScales(Epetra_Vector& rowScales) const
{

#ifdef DEBUG_MODE
    // Make sure the lengths are compatible
    int rlength = rowScales.MyLength();
    if (rlength != NumOwnedEqns) {
        throw m1d_Error("LocalRowNodeVBRIndices::calcRowSumScale", "wrong size of rowScales");
    }
#endif

    // zero the row scale vector
    rowScales.PutScalar(0.0);
    int iLc_Start = 0;
    for (int iBlockRow = 0; iBlockRow < NumLcRowNodes; iBlockRow++) {
        int numBlocks = NumColBlocks_LcRNodes[iBlockRow];
        //BlockMatrices[i] = new Epetra_SerialDenseMatrix *[numBlocks];
        // Get the column size of each block column in the current block row.
        int* ColSizeColBlock = ColSizeColBlock_LcRNodes[iBlockRow];
        /*
         * Find the row size of the copied block and see that it is the same
         */
        int numRowEqns = RowSize_LcRNodes[iBlockRow];
        Epetra_SerialDenseMatrix** rowBlocks = BlockMatrices[iBlockRow];
        for (int jColIndex = 0; jColIndex < numBlocks; jColIndex++) {
            int colSize = ColSizeColBlock[jColIndex];
            const Epetra_SerialDenseMatrix& rowColBlock = *(rowBlocks[jColIndex]);
            for (int je = 0; je < colSize; je++) {
                for (int ie = 0; ie < numRowEqns; ie++) {
                    rowScales[iLc_Start + ie] += fabs(rowColBlock(ie, je));
                }
            }
        }
        iLc_Start += numRowEqns;
    }
    for (int i = 0; i < NumOwnedEqns; i++) {
#ifdef DEBUG_MODE
        if (fabs(rowScales[i]) < 1.0E-100) {
            int iLcNode;
            int iGbNode;
            int iNodeEqnNum;
            EqnType var;
            EQ_TYPE_SUBNUM vtsub;
            std::string ename = Func_ptr_->equationID(i, iLcNode, iGbNode, iNodeEqnNum, var, vtsub);

            throw m1d_Error("LocalRowNodeVBRIndices::calcRowSumScale",
                            "row is nearly zero:" + ZZCantera::int2str(i) + ", value =  " +
                            ZZCantera::fp2str(rowScales[i]) + "\n" + ename + "\n");
        }
#endif
        rowScales[i] = 1.0 /  rowScales[i];
    }
}
//===========================================================================
// Print out the matrix structure of the VBR Matrix
/*
 *   Using this routine to learn more about Epetra
 * @param A  Epetra A
 * @param oo  ostream
 */
void
LocalRowNodeVBRIndices::queryMatrixStructure(const Epetra_VbrMatrix* const A, std::ostream& oo)
{
    oo << "queryMatrix for processor " << Comm_ptr_->MyPID() << endl;
    int numMyBlockRows = A->NumMyBlockRows();
    oo << "number of rows = " << numMyBlockRows << endl;
    int row_dim;
    int num_block_entries;
    int* block_indices;

    // Check total rows
    int numLcRowNodes = NumLcRowNodes;
    AssertTrace(numLcRowNodes == numMyBlockRows);

    for (int iBlockRow = 0; iBlockRow < numMyBlockRows; iBlockRow++) {

        A->BeginExtractMyBlockRowView(iBlockRow, row_dim, num_block_entries, block_indices);

        int numBlocks = NumColBlocks_LcRNodes[iBlockRow];
        AssertTrace(numBlocks == num_block_entries);

        int* bIndexes = IndexLcNodeColBlock_LcRNodes[iBlockRow];
        int* bColSizes = ColSizeColBlock_LcRNodes[iBlockRow];

        std::vector<Epetra_SerialDenseMatrix*> matrices(num_block_entries);
        std::vector<int> gnc(num_block_entries);
        oo.width(10);
        int grn = A->GRID(iBlockRow);

        // Check the global node id agreement
        int GbNode = LI_ptr_->IndexGbNode_LcNode[iBlockRow];
        AssertTrace(GbNode == grn);

        oo << "row " << iBlockRow << ":" << grn << " with neqs=";
        oo << row_dim;
        oo << " ( ";
        for (int j = 0; j < num_block_entries; j++) {

            if (bIndexes[j] != block_indices[j]) {
                oo << "j = " << j << ", bIndexes[j] = " << bIndexes[j];
                oo << ", block_indices[j] = " << block_indices[j] << "\n";
            }

            gnc[j] = A->GCID(block_indices[j]);
            oo << block_indices[j];
            oo << ":" << gnc[j] << " colSize=" << bColSizes[j];
            if (j < num_block_entries - 1) {
                oo << ", ";
            }
        }

        Epetra_SerialDenseMatrix** rowBlock = BlockMatrices[iBlockRow];
        oo << " )\n";
        oo << "                          ";
        bool doMorePrinting = true;
        for (int j = 0; j < num_block_entries; j++) {
            // Epetra_SerialDenseMatrix *rowColBlock = rowBlock[j];
            // Go get the pointer to the matrix block
            A->ExtractEntryView(matrices[j]);
            oo << block_indices[j];
            oo << ":" << (void*)(matrices[j]->A());
            oo << "   , ";
            if (rowBlock[j]->A() != matrices[j]->A()) {
                doMorePrinting = true;
            }
        }
        oo << " )\n";
        if (doMorePrinting) {
            oo << "          LocalVBRBlock:  ";
            for (int j = 0; j < num_block_entries; j++) {
                // Epetra_SerialDenseMatrix *rowColBlock = rowBlock[j];
                // Go get the pointer to the matrix block
                oo << block_indices[j];
                oo << ":" << (void*)(rowBlock[j]->A());
                oo << "   , ";
            }
            oo << " )\n";
        }
    }
}//===========================================================================
void
LocalRowNodeVBRIndices::informProblem(ProblemResidEval* r_ptr)
{
    Func_ptr_ = r_ptr;
}
//====================================================================================================================

//===========================================================================
} // end of namespace m1d
//===========================================================================
