/**
 * @file m1d_VBRIndices.h Header for the class LocalRowNodeVBRIndices
 * that manages the interaction with the VBR matrix and the storage
 * for the blocks of the VBR matrix (see \ref matrixInteraction
 * and class \link m1d::LocalRowNodeVBRIndices LocalRowNodeVBRIndices\endlink).
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef _M1D_VBRINDICES_H
#define _M1D_VBRINDICES_H

#include "m1d_defs.h"
#include "m1d_LocalNodeIndices.h"

#include "Epetra_VbrRowMatrix.h"

#include <vector>
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

class GlobalIndices;
class ProblemResidEval;

//==================================================================================================================================
/**
 * @defgroup matrixInteraction How to interact with the VBR Matrix
 *
 *  We describe how to interact with the VBR Matrix stored within
 *  an Epetra_VbrMatrix object.
 *
 *  We describe how to fill the matrix.
 *
 *  We describe the difference between copy and view.
 *
 *
 *
 */

//==================================================================================================================================
//! This class stores and has facilities for accessing the VBR matrix on this processor
/*!
 *   All numbers are in local Row Node Format. Local row node format is a
 *   sequential local listing of the nodes on the current processor. An index
 *   is then used to associate the local node number with the global node number.
 *   The local row node format lists all of the locally owned nodes on the processor
 *   first. Then, the external nodes are listed next.
 *
 *   The first numLRNodes are nodes that are owned by this processor.
 *
 *     0
 *     . . .
 *     numLRNodes
 *
 *
 *
 */
class LocalRowNodeVBRIndices
{
public:

    //! Default constructor
    /*!
     *  @param[in]           comm_ptr             Communications pointer
     *  @param[in]           copyMode             Mode of the interaction
     *  @param[in]           gi_ptr               Global indices
     *  @param[in]           li_ptr               Local node indices
     */
    LocalRowNodeVBRIndices(Epetra_Comm* comm_ptr, bool copyMode, GlobalIndices* gi_ptr, LocalNodeIndices* li_ptr);

    //! Default destructor
    ~LocalRowNodeVBRIndices();

    //! Copy constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    LocalRowNodeVBRIndices(const LocalRowNodeVBRIndices& r);

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied
     *  @return                                  returns the copied object
     */
    LocalRowNodeVBRIndices& operator=(const LocalRowNodeVBRIndices& r);

    //! Resize the arrays when we know the value of
    //! the number of block rows on this processor.
    /*!
     *  @param               numLcRNodes         Number of block rows on each processor
     */
    void initSize(int numLcRNodes);

    //! Determine the matrix stencil needed to fully solve the system of equations
    /*!
     *  We do this by calculating the block column for each block row.
     *
     *  This in effect determines the matrix stencil
     */
    void determineBlockCol();

    //! Allocate the blocks that make up the matrix structure.
    /*!
     *
     */
    void mallocBlockMatrices();

    //! Allocate the matrix, creating and fixing the matrix stencil
    /*!
     *  @return                                  returns the matrix object as a pointer to an Epetra_VbrMatrix object.
     */
    Epetra_VbrMatrix* alloc_VbrMatrix();

    //! Allocate an Epetra_VbrRowMatrix matrix given a pointer to an Epetra_VbrMatrix
    /*!
     *  The VbrRowMatrix is suppose to fix up the nasty bug in the VbrMatix class. This will
     *  be an attempt to extend the matrix format to see if this is the case.
     *
     *  @param[in]           A                   Pointer to a previously allocated Epetra_VbrMatrix
     *  @return                                  returns the matrix object as a pointer to an Epetra_VbrRowMatrix object.
     */
    Epetra_VbrRowMatrix* alloc_VbrRowMatrix(Epetra_VbrMatrix* A);

    //! Calculate the row sum scales
    /*!
     *  This calculates the sum of the abs value of all entries on a row and returns its inverse in a vector.
     *
     *  @param[out]          rowScales           returns the row sum scales
     */
    void calcInvRowSumScales(Epetra_Vector& rowScales) const;

    //! Set the copy mode for access to the VBR matrix
    /*!
     *  We must set the copy mode before calling the fillComplete Epetra function.
     *
     *  @param[in]           copyMode            boolean indicating what the mode is
     */
    void setCopyMode(bool copyMode);

    //! Report the copy mode
    /*!
     *  A true value means we are in copy mode. A false value means we are in view mode.
     *
     *  @return                                  Returns whether we are in copy mode or view mode
     */
    bool copyMode() const;

    //! Print out the matrix structure of the VBR Matrix
    /*!
     *  @param[in]           A                   Pointer to an Epetra_VbrMatrix object 
     *  @param[out]          oo                  ostream
     */
    void queryMatrixStructure(const Epetra_VbrMatrix* const A, std::ostream& oo);


    //! Store the ProblemResidEval pointer so that we can access names of rows and columns
    /*!
     *  @param[in]           r_ptr               Pointer to the current ProblemResidEval
     */
    void informProblem(ProblemResidEval* r_ptr);

private:

    //!  Free the Block matrices
    /*!
     * This utility routine frees the block matrices
     */
    void freeBlockMatrices();

    //! Copy the block matrices into this structure
    /*!
     *  The supporting structures must already be prepped within the current object
     *
     *  @param[in]           rBlockMatrices      Vector of pointers to pointers of EpetraSerialDenseMatrix objects
     */
    void copyBlockMatrices(const std::vector<Epetra_SerialDenseMatrix**>& rBlockMatrices);

    /* ---------------------------------------------------------- D A T A ------------------------------------------------ */
public:

    //! Local copy of the Epetra_Comm ptr object
    /*!
     *  This is a shallow copy. We do not own this object
     */
    Epetra_Comm* Comm_ptr_;

    //! If true, we are using VBR_copy mode. If not, we are using view mode.
    bool CopyMode;

    //! Number of block rows on this processor in the VBR matrix. This processor owns the equations on these rows.
    /*!
     *  This is usually equal to the number of owned nodes on the processor.
     *  This is equal to GI_ptr->numLcNodes_Proc[MyProcID];
     */
    int NumLcRowNodes;

    //! Number of local equations on this processor in the VBR matrix. This processor owns the equations on these rows.
    int NumOwnedEqns;

    //! Vector containing the number of block columns for each block row within the VBR matrix
    /*!
     *  NumColBlocks_LcRNodes[iBlockRow] is the number of block columns at each iBlockRow.
     *  Note this is equivalent to the bandwidth for LAPACK banded matrices.
     *
     *  Length: NumLcRNodes
     */
    std::vector<int> NumColBlocks_LcRNodes;

    //! Vector containing the row size of each block matrix within the VBR matrix on the processor for each block row
    /*!
     *  RowSize_LcRNodes[iBlockRow] is the number of rows in each block matrix within the block row, iBlockRow,
     *  on the current processor.
     *
     *  Length: NumLcRNodes
     */
    std::vector<int> RowSize_LcRNodes;

    //! Vector over each block row containing a vector of the sizes of each block column within the VBR matrix
    /*!
     *  To get the column size for all of the blocks in a particular block row:
     *
     *      int *ColSize_ColBlock = ColSizeColBlock_LcRNodes[iBlockRow]
     *
     *  ColSizeColBlock_LcRNodes[iBlockRow] is a vector containing a vector, ColSize_ColBlock[], of length NumColblocks[iBlockRow].
     *
     *      ColSize_ColBlock[iBlockCol]
     *
     *  is the number of  columns in the block matrix located at the block position: (iBlockRow, iBlockCol)
     *
     *  The vector has a length equal to NumLcRNodes, the number of owned nodes on the current processor.
     *  Length: NumLcRNodes
     */
    std::vector<int*> ColSizeColBlock_LcRNodes;

    //! Each entry is a vector containing the block column indices for each block row.
    /*!
     *  To get the block  column indices for the block columns in a block row:
     *
     *      int* LnNodeBlockCols = IndexLcNodeColBlock_LcRNodes[iBlockRow]
     *
     *  IndexLcNodeColBlock_LcRNodes[iBlockRow] is a vector containing a vector, LnNodeBlockCols[] of 
     *  length NumColblocks[iBlockRow].
     *
     *      LnNodeBlockCols[iBlockCol]
     *
     *  is the local node id of the iBlockCol'th block column of the iBlockRow block row.
     *
     *  The vector has a length equal to NumLcRNodes, the number of owned nodes on the current processor.
     *  Length: NumLcRNodes
     */
    std::vector<int*> IndexLcNodeColBlock_LcRNodes;

    //! Copy of the Vector containing the block column indices for each block row.
    /*!
     *     int* LnNodeBlockCols  =      IndexLcNodeColBlock_LcRNodes[iBlockRow]
     *
     *     LnNodeBlockCols[iBlockCol] is the local node id of the iBlockCol'th block column of the iBlockRow  block row.
     */
    std::vector<int*> Copy_IndexLcNodeColBlock_LcRNodes;

    //! Vector of Pointers into the value storage of each Block Matrices
    /*!
     *  This contains pointers into each block matrix.
     *
     *  The outer loop is over block rows owned by the processor.
     *         Epetra_SerialDenseMatrix** rowBlock = BlockMatrices[iBlockRow];
     *
     *  The next loop is over malloced  pointers to Epetra_SerialDenseMatrices of length numBlocks:
     *          BlockMatrices[i] = new Epetra_SerialDenseMatrix* [numBlocks]
     *
     *  The inner loop is malloced Epetra_SerialDenseMatrix objects of size (rowSize, colSize)
     */
    std::vector<Epetra_SerialDenseMatrix**> BlockMatrices;

    //! Global Indices for the current problem
    /*!
     *  This is a shallow copy.
     *  We do not own this pointer.
     */
    GlobalIndices* GI_ptr_;

    //! Local Indices for the current problem
    /*!
     *  This is a shallow copy.
     *  We do not own this pointer.
     */
    LocalNodeIndices* LI_ptr_;

    //! Pointer to the ProblemResidEval object
    /*!
     *  This is a shallow pointer, and we do not own this pointer.
     *  We use this to get the equation and variable names of rows and columns.
     */
    ProblemResidEval* Func_ptr_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
