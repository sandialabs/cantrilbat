/**
 *  @file m1d_EpetraJac.h This class stores an Epetra block matrix
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#ifndef M1D_EPETRAJAC_H
#define M1D_EPETRAJAC_H

#include "m1d_ProblemResidEval.h"
#include "m1d_VBRIndices.h"

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_VbrMatrix.h"
#include "BlockEntryGlobal.h"

#include "m1d_RecordTree_base.h"

//----------------------------------------------------------------------------------------------------------------------------------
namespace BEInput
{
class BlockEntry;
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
enum SolverType {
    Direct = 0,
    Iterative
};
//==================================================================================================================================
//! Class EpetraJac evaluates the Jacobian of a system of equations
/*!
 *  Class EpetraJac evaluates the Jacobian of a system of equations defined by a residual function supplied
 *  by an instance of class 'ProblemResidEval.' The residual function may consist of several linked
 *  1D domains, with different variables in each domain.
 */
class EpetraJac
{

public:

    //! Record tree format for input to the solvers
    /*!
     *  Additions to the record tree for the solver.
     */
    class BEinput_EpetraJac : public m1d::RecordTree_base
    {
    public:
        //! Default Constructor
        BEinput_EpetraJac();

        //! virtual destructor
        ~BEinput_EpetraJac();

        //! Solver type
        m1d::SolverType I_solverType;

        //! Direct solver name
        std::string S_directSolverName;
    };

    //! Static function to create a parsing block for the jacobian object to read
    static RecordTree_base* setupMDinput_pass1(BEInput::BlockEntry* Parent, RecordTree_base* dbb = 0);

    //! Constructor.
    /*!
     *  @param[in]             r                    function which calculates the residual
     */
    EpetraJac(ProblemResidEval& r);

    /// Destructor. Does nothing.
    virtual ~EpetraJac();

    //! Copy Constructor
    /*!
     *  (NOT OPERATIONAL)
     *  @param[in]           r                      Object to be copied
     */
    EpetraJac(const EpetraJac& r);

    //! Assigment operator
    /*!
     *  (NOT OPERATIONAL)
     *  @param[in]           r                      Object to be copied
     *
     *  @return                                     Returns a reference to the newly copied matrix
     */
     EpetraJac& operator=(const EpetraJac& r);

    //! Process the input file
    /*!
     *  @param[in]             dbb                 Database of entries from input file processing
     */
    void process_BEinput(const RecordTree_base* const dbb);

    //! Process the input file BlockEntry structure directly, extracting needed input for setting up the linear solver
    /*!
     *  @param[in]             cf                  Pointer to the BlockEntry database
     */
    void process_input(BEInput::BlockEntry* const cf);

    //! Allocate the matrix and other vectors
    void allocateMatrix();

    //! Calculate a Jacobian -> internal
    /*!
     *
     *  @param[in]           doTimeDependentResid Boolean indicating whether we should formulate the time dependent residual
     *  @param[in]           solnBase             Solution vector containing the base solution
     *  @param[in]           solnDotBase          Solution derivative vector containing the base solution
     *  @param[in]           t                    current value of the time
     *  @param[in]           rdelta_t             Inverse of the deltaT value
     *  @param[in]           solveType            Type of the solution solve
     */
    void matrixEval(const bool doTimeDependentResid, const Epetra_Vector* const solnBase,
                    const Epetra_Vector* const solnDotBase, const double t, const double rdelta_t,
                    const m1d::Solve_Type_Enum solveType);

    //! Calculate a Jacobian and a residual
    /*!
     *
     *  @param[in]           doTimeDependentResid Boolean indicating whether we should formulate the time dependent residual
     *  @param[in]           solnBase             Solution vector containing the base solution
     *  @param[in]           solnDotBase          Solution derivative vector containing the base solution
     *  @param[out]          resBase              Base value of the residual
     *  @param[in]           t                    current value of the time
     *  @param[in]           rdelta_t             Inverse of the deltaT value
     *  @param[in]           solveType            Type of the solution solve
     */
    void matrixResEval(const bool doTimeDependentResid, const Epetra_Vector* const solnBase,
                       const Epetra_Vector* const solnDotBase, Epetra_Vector* const resBase,
                       const double t, const double rdelta_t, const m1d::Solve_Type_Enum solveType);

    //! Calculate a Jacobian -> internal
    /*!
     *
     *  @param[in]           doTimeDependentResid Boolean indicating whether we should formulate the time dependent residual
     *  @param[in]           solnBase             Solution vector containing the base solution
     *  @param[in]           solnDotBase          Solution derivative vector containing the base solution
     *  @param[out]          resBase              Base value of the residual
     *  @param[in]           t                    current value of the time
     *  @param[in]           rdelta_t             Inverse of the deltaT value
     *  @param[in]           solveType            Type of the solution solve
     */
    void matrixEval1(const bool doTimeDependentResid, const Epetra_Vector* const solnBase, const Epetra_Vector* const solnDotBase,
                     Epetra_Vector* const resBase, const double t, const double rdelta_t, const m1d::Solve_Type_Enum solveType);

private:
    //! Evaluate the Jacobian at x0.
    /*!
     * The unperturbed residual function is resid0, which must be supplied on input. The
     * third parameter 'rdt' is the reciprocal of the time step. If zero, the steady-state Jacobian is evaluated.
     *
     *  @param[in]           doTimeDependentResid Boolean indicating whether we should formulate the time dependent residual
     *  @param[in]           solnBase             Solution vector containing the base solution
     *  @param[in]           solnDotBase          Solution derivative vector containing the base solution
     *  @param[out]          resBase              Base value of the residual
     *  @param[in]           t                    current value of the time
     *  @param[in]           rdelta_t             Inverse of the deltaT value
     */
    void eval(const bool doTimeDependentResid, const Epetra_Vector* const solnBase, const Epetra_Vector* const solnDotBase,
              const Epetra_Vector& residBase, const double t, const double rdelta_t);

public:

    //! Fill in the Block VBR matrix
    /*!
     *   We move over each blockRow and then each block col for that block row.
     *   We then use the SubmitBlockEntry() option to fill in the VBR matrix from the  LRN_VBR_ptr_ global structure
     *
     *    A_->SubmitBlockEntry(*rowColBlock)
     */
    void fillVbr();

    //! Get the row scales for the matrix
    /*!
     *  In this calculation we sum up the absolute value of all elements on a row. Then take the inverse.
     *
     *  @param[out]          rowScales           Epetra_Vector of length num local  equations on the processor.
     */
    void getRowScales(Epetra_Vector* const rowScales) const;

    //! Scale the columns of the matrix with a vector
    /*!
     *  @param[in]           colScales           This is a vector with map columnMap() that scales the columns of the matrix
     */
    void columnScale(const Epetra_Vector* const colScales);

    //! Scale the columns of the matrix with a vector
    /*!
     * @param rowScales  This is a vector with map rangeMap()
     *                   that scales the rows of the matrix
     */
    void rowScale(const Epetra_Vector* const rowScales);

    //! Print out the matrix structure of the VBR Matrix
    /*!
     * @param oo  ostream
     */
    void queryMatrixStructure(std::ostream& oo);

    //! Elapsed CPU time spent computing the Jacobian.
    double elapsedTime() const;

    //! Number of Jacobian evaluations.
    /*!
     *  @return                                  Returns the number of jacobian evaluations done in the entire calculation
     */
    int nEvals() const;

    //! Return the number of times 'incrementAge' has been called since the
    /*!
     *  @return                                  Returns the age of the jacobian
     */
    int age() const;

    //! Increment the Jacobian age.
    void incrementAge();

    //! Add the transient jacobian term back into the diagonal
    /*!
     *  (INOPERABLE)
     *
     *  @param[in]           rdt                 inverse of deltaT
     *  @param[in]           mask                mask declaring which terms get a deltat term
     */
    void updateTransient(double rdt, int* mask);

    //! Set the age of the jacobian
    /*!
     *  @param[in]           age                 Value to set the age
     */
    void setAge(int age);

    //! Accessor routine for the mask variable
    /*!
     *  @return                                  Returns a changeable reference to the mask variable
     */
    Epetra_IntVector& transientMask();

    //! Increment the diagonal of the matrix
    /*!
     *  (INOPERABLE)
     *
     *  @param[in]           j                   local row number
     *  @param[in]           d                   Value of the diagonal of the matrix to replace the entry with
     */
    void incrementDiagonal(int j, double d);

    //! Return a reference to the delta for the solution variables used to calculate the one-sided numerical jacobian
    /*!
     *  @return                                  Returns an Epetra_Vector with ghost nodes containing the delta variables
     */ 
    const Epetra_Vector& deltaSolnJac() const;

    //! Zero the matrix
    void zeroMatrix();

    //! Return a changeable pointer into the matrix given Global Block Row indices
    /*!
     *   I don't think this will ever be used in practice, because of the speed. But,
     *   here it is.
     * @param gbRow          Global row number
     * @param lcRowIndex     local row index of the equation on that row
     * @param gbCol          Global col number
     * @param lcColIndex     local col index of the variable on that col
     * @return    Returns a pointer to the position in the block matrix corresponding
     *            to that row and column. Will return 0 if this number is not on the current
     *            processor. Will return 0 if this number isn't in the sparce matrix stencil.
     */
    double* GbBlkValue(int gbRow, int lcRowIndex, int gbCol, int lcColIndex) const;

    //! Return a changeable reference into the matrix given Global equation numbers
    /*!
     *   I don't think this will ever be used in practice, because of the speed. But,
     *   here it is. This is overloaded with the (i,j) operator
     *
     * @param gbRow          Global row number
     * @param lcRowIndex     local row index of the equation on that row
     * @param gbCol          Global col number
     * @param lcColIndex     local col index of the variable on that col
     * @return    Returns a pointer to the position in the block matrix corresponding
     *            to that row and column. Will return 0 if this number is not on the current
     *            processor. Will return 0 if this number isn't in the sparce matrix stencil.
     */
    double& operator()(const int iGlobalEqn, const int jGlobalEqn);

    //! Return a changeable reference into the matrix given the global row number and the global column number
    /*!
     *  I don't think this will ever be used in practice, because of the speed. But,
     *  here it is. This is overloaded with the (i,j) operator.
     *
     *  @param[in]           iGlobalEqn          Global row number
     *  @param[in]           jGlobalEqn          Global col number
     * 
     *  @return                                  Returns a changeable reference to the position in the block matrix corresponding
     *                                           to that row and column. 
     *                                           Will throw an error if this number is not on the current
     *                                           processor. Will throw an error if this number isn't in the sparce matrix stencil.
     */
    double& value(const int iGlobalEqn, const int jGlobalEqn);

    //! Return the value of the matrix given by Global row numbers and global column numbers
    /*!
     *  I don't think this will ever be used in practice, because of the speed. But,
     *  here it is. This is overloaded with the (i,j) operator
     *
     *  @param[in]           iGlobalEqn          Global row equation number
     *  @param[in]           jGlobalEqn          Global col equation number
     *
     *  @return                                  Returns the value of the matrix entry
     *                                           Returns zer if this number is not on the current processor, or if the number isn't
     *                                           in the sparse matrix stencil.
     */
    double value(const int iGlobalEqn, const int jGlobalEqn) const;

    //! Returns the number of global equations in the matrix
    /*!
     *  @return                                  Returns the number of global equations in the matrix
     */
    int nRows() const;

    //! Number of columns
    int nColumns() const;

    //! Number of subdiagonals
    int nSubDiagonals() const;

    //! Number of superdiagonals
    /*!
     *   @return                                 Returns the number of superdiagonals
     */
    int nSuperDiagonals() const;

    //! Returns the bandwidth of the matrix
    /*!
     *   @return                                 Returns the value of  2 * m_kl + m_ku + 1
     */
    int ldim() const;

    //! Multiply A*b and write result to prod.
    /*!
     *  @param[in]           b                   The const value of b
     *  @param[out]          prod                The vector, A * b , which is the result
     */
    void mult(const Epetra_Vector& b, Epetra_Vector& prod) const;

    //! Multiply b*A and write result to prod.
    /*!
     *  @param[in]           b                   The const value of b
     *  @param[out]          prod                The vector, b * A , which is the result
     */
    void leftMult(const Epetra_Vector& b, Epetra_Vector& prod) const;
  
    //! Factor the matrix
    /*!
     *  (NOT USED)
     *  @return                                  Returns 0 if the factorization was successful
     */
    int factor();

    //! Main routine to solve the linear system
    /*!
     *
     *  @param[in]             b                   This is the input rhs
     *  @param[out]             x                   Returns the solved system in x
     *  @param[out]            its                 Returns the number of iterations to solve the linear system
     *  @param[out]            norm                The norm of the value of (Ax-b) on return if doRes is true
     *  @param[in]             doRes               We check the Ax-b calc if true
     *
     *  @return                                    0 for a successful routine.
     *                                             Any other number indicates a failure
     */
    int solve(Epetra_Vector* const b, Epetra_Vector* const x, int& its, double& norm, bool doRes=false);

    // ------------------------------------------- D A T A -------------------------------------------------------------------

public:
    //! This is true if we have factored the matrix
    /*!
     *  Factoring occurs after column and row scaling operations.
     */
    bool m_factored;

    //! This is true if we have column Scaled the matrix
    /*!
     *  Column Scaling always occurs first before row scaling and/or factoring.
     */
    bool m_columnScaled;

    //! This is true if we have row-scaled the matrix
    /*!
     *  Row Scaling always occurs after column scaling and before factoring
     */
    bool m_rowScaled;

protected:

    //! Number of rows
    int m_n;

    //! Number of lower rows
    int m_kl;

    //! Number of upper rows
    int m_ku;

public:
    //! pointer to the Epetra Vbr matrix
    Epetra_VbrMatrix* A_;

    //! Pointer to the Epetra Vbr row matrix representation
    Epetra_VbrRowMatrix* Arow_;

    //! local copy of the Epetra_Comm ptr object
    Epetra_Comm* Comm_ptr_;

    //! Storage for the B - AX values
    /*!
     *  This is related to the analysis of the overall solve of the system
     */
    Epetra_Vector* m_BmAX;

    //! Storage for the AX product
    Epetra_Vector* m_AX;

    //! Pointer to the delta solution vector
    /*!
     *  This is used to store the deltas for the numerical jacobian
     */
    Epetra_Vector* deltaSoln_;

    //! Pointer to the residual evaluator.
    /*!
     *  This is not owned by this object
     */
    ProblemResidEval* m_resid;

    //! Relative tolerance
    double m_rtol;

    //! abslute tolerance
    double m_atol;

    //! Elapsed time
    double m_elapsed;

    //! Vector containing the diagonal entries of the jacobian
    Epetra_Vector* m_ssdiag;

    //! Matrix as represented by an CrsMatrix format
    Epetra_CrsMatrix* Acrs_;

    //! Mask indicating which variable is an algebraic constraint and which
    //! variable has a time derivative
    /*!
     *  1 is algebraic
     *  0 has a time derivative
     */
    Epetra_IntVector* m_isAlgebraic;


    //! Number of Jacobian formation evaluations
    int m_nevals;

    //! Age of the jacobian
    /*!
     * Note this variable is entirely user driven. Right now it is not connected to anything.
     */
    int m_age;

    //! Number of global equations
    /*!
     * This is the number of equations that are on all processors
     * combined
     */
    int m_NumGbEqns;

    //! Type of the solver to use: either Direct or Iterative
    SolverType solverType_;

    //! Name of the direct solver to use if using a direct solver
    std::string directSolverName_;

    //! Number of global block rows defined in the problem
    /*!
     * This is frequently equal to the number of global mesh points
     * in the problem
     */
    int m_GbBlockRows;

    //! Global indices class is the same on all processors
    /*!
     *  This is a global structure, containing global information.
     *  All information in this structure, except for my_procID, numLnNodes,
     *  and numLnEqns, is the same on all processors.
     *
     */
    GlobalIndices* GI_ptr_;

    //! Class setting up the local node indices on the current processor
    /*!
     *   These indices have to do with accessing the local nodes on the processor
     *
     *   All numbers are in local Row Node Format. Local Row node format is
     *   defined as the following. All owned nodes come first. They are
     *   ordered in terms of increasing global node number.
     *   Then the right ghost node is listed.
     *   Then the left ghost node is listed
     *   Then, the "globally-all-connected node is listed, if available.
     */
    LocalNodeIndices* LI_ptr_;

    //! This class stores and has facilities for accessing the VBR matrix
    /*!
     *  This object is owned and created by this EpetraJac object.
     */
    LocalRowNodeVBRIndices* LRN_VBR_ptr_;

    //! This vector stores a string describing the current variable changes used
    //! to calculate the jacobian by numerical differencing
    std::vector<std::string> varDiffString_LcNode;

public:
    //! jacobian type
    /*!
     *      0    SteadyState_Solve
     *      1    TimeDependentAccurate_Solve,
     *      2    TimeDependentRelax_Solve,
     *      3    DAESystemInitial_Solve
     */
    m1d::Solve_Type_Enum solveType_;
};
//==================================================================================================================================

//! Write the contents of a matrix to an ostream
/*!
 *
 *  @param[in]               os                  ostream
 *  @param[in]               m                   reference to an EpetraJac object
 *  @return                                      Returns a reference to the ostream
 */
std::ostream& operator<<(std::ostream& os, const EpetraJac& m);

//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------

#endif

