/**
 *  @file m1d_EpetraJac.h
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2012-02-23 14:34:18 -0700 (Thu, 23 Feb 2012) $
 *  $Revision: 5 $
 *
 */

#ifndef M1D_EPETRAJAC_H
#define M1D_EPETRAJAC_H

#include "m1d_ProblemResidEval.h"

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_VbrMatrix.h"
#include "BlockEntryGlobal.h"

#include "m1d_RecordTree_base.h"

namespace BEInput
{
  class BlockEntry;
}

namespace m1d
{

enum SolverType {
  Direct = 0, Iterative
};



/**
 * Class EpetraJac evaluates the Jacobian of a system of equations
 * defined by a residual function supplied by an instance of class
 * 'ProblemResidEval.' The residual function may consist of several linked
 * 1D domains, with different variables in each domain.
 */
class EpetraJac {

public:

  class BEinput_EpetraJac : public m1d::RecordTree_base  {
  public:
    BEinput_EpetraJac();
    ~BEinput_EpetraJac();
    m1d::SolverType I_solverType;
  };
  
  static RecordTree_base *
  setupMDinput_pass1(BEInput::BlockEntry * Parent, RecordTree_base *dbb = 0);

  void process_BEinput(RecordTree_base *dbb);

  //!Constructor.
  /*!
   * @param r  function which calculates the residual
   */
  EpetraJac(ProblemResidEval& r);

  /// Destructor. Does nothing.
  virtual
  ~EpetraJac();

  //! Allocate the matrix
  void
  allocateMatrix();

  //! Calculate a Jacobian -> internal
  /*!
   *
   *
   * @param doTimeDependentResid    Boolean indicating whether we should
   *                                formulate the time dependent residual
   * @param soln                    Solution vector
   */
  void
  matrixEval(const bool doTimeDependentResid,
             const Epetra_Vector * const solnBase,
             const Epetra_Vector * const solnDotBase,
             const double t,
             const double rdelta_t,
	     const Solve_Type_Enum solveType);

  //! Calculate a Jacobian and a residual
  /*!
   *
   *
   * @param doTimeDependentResid    Boolean indicating whether we should
   *                                formulate the time dependent residual
   * @param soln                    Solution vector
   */
  void
  matrixResEval(const bool doTimeDependentResid,
                const Epetra_Vector * const solnBase,
                const Epetra_Vector * const solnDotBase,
                Epetra_Vector * const resBase,
                const double t,
                const double rdelta_t,
		const Solve_Type_Enum solveType);

  //! Calculate a Jacobian -> internal
  /*!
   *
   *
   * @param doTimeDependentResid    Boolean indicating whether we should
   *                                formulate the time dependent residual
   * @param soln                    Solution vector
   */
  void
  matrixEval1(const bool doTimeDependentResid,
              const Epetra_Vector * const solnBase,
              const Epetra_Vector * const solnDotBase,
              Epetra_Vector * const resBase,
              const double t,
              const double rdelta_t,
	      const Solve_Type_Enum solveType);

private:
  //! Evaluate the Jacobian at x0.
  /*!
   * The unperturbed residual function is resid0, which must be supplied on input. The
   * third parameter 'rdt' is the reciprocal of the time
   * step. If zero, the steady-state Jacobian is evaluated.
   */
  void
  eval(const bool doTimeDependentResid,
       const Epetra_Vector * const solnBase,
       const Epetra_Vector * const solnDotBase,
       const Epetra_Vector & residBase,
       const double t,
       const double rdelta_t);

public:

  void
  fillVbr();

  //! Get the row scales for the matrix
  /*!
   * In this calculation we sum up the absolute value of
   * all elements on a row. Then take the inverse.
   *
   * @param rowScales  Epetra_Vector of length num local
   *                   equations on the processor.
   */
  void
  getRowScales(Epetra_Vector * const rowScales) const;

  //! Scale the columns of the matrix with a vector
  /*!
   * @param colScales  This is a vector with map columnMap()
   *                   that scales the columns of the matrix
   */
  void
  columnScale(const Epetra_Vector * const colScales);

  //! Scale the columns of the matrix with a vector
  /*!
   * @param rowScales  This is a vector with map rangeMap()
   *                   that scales the rows of the matrix
   */
  void
  rowScale(const Epetra_Vector * const rowScales);

  //! Print out the matrix structure of the VBR Matrix
  /*!
   * @param oo  ostream
   */
  void
  queryMatrixStructure(std::ostream &oo);

  //! Elapsed CPU time spent computing the Jacobian.
  double
  elapsedTime() const;

  //! Number of Jacobian evaluations.
  int
  nEvals() const;

  //! Return the number of times 'incrementAge' has been called since the
  int
  age() const;

  //! Increment the Jacobian age.
  void
  incrementAge();

  void
  updateTransient(double rdt, int* mask);

  //! Set the age.
  /*!
   *
   * @param age  Value to set the age
   */
  void
  setAge(int age);

  //! Accessor routine for the mask variable
  /*!
   *
   * @return Returns a changeable reference to the mask variable
   */
  Epetra_IntVector&
  transientMask();

  //! Increment the diagonal of the matrix
  /*!
   *
   */
  void
  incrementDiagonal(int j, double d);


  //! Zero the matrix
  void
  zeroMatrix();

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
  double *
  GbBlkValue(int gbRow, int lcRowIndex, int gbCol, int lcColIndex) const;

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
  double&
  operator()(const int iGlobalEqn, const int jGlobalEqn);

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
  double&
  value(const int iGlobalEqn, const int jGlobalEqn);

  //! Return the value of the matrix given by Global equation numbers
  /*!
   *   I don't think this will ever be used in practice, because of the speed. But,
   *   here it is. This is overloaded with the (i,j) operator
   *
   * @param iGlobalEqn          Global row equation number
   * @param jGlobalEqn          Global col equation number
   * @return    Returns the value of the matrix entry
   */
  double
  value(const int iGlobalEqn, const int jGlobalEqn) const;

  //! Return the value of the matrix given by Global equation numbers
  /*!
   *   I don't think this will ever be used in practice, because of the speed. But,
   *   here it is. This is overloaded with the (i,j) operator
   *
   * @param iGlobalEqn          Global row equation number
   * @param jGlobalEqn          Global col equation number
   * @return   Returns the value of the matrix entry
   */
  double
  _value(const int iGlobalEqn, const int jGlobalEqn) const;

  //! Returns the number of global equations in the matrix
  int
  nRows() const;

  //! Number of columns
  int
  nColumns() const;

  //! Number of subdiagonals
  int
  nSubDiagonals() const;

  //! Number of superdiagonals
  int
  nSuperDiagonals() const;

  int
  ldim() const;

  //! Multiply A*b and write result to prod.
  void
  mult(const Epetra_Vector &b, Epetra_Vector &prod) const;

  /// Multiply b*A and write result to prod.
  void
  leftMult(const Epetra_Vector &b, Epetra_Vector &prod) const;

  int
  factor();

  //! Main routine to solve the linear system
  /*!
   *
   * @param b  This is the input rhs
   * @param x  Returns the solved system in x
   * @param its Returns the number of its to solve the linear system
   * @param norm of the value Ax-b on return if doRes is postivite
   * @param doRes We check the Ax-b calc if true
   * @return  0 for a successful routine.
   */
  int
  solve(Epetra_Vector *b, Epetra_Vector *x, int &its,  double &norm, bool doRes=false);

public:
  //! This is true if we have factored the matrix.
  /*!
   *  Factoring occurs after column and row scaling operations.
   */
  bool m_factored;

  //! This is true if we have column Scaled the matrix
  /*!
   *  Column Scaling always occurs first before row scaling and/or factoring.
   */
  bool m_columnScaled;

  //! This is true if we have row Scaled the matrix
  /*!
   * row Scaling always occurs after column scaling and before factoring
   */
  bool m_rowScaled;

protected:

  int m_n, m_kl, m_ku;
public:
  //! pointer to the Epetra Vbr matrix
  Epetra_VbrMatrix * A_;

  //! local copy of the Epetra_Comm ptr object
  Epetra_Comm *Comm_ptr_;

  Epetra_Vector *m_BmAX;
  Epetra_Vector *m_AX;
  //! pointer to the residual evaluator.
  /*!
   * This is not owned by this object
   */
  ProblemResidEval* m_resid;
  double m_rtol, m_atol;
  double m_elapsed;
  Epetra_Vector *m_ssdiag;

  Epetra_CrsMatrix * Acrs_;


  //! Mask indicating which variable is an algebraic constraint and which
  //! variable has a time derivative
  /*!
   *  1 is algebraic
   *  0 has a time derivative
   */
  Epetra_IntVector *m_isAlgebraic;


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

  SolverType solverType_;

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
  GlobalIndices *GI_ptr_;

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
  LocalNodeIndices *LI_ptr_;

  //! This class stores and has facilities for accessing the VBR matrix
  /*!
   *  This object is owned and created by this EpetraJac object.
   */
  LocalRowNodeVBRIndices *LRN_VBR_ptr_;


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
  Solve_Type_Enum solveType_;
};

//! Write the contents of a matrix to an ostream
/*!
 *
 * @param os   ostream
 * @param m    reference to an EpetraJac object
 * @return     Returns a reference to the ostream
 */
std::ostream&
operator<<(std::ostream& os, const EpetraJac& m);

//======================================================================================
} // end of m1d namespace
//======================================================================================

#endif

