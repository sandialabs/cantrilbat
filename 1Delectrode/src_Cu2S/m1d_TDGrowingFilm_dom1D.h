/*
 * m1d_SimpleTDDiff_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#ifndef M1D_TDGROWINGFILM_DOM1D_H_
#define M1D_TDGROWINGFILM_DOM1D_H_

//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_DomainDescription.h"
#include "m1d_BulkDomain1D.h"

namespace m1d
{
class LocalNodeIndices;

//! Class that does a simple diffusion operator on the
//! of the bulk domain class
class TDGrowingFilm_dom1D : public BulkDomain1D {

public:

  //! Constructor
  /*!
   * @param bdd   Contains the bulk domain description.
   */
  TDGrowingFilm_dom1D(m1d::BulkDomainDescription &bdd);

  //! Copy constructor
  /*!
   * @param r      Object to be copied into the current object
   */
  TDGrowingFilm_dom1D(const TDGrowingFilm_dom1D &r);

  //! Destructor
  virtual
  ~TDGrowingFilm_dom1D();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  TDGrowingFilm_dom1D &
  operator=(const TDGrowingFilm_dom1D &r);

  //! Prepare all of the indices for fast calculation of the residual
  /*!
   *  Ok, at this point, we will have figured out the number of equations
   *  to be calculated at each node point. The object NodalVars will have
   *  been fully formed.
   *
   *  We use this to figure out what local node numbers/ cell numbers are
   *  needed and to set up indices for their efficient calling.
   *
   *  Child objects of this one will normally call this routine in a
   *  recursive fashion.
   */
  virtual void
  domain_prep(LocalNodeIndices *li_ptr);

  //! Basic function to calculate the residual for the domain.
  /*!
   *  All residual terms are written with the following sign convention
   *  based on keeping the time derivative term positive.
   *
   *       res = dcdt - dc2 /dx2 - src = 0
   *
   * @param res  Output vector containing the residual
   * @param doTimeDependentResid  boolean indicating whether the time
   *                         dependent residual is requested
   * @param soln_ptr     solution vector at which the residual should be
   *                     evaluated
   * @param solnDot_ptr  solution dot vector at which the residual should
   *                     be evaluated.
   *  @param t           time
   *  @param rdelta_t    inverse of delta_t
   */
  virtual void
  residEval(Epetra_Vector &res,
            const bool doTimeDependentResid,
            const Epetra_Vector *soln_ptr,
            const Epetra_Vector *solnDot_ptr,
            const Epetra_Vector *solnOld_ptr,
            const double t,
            const double rdelta_t,
            const ResidEval_Type_Enum residType = Base_ResidEval,
	    const Solve_Type_Enum solveType = TimeDependentAccurate_Solve);


  //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
  /*!
   * @param atolDefault             Default atol value
   * @param soln                    Solution vector. This is a constant
   *                                the residual calculation.
   * @param atolVector              (OUTPUT) Reference for the atol vector to fill up
   */
  virtual void setAtolVector(double atolDefault, const Epetra_Vector_Ghosted & soln, 
			     Epetra_Vector_Ghosted & atolVector,
			     const Epetra_Vector_Ghosted * const atolV = 0);

  //!  Fill the vector atolVector with the values from the DomainDescription for abs tol
  /*!
   * @param atolDefault             Default atol value
   * @param soln                    Solution vector. This is a constant
   *                                the residual calculation.
   * @param atolVector              (OUTPUT) Reference for the atol vector to fill up
   */
  virtual void setAtolVector_DAEInit(double atolDefault, const Epetra_Vector_Ghosted & soln, 
				     const Epetra_Vector_Ghosted & solnDot,
				     Epetra_Vector_Ghosted & atolVector_DAEInit,
				     const Epetra_Vector_Ghosted * const atolV = 0);

  //! Evaluates the atol vector used in the delta damping process.
  /*!
   *   @param relcoeff     Relative constant to multiply all terms by
   *   @param soln         current solution vector.
   *   @param atolDeltaDamping      If non-zero, this copies the vector into the object as input
   *                      The default is zero.
   */
  virtual void
  setAtolDeltaDamping(double atolDefault, double relcoeff, 
		      const Epetra_Vector_Ghosted & soln, 
		      Epetra_Vector_Ghosted & atolDeltaDamping,
		      const Epetra_Vector_Ghosted * const atolV = 0);
		    

 //! Evaluates the atol vector used in the delta damping process for the DAE problem
  /*!
   *   @param relcoeff     Relative constant to multiply all terms by
   *   @param soln         current solution vector.
   *   @param solnDot      Current solutionDot vector.
   *   @param atolDeltaDamping       If non-zero, this copies the vector into the object as input
   *                       The default is zero.
   */
  virtual void
  setAtolDeltaDamping_DAEInit(double atolDefault, double relcoeff, 
			      const Epetra_Vector_Ghosted & soln,
			      const Epetra_Vector_Ghosted & solnDot,
			      Epetra_Vector_Ghosted & atolDeltaDamping,
			      const Epetra_Vector_Ghosted * const atolV = 0);
  

private:
  void
  err(const char *msg);

public:
  std::vector<double> cellBoundVeloc_left;
  std::vector<double> cellBoundVeloc_right;

};

}

#endif /* M1D_TDGrowingFilm_dom1D_H_ */
