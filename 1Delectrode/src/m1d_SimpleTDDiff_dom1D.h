/*
 * m1d_SimpleTDDiff_dom1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_SIMPLETDDIFF_DOM1D_H_
#define M1D_SIMPLETDDIFF_DOM1D_H_

//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_DomainDescription.h"
#include "m1d_BulkDomain1D.h"

namespace m1d
{
class LocalNodeIndices;

//! Class that does a simple diffusion operator on the
//! of the bulk domain class
class SimpleTDDiff_dom1D : public BulkDomain1D {

public:

  //! Constructor
  /*!
   * @param bdd   Contains the bulk domain description.
   */
  SimpleTDDiff_dom1D(m1d::BulkDomainDescription &bdd);

  //! Copy constructor
  /*!
   * @param r      Object to be copied into the current object
   */
  SimpleTDDiff_dom1D(const SimpleTDDiff_dom1D &r);

  //! Destructor
  virtual
  ~SimpleTDDiff_dom1D();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  SimpleTDDiff_dom1D &
  operator=(const SimpleTDDiff_dom1D &r);

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

  //! Base class for saving the solution on the domain in an xml node.
  /*!
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
  virtual void
  saveDomain(Cantera::XML_Node& oNode,
             const Epetra_Vector *soln_GlAll_ptr,
             const Epetra_Vector *solnDot_GlAll_ptr,
             const double t,
             bool duplicateOnAllProcs = false);

private:
  void
  err(const char *msg);

};

}

#endif /* M1D_DOMAIN1D_H_ */
