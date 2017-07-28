/*
 * m1d_Domain1D.h
 *
 *  Created on: May 19, 2009
 *      Author: hkmoffa
 */

#ifndef M1D_SIMPLEDIFF_DOM1D_H_
#define M1D_SIMPLEDIFF_DOM1D_H_

//! This is a heavyweight base class that provides the function
//! evaluation for a single bulk domain.

#include "m1d_DomainDescription.h"
#include "m1d_BulkDomain1D.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
class LocalNodeIndices;

//==================================================================================================================================
//! Class that does a simple diffusion operator on the of the bulk domain class
class SimpleDiff_dom1D : public BulkDomain1D {

public:

  //! Constructor
  /*!
   * @param[in]              bdd_ptr             Contains the bulk domain description.
   */
  SimpleDiff_dom1D(m1d::BulkDomainDescription* bdd_ptr);

  //! Copy constructor
  /*!
   *  @param r      Object to be copied into the current object
   */
  SimpleDiff_dom1D(const SimpleDiff_dom1D &r);

  //! Destructor
  virtual ~SimpleDiff_dom1D();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  SimpleDiff_dom1D& operator=(const SimpleDiff_dom1D &r);

  //! Prepare all of the indices for fast calculation of the residual
  /*!
   *  Ok, at this point, we will have figured out the number of equations
   *  to be calculated at each node point. The object NodalVars will have
   *  been fully formed.
   *
   *  We use this to figure out what local node numbers/ cell numbers are
   *  needed and to set up indices for their efficient calling.
   *
   *  Child objects of this one will normally call this routine in a  recursive fashion.
   *
   *  @param[in]             li_ptr              Pointer to the LocalNodeIndices
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
   *  @param[in]             res                 Output vector containing the residual
   *  @param[in]             doTimeDependentResid  Boolean indicating whether the time
   *                                               dependent residual is requested
   *  @param[in]             soln_ptr            Solution vector at which the residual should be evaluated
   *  @param[in]             solnDot_ptr         Solution dot vector at which the residual should be evaluated.
   *  @param[in]             solnOld_ptr         Solution vector at the old time
   *  @param[in]             t                   time
   *  @param[in]             rdelta_t            inverse of delta_t
   *  @param[in]             residType           Enum ResidEval_Type_Enum describing the type of residual
   *                                                Defaults to the Base_ResidEval type
   *  @param[in]             solveType           Enum Solven_Type_Enum descring the type of solve.
   *                                                Defaults to the TimeDependentAccurate_Solve
   */
  virtual void
  residEval(Epetra_Vector &res, const bool doTimeDependentResid,
            const Epetra_Vector* const soln_ptr, const Epetra_Vector* const solnDot_ptr,
            const Epetra_Vector* const solnOld_ptr, const double t, const double rdelta_t,
            const ResidEval_Type_Enum residType = Base_ResidEval,
	    const Solve_Type_Enum solveType = TimeDependentAccurate_Solve) override;

  //! Base class for saving the solution on the domain in an xml node.
  /*!
   *
   *  @param[in]             oNode               Reference to the XML_Node
   *  @param[in]             soln_GlAll_ptr      Pointer to the Global-All solution vector
   *  @param[in]             solnDot_GlAll_ptr   Pointer to the time derivative of the Global-All solution vector
   *  @param[in]             t                   time
   *  @param[in]             duplicateOnAllProcs  If this is true, all processors will include the same XML_Node information
   *                                              as proc 0. If false, the xml_node info will only exist on proc 0.
   */
  virtual void
  saveDomain(ZZCantera::XML_Node& oNode, const Epetra_Vector* const soln_GlAll_ptr,
             const Epetra_Vector* const solnDot_GlAll_ptr, const double t, bool duplicateOnAllProcs = false) override;

private:
  //! Error message handler
  /*!
   *  @param[in]             msg                 Message
   */
  void err(const char *msg);

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif /* M1D_DOMAIN1D_H_ */
