/**
 * @file m1d_porousFlow_dom1D.h
 */

/*
 *   $Id: m1d_porousFlow_dom1D.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_POROUSFLOW_DOM1D_H_
#define M1D_POROUSFLOW_DOM1D_H_

#include "m1d_BulkDomain1D.h"

#include <cantera/transport.h>    



namespace m1d
{
class LocalNodeIndices;

//======================================================================================================================
//! This is derived class  provides the function
//! evaluation for a porous electrode.
/*!
 * The porous electrolyte domain is characterized by a 
 * current conservation equation and several species 
 * conservation equations describing the electrolyte.
 * A porosity/tortuosity is also associated with the domain.
 *
 * There is a 1 to 1 mapping between the local control volume indexing and
 * the Local Consecutive Ordering indexing
 *
 * There is a 1 to 1 mapping between the global control volume indexing
 * and the Global node number indexing that is given by a single offset.
 */
class porousFlow_dom1D : public BulkDomain1D
{

public:

  //! Constructor
  /*!
   * @param bdd   Contains the bulk domain description.
   */
  porousFlow_dom1D(m1d::BulkDomainDescription &bdd);

  //! Copy constructor
  /*!
   * @param r      Object to be copied into the current object
   */
  porousFlow_dom1D(const porousFlow_dom1D &r);

  //! Destructor
  virtual  ~porousFlow_dom1D();

  //! Assignment operator
  /*!
   * @param r      Object to be copied into the current object
   * @return       Returns a changeable reference to the current object
   */
  porousFlow_dom1D&
  operator=(const porousFlow_dom1D&r);

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


  // -------------------------------------------------------------------------------------------------------------------
  // --- DATA 
  // -------------------------------------------------------------------------------------------------------------------

  // ------------------- Thermodynamics quantities on the domain -------------------------------------------------------

  //!  Partial molar Heat Capacity  of the electrolyte species located in all of the cells
  /*!
   *   Vector of partial molar heat capacity const press (KRSpecies, iCell)
   *   Units of Joules/(kmol K)
   */
  std::vector<doublereal> CpPM_lyte_Cell_;

  //!  Partial molar Enthalpy  of the electrolyte species located in all of the cells
  /*!
   *   Vector of partial molar enthalpy  (KRSpecies, iCell)
   *   Units of Joules/(kmol)
   */
  std::vector<doublereal> EnthPM_lyte_Cell_;




  // ------------------- Porosity of the Domain -----------------------------------------------------------------------

protected:
  //! Volume Fraction of the electrolyte within each control volume
  /*!
   * (change to CV)
   *  Length is number of cells on the processor.
   */
  std::vector<double> porosity_Cell_;

  //! Volume Fraction of the electrolyte within the cell at the previous time step
  /*!
   *  Length is number of cells on the processor.
   */
  std::vector<double> porosity_Cell_old_;

  

};
//======================================================================================================================
}
//======================================================================================================================
#endif 
