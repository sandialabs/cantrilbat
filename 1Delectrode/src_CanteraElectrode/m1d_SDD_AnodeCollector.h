/**
 * @file m1d_SDD_AnodeCollector.h  Definitions for the Surface demain description class that handles the boundary
 *           conditions on the anode current collector
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#ifndef M1D_SDD_ANODECOLLECTOR_H_
#define M1D_SDD_ANODECOLLECTOR_H_

#include "m1d_SDD_Mixed.h"

namespace Zuzax
{
class ELECTRODE_MODEL;
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

//==================================================================================================================================
//! This class specifies that all equations are handled by a simple Dirichlet condition
/*!
 *
 */
class SDD_AnodeCollector : public SDD_Mixed
{
public:

  //! Constructor
  /*!
   *  We construct the object but don't actually specify any Dirichlet conditions.
   *  Later we can add dirichlet conditions into the object.
   *
   *  In the constructor, we have typically been laying out what the unknowns are
   *  and what the equations are, that are solved within the domain.
   *
   *  @param[in]             dl_ptr              Domain Layout object that owns this description.
   *  @param[in]             position            Position within the domain
   *  @param[in]             domainName          String Name of the domain
   */
  SDD_AnodeCollector(DomainLayout *dl_ptr, int position, const std::string& domainName = "");

  //! Destructor
  virtual ~SDD_AnodeCollector();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  SDD_AnodeCollector(const SDD_AnodeCollector &r);

  //! Assignment operator
  /*!
   *  @param[in]             r                   Object to be copied
   *  @return     Returns a changeable reference to the current object
   */
  SDD_AnodeCollector& operator=(const SDD_AnodeCollector &r);

  //! Set the equation description
  /*!
   *  This routine is responsible for setting the variables:
   *    - NumEquationsPerNode
   *    - VariableNameList
   *    - EquationNameList
   *    - EquationIndexStart_EqName
   */
  virtual void SetEquationDescription() override;

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   *
   * @return                                     Returns a pointer to the object that will calculate the residual
   *                                             efficiently
   */
  virtual SurDomain1D * mallocDomain1D() override;

  // --------------------------------------------------------------------------------------------------------------
  //                                   DATA
  // --------------------------------------------------------------------------------------------------------------

  //! top or bottom of the domain
  /*!
   *   0 - top, right
   *   1 - bottom, left
   */
  int m_position;

  //!  Type of the voltage boundary condition
  /*!
   *   There are two possibilities
   *
   *     0  Set the anode voltage to zero
   *    10  Set a Robin boundary condition ---   current = R_anodeCC (v_acc - vanode)
   */
  int voltageVarBCType_;

  //!  Type of the temperature boundary condition
  /*!
   *  There are three possibilities
   *     0 set the cathode temperature to a constant
   *     1 set the flux to a constant
   *    10 Set a Robin boundary condition - heatflux = h ( T - T_cath) 
   */
  int anodeTempBCType_;

  //! Temperature reference of the anode, as set by the input file
  /*!
   *  This is the temperatue to be used in temperature boundary conditions, not necessarily the actual temperature
   *  Units = Kelvin
   */
  double anodeTempRef_;

  //! Heat transfer coefficient for the anode collector per cross-sectional area
  /*!
   *   Units:   Watts m-2 K-1
   */
  double anodeHeatTranCoeff_;

  //!  Thickness of the anode current collector
  double anodeCCThickness_;

  //! Make the SurDomain1D class a friend so that it can access all of the stuff in this class
  friend class SurDomain_AnodeCollector;
};

}

#endif /*  */
