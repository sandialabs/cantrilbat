/**
 * @file  m1d_SDD_CathodeCollector.h   Definitions for the surface domain description of the cathode collector
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_SDD_CATHODECOLLECTOR_H_
#define M1D_SDD_CATHODECOLLECTOR_H_

#include "m1d_SDD_Mixed.h"

//----------------------------------------------------------------------------------------------------------------------------------
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
{
class ELECTRODE_MODEL;
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Class for the Surface domain description of the cathode collector
/*!
 *   This class specifies that all equations are handled
 *   by a simple Dirichlet condition or simple flux conditions of the third type
 */
class SDD_CathodeCollector : public SDD_Mixed
{
public:

  //! Default Constructor
  /*!
   *  We construct the object but don't actually specify any Dirichlet conditions.
   *  Later we can add dirichlet conditions into the object.
   *
   *  In the constructor, we have typically been laying out what the unknowns are
   *  and what the equations are, that are solved within the domain.
   *
   *  @param[in]             dl_ptr              Domain Layout object that owns this description.
   *  @param[in]             position            Position within the list of SurfaceDomains that comprise the problem statement 
   *                                               (this is currently unused and maybe not created correctly)
   *                                                -> used to signify the top or bottom of the domain as well
   *  @param[in]             domainName          Name of the surface domain
   *                                                 Default = ""
   */
  SDD_CathodeCollector(DomainLayout *dl_ptr, int position, const char *domainName = "");

  //! Destructor
  virtual ~SDD_CathodeCollector();

  //! Copy Constructor
  /*!
   *  @param[in]             r                   Object to be copied
   */
  SDD_CathodeCollector(const SDD_CathodeCollector &r);

  //! Assignment operator
  /*!
   *  @param[in]             r                   Object to be copied
   *  @return                                    Returns a changeable reference to the current object
   */
  SDD_CathodeCollector& operator=(const SDD_CathodeCollector &r);

  //! Set the equation description
  /*!
   *  (virtual from DomainDescription)
   *
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
   *  @return                                    Returns a pointer to the object that will calculate the residual
   *                                             efficiently
   */
  virtual SurDomain1D* mallocDomain1D() override;

  // --------------------------------------------------------------------------------------------------------------
  //                                   DATA
  // --------------------------------------------------------------------------------------------------------------

  //! Top or bottom of the domain
  /*!
   *   0 - top, right
   *   1 - bottom, left
   */
  int m_position;

  //! Type of the boundary condition specified on the cathode
  /*!
   *   0 specify the voltage
   *   1 specify the current
   *
   *  10 Set a Robin boundary condition ---   current = R_anodeCC (v_acc - vanode)
   */
  int voltageVarBCType_;

  //! Specified current in the cathode
  /*!
   *  This is actually the current from the cathode into the electrolyte.
   *  Therefore, during a normal discharge operation of the battery, this will be a
   *  negative quantity.
   *
   *   Note, this is only relevant when voltageVarBCType_ = 1, 3, 5, 7, 9
   */
  double icurrCathodeSpecified_;

  //! Specified voltage at the cathode and sometimes at the cathode current collector
  double voltageCathodeSpecified_;

  //!  Thickness of the cathode current collector
  /*!
   *   This is the thickness of the cathode current collector, when the collector has the cross-sectional area
   *   assumed for the rest of the battery. Another part of the collector will not be proportional to the cross-sectional
   *   area, but will be assued to be a "wire" attached between this part and load.
   *
   *   For some boundary conditions, there will be a voltage drop and a Joule heating term within the current 
   *   collector due to this finite thickness. It will depend on the resistivity of the current collector metal.
   *
   *   -> This is read in from the input file
   *
   *   Units:  m
   */
  double cathodeCCThickness_;

  //! Extra resistance attached to the entire battery that is still included in heat buildup
  /*!
   *    Note this is not cross section based. This is for the entire battery. Therefore, the
   *    cross section must be multiplied in.
   *
   *  -> This is read in from the input file
   *
   *    Units: ohms
   */
  double extraResistanceCathode_;

  //! Load resistance attached to the entire battery 
  /*!
   *    Note this is not cross-section based. This is for the entire battery. Therefore, the
   *    cross section must be multiplied in to connect with the rest of the problem.
   *
   *  -> This is read in from the input file
   *
   *    Units:   ohm
   */
  double ResistanceLoad_;

  //! Voltage load attached to the entire battery 
  /*!
   *  Voltage that is attached to the load. Acts as a battery element in the global circuit that
   *  the battery is attached to.
   *
   *  -> This is read in from the input file
   *
   *    Note this is not cross-section based. This is for the entire battery. Therefore, the
   *    cross section must be multiplied in to connect with the rest of the problem.
   *
   *    Units:   volt
   */
  double VoltageLoad_;

  //!  Type of the temperature boundary condition
  /*!
   *  There are three possibilities
   *     0 set the cathode temperature to a constant
   *     1 set the flux to a constant
   *    10 Set a Robin boundary condition - heatflux = h ( T - T_cath)
   *
   *  -> This is read in from the input file
   */
  int cathodeTempBCType_;
  double cathodeTempCollector_;
  double cathodeHeatTranCoeff_;


  //! Make the SurDomain1D class a friend so that it can access all of the stuff in this class
  friend class SurDomain_CathodeCollector;
};

}

#endif /*  */
