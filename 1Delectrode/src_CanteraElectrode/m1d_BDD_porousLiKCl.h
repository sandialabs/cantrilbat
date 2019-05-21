/**
 * @file m1d_BDD_porousLiKCl.h
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_BDD_POROUSLIKCL_H_
#define M1D_BDD_POROUSLIKCL_H_

#include "m1d_BulkDomainTypes.h"
#include "m1d_BDD_porousFlow.h"

#include "zuzax/transport.h"      // transport properties
#include "zuzax/thermo.h"      // transport properties
#include "zuzax/thermo/IonsFromNeutralVPSSTP.h"  // ion properties
#include "zuzax/thermo/StoichSubstance.h"  // separator

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This class consists of multiple species diffusing in a time
//! dependent manner.  There is a net flow and a net electric current.
/*!
 *  This class is used to test the implementation
 */
class BDD_porousLiKCl : public BDD_porousFlow
{
public:

  //! Constructor
  /*!
   * This constructor constructs the bulk domain from a MultiPhase object.
   *
   * In the constructor, we have typically been laying out what the unknowns are
   * and what the equations are, that are solved within the domain.
   *
   */
  BDD_porousLiKCl(DomainLayout *dl_ptr, std::string domainName = ""); 

  //! Destructor
  virtual
  ~BDD_porousLiKCl();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  BDD_porousLiKCl(const BDD_porousLiKCl &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  BDD_porousLiKCl& operator=(const BDD_porousLiKCl &r);

  //! Read in the possible models for each domain
  /*!
   *  This procedure is done before the Equations anv variable list are set up.
   *  Needed information about what is possible is input here.
   *  We read the Zuzax ThermoPhase and transport object into DomainDescriptions here.
   *
   *   We loop over volume and then surface domains.
   */
  virtual void ReadModelDescriptions();

  //! Determine the list of Equations and Variables
  /*!
   *  This routine is responsible for setting the variables:
   *    - VariableNameList
   *    - EquationNameList
   */
  virtual void SetEquationsVariablesList();

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual BulkDomain1D *
  mallocDomain1D();

  //! This is done after the equations are set up
  /*!
   *  We loop over volume and then surface domains here.
   */
  virtual void
  DetermineConstitutiveModels();

  // --------------------------------------------------------------------------------------------

  //! Equation type and var type to apply them
  std::vector<VarType> EquationID;

  //! Pointer to the thermo object for the molten salt
  /*!
   *   This is a shallow pointer
   */
  Zuzax::IonsFromNeutralVPSSTP *ionicLiquidIFN_;

};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
