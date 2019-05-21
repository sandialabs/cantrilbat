/**
 * @file m1d_BDD_porousLiKCl.h
 */
/*
 *   $Id: m1d_BDD_porousLiKCl.h 506 2013-01-07 22:43:59Z hkmoffa $
 */

#ifndef M1D_BDD_POROUSLIKCL_H_
#define M1D_BDD_POROUSLIKCL_H_

#include "m1d_BDD_porousFlow.h"

#include "zuzax/transport.h"      // transport properties
#include "zuzax/thermo.h"      // transport properties
#include "zuzax/thermo/IonsFromNeutralVPSSTP.h"  // ion properties
#include "zuzax/thermo/StoichSubstance.h"  // separator

namespace m1d
{

//==================================================================
//==================================================================
//==================================================================

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
  BDD_porousLiKCl(DomainLayout *dl_ptr);

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
  BDD_porousLiKCl &
  operator=(const BDD_porousLiKCl &r);

  virtual void SetEquationsVariablesList();

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual BulkDomain1D *
  mallocDomain1D();

  virtual void setupTransport();

  //! This is done after the equations are set up
  /*!
   *  We loop over volume and then surface domains at that point
   */
  virtual void
  DetermineConstitutiveModels();


  // --------------------------------------------------------------------------------------------

  //! Equation type and var type to apply them
  std::vector<VarType> EquationID;

  //! Pointer to the thermo object for the molten salt
  /*!
   *   We own this object
   */
  Zuzax::IonsFromNeutralVPSSTP *ionicLiquidIFN_;

  //! Pointer to the transport object for the molten salt
  /*!
   * We own this object
   */
  //Zuzax::Transport* trans_;

};

}

#endif /* M1D_BULKDOMAINTYPES_H_ */
