/*
 * $Id: m1d_PS_InfPlateBat.h 506 2013-01-07 22:43:59Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _M1D_PS_INFPLATEBAT_H
#define _M1D_PS_INFPLATEBAT_H

#include "zuzax/equilibrium.h"
#include "tok_input_util.h"
#include <string>
#include <vector>

#include "m1d_BEinput_base.h"
#include "m1d_ProblemStatement.h"

namespace BEInput
{
class BlockEntry;
}
namespace m1d
{

//! storage for Command file input
/*!
 * Complete problem statement
 *
 * This is the current command file specification
 *                       of the problem statement.
 */

class PS_InfPlateBat : public ProblemStatement {
public:
  //! Constructor
  PS_InfPlateBat();

  //! Destructor
  virtual
  ~PS_InfPlateBat();


  virtual void setup_input_pass3(BEInput::BlockEntry *cf);

  //!        DATA

  //! Type of the boundary condition specified on the cathode
  /*!
   *   0 specify the voltage
   *   1 specify the current
   */
  int cathodeBCType_;

  //! Specified current in the cathode
  /*!
   *  This is actually the current from the cathode into the electrolyte.
   *  Therefore, during a normal discharge operation of the battery, this will be a
   *  negative quantity.
   *
   *   Note, this is only relevant when voltageVarBCType_ = 1
   */
  double icurrDischargeSpecified_;

  //!
  double CathodeVoltageSpecified_;

};
//=====================================================================================================================
}
//=====================================================================================================================
#endif
//=====================================================================================================================

