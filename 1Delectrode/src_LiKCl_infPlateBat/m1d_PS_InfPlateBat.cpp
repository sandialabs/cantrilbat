/*
 * $Id: m1d_PS_InfPlateBat.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"

#include "m1d_globals.h"
#include "m1d_PS_InfPlateBat.h"

#include "Epetra_Comm.h"

#include "BlockEntryGlobal.h"

using namespace std;
using namespace BEInput;
using namespace TKInput;

namespace m1d
{
//=====================================================================================================================
/**************************************************************************
 *
 *  MPEQUIL_INPUT: constructor
 *
 *  We initialize the arrays in the structure to the appropriate sizes.
 *  And, we initialize all of the elements of the arrays to defaults.  
 */

PS_InfPlateBat::PS_InfPlateBat() :
  ProblemStatement(), cathodeBCType_(0), icurrDischargeSpecified_(1.0), CathodeVoltageSpecified_(1.9)
{

}
//=====================================================================================================================
PS_InfPlateBat::~PS_InfPlateBat()
{

}
//=====================================================================================================================
void
PS_InfPlateBat::setup_input_pass3(BlockEntry *cf)
{

  ProblemStatement::setup_input_pass3(cf);

  int reqd = 0;

  /* -------------------------------------------------------------------------
   *
   * Cathode BC Type - int (required)
   *     0 - Specify a fixed voltage at the cathode
   *     1 - Specify a fixed discharge current through the battery
   */
  reqd = 1;
  LE_OneInt *i2 = new LE_OneInt("Cathode BC Type", &(cathodeBCType_), reqd, "cathode_bc_type");
  i2->set_default(0);
  i2->set_limits(1, 0);
  cf->addLineEntry(i2);

  /* -------------------------------------------------------------------------
   *
   * Discharge Current - double [Amps / cm2] (no default)
   * (conditionally required if Problem Type = 1)
   *
   * Specify the current in amps per cm2 for the battery.
   * This is only needed when the Problem Type is 1
   */
  LE_OneDbl *d1 = new LE_OneDbl("Discharge Current", &(icurrDischargeSpecified_), 0, "icurrDischargeSpecified");
  d1->set_default(1.0);
  d1->set_limits(1.0E3, 0.0);
  cf->addLineEntry(d1);

  /* -------------------------------------------------------------------------
   *
   * Cathode Voltage - double [volts] (no default) (required)
   *
   * Specify the initial cathode voltage. Note the anode voltage is
   * always set to 0.0.
   * For a problem type of 0, this is the Dirichlet condition on the
   * cathode voltage. For a problem type of 1, this value is used to
   * set the initial condition of the cathode voltage.
   */
  LE_OneDbl *d2 = new LE_OneDbl("Cathode Voltage", &(CathodeVoltageSpecified_), 0, "icurrCathodeSpecified");
  d2->set_default(1.9);
  d2->set_limits(2.7, 0.1);
  cf->addLineEntry(d2);


  BaseEntry::set_SkipUnknownEntries(false);
}
//=====================================================================================================================
}
//=====================================================================================================================
