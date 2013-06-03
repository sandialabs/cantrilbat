/*
 * $Id: cell_input.h 504 2013-01-07 22:32:48Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _CELL_INPUT_H
#define _CELL_INPUT_H

#include "cantera/equilibrium.h"
#include "tok_input_util.h"
#include <string>
#include <vector>
/*
 *-----------------------------------------------------------------------------
 *
 * Include file containing constant declarations for inputs to 
 * mpequil
 *
 *-----------------------------------------------------------------------------
 */
#define MPEQUIL_MAX_NAME_LEN_P1 81
#define MPEQUIL_MAX_NAME_LEN    80

#define MPEQUIL_SUCCESS 0


//! storage for Command file input
/*!
 * This is the current command file specification
 *                       of the problem statement.
 */
class CELL_KEY_INPUT {
public:
  CELL_KEY_INPUT() ;
  ~CELL_KEY_INPUT();
  
  std::string anodeElectrodeFileName_;
  std::string cathodeElectrodeFileName_;
  bool doCathode_;
  double netCurrent_;
  double AnodeVoltage_;
  double CathodeVoltage_;
  double NetVoltage_;

  std::string Title_;
  void InitForInput();
  
};
extern CELL_KEY_INPUT CellO;


extern int cell_input(std::string commandFile);


#endif 
/*****************************************************************************/
