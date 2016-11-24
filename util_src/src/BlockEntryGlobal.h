/**
 * @file BlockEntryGlobal.h
 *    Entire list of external declarations for the BlockEntry system
 *    and also contains the blockentryModule module (see \ref blockentryModule).
 */
/*
 * $Id: BlockEntryGlobal.h 244 2012-07-03 17:54:30Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#ifndef BLOCKENTRYGLOBAL_H
#define BLOCKENTRYGLOBAL_H

/**
  @defgroup blockentryModule Block Entry Input
 
  These classes implement free form block
  entry of ascii input from an input deck
 
  Typically one can use multiple passes to process the
  input deck. This is needed for constitutive modeling packages, typically. 
 
  In order to use multiple passes, the input deck
  must be read into a file, and not just read in
  as a stream from stdin. A standard C file pointer
  is used to step through the input deck.
 
  A typical input deck might look like this
 
\verbatim
Number of Cantera Files = 1
Cantera File Name = airNASA9.xml

DebugPrinting = true
Bath Temperature = 6500.
Bath Pressure = 1.0 atm

START BLOCK Temperature Table Format
  Number of Points = 22
  Delta Temperature = 400.

  !----------------------------------------------------------
  ! Low Temperature = [double]
  !    (optional)
  !    (default = 300.)
  !    Low temperature for the temperature table
  !----------------------------------------------------------
  Low Temperature = 300.

END BLOCK Temperature Table Format

start block Bath specification for Phase airNASA9
  Bath Species ID = N2
  start block Bath Species Mole Fraction
    N2 = 0.7
    O2 = 0.3
  End block Bath Species Mole Fraction
end block Bath specification for Phase airNASA9
\endverbatim
 
  Dependencies amongst the cards are also handled. This allows for a fairly tight error checking 
  on the  setup of problems. We require a conversation between the user and the program to occur.
 
  Identical Multiple blocks are allowed to occur. For
  example if the particle initial distribution 
  can be described by 2 gaussians of varying distributions
  it's possible to specify each of the gaussians 
  using one block of input each. 

  The Current list of \link BEInput::LineEntry LE_LineEntry\endlink 
  derived types and their functionality is:

     - \link BEInput::LE_OneInt LE_OneInt\endlink Entry of a single integer value
     - \link BEInput::LE_OneDbl LE_OneDbl\endlink Entry
                of a single double value
     - \link BEInput::LE_OneDblUnits LE_OneDblUnits\endlink Entry
                of a single double value with an
                optional units converter
     - \link BEInput::LE_OneBool LE_OneBool\endlink Entry
                of a single boolean value stored as a bool
     - \link BEInput::LE_OneBoolInt LE_OneBoolInt\endlink Entry
                of a single boolean value stored as an int
     - \link BEInput::LE_OneStr LE_OneStr\endlink Entry
                of a single C++ string
     - \link BEInput::LE_OneCStr LE_OneCStr\endlink Entry
                of a single C-style string
     - \link BEInput::LE_PickList LE_PickList\endlink Entry
                of a List choice. The choice is entered
                as a string selection from a list. The
                choice is saved as an int value.
     - \link BEInput::LE_MultiCStr LE_MultiCStr\endlink Entry
                of multiple C-style string, saved as a
                vector of C-style strings.
     - \link BEInput::LE_StrListDbl LE_StrListDbl\endlink Entry
               of multiple related doubles into a single 
               vector of doubles. The rhs consists of a
               list token, which determines the vector index,
               and the double.
     - \link BEInput::LE_VecDbl LE_VecDbl\endlink Entry
               of a vector of doubles. The length of the
               vector is specified at construction.
     - \link BEInput::LE_VecDblVarLength LE_VecDblVarLength\endlink 
               Entry of a vector of doubles. The length of the
               vector is a variable.
     - \link BEInput::LE_StdVecDbl LE_StdVecDbl\endlink Entry
               of a std vector of doubles. The length of the
               vector is specified at construction.
     - \link BEInput::LE_StdVecDblVarLength LE_StdVecDblVarLength\endlink 
               Entry of a std vector of doubles. The length of the
               vector is a variable.


  The current list of \link BEInput::BlockEntry BlockEntry\endlink 
  derived types and their functionality is:

     - \link BEInput::BE_StrDbl BE_StrDbl\endlink 
             Specification of the properties of a single
             item all within a block format. This is essentially
             the same as LE_StrListDbl, as the results are
             written to a single vector of doubles, but the format is
             a block format.
     - \link BEInput::BE_MoleComp BE_MoleComp\endlink 
             Specification of the mole fractions of a phase
             all within a block format.
     - \link BEInput::BE_MolalityComp BE_MolalityComp\endlink 
             Specification of the molality of solute in a phase
             all within a block format.
     - \link BEInput::BE_MultiBlock BE_MultiBlock\endlink 
             Specification of a block in the input file that
             can be repeated an arbitray number of times
             in the input file. Special requirements for
             the external writing of information are specified.

   The current list of  \link BEInput::BE_UnitConversion BE_UnitConversion\endlink 
  derived types and their functionality is:

     - \link BEInput::BE_UnitConversionConcentration BE_UnitConversionConcentration\endlink 
             Converts units to the MKS value of kmol m-3. It checks
             to see that the units are appropriate for specification
             of the concentration
     - \link BEInput::BE_UnitConversionLength BE_UnitConversionLength\endlink 
             Converts units to the MKS value of m. It checks
             to see that the units are appropriate for specification
             of the length.
     - \link BEInput::BE_UnitConversionEnergy BE_UnitConversionEnergy\endlink 
             Converts units to the MKS value of Joules or kg m2 s-2. It checks
             to see that the units are appropriate for specification
             of the energy
     - \link BEInput::BE_UnitConversionPressure BE_UnitConversionPressure\endlink 
             Converts units to the MKS value of pascals. It checks
             to see that the units are appropriate for specification
             of the pressure.
     - \link BEInput::BE_UnitConversion BE_UnitConversion\endlink 
             Converts units to the MKS value. No unit checks
             are carried out.
    
   The user may gather the parsed information via 
   multiple ways. One way is that each of the cards
   allows for specification of where to write its
   output via an external address. 
   The other way allows for the code to read the 
   contents of what was input from the objects
   themselves, in a postprocessing fashion. 
  
 */

//
// List of individual Line Entry classes
//
#include "LE_MultiCStr.h"
#include "LE_OneBool.h"
#include "LE_OneBoolInt.h"
#include "LE_OneCStr.h"
#include "LE_OneDbl.h"
#include "LE_OneDblUnits.h"
#include "LE_OneInt.h"
#include "LE_OneSizet.h"
#include "LE_OneStr.h"
#include "LE_PickList.h"
#include "LE_StrListDbl.h"
#include "LE_VecDbl.h"
#include "LE_VecDblVarLength.h"
#include "LE_StdVecDbl.h"
#include "LE_StdVecDblVarLength.h"
//
//  List of individual Block Entry Classes
//
#include "BE_MoleComp.h"
#include "BE_MultiBlock.h"
#include "BE_MultiBlockNested.h"
#include "BE_StrDbl.h"
#include "BE_MolalityComp.h"
#include "BE_StrVecDbl.h"
#include "BE_MoleComp_VecDbl.h"
#include "BE_MultiBlockVec.h"
//
// List of individual Units conversion classes
//
#include "BE_UnitConversionConcentration.h"
#include "BE_UnitConversionLength.h"
#include "BE_UnitConversionPressure.h"
//
// List of Dependency Classes
//
#include "BI_DepIntMaxMin.h"
//
//  Definitions for error class input
//
#include "BI_InputError.h"

#endif
