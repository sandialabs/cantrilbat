/*
 * $Id: m1d_ProblemStatementCell.h 593 2013-05-13 21:25:47Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _M1D_PROBLEMSTATEMENTCELL_H
#define _M1D_PROBLEMSTATEMENTCELL_H

#include "cantera/equilibrium.h"
#include "cantera/thermo.h"
#include "PhaseList.h"
#include "importPL.h"
#include "tok_input_util.h"

#include "BlockEntryGlobal.h"

#include <string>
#include <vector>

#include "m1d_ProblemStatement.h"
#include "m1d_BoundaryCondition.h"
#include "m1d_exception.h"

#include "Electrode_input.h"
#include "Electrode.h"
#include "Electrode_Factory.h"

#define MPEQUIL_MAX_NAME_LEN_P1 81
#define MPEQUIL_MAX_NAME_LEN    80

namespace BEInput
{
class BlockEntry;
}
namespace m1d
{

//! Storage for Command file input
/*!
 * Complete problem statement
 *
 * This is the current command file specification
 *                       of the problem statement.
 */

class ProblemStatementCell : public ProblemStatement {
public:
  //! Constructor
  ProblemStatementCell();

  //! Destructor
  virtual
  ~ProblemStatementCell();

  /** 
   * The first pass through the input file determines
   *  the number of Cantera files.
   */
  virtual void setup_input_pass1(BEInput::BlockEntry *cf);

  /** 
   * The second pass through the input file parses the Cantera files.
   */
  virtual void setup_input_pass2(BEInput::BlockEntry *cf);

  /** 
   * Provide subclass specific parsing information.
   * This primarily fills the data structures defined in 
   * this subclass.
   */
  virtual void setup_input_pass3(BEInput::BlockEntry *cf);

   //! other preparation steps
  /**
   * This processes the phases in the Cantera input files,
   * fills the PhaseList_ object and other auxiliary data like
   * the numbers and names of species, elements and phases.
   */
  virtual void InitForInput();

  /**
   * Do any post processing required.
   * This might include unit conversions, opening files, etc.
   */
  virtual void post_process_input();


  virtual void readAnodeInputFile(Electrode_Factory *f = 0);
  virtual void readCathodeInputFile(Electrode_Factory *f = 0);

  //!  Test whether the anode and the cathode are compatible
  /*!
   *    The only test so far is to determine if the electrode areas input from
   *    the anode and the cathode are the same.
   *
   *    @return true if they are. False if they are not compatible.
   */
  bool AnodeCathodeCompatibility();

  //!        DATA

  //! Number of cantera files that will be used
  int NumberCanteraFiles;

  //! Vector containing the names of the cantera files
  char ** CanteraFileNames;

  //! type of the boundary condition specified on the current of the anode
  /*!
   *
   *   0 - Specify a fixed voltage of zero at the anode
   *  10 - anode Collector - robin boundary condition on the current
   */
  int anodeBCType_;

  //! Type of the boundary condition specified on the cathode
  /*!
   *     0 - Specify a fixed voltage at the cathode
   *     1 - Specify a fixed discharge current through the battery
   *     2 - Time dependent voltage at the cathode
   *     3 - Time dependent current at the cathode
   *     4 - specify time dependent voltage BoundaryCondition BCconstant
   *     5 - specify time dependent current BoundaryCondition BCconstant
   *     6 - specify time dependent voltage BoundaryCondition BCsteptable
   *     7 - specify time dependent current BoundaryCondition BCsteptable
   *     8 - specify time dependent voltage BoundaryCondition BClineartable
   *     9 - specify time dependent current BoundaryCondition BClineartable
   *    10 - Cathode Collector - robin boundary condition
   */
  int cathodeBCType_;

  //! Provides XML formatted input for the boundary condition specified on the cathode.  Only relevant for cathodeBCType_ = 6, 7, 8 or 9. 
  /*!
   * For cathodeBCType_ = 6 or 7 the voltage or current 
   * is specified as a step function with value given 
   * held until the next time value.  See class BCsteptable.
   *
   * For cathodeBCType_ = 8 or 9 the voltage or current 
   * is linearly interpolated between time-value pairs.  
   * See class BClineartable.
   */
  std::string cathodeBCFile_;

  //! Specified current in the cathode
  /*!
   *  This is actually the current from the cathode into the electrolyte.
   *  Therefore, during a normal discharge operation of the battery, this will be a
   *
   *   Note, this is only relevant when voltageVarBCType_ = 1
   */
  double icurrDischargeSpecified_;

  //! Specify the voltage at the cathode
  double CathodeVoltageSpecified_;

  //!  Treatment of the constant current condition
  /*!
   *   0  Default treatment is to assign a flux boundary condition to the current collector
   *   1  Constant current is calculated as a root finder algorithm using a constant voltage inner
   *      loop
   *   2  Hybrid between the two.  we sense when we need to switch ?!?
   */
  int rootFinderForConstantCurrent_;

  //! Function Pointers for time dependent BC for BC_Type = 2, 3
  //! NOTE THAT THIS MAY CHANGE WHEN WE IMPLEMENT MORE GENERAL BC
  /*!
   *   Vector has length equal to the number of equations defined at the node
   */
  double (*TimeDepFunction_)(double time);

  //! BoundaryCondition Function Pointers for time dependent BC for BC_Type = 4-9
  //! NOTE THAT THIS MAY CHANGE WHEN WE IMPLEMENT MORE GENERAL BC
  /*!
   *   Vector has length equal to the number of equations defined at the node
   */
  BoundaryCondition *BC_TimeDep_;

  //! Name of Anode input file
  std::string anodeFile_;

  //! Name of Cathode input file
  std::string cathodeFile_;

  //! Name of Electrolyte Cantera phase
  std::string electrolytePhase_;

  //! Mole fractions of electrolyte
  double* electrolyteMoleFracs_;

  //! Name of Separator species
  std::string separatorPhase_;

  //! Separator species mass 
  double separatorMass_;

  //! Separator area (m2)
  double separatorArea_;

  //! Separator thickness (m)
  double separatorThickness_;

  //! Separator diameter (m)
  double separatorDiameter_;

  //! Electrical conductivity of the anode
  /*!
   *    units are S/m
   */
  double conductivityAnode_;

  //! Electrical conductivity of the cathode
  /*!
   *    units are S/m
   */
  double conductivityCathode_;

  //! Thickness of anode current collector (m)
  double anodeCCThickness_;

  //! Thickness of cathode current collector (m)
  double cathodeCCThickness_;

  //! Extra resistance put in series with the cathode (ohms)
  /*!
   *  Note the effective resistance for the battery stack is
   *
   *     extraCathodeResistance_ * crossSectionalArea_
   */
  double extraCathodeResistance_;

  //! Flag to use Dakota I/O
  bool useDakota_;

  //! Maximum number of subgrid time steps per global time step
  int  maxSubgridTimeStepsPerGlobalTimeStep_;
  
  //! file name to receive parameters from Dakota
  //! (default =  'params.in')
  std::string fromDakotaFileName_;

  //! file name to send results to Dakota
  //! (default =  'results.out')
  std::string toDakotaFileName_;

  ////////////////////////////////////////////////////
  //    Member data derived from other members
  //    (not directly specified in input file)

  //! Master PhaseList object 
  Cantera::PhaseList *PhaseList_;

  //! total number of phases in PhaseList_
  int nTotPhases_;

  //! total number of species in PhaseList_
  int nTotSpecies_;
  
  //! Total number of elements
  int nTotElements_;

  //!

  char **SpeciesNames_;
  char **PhaseNames_;
  char **ElementNames_;


  //! Pointer to input parameters for anode and cathode model
  /**
   * This object will carry some data that is needed for the
   * electrolyte transport algorithms including information
   * required to compute the porosity
   */
  ELECTRODE_KEY_INPUT *anode_input_;
  ELECTRODE_KEY_INPUT *cathode_input_;

  //! Initial number of nodes
  int initDefaultNumCVsAnode_;

  //! Initial number of nodes
  int initDefaultNumCVsCathode_;

  //! Initial number of cells in the cathode
  int initDefaultNumCVsSeparator_;

  int doHeatSourceTracking_;

};
//=====================================================================================================================
}
//=====================================================================================================================
#endif //_M1D_PROBLEMSTATEMENTCELL_H
//=====================================================================================================================

