/*
 * $Id: InterfacialMassTransfer_input.h 507 2013-01-07 22:48:29Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _INTERFACIALMASSTRANSFER_INPUT_H
#define _INTERFACIALMASSTRANSFER_INPUT_H

#include "cantera/equilibrium.h"
#include "tok_input_util.h"


#include "PhaseList.h"
#include "ReactingSurDomain.h"


#include "BlockEntryGlobal.h"




#include "mdp_allo.h"
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

namespace Cantera {
  class Electrode;
}

//!  HKM -> getting rid of the Bath class
class ElectrodeBath {
public:
  //! species mole fractions in one phase
  double * XmolPLSpecVec;
  //! species mole fractions in each phase of phase list
  double **XmolPLPhases;
  //! species molalities in one phase
  double * MolalitiesPLSpecVec;
  //! species molalities in each phase of phase list
  double **MolalitiesPLPhases;

  //! total moles of given phase
  std::vector<double> PhaseMoles;
  //! total mass of given phase
  std::vector<double> PhaseMass;


  ElectrodeBath() :
    XmolPLSpecVec(0),
    XmolPLPhases(0),
    MolalitiesPLSpecVec(0),
    MolalitiesPLPhases(0)
  {
  }
  ~ElectrodeBath() {
    //  mdpUtil::mdp_safe_free((void **) &Xmol);
    mdpUtil::mdp_safe_free((void **) &XmolPLSpecVec);
    mdpUtil::mdp_safe_free((void **) &XmolPLPhases);
    mdpUtil::mdp_safe_free((void **) &MolalitiesPLSpecVec);
    mdpUtil::mdp_safe_free((void **) &MolalitiesPLPhases);
  
    //   mdpUtil::mdp_safe_free((void **) &BathSpeciesIDVec);
  }
} ;

namespace Cantera {


//! storage for Command file input
/*!
 * This is the current command file specification
 *                       of the problem statement.
 */
class IMT_KEY_INPUT {
public:
  //! Constructor
  IMT_KEY_INPUT();
  IMT_KEY_INPUT(const IMT_KEY_INPUT &right);
  ~IMT_KEY_INPUT();

  IMT_KEY_INPUT & operator=(const IMT_KEY_INPUT &right);

  //! Initialize the fields for reading
  /*!
   *
   */
  void InitForInput(const Cantera::PhaseList  * const pl);

  //! This is assigned to the first surface phase found in the PhaseList
  std::string CanteraFNSurface;

  int NumberCanteraFiles;
  char **CanteraFileNames;
  double Temperature;
  double PressureA;
 
  ElectrodeBath *m_BG;
  double *MoleNumber;
  double *MoleFraction;
 
  int    *PhaseInclude;
  int    ProblemType;
  char **SpeciesNames;
  char **PhaseNames;
  char **ElementNames;
  double *ElementAbundances;
  bool  specifiedElementAbundances;
  int   specifiedBlockKmolSpecies;

  //! Name of the phase that is on the left side of the interface
  std::string PhaseAName;

 //! Name of the phase that is on the right side of the interface
  std::string PhaseBName;

  int solnAIndex_;

  int solnBIndex_;

  //! Title
  std::string Title;

  //! IMT Model Name
  /*!
   * This is used in the factory routine
   */
  std::string IMT_ModelName;



  //! Electrode Area
  /*! 
   *   Area of the electrode to consider. This will be combined with the width
   *   and the porosity information to calculate a net amount of electrode material
   *
   * units = m**2
   */
  double SurfaceArea;
 

 
  //! Boundary Layer Thickness on the A side 
  /*! 
   *  
   *
   * units = m
   */
  double BLThickness_A;

  //! Boundary Layer Thickness on the B side 
  /*! 
   *  
   *
   * units = m
   */
  double BLThickness_B;


  //! species mole fractions in each phase of phase list
  std::vector<double> XmfPhaseA_;


  //! species mole fractions in each phase of phase list
  std::vector<double> XmfPhaseB_;

  int nTotPhases;


  //! Total number of species in the kinetic species list
  int nTotSpecies;

  //! Total number of elements
  int nTotElements;


  //! PhaseList object 
  /*!
   * this includes all of the phases, "period".
   *  In particular this includes the surface phases
   */
  Cantera::PhaseList *m_pl;

	
};


int imt_model_init(IMT_KEY_INPUT *ei, BEInput::BlockEntry *cf);


int imt_input(IMT_KEY_INPUT *input,  std::string commandFile,  BEInput::BlockEntry *cf);

bool process_imt_input(BEInput::BlockEntry *cf, std::string fileName, int printFlag);
}
#endif 
/*****************************************************************************/

