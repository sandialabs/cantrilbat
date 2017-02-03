/**
 * @file m1d_Eqn_Names.h
 *
 */


/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_EQN_NAMES_H
#define M1D_EQN_NAMES_H

#include "m1d_defs.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! Enum representing the types of variables
/*!
 *    Right now, there is a 1-1 correspondence between the variable number and the
 *    usual equation number.
 */
enum Var_Name_Enum
{
  Variable_Type_Unspecified = -2,   //!< Variable_Type_Unspecified
  Variable_Type_Any,                //!< Variable_Type_Any
  Displacement_Axial,       //0     //!< Displacement_axial
  Velocity_Axial,          // 1     //!< Axial velocity
  Velocity_Radial,                  //!< Velocity_radial
  Velocity_Swirl,                   //!< Velocity_swirl
  Pressure_Axial,           // 4    //!< Pressure, in the axial direction, Thermodynamic pressure has been subtracted out
  PressureGradient_Radial,  //5     //!< PressureGradient_radial
  Concentration_Species,            //!< Concentration of species (kmol m-3)
  MoleFraction_Species,     // 7    //!< Mole fraction of species
  VolumeFraction_Species,           //!< VolumeFraction_species
  VolumeFraction_Phase,             //!< VolumeFraction_phase
  Voltage,                 // 10    //!< Voltage
  Temperature,             // 11    //!< Temperature
  // re-use Displacement_Axial
  Solid_Stress_Axial, // scalar, for 1 d
  Max_Var_Name                      //!< Max_Var_Name, must be last in the list
};

//==================================================================================================================================
//! Enum representing the types of equations
/*!
 *  Below are the types of equations that are defined in all of the child classes
 */
enum EQ_Name_Enum
{
  Equation_Type_Unspecified = -2, 
  Equation_Type_Any,
  Momentum_Axial,            // 0
  Momentum_Radial,
  Momentum_Swirl,
  MeshDisplacement_Axial,
  Enthalpy_Conservation,
  Thermal_Conservation,
  Thermal_Dirichilet,
  Continuity,                // 5
  Continuity_Global,
  Species_Conservation,        // 7
  VolumeFraction_Conservation,
  PhaseVolumeFraction_Conservation,
  Current_Conservation,        // 10
  MoleFraction_Summation,      // 11
  ChargeNeutrality_Summation,   // 12
  Current_Specification,
  Voltage_Specification,
  Dirichlet_Specification,
  Mechanical_Model_Axial,
  Mechanical_Stress_Axial,
  Species_Eqn_Offset,   //! Special equation name representing the first species equation offset
  Max_Eqn_Name
};
//==================================================================================================================================
//! Typedef so that you don't need to type enum all the time for the Variable type
typedef enum Var_Name_Enum VAR_TYPE;

//! Typedef so that you don't need to type enum all the time for the equation type
typedef enum EQ_Name_Enum EQ_TYPE;
//==================================================================================================================================
//! Subtype of the equation type is defined as an integer
/*!
 *  This is used to differentiate multiple instances of the same equation type.
 *  For example the Species_Conservation equation type is differentiated by the species id.
 */
typedef int EQ_TYPE_SUBNUM;
//==================================================================================================================================
//! Subtype of the variable type is defined as an integer
/*!
 *  This is used to differentiate multiple instances of the same variable type.
 *  For example the MoleFraction_Species type is differentiated by the species id.
 */
typedef int VAR_TYPE_SUBNUM;
//==================================================================================================================================
//! Increment overload definition for the VAR_TYPE enum
/*!
 *  This is needed so that we use the VAR_TYPE variable in for loops
 *
 *  @param[in,out]           varT                Reference to a VAR_TYPE variable that will be incremented by one.
 */
inline void operator++(VAR_TYPE& varT) 
{
  varT = VAR_TYPE(varT + 1);
}
//==================================================================================================================================
//! Increment overload definition for the EQ_TYPE enum
/*!
 *  This is needed so that we use the EQ_TYPE variable in for loops
 *
 *  @param[in,out]           eqnT                Reference to a EQ_TYPE variable that will be incremented by one.
 */
inline void operator++(EQ_TYPE & eqnT)
{
  eqnT = EQ_TYPE(eqnT + 1);
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

