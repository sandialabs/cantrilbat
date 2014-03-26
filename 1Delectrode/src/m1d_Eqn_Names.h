/**
 * @file m1d_Eqn_Names.h
 *
 */

/*
 *  $Id: m1d_Eqn_Names.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_EQN_NAMES_H
#define M1D_EQN_NAMES_H

namespace m1d
{

//! Enum representing the types of variables
/*!
 *    Right now, there is a 1-1 correspondence between the variable number and the
 *    usual equation number.
 */
enum Var_Name_Enum
{
  Variable_Type_Unspecified = -2,//!< Variable_Type_Unspecified
  Variable_Type_Any,             //!< Variable_Type_Any
  Velocity_Axial,          // 0
  Velocity_Radial,               //!< Velocity_radial
  Velocity_Swirl,                //!< Velocity_swirl
  Displacement_Axial,            //!< Displacement_axial
  Temperature,                   //!< Temperature
  Pressure_Axial,           // 5
  PressureGradient_Radial,       //!< PressureGradient_radial
  MoleFraction_Species,     // 7
  VolumeFraction_Species,        //!< VolumeFraction_species
  VolumeFraction_Phase,          //!< VolumeFraction_phase
  Voltage,                 // 10
  Concentration_Species,     // adding in extra variable types
  Max_Var_Name                   //!< Max_Var_Name
//! must be last in the list
};

//! Enum representing the types of equations
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
  Max_Eq_Name
//! must be last in the list
};

//! Typedef so that you don't need to type enum all the time for the Variable type
typedef enum Var_Name_Enum VAR_TYPE;

//! Typedef so that you don't need to type enum all the time for the equation type
typedef enum EQ_Name_Enum EQ_TYPE;

typedef int EQ_TYPE_SUBNUM;
typedef int VAR_TYPE_SUBNUM;

inline void operator++(VAR_TYPE & eVal) {
  eVal = VAR_TYPE(eVal+1);
};

}
#endif

