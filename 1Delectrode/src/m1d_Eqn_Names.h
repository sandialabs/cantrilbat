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
  Displacement_Axial,       //0     //!< Displacement_axial
  Velocity_Axial,          // 1
  Velocity_Radial,               //!< Velocity_radial
  Velocity_Swirl,                //!< Velocity_swirl
  Pressure_Axial,           // 4
  PressureGradient_Radial,  //5     //!< PressureGradient_radial
  Concentration_Species,     // adding in extra variable types
  MoleFraction_Species,     // 7
  VolumeFraction_Species,        //!< VolumeFraction_species
  VolumeFraction_Phase,          //!< VolumeFraction_phase
  Voltage,                 // 10
  Temperature,             // 11      //!< Temperature
#ifdef MECH_MODEL
  // \todo CBL I think these can be computed from VolumeFraction_Species
  //  Solid_Density, 
  //  Electrolyte_Density,
  //  Gas_Density,
  Solid_Stress_Axial, // scalar, for 1 d
  // CBL this is not something this solved for, it's an initial state
  //  IStress_Free_Strain_Axial, // this initial state may _not_ be uniform at t=0.
  // things like Bulk Modulus, Poisson's Ratio, Cv/Cp, thermal expansion coeffs
  //  are stored in hard coded functions, for anode/spacer/cathode materials   
  Solid_Stress_Transverse,
#endif // MECH_MODEL
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
#ifdef MECH_MODEL
  Mechanical_Model_Axial,
  Mechanical_Model_Transverse,
#endif
  Species_Eqn_Offset,   //! Special equation name representing the first species equation offset
  Max_Eqn_Name
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

inline void operator++(EQ_TYPE & eVal) {
  eVal = EQ_TYPE(eVal+1);
};

}
#endif

