/**
 * @file m1d_Eqn_Names.h   File that contains the Variable types and the equation type enums for the m1d problem class
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
enum Var_Name_Enum {
    Variable_Type_Unspecified = -2,   //!< Variable_Type_Unspecified
    Variable_Type_Any,                //!< Variable_Type_Any
    Displacement_Axial,      //0      //!< Displacement_axial
    Velocity_Axial,          //1      //!< Axial velocity
    Velocity_Radial,                  //!< Velocity_radial
    Velocity_Swirl,                   //!< Velocity_swirl
    Pressure_Axial,          //4      //!< Pressure, in the axial direction, Thermodynamic pressure has been subtracted out
    PressureGradient_Radial, //5      //!< PressureGradient_radial
    Concentration_Species,            //!< Concentration of species (kmol m-3)
    MoleFraction_Species,    //7      //!< Mole fraction of species
    VolumeFraction_Species,           //!< VolumeFraction_species
    VolumeFraction_Phase,             //!< VolumeFraction_phase
    Voltage,                 //10     //!< Voltage
    Temperature,             //11     //!< Temperature
    Solid_Stress_Axial,      //12     //!< Value of the stress in the axial direction, scalar, for 1 d
    Max_Var_Name                      //!< Max_Var_Name, must be last in the list
};

//==================================================================================================================================
//! Enum representing the types of equations
/*!
 *  Below are the types of equations that are defined in all of the child classes
 */
enum EQ_Name_Enum {
    Equation_Type_Unspecified = -2,  //!< Equation type is unspecified -> used to initialize before setting
    Equation_Type_Any,               //!< Equation type matches any equation type
    Momentum_Axial,           //0    //!< Axial momentum equation
    Momentum_Radial,                 //!< Radial momentum equation
    Momentum_Swirl,                  //!< Swirl equation
    MeshDisplacement_Axial,          //!< Mesh displacement equation in the axial direction
    Enthalpy_Conservation,           //!< Enthalpy conservation equation
    Thermal_Conservation,            //!< Temperature equation
    Thermal_Dirichilet,              //!< Temperature dirichlet condition
    Continuity,               //5    //!< Total continuity equation for mass or moles
    Continuity_Global,               //!< Global total continuity equation of mass or moles
    Species_Conservation,     //7    //!< Conservation equation for a single species, mass or moles
    VolumeFraction_Conservation,     //!< Conservation equation for the volume fraction of a single species
    PhaseVolumeFraction_Conservation, //!< Conservation equation for the volume fractino of a single phase
    Current_Conservation,       //10  //!< Conservation equation of the current
    MoleFraction_Summation,     //11  //!< Mole fraction summation equals one equation
    ChargeNeutrality_Summation, //12 //!< Charge neutrality equation
    Current_Specification,           //!< Specification of the current in a boundary condition or in the volume
    Voltage_Specification,           //!< Specification of the voltage in a boundary condition or in the volume
    Dirichlet_Specification,         //!< Specification of a Dirichlet condition
    Mechanical_Model_Axial,          //!< Mechanics model in the axial direction
    Mechanical_Stress_Axial,         //!< Stress in the axial direciotn
    Species_Eqn_Offset,       //18   //!< Special equation name representing the first species equation offset
    Max_Eqn_Name                     //!< Max value of list of equation types
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
inline void operator++(EQ_TYPE& eqnT)
{
    eqnT = EQ_TYPE(eqnT + 1);
}
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

