/**
 * @file m1d_EqnVarTypes.cpp
 *
 */

/*
 *  $Id: m1d_EqnVarTypes.cpp 5 2012-02-23 21:34:18Z hkmoffa $
 */

#include "m1d_EqnVarTypes.h"

#include "m1d_exception.h"

#include <cstring>
#include <cstdio>
#include <sstream>
#include <ostream>

using namespace std;

namespace m1d
{
//=====================================================================================================================
std::string
VarType::VarMainName(const VAR_TYPE variableType)
{
  switch (variableType)
  {
    case Displacement_Axial:
      return string("Displacement_Axial");
      break;
    case Velocity_Axial:
      return string("Vel_axial");
      break;
    case Velocity_Radial:
      return string("Vel_radial");
      break;
    case Velocity_Swirl:
      return string("Vel_swirl");
      break;
    case Temperature:
      return string("Temperature");
      break;
    case Pressure_Axial:
      return string("Pres_axial");
      break;
    case PressureGradient_Radial:
      return string("PresGrad_radial");
      break;
    case MoleFraction_Species:
      return string("MF_sp");
      break;
    case Concentration_Species:
      return string("Conc_sp");
      break;
    case VolumeFraction_Species:
      return string("VF_sp");
      break;
    case VolumeFraction_Phase:
      return string("VF_ph");
      break;
    case Voltage:
      return string("Volt");
      break;
#ifdef MECH_MODEL
    case Solid_Stress_Axial:
      return string("Solid_Stress_Axial");
      break;
#endif
    case Max_Var_Name:
      return string("Max_Variable");
    default:
      throw m1d_Error("VarMainName", "unknown");
      break;
  }
  return string("");
}
//=====================================================================================================================
std::string
EqnType::EqnMainName(const EQ_TYPE equationType)
{
  switch (equationType)
  {
    case Momentum_Axial:
      return std::string("Mom_axial");
      break;
    case Momentum_Radial:
      return std::string("Mom_radial");
      break;
    case Momentum_Swirl:
      return std::string("Mom_swirl");
      break;
    case MeshDisplacement_Axial:
      return std::string("MeshDisplacement_axial");
      break;
    case Enthalpy_Conservation:
      return std::string("Enthalpy_conservation");
      break;
    case Thermal_Conservation:
      return std::string("Thermal_conservation");
      break;
    case Thermal_Dirichilet:
      return std::string("Thermal_Dirichilet");
      break;
    case Continuity:
      return std::string("Continuity");
      break;
    case Continuity_Global:
      return std::string("Continuity_global");
      break;
    case Species_Conservation:
      return std::string("Species_conservation");
      break;
    case VolumeFraction_Conservation:
      return std::string("VF_sp");
      break;
    case PhaseVolumeFraction_Conservation:
      return string("VF_ph");
      break;
    case Current_Conservation:
      return std::string("Current_conservation");
      break;
    case MoleFraction_Summation:
      return std::string("MF_summation");
      break;
    case ChargeNeutrality_Summation:
      return std::string("ChargeNeutrality_summation");
      break;
    case Current_Specification:
      return std::string("Current_specification");
      break;
    case Voltage_Specification:
      return std::string("Voltage_specification");
      break;
    case Species_Eqn_Offset:
      return std::string("Species_Eqn_Offset");
      break;
#ifdef MECH_MODEL
  case Mechanical_Model_Axial:
    return std::string("Mechanical_Model_Axial");
    break;
  case Mechanical_Stress_Axial:
    return std::string("Mechanical_Stress_Axial");
    break;
#endif 
    case Max_Eqn_Name:
      return std::string("Max_Eqn_Name");
      break;
    default:
      {
	ostringstream oss;
	oss << " eqn num "<<equationType;

	throw m1d_Error("EqnMainName", oss.str() );
	break;
      }
  }
  return string("");
}
//=====================================================================================================================
// Constructor
VarType::VarType() :
  VariableType(Variable_Type_Unspecified), VariableSubType(-1)
{
  VariableMainName[0] = '\0';
  VariableSubTypeName[0] = '\0';
}
//=====================================================================================================================
VarType::VarType(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType, const char *subName) :
  VariableType(variableType), VariableSubType(variableSubType)
{
  VariableMainName[0] = '\0';
  std::string h = VarMainName(VariableType);
  strncpy(VariableMainName, h.c_str(), 23);
  VariableSubTypeName[0] = '\0';
  if (subName) {
    strncpy(VariableSubTypeName, subName, 23);
  }
  VariableSubTypeName[23] = '\0';
}
//=====================================================================================================================
//! Destructor
VarType::~VarType()
{
}
//=====================================================================================================================
// Copy Constructor
/*
 * @param r Object to be copied
 */
VarType::VarType(const VarType &r) :
  VariableType(Variable_Type_Unspecified), VariableSubType(-1)
{
  VariableMainName[0] = '\0';
  VariableSubTypeName[0] = '\0';
  operator=(r);
}
//=====================================================================================================================
// Assignment Operator
/*
 * @param r Object to be copied.
 * @return Returns a variable reference to the current object
 */
VarType &
VarType::operator=(const VarType &r)
{
  if (this == &r) {
    return *this;
  }
  VariableType = r.VariableType;
  VariableSubType = r.VariableSubType;
  strncpy(VariableMainName, r.VariableMainName, 23);
  strncpy(VariableSubTypeName, r.VariableSubTypeName, 23);
  return *this;
}

//=====================================================================================================================
void
VarType::setID(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType, const char *subName)
{
  VariableType = variableType;
  VariableSubType = variableSubType;
  if (subName) {
    strncpy(VariableSubTypeName, subName, 23);
  } else {
    VariableSubTypeName[0] = '\0';
  }
  VariableMainName[0] = '\0';
  VariableSubTypeName[23] = '\0';
}
//=====================================================================================================================
std::string
VarType::VariableName(const int len) const
{
  char buf[128];
  if (strlen(VariableMainName) > 0) {
    strncpy(buf, VariableMainName, 63);
  } else {
    sprintf(buf, "Var_%d ", (int) VariableSubType);
  }
  int ll = strlen(buf);
  if (len > (ll + 3)) {
    int left = len - (ll);
    if (strlen(VariableSubTypeName) > 0) {
      snprintf(buf + ll, left, "(%s", VariableSubTypeName);
      int ll = strlen(buf);
      sprintf(buf + ll, ")");
    } else {
      if ((int) VariableSubType > 0) {
        snprintf(buf + ll, left, "(%d", (int) VariableSubType);
        int ll = strlen(buf);
        sprintf(buf + ll, ")");
      }
    } 
  }
  return std::string(buf);
}
//=====================================================================================================================
// Constructor
EqnType::EqnType() :
  EquationType(Equation_Type_Unspecified), EquationSubType(-1)
{
  EquationMainName[0] = '\0';
  EquationSubTypeName[0] = '\0';
}
//=====================================================================================================================
EqnType::EqnType(const EQ_TYPE equationType, const EQ_TYPE_SUBNUM equationSubType, const char *subName) :
  EquationType(equationType), EquationSubType(equationSubType)
{
  EquationMainName[0] = '\0';
  string h = EqnMainName(equationType);
  strncpy(EquationMainName, h.c_str(), 23);
  EquationSubTypeName[0] = '\0';
  if (subName) {
    strncpy(EquationSubTypeName, subName, 23);
  }
  EquationMainName[23] = '\0';
  EquationSubTypeName[23] = '\0';
}
//=====================================================================================================================
//! Destructor
EqnType::~EqnType()
{
}
//=====================================================================================================================
// Copy Constructor
/*
 * @param r Object to be copied
 */
EqnType::EqnType(const EqnType &r) :
  EquationType(Equation_Type_Unspecified), EquationSubType(-1)
{
  EquationMainName[0] = '\0';
  EquationSubTypeName[0] = '\0';
  operator=(r);
}
//=====================================================================================================================
// Assignment Operator
/*
 * @param r Object to be copied.
 * @return Returns a variable reference to the current object
 */
EqnType &
EqnType::operator=(const EqnType &r)
{
  if (this == &r) {
    return *this;
  }
  EquationType = r.EquationType;
  EquationSubType = r.EquationSubType;
  strcpy(EquationMainName, r.EquationMainName);
  strcpy(EquationSubTypeName, r.EquationSubTypeName);
  return *this;
}
//=====================================================================================================================
void
EqnType::setID(const EQ_TYPE equationType, const EQ_TYPE_SUBNUM equationSubType, const char *subName)
{
  EquationType = equationType;
  EquationSubType = equationSubType;
  if (subName) {
    strncpy(EquationSubTypeName, subName, 23);
  } else {
    EquationSubTypeName[0] = '\0';
  }
  EquationMainName[0] = '\0';
  EquationSubTypeName[23] = '\0';
}
//=====================================================================================================================
std::string
EqnType::EquationName(const int len) const
{
  char buf[128];
  if (strlen(EquationMainName) > 0) {
    strncpy(buf, EquationMainName, 63);
  } else {
    sprintf(buf, "Eqn_%d ", (int) EquationSubType);
  }
  int ll = strlen(buf);
  if (len > (ll + 3)) {
    int left = len - (ll);
    if (strlen(EquationSubTypeName) > 0) {
      snprintf(buf + ll, left, "(%s", EquationSubTypeName);
      ll = strlen(buf);
      sprintf(buf + ll, ")");
    } else {
      if ((int) EquationSubType > 0) {
        snprintf(buf + ll, left, "(%d", (int) EquationSubType);
        int ll = strlen(buf);
        sprintf(buf + ll, ")");
      }
    }
  }
  return std::string(buf);
}
//=====================================================================================================================
bool
operator==(const VarType &a, const VarType &b)
{
  if (a.VariableType == b.VariableType) {
    if (a.VariableSubType == b.VariableSubType) {
      return true;
    }
  }
  return false;
}
//=====================================================================================================================
bool
operator!=(const VarType &a, const VarType &b)
{
  if (a.VariableType == b.VariableType) {
    if (a.VariableSubType == b.VariableSubType) {
      return false;
    }
  }
  return true;
}
//=====================================================================================================================
bool
operator>(const VarType &a, const VarType &b)
{
  if (a.VariableType > b.VariableType) {
    return true;
  } else if (a.VariableType < b.VariableType) {
    return false;
  } else {
    if (a.VariableSubType > b.VariableSubType) {
      return true;
    } else if (a.VariableSubType < b.VariableSubType) {
      return false;
    }
  }
  return false;
}
//=====================================================================================================================
bool
operator==(const EqnType &a, const EqnType &b)
{
  if (a.EquationType == b.EquationType) {
    if (a.EquationSubType == b.EquationSubType) {
      return true;
    }
  }
  return false;
}
//=====================================================================================================================
bool
operator!=(const EqnType &a, const EqnType &b)
{
  if (a.EquationType == b.EquationType) {
    if (a.EquationSubType == b.EquationSubType) {
      return false;
    }
  }
  return true;
}
//=====================================================================================================================
VAR_TYPE
EqnToVarEnum(EQ_TYPE eqType)
{
  if (eqType == Momentum_Axial) {
    return Velocity_Axial;
  } else if (eqType == Momentum_Radial) {
    return Velocity_Axial;
  } else if (eqType == Momentum_Swirl) {
    return Velocity_Swirl;
  } else if (eqType == MeshDisplacement_Axial ) {
    return Displacement_Axial; 
  } else if (eqType == Enthalpy_Conservation) {
    return Temperature;
  } else if (eqType == Thermal_Conservation) {
    return Temperature;
  } else if (eqType == Thermal_Dirichilet) {
    return Temperature;
  } else if (eqType == Continuity) {
    return Pressure_Axial;
  } else if (eqType == Continuity_Global) {
    return PressureGradient_Radial;
  } else if (eqType == Species_Conservation) {
    return Concentration_Species;
  } else if (eqType == VolumeFraction_Conservation) {
    return VolumeFraction_Phase;
  } else if (eqType == PhaseVolumeFraction_Conservation) {
    return VolumeFraction_Species;
  } else if (eqType == Current_Conservation) {
    return Voltage;
  } else if (eqType == Current_Specification) {
    return Voltage;
  } else if (eqType == Voltage_Specification) {
    return Voltage;
  } else if (eqType == Equation_Type_Any) {
    return Variable_Type_Any;
  } else if (eqType == Equation_Type_Unspecified) {
    return Variable_Type_Unspecified;
  } else if (eqType == Species_Eqn_Offset) {
    return Concentration_Species;
#ifdef MECH_MODEL
  } else if (eqType == Mechanical_Model_Axial) {
    return Displacement_Axial;
 } else if (eqType == Mechanical_Stress_Axial) {
    return Solid_Stress_Axial;
#endif
  } else {
    m1d_Error("EqnToVarEnum", "Unknown Conversion");
  }
  return Variable_Type_Unspecified;
}
//=====================================================================================================================
VarType
EqnTypeToVarType(EqnType et)
{
  EQ_TYPE eqType = et.EquationType;
  VAR_TYPE vtType = EqnToVarEnum(eqType);
  return VarType(vtType, et.EquationSubType, et.EquationSubTypeName);
}
//=====================================================================================================================
}
/* end of namespace */
//=====================================================================================================================
