/**
 * @file m1d_EqnVarTypes.h
 *
 */

/*
 *  $Id: m1d_EqnVarTypes.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_EQNVARTYPES_H
#define M1D_EQNVARTYPES_H

#include "m1d_defs.h"
#include "m1d_Eqn_Names.h"

#include <string>

//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This class contains a name of a variable
/*!
 *   The class is used to identify variables.
 *
 *              Identification of Variables within the code.
 * 
 *      VAR_TYPE and VAR_TYPE_SUBNUM  pair completely identify a variable within the program.m1d_Eqn_Names.h
 *      What that means is that if the axial velocity has a step jump located at an interface,
 *      then the axial velocity variable in the two adjacent bulk regions must have different
 *      VAR_TYPE_SUBNUM values.
 *
 *      Also, if mole fraction unknowns for adjacent domains refer to different species in 
 *      different phases, then they must have different VAR_TYPE_SUBNUM values. Moreover, within
 *      these phases they have to have continuous values of the index. This means that the
 *      application is responsible for assigning indexes from 0 to (nsp_i - 1) in phase i
 *      and then (nsp_i) to (nsp_i + nsp_j - 1) in the adjacent phase j.
 *
 */
class VarType {
public:
  //! Constructor
  VarType();

  //! Constructor with specification
  /*!
   *  @param[in]             variableType        VAR_TYPE enum specifying what type of variable this is
   *  @param[in]             variableSubType     Integer value of the subvariable type. Defaults to 0.
   *  @param[in]             subName             Name of the subtype. Defaults to nullptr
   */
  VarType(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType = 0, const char *subName = nullptr);

  //! Destructor
  ~VarType();

  //! Returns the variable main name given the variable type (Static function)
  /*!
   *  @return                                    Returns a string representation of the main VAR_TYPE enum
   */
  static std::string VarMainName(const VAR_TYPE variableType);

  //! Copy Constructor
  /*!
   *  @param[in]             r                   Object to be copied
   */
  VarType(const VarType &r);

  //! Assignment Operator
  /*!
   *  @param[in]             r                   Object to be copied.
   *
   *  @return                                    Returns a variable reference to the current object
   */
  VarType& operator=(const VarType& r);

  //! Set the variable type of the object
  /*!
   *
   *  @param[in]             variableType        VAR_TYPE enum specifying what type of variable this is
   *  @param[in]             variableSubType     Integer value of the subvariable type. Defaults to 0.
   *  @param[in]             subName             Name of the subtype. Defaults to nullptr
   */
  void setID(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType = 0, const char *subName = nullptr);

  std::string
  VariableName(const int len = 128) const;

  //! Variable type
  m1d::VAR_TYPE VariableType;

  //! Variable subtype
  m1d::VAR_TYPE_SUBNUM VariableSubType;

  //!  Main Part of the variable Name
  char VariableMainName[24];

  //! Sub type name of the variable
  char VariableSubTypeName[24];
};
//==================================================================================================================================
//! This class contains a name of a variable
/*!
 *   The class is used to identify variables.
 */
class EqnType {
public:
  //! Construct
  EqnType();

  //! Constructor with specification
  /*!
   *
   * @param equationType
   * @param equationSubType
   * @param subName
   */
  EqnType(const EQ_TYPE equationType, const EQ_TYPE_SUBNUM equationSubType = 0, const char *subName = 0);

  //! Destructor
  ~EqnType();

  //! Copy Constructor
  /*!
   *  @param[in]             r                   Object to be copied
   */
  EqnType(const EqnType &r);

  //! Assignment Operator
  /*!
   *  @param[in]             r                   Object to be copied.
   *
   *  @return                                    Returns a variable reference to the current object
   */
  EqnType & operator=(const EqnType &r);

  //! Returns the equation main name given the equation type
  /*!
   *  @param[in]             equationType        Equation type as input
   *  @return                                    Returns a string
   */
  static std::string EqnMainName(const EQ_TYPE equationType);

  //! Set the equation type of the object
  /*!
   *
   *  @param[in]             equationType        EQ_TYPE enum for the equation type
   *  @param[in]             equationSubType     Integer like variable for the sub type of the equation 
   *  @param[in]             subName             Character Name of the equation. Will use a default value if this nullptr.
   *                                             Defaults to nullptr
   */
  void setID(const EQ_TYPE equationType, const EQ_TYPE_SUBNUM equationSubType = 0, const char *subName = nullptr);

  //! Return the name of the equation within a certain character length string
  /*!
   *  @param[in]             len                 Maximum length of the string
   *  @return                                    returns a C++ string
   */
  std::string EquationName(const int len = 128) const;

  // --------------------------------------------- D A T A ------------------------------------------------------------

  //! Equation type
  m1d::EQ_TYPE EquationType;

  //! Equation subtype
  m1d::EQ_TYPE_SUBNUM EquationSubType;

  //! Equation Main Name
  char EquationMainName[24];

  //! Equation sub type
  char EquationSubTypeName[24];
};

//==================================================================================================================================
//! Equality boolean operator for VarType Objects
/*!
 *  @param[in]               a                   Object 1
 *  @param[in]               b                   Object 2
 *  @return                                      Returns whether the two are the same
 */
bool operator==(const VarType &a, const VarType &b);

//==================================================================================================================================
//! Inequality boolean operator for VarType Objects
/*!
 *  @param[in]               a                   Object 1
 *  @param[in]               b                   Object 2
 *  @return                                      Returns whether the two are different
 */
bool operator!=(const VarType &a, const VarType &b);

//==================================================================================================================================
//! Greater than operator for VarType Objects
/*!
 *  VarTypes are first compared numerically via their enum values. Then, they are compared via their
 *  subindex numbers. This provides an ordering.
 *
 *  @param[in]               a                   Object 1
 *  @param[in]               b                   Object 2
 *  @return                                      Returns whether object a is greater than object b.
 */
bool operator>(const VarType &a, const VarType &b);

//==================================================================================================================================
//! Equality boolean operator for EqnType Objects
/*!
 *  @param[in]               a                   Object 1
 *  @param[in]               b                   Object 2
 *  @return                                      Returns whether the two are the same
 */
bool operator==(const EqnType &a, const EqnType &b);

//==================================================================================================================================
//! Inequality boolean operator for VarType Objects
/*!
 *  @param[in]               a                   Object 1
 *  @param[in]               b                   Object 2
 *  @return                                      Returns whether the two are different
 */
bool operator!=(const EqnType &a, const EqnType &b);

//==================================================================================================================================
//! Returns a mapping between the equation type and the variable type
/*!
 *  @param[in]               eqType              EQ_TYPE enum input describing the equation type
 *
 *  @return                                      Returns the corresponding variable type
 */
VAR_TYPE EqnToVarEnum(EQ_TYPE eqType);

//==================================================================================================================================
//! Returns a mapping between the equation type and the variable type
/*!
 *  @param[in]               et                  EqnType input describing everything about the equation
 *
 *  @return                                      Returns the corresponding VarType value
 */
VarType EqnTypeToVarType(EqnType et);

//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif

