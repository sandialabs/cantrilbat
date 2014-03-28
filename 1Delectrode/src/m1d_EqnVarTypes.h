/**
 * @file m1d_EqnVarTypes.h
 *
 */

/*
 *  $Id: m1d_EqnVarTypes.h 5 2012-02-23 21:34:18Z hkmoffa $
 */
#ifndef M1D_EQNVARTYPES_H
#define M1D_EQNVARTYPES_H
#include "m1d_Eqn_Names.h"

#include <string>

namespace m1d
{

//! This class contains a name of a variable
/*!
 *   The class is used to identify variables.
 */
class VarType {
public:
  //! Constructor
  VarType();

  //! Constructor with specification
  /*!
   * @param variableType
   * @param variableSubType
   * @param subName
   */
  VarType(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType = 0, const char *subName = 0);

  //! Destructor
  ~VarType();

  //! Returns the variable main name given the variable type
  static std::string
  VarMainName(const VAR_TYPE variableType);

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  VarType(const VarType &r);

  //! Assignment Operator
  /*!
   * @param r Object to be copied.
   * @return Returns a variable reference to the current object
   */
  VarType &
  operator=(const VarType &r);

  //! Set the variable type of the object
  /*!
   *
   * @param variableType
   * @param variableSubType
   * @param subName
   */
  void
  setID(const VAR_TYPE variableType, const VAR_TYPE_SUBNUM variableSubType = 0, const char *subName = 0);

  std::string
  VariableName(const int len) const;

  //! Variable type
  m1d::VAR_TYPE VariableType;

  //! Variable subtype
  m1d::VAR_TYPE_SUBNUM VariableSubType;

  //!  Main Part of the variable Name
  char VariableMainName[24];

  //! Sub type name of the variable
  char VariableSubTypeName[24];
};

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
   * @param r Object to be copied
   */
  EqnType(const EqnType &r);

  //! Assignment Operator
  /*!
   * @param r Object to be copied.
   * @return Returns a variable reference to the current object
   */
  EqnType &
  operator=(const EqnType &r);

  //! Returns the equation main name given the equation type
  /*!
   * @param equationType   Equation type as input
   * @return   Returns a string
   */
  static std::string
  EqnMainName(const EQ_TYPE equationType);

  //! Set the equation type of the object
  /*!
   *
   * @param equationType
   * @param equationSubType
   * @param subName
   */
  void
  setID(const EQ_TYPE equationType, const EQ_TYPE_SUBNUM equationSubType = 0, const char *subName = 0);

  //! return the name of the equation
  /*!
   *
   * @param len  Maximum length of the string
   * @return  returns a malloced C++ string
   */
  std::string
  EquationName(const int len) const;

  //! Equation type
  m1d::EQ_TYPE EquationType;

  //! Equation subtype
  m1d::EQ_TYPE_SUBNUM EquationSubType;

  //! Equation Main Name
  char EquationMainName[24];

  //! Equation sub type
  char EquationSubTypeName[24];

};

//! Equality boolean operator for VarType Objects
/*!
 *
 * @param a Object 1
 * @param b Object 2
 * @return Returns whether the two are the same
 */
bool
operator==(const VarType &a, const VarType &b);

bool
operator!=(const VarType &a, const VarType &b);

//! Greater than  operator for VarType Objects
/*!
 * VarTypes are first compared numerically via their
 * enum values. Then, they are compared via their
 * subindex numbers. This provides an ordering.
 *
 * @param a Object 1
 * @param b Object 2
 * @return Returns whether object a is greater
 *         than object b.
 */
bool
operator>(const VarType &a, const VarType &b);



bool
operator==(const EqnType &a, const EqnType &b);
bool
operator!=(const EqnType &a, const EqnType &b);

VAR_TYPE
EqnToVarEnum(EQ_TYPE eqType);

VarType
EqnTypeToVarType(EqnType et);

}
#endif

