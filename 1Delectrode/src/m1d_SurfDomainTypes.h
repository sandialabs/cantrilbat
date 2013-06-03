/**
 * @file m1d_SurfDomainTypes.h
 */
/*
 * $Id: m1d_SurfDomainTypes.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_SURFDOMAINTYPES_H_
#define M1D_SURFDOMAINTYPES_H_

#include "m1d_SurfDomainDescription.h"
#include "m1d_BoundaryCondition.h"

namespace m1d
{

//! This class specifies that all equations are handled
//! by a simple Dirichlet condition
/*!
 *
 */
class SDT_Dirichlet : public SurfDomainDescription
{
public:

  //! Constructor
  /*!
   *   We construct the object but don't actually specify any Dirichlet conditions.
   *   Later we can add dirichlet conditions into the object.
   *
   * @param dl_ptr  Domain Layout object that owns this description.
   */
  SDT_Dirichlet(DomainLayout *dl_ptr,  std::string domainName = "");

  //! In this constructor, we set all variables in
  //! adjoining bulk domains to a single constant value
  /*!
   *  Note this specification isn't for all cases. Setting all values of all variables
   *  to a single number must be considered to be a very, very special case.
   *
   * @param value    Value to set all boundary conditions to.
   */
  SDT_Dirichlet(DomainLayout *dl_ptr, double value, std::string domainName = "");

  //! Destructor
  virtual
  ~SDT_Dirichlet();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  SDT_Dirichlet(const SDT_Dirichlet &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  SDT_Dirichlet &
  operator=(const SDT_Dirichlet &r);

  //! Add a Dirichlet Condition
  /*!
   *
   * @param  equationID  Equation ID to apply the Dirichlet condition to
   * @param  variableID  VariableID to apply the Dirichlet condition to
   * @param  value  Value to apply
   */
  void
  addDirichletCondition(EqnType equationID, VarType VariableID, double value);

  //! Add a Dirichlet Condition, assuming the default mapping between variable
  //! and equation ID.
  /*!
   * @param  variableID  VariableID to apply the Dirichlet condition to
   * @param  value  Value to apply
   */
  void
  addDirichletCondition(VarType VariableID, double value);

  //! Add a Dirichlet Condition with time dependent function pointer
  /*!
   *
   * @param  equationID  Equation ID to apply the flux condition to
   * @param  variableID  VariableID to apply the flux condition to
   * @param  value  Value to apply
   * @param  timeDep function pointer
   */
  void
  addDirichletCondition(EqnType equationID, VarType VariableID, double value, double (*timeDep)(double));

  //! Add a Dirichlet Condition with BoundaryCondition class
  /*!
   *
   * @param  equationID  Equation ID to apply the flux condition to
   * @param  variableID  VariableID to apply the flux condition to
   * @param  BC_Type boundary condition type 3, 4, 5  for const, stepTable, linearTable
   * @param  BC_timeDep time dependent boundary condition pointer
   */
  void
  addDirichletCondition(EqnType equationID, VarType VariableID, int BC_Type, BoundaryCondition * BC_timeDep);

  //! Set the equation description
  /*!
   *  This routine is responsible for setting the variables:
   *    - NumEquationsPerNode
   *    - VariableNameList
   *    - EquationNameList
   *    - EquationIndexStart_EqName
   */
  virtual void
  SetEquationDescription();

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   *
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual SurDomain1D *
  mallocDomain1D();

  //! Number of Dirichlet conditions to apply
  /*!
   *   This may be applied to multiple subtypes within an equation type.
   *   by definitions with the EqnType and VarType objects.
   */
  int NumConditions;

  //! Equation type to apply them
  /*!
   *  Length is equal to NumConditions
   */
  std::vector<EqnType> EquationID;

  //! Variable type to apply them on
  /*!
   *  Length is equal to NumConditions
   */
  std::vector<VarType> VariableID;

  //! Value of the variable
  /*!
   *  Length is equal to NumConditions
   */
  std::vector<double> Value;

  //! Function Pointers for time dependent BC for BC_Type = 2, 3
  /*!
   *   Vector has length equal to the number of equations defined at the node
   */
  typedef double (*TimeDepFunction)(double time);
  std::vector<TimeDepFunction> TimeDep;

  //! BoundaryCondition Pointers for time dependent BC for BC_Type_NE = 4 - 9
  /*!
   *   Vector has length equal to the number of equations defined at the node
   */
  std::vector<BoundaryCondition*> BC_TimeDep_;

  //!  Type of the boundary condition
   /*!
    *  0 Dirichlet
    *  1 pure flux condition
    *  2 Time Dependent Dirichlet
    *  3 Time Dependent pure flux
    *  4 Time Dependent Dirichlet using BCconstant
    *  5 Time Dependent pure flux using BCconstant
    *  6 Time Dependent Dirichlet using BCsteptable
    *  7 Time Dependent pure flux using BCsteptable
    *  8 Time Dependent Dirichlet using BClineartable
    *  9 Time Dependent pure flux using BClineartable
    */
   std::vector<int> BC_Type_;
};

//! This class specifies that some of the equations are handled by Dirichlet conditions
/*!
 *  However other equations are handled by flux and reaction conditions
 */
class SDT_Mixed : public SDT_Dirichlet
{
public:

  //! Constructor
  /*!
   *   We construct the object but don't actually specify any Dirichlet conditions.
   *   Later we can add dirichlet conditions into the object.
   *
   * @param dl_ptr  Domain Layout object that owns this description.
   */
  SDT_Mixed(DomainLayout *dl_ptr, std::string domainName = "");

  //! Destructor
  virtual
  ~SDT_Mixed();

  //! Copy Constructor
  /*!
   * @param r Object to be copied
   */
  SDT_Mixed(const SDT_Mixed &r);

  //! Assignment operator
  /*!
   * @param r    Object to be copied
   * @return     Returns a changeable reference to the current object
   */
  SDT_Mixed &
  operator=(const SDT_Mixed &r);

  //! Add a flux Condition
  /*!
   *
   * @param  equationID  Equation ID to apply the flux condition to
   * @param  variableID  VariableID to apply the flux condition to
   * @param  value  Value to apply
   */
  void
    addFluxCondition(EqnType equationID, VarType VariableID, double value);

  //! Add a flux Condition using time dependent continuous function
  /*!
   *
   * @param  equationID  Equation ID to apply the flux condition to
   * @param  variableID  VariableID to apply the flux condition to
   * @param  value  Value to apply
   */
  void
    addFluxCondition(EqnType equationID, VarType variableID, double value, double (*timeDep)(double));

  //! Add a flux Condition using Boundary Condition class
  /*!
   *
   * @param  equationID  Equation ID to apply the flux condition to
   * @param  variableID  VariableID to apply the flux condition to
   * @param  value  Value to apply
   */
  void
    addFluxCondition(EqnType equationID, VarType variableID, int BC_Type, BoundaryCondition *BC_timeDep);

  //! Malloc and Return the object that will calculate the residual efficiently
  /*!
   *
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual SurDomain1D *
  mallocDomain1D();

  //! SBC type
  /*!
   *   value of zero here indicates that an inherited class handles all setup functions
   */
  int SBC_Type_;

};

}

#endif /* M1D_SURFDOMAINTYPES_H_ */
