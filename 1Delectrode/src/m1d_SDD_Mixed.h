/**
 * @file m1d_SDD_Mixed.h
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_SURFDOMAINTYPES_H_
#define M1D_SURFDOMAINTYPES_H_

#include "m1d_SurfDomainDescription.h"
#include "m1d_BoundaryCondition.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This class specifies that all equations are handled by a simple Dirichlet condition
/*!
 *  We provide hooks for adding Dirichlet conditions of all types onto arbitrary equation - variable combinations.
 *  This class also adds on to the Dirichlet condition surface domain description, by adding flux specification conditions.
 *  This is different than Dirichlet conditions, because the main continuity equation is not thrown out in these cases.
 */
class SDD_Mixed : public SurfDomainDescription
{
public:

    //! Constructor
    /*!
     *  We construct the object but don't actually specify any Dirichlet conditions.
     *  Later we can add dirichlet conditions into the object.
     *
     *  @param[in]           dl_ptr              Domain Layout object that owns this description.
     *  @param[in]           domainFunctionName  Functional name of domain. Defaults to "".
     *  @param[in]           domainName          name of domain. Defaults to "".
     */
    SDD_Mixed(DomainLayout* dl_ptr, std::string domainFunctionName = "", std::string domainName = "");

    //! In this constructor, we set all variables in adjoining bulk domains to a single constant value
    /*!
     *  Note this specification isn't for all cases. Setting all values of all variables
     *  to a single number must be considered to be a very, very special case.
     *
     *  @param[in]           dl_ptr              Domain Layout object that owns this description.
     *  @param[in]           value               Value to set all boundary conditions to.
     *  @param[in]           domainFunctionName  Functional name of domain. Defaults to "".
     *  @param[in]           domainName          name of domain. Defaults to "".
     */
    SDD_Mixed(DomainLayout* dl_ptr, double value, std::string domainFunctionName = "", std::string domainName = "");

    //! Virtual Destructor
    virtual ~SDD_Mixed();

    //! Copy Constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    SDD_Mixed(const SDD_Mixed& r);

    //! Assignment operator
    /*!
     *  @param[in]           r                   Object to be copied
     *  @return                                  Returns a changeable reference to the current object
     */
    SDD_Mixed& operator=(const SDD_Mixed& r);

    //! Add a constant Dirichlet Condition
    /*!
     *  This sets up a Dirichlet condition on an equation, setting a variable equal to a constant.
     *  It doesn't have to be the variable associated with the equation.
     *
     *  @param[in]           equationID          Equation ID to apply the Dirichlet condition to
     *  @param[in]           variableID          VariableID to apply the Dirichlet condition to
     *  @param[in]           value               Value to apply
     */
    void
    addDirichletCondition(EqnType equationID, VarType VariableID, double value);

    //! Add a Dirichlet Condition, assuming the default mapping between variable and equation ID.
    /*!
     *  This sets up a Dirichlet condition on an equation, setting a variable equal to a constant.
     *  The variable is defined to be the one associated with the equation.
     *
     *  @param[in]           variableID          VariableID to apply the Dirichlet condition to
     *  @param[in]           value               Value to apply
     */
    void
    addDirichletCondition(VarType VariableID, double value);

    //! Add a Dirichlet Condition with time dependent function pointer
    /*!
     *  This sets up a Dirichlet condition on an equation, setting a variable equal to the results of a
     *  function of time.
     *  The variable  doesn't have to be the variable associated with the equation.
     *
     * @param[in]            equationID          Equation ID to apply the flux condition to
     * @param[in]            variableID          VariableID to apply the flux condition to
     * @param[in]            value               Value to apply
     * @param[in]            timeDep             function pointer
     */
    void
    addDirichletCondition(EqnType equationID, VarType VariableID, double value, double (*timeDep)(double));

    //! Add a Dirichlet Condition with BoundaryCondition class
    /*!
     *  This sets up a Dirichlet condition on an equation, setting a variable equal to the results of a
     *  function of time.
     *  The variable  doesn't have to be the variable associated with the equation.
     *
     *  @param[in]           equationID          Equation ID to apply the flux condition to
     *  @param[in]           variableID          VariableID to apply the flux condition to
     *  @param[in]           BC_Type             boundary condition type 3, 4, 5  for const, stepTable, linearTable
     *  @param[in]           BC_timeDep          time dependent boundary condition pointer
     */
    void
    addDirichletCondition(EqnType equationID, VarType VariableID, int BC_Type, BoundaryCondition* BC_timeDep);

    //! Add a flux Condition to an equation
    /*!
     *  We reuse the Dirichlet condition structure to add flux conditions onto equations, by adding boundary condition types.
     *
     *  We add a flux condition which sets the derivative of the flux from the continuity equation equal to a constant value
     *  specified in the parameter list. The specification of the flux is up to the equation system. However, it's usually
     *  the pertinent one needed for the closure of the continuity equation.
     *
     *  \todo The specification of the variable seems superfluous. We are setting the specification of the flux in the
     *        equation system. The variable can not be varied. Check this out and get rid of the variable id.
     *
     *  @param[in]           equationID          Equation ID to apply the flux condition to
     *  @param[in]           variableID          VariableID to apply the flux condition to
     *  @param[in]           value               Value to apply
     */
    void addFluxCondition(const EqnType& equationID, const VarType& VariableID, double value);

    //! Add a flux Condition using time dependent continuous function
    /*!
     *
     * @param  equationID  Equation ID to apply the flux condition to
     * @param  variableID  VariableID to apply the flux condition to
     * @param  value  Value to apply
     */
    void addFluxCondition(const EqnType& equationID, const VarType& variableID, double value, double (*timeDep)(double));

    //! Add a flux Condition using Boundary Condition class
    /*!
     *
     * @param  equationID  Equation ID to apply the flux condition to
     * @param  variableID  VariableID to apply the flux condition to
     * @param  value  Value to apply
     */
    void
    addFluxCondition(const EqnType& equationID, const VarType& variableID, int BC_Type, BoundaryCondition* BC_timeDep);

    //! Add a Robin Mixed boundary Condition
    /*!
     *
     * @param  equationID  Equation ID to apply the flux condition to
     * @param  variableID  VariableID to apply the flux condition to
     * @param  value  Value to apply
     */
    void addRobinCondition(EqnType equationID, VarType variableID, BoundaryCondition* BC_timeDep, int bc_type=10);

    //! Set the equation description
    /*!
     *  (virtual from SurfDomainDescription)
     *  This routine is responsible for setting the variables:
     *    - NumEquationsPerNode
     *    - VariableNameList
     *    - EquationNameList
     *    - EquationIndexStart_EqName
     */
    virtual void
    SetEquationDescription() override;

    //! Malloc and Return the object that will calculate the residual efficiently
    /*!
     *  (virtual from SurfDomainDescription)
     *
     *  @return                                  Returns a pointer to the object that will calculate the residual  efficiently
     */
    virtual SurDomain1D*
    mallocDomain1D() override;

    // ----------------------------------------- D A T A --------------------------------------------------------------

    //! Number of Dirichlet conditions to apply
    /*!
     *   This may be applied to multiple subtypes within an equation type.
     *   by definitions with the EqnType and VarType objects.
     */
    int NumConditions;

    //! Equation type to apply them
    /*!
     *  Length is equal to NumConditions
     *  All vectors in this class have length of NumConditions, and the ith entry refers to the ith Dirichlet Condition
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

    //! Vector of pointers to Time dependent functions   
    /*!
     *   Vector has length equal to the number of Dirichlet conditions, NumConditions, defined at the node
     */
    std::vector<TimeDepFunction> TimeDep;

    //! BoundaryCondition Pointers for time dependent BC for BC_Type_NE = 4 - 9
    /*!
     *   Vector has length equal to the number of Dirichlet conditions, NumConditions, defined at the node
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
     * 10 Robin boundary condition (time dependent boundary condition [ flux = h (T - T0) ]
     *
     *  Vector has length equal to the number of Dirichlet conditions, NumConditions, defined at the node
     */
    std::vector<int> BC_Type_;

    //! SBC type
    /*!
     *   value of zero here indicates that an inherited class handles all setup functions
     */
    int SBC_Type_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif /* M1D_SURFDOMAINTYPES_H_ */
