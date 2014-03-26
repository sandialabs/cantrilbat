/**
 * @file m1d_DomainDescription.h
 *
 */

/*
 *  $Id: m1d_DomainDescription.h 567 2013-03-21 23:03:11Z hkmoffa $
 */

#ifndef M1D_DOMAINDESCRIPTION_H
#define M1D_DOMAINDESCRIPTION_H

#include "m1d_EqnVarTypes.h"
//#include "cantera/equil/MultiPhase.h"

#include <vector>
#include <string>

namespace m1d
{
//forward declarations
class SurfDomainDescription;
class SurDomain1D;
class BulkDomain1D;
class GlobalIndices;
class DomainLayout;

//! This is a light weight base class that describes a surface
//! or bulk domain using global indexing, which will be the same for
//! all processors.
/*!
 *  The base class for the heavyweight objects that calculate the actual residuals
 *  is called Domain1D. 
 *
 */
class DomainDescription {

public:

  //! Constructor
  DomainDescription(DomainLayout *dl_ptr, std::string domainName = "");

  //! Copy Constructor
  /*!
   * @param r      Object to be copied
   */
  DomainDescription(const DomainDescription &r);

  virtual
  ~DomainDescription();

  DomainDescription &
  operator=(const DomainDescription &r);

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

  // -------------------------------------------------------------------------------
  //               DATA
  // -------------------------------------------------------------------------------

  //! Number of equations
  /*!
   *   This is the number of equations per grid point
   *   within the domain. For surface domains, this is the number of
   *   unique additional equations at the node, over and above those
   *   represented by bulk domains.
   */
  int NumEquationsPerNode;

  //! unique id number for each domain, whether it's a surface
  //! or whether it's a volume.
  // int UniqueID;

  //! Vector containing the variable names as they appear in the
  //! solution vector
  /*!
   * This vector defines the ordering of the equations on the domain
   *
   * Length = number of equations defined on the domain
   */
  std::vector<VarType> VariableNameList;

  //! Listing of the equations types as a vector
  /*!
   * This vector defines the ordering of the equations on the domain.
   *
   * Length = number of equations defined on the domain
   */
  std::vector<EqnType> EquationNameList;

  //! Listing of the index of the first equation
  //! of a particular type within the domain
  /*!
   *  Length is the length of the EQ_Name_Enum structure
   */
  std::vector<int> EquationIndexStart_EqName;

  //! Listing of the index of the first equation
  //! of a particular type within the domain
  /*!
   *  Length is the length of the EQ_Name_Enum structure
   */
  std::vector<int> VariableIndexStart_VarName;

  //! Vector containing the info on whether these variables are algebraic constraints
  //! or whether they have a time derivative
  /*!
   * Length = number of equations defined on the domain
   */
  std::vector<int> IsAlgebraic_NE;

  //! Vector containing the info on whether these variables are arithmetically scaled variables
  /*!
   * Length = number of equations defined on the domain
   */
  std::vector<int> IsArithmeticScaled_NE;

  //! Name of the domain.
  std::string DomainName;

  //! Shallow pointer to the domain layout for this Domain Description
  DomainLayout *DL_ptr_;


  //! Print level that is set through the input file.
  /*!
   *   0 -> Don't print anything
   *   1 -> Print only about significant issues going on
   *   2 -> Print status information at regular intervals.
   *   3 -> Print ShowSolution at regular intervals
   *   4 -> Print ShowSolution at all successful time steps
   *   5 -> Print additional information at each time step
   *   6 -> Print some information about each electrode object at each time step
   *   7 -> Print a lot of information about each electrode object at each time step
   */
  int SolutionBehavior_printLvl_;

  //! Level of residual information printing done to stdout
  /*!
   *   0 -> Don't print anything
   *   1 -> Print only about significant issues going on
   *   2 -> Print status information at regular intervals.
   *   3 -> Print ShowResidual at regular intervals
   *   4 -> Print ShowResidual at all successful time steps
   *   5 -> Print additional information when ShowSolution is called.
   *   6 -> Print additional information when any residual is called.
   *   7 -> Print a lot of information about each when ShowSolution is called.
   *   8 -> Print a lot of information when any base or show residual is called
   *   9 -> Print a lot of information when any residual is called
   */
   int Residual_printLvl_;


};

// ==================================================================================



}

#endif
