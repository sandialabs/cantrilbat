/**
 * @file m1d_SurfDomainDescription.h
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_SURFDOMAINDESCRIPTION_H
#define M1D_SURFDOMAINDESCRIPTION_H

#include "m1d_BulkDomainDescription.h"

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

// ====================================================================================================================

//! The light-weight base class for boundaries between one-dimensional spatial
//! domains. The boundary may have its own internal variables, such as surface species coverages.
/*!
 *
 * The SurfDomainDescription boundary types are an inlet, an outlet, a symmetry plane,
 * a fluxMatching plane, and a surface.
 * *
 *  The variable NumEquationsPerNode represents the number of extra equations that
 *  are defined by this surface domain description. It is not necessarily equal
 *  to the number of surface boundary conditions that are defined on this node.
 */
class SurfDomainDescription : public DomainDescription {
public:

  //! Constructor
  /*!
   * In the constructor, we have typically been laying out what the unknowns are
   * and what the equations are, that are solved within the domain.
   *
   */
  SurfDomainDescription(DomainLayout *dl_ptr, std::string domainFunctionName = "", std::string domainName = "");

  //! Copy Constructor
  /*!
   * @param  r      Object to be copied
   */
  SurfDomainDescription(const SurfDomainDescription &r);

  //! Destructor
  virtual
  ~SurfDomainDescription();

  //! Assignment operator
  /*!
   * @param r       Object to be copied
   * @return        Return a changeable reference to the current object
   */
  SurfDomainDescription &
  operator=(const SurfDomainDescription &r);

  //! sets the id of the surface domain
  void
  setID(int id);

  //! Reports the id of the surface domain
  int
  ID() const;

  //! set the adjacent domains
  void
  setAdjBulkDomains(BulkDomainDescription *leftBulk, BulkDomainDescription *rightBulk);

  //! Set the Global node ID of the surface
  /*!
   *
   * @param locGbNode Value of the global node
   */
  void
  setGbNode(const int locGbNode);

  //! Determine the list of Equations and Variables
  /*!
   *  This routine is responsible for setting the variables:
   *    - VariableNameList
   *    - EquationNameList
   *  Here we clear these variables and assume that there are no equations and/or variables 
   *  assigned on the surface.
   */
  virtual void
  SetEquationsVariablesList();

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

  //! Malloc and return the object that will calculate the residual for the surface domain efficiently
  /*!
   * @return  Returns a pointer to the object that will calculate the residual
   *          efficiently
   */
  virtual SurDomain1D *
  mallocDomain1D();

  //! Sets the mapping between right and left domains to either
  //! of two extremes
  /*!
   * This sets two common mapping extremes. A different function
   * sets the inbetween cases.
   *
   * @param mapType     0 If the map type is 0, this means that
   *                      the right equations are mapped into the
   *                      left equations to the full extent possible
   *                    1 If the map type is 1, his means that
   *                      the right equations are never mapped into the
   *                      left equations at all. They are separate equations
   */
  void
  setRLMapping(int mapType = 0);

  //! Sets the mapping between right and left domains to the arbritrary case
  /*!
   *  Length is equal to the number of equations in the right domain
   *  If an entry is 0, this means that that equation maps into the left.
   *  If an entry is 1, this means that the right is a separate degree of
   *  freedom
   */
  void
  setRLMapping(const int * const rightConnectivity);

  // **********************************************************************
protected:
  //! ID of this surface domain.
  /*!
   *  this is a unique ID of the surface domain. It is used as an index to look up
   *  stuff in the DomainLayout object
   */
  int IDSurfDomain;

public:
  //! Global node number where the surface domain is located
  int LocGbNode;

  //! Left domain surface or bulk
  /*!
   *  Note if there is no left domain, a NULL value is used.
   *  frequently his is equal to LeftBulk
   */
  DomainDescription *LeftDomain;

  //! Left domain surface or bulk
  /*!
   *  Note if there is no left domain, a NULL value is used.
   */
  BulkDomainDescription *LeftBulk;

  //! Right bulk domain
  /*!
   *  Note if there is no left domain, a NULL value is used.
   */
  DomainDescription *RightDomain;

  //! Right bulk domain
  /*!
   *  Note if there is no left domain, a NULL value is used.
   */
  BulkDomainDescription *RightBulk;

  //! Right to left mapping of the bulk equations
  std::vector<int> RightToLeftBulkEqnMapping;

  SurDomain1D* SurDomain1DPtr_;
};

}

#endif
