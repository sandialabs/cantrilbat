/**
 * @file m1d_BulkDomainDescription.h
 *
 */

/*
 *  $Id: m1d_BulkDomainDescription.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_BULKDOMAINDESCRIPTION_H
#define M1D_BULKDOMAINDESCRIPTION_H

#include "m1d_DomainDescription.h"

namespace m1d
{

//! This is a light weight class that describes the domain
//! using global indexing, which will be the same for
//! all processors.
class BulkDomainDescription : public DomainDescription {
public:

  //! Constructor
  BulkDomainDescription(DomainLayout *dl_ptr, std::string domainName = "");

  //! Copy Constructor
  /*!
   * @param r   Object to be copied
   */
  BulkDomainDescription(const BulkDomainDescription &r);

  //! Destructor
  virtual
  ~BulkDomainDescription();

  //! Assignment operator
  /*!
   * There are some issues to be resolved here.
   *
   * @param r  Object to be copied
   * @return Returns a changeable reference to the current object
   */
  BulkDomainDescription &
  operator=(const BulkDomainDescription &r);

  //! Sets the id of the bulk domain
  void
  setID(const int id);

  //! Reports the id of the bulk domain
  int
  ID() const;

  //! Specify the left and right surface domains for this
  //! bulk domain
  /*!
   *  Note, there must be viable surface domains on all sides
   *  of a bulk domain
   *
   * @param leftSurf   Pointer to the surface domain on the left
   * @param rightSurf  Pointer to the surface domain on the right
   */
  void
  setAdjSurfDomains(SurfDomainDescription *leftSurf, SurfDomainDescription *rightSurf);

  //! Specify the left and right global nodes for this bulk domain
  /*!
   *
   * @param leftGbNode   Left value of the node
   * @param rightGbNode  Right value of the node
   */
  void
  setGbNodeBounds(const int leftGbNode, const int rightGbNode);

  //! Specify the left and right X Positions of the domain
  /*!
   * @param xleft   Left value of the domain boundary
   * @param xright  right value of the domain boundary
   */
  void
  setXposBounds(const double xleft, const double xright);

  //! Initialize the positions of the nodes in the domain
  void
  InitializeXposNodes(GlobalIndices *gi_ptr);

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
  virtual BulkDomain1D *
  mallocDomain1D();

  // --------------------------------------------------------------------------------------------
  //                       DATA
  // --------------------------------------------------------------------------------------------
protected:
  //! ID of this bulk domain.
  /*!
   *  This is a unique ID of the bulk domain. It is used as an index to look up
   *  stuff in the DomainLayout object. 
   */
  int IDBulkDomain;

public:
  //! First global node number of this bulk domain
  /*!
   *  There is an implicit assumption that all global nodes in this domain are
   *  contiguously numbered
   */
  int FirstGbNode;

  //! Last global node number of this bulk domain
  /*!
   *  There is an implicit assumption that all global nodes in this domain are
   *  contiguously numbered
   */
  int LastGbNode;

  //!  Axial position of the first global node
  double Xpos_start;

  //!  Axial position of the last global node
  double Xpos_end;

  //! Left surface domain
  /*!
   *  Note if there is no surf domain, a NULL value is used.
   */
  SurfDomainDescription *LeftSurf;

  //! Right surf domain
  /*!
   *  Note if there is no right surf domain, a NULL value is used.
   */
  SurfDomainDescription *RightSurf;

  //! Pointer to the bulk domain object corresponding to the Bulk Domain Description
  BulkDomain1D *BulkDomainPtr_;
};

}

#endif
