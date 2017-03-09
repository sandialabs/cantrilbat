/**
 * @file m1d_BulkDomainDescription.h
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef M1D_BULKDOMAINDESCRIPTION_H
#define M1D_BULKDOMAINDESCRIPTION_H

#include "m1d_DomainDescription.h"
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This is a light weight class that describes the domain
//! using global indexing, which will be the same for all processors.
/*!
 *
 */
class BulkDomainDescription : public DomainDescription
{
public:

    //! Constructor
    /*!
     *  Functional names of domain may be used to define objects that may be mapped to input file parameters
     *  For example, functional names for bulk domains are "anode", "separator", and "cathode".
     *
     *  @param[in]           dl_ptr              Domain layout pointer
     *  @param[in]           domainFunctionName  Functional name of domain
     *  @param[in]           domainName          name of domain
     */
    BulkDomainDescription(DomainLayout* dl_ptr, std::string domainFunctionName = "", std::string domainName = "");

    //! Copy Constructor
    /*!
     *  @param[in]           r                   Object to be copied
     */
    BulkDomainDescription(const BulkDomainDescription& r);

    //! Destructor
    virtual ~BulkDomainDescription();

    //! Assignment operator
    /*!
     * There are some issues to be resolved here.
     *
     * @param[in]            r                   Object to be copied
     * @return                                   Returns a changeable reference to the current object
     */
    BulkDomainDescription& operator=(const BulkDomainDescription& r);

    //! Sets the id of the bulk domain
    /*!
     *  @param[in]           id                 Index used within lists of bulk domains
     */
    void setID(const int id);

    //! Reports the id index of the bulk domain
    /*!
     *  @return                                 Return the index used within lists of bulk domains.
     */
    int ID() const;

    //! Specify the left and right surface domains for this bulk domain
    /*!
     *  (virtual from BulkDomainDescription)
     *  Note, there must be viable surface domains on all sides  of a bulk domain
     *
     *  @param[in]           leftSurf            Pointer to the surface domain on the left
     *  @param[in]           rightSurf           Pointer to the surface domain on the right
     */
    virtual void 
    setAdjSurfDomains(SurfDomainDescription* const leftSurf, SurfDomainDescription* const rightSurf);

    //! Specify the left and right global nodes for this bulk domain
    /*!
     *  (virtual from BulkDomainDescription)
     *  @param[in]           leftGbNode          Left value of the node
     *  @param[in]           rightGbNode         Right value of the node
     */
    virtual void setGbNodeBounds(const int leftGbNode, const int rightGbNode);

    //! Specify the left and right X positions of the domain
    /*!
     *  (virtual from BulkDomainDescription)
     *  @param[in]           xleft               Left value of the domain boundary
     *  @param[in]           xright              right value of the domain boundary
     */
    virtual void setXposBounds(const double xleft, const double xright);

    //! Initialize the positions of the nodes in the domain
    /*!
     *  (virtual from BulkDomainDescription)
     *  @param[in]           gi_ptr              Pointer to the GlobalIndices structure
     */
    virtual void InitializeXposNodes(GlobalIndices* const gi_ptr);

    //! Set the equation description
    /*!
     *  (virtual from DomainDescription)
     *  This routine is responsible for setting the variables:
     *    - NumEquationsPerNode
     *    - VariableNameList
     *    - EquationNameList
     *    - EquationIndexStart_EqName
     */
    virtual void SetEquationDescription() override;

    //! Malloc and Return the object that will calculate the residual efficiently
    /*!
     *  (virtual from BulkDomainDescription)
     *  @return                                  Returns a pointer to the object that will calculate the residual
     *                                           efficiently
     */
    virtual BulkDomain1D* mallocDomain1D();

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
    SurfDomainDescription* LeftSurf;

    //! Right surf domain
    /*!
     *  Note if there is no right surf domain, a NULL value is used.
     */
    SurfDomainDescription* RightSurf;

    //! Pointer to the bulk domain object corresponding to the Bulk Domain Description
    BulkDomain1D* BulkDomainPtr_;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif
