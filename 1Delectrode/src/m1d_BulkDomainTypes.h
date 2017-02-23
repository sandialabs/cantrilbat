/**
 * @file m1d_BulkDomainTypes.h
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */


#ifndef M1D_BULKDOMAINTYPES_H_
#define M1D_BULKDOMAINTYPES_H_

#include "m1d_BulkDomainDescription.h"

namespace m1d
{

//! This class consists of one species diffusing
/*!
 *  This class is used to test the implementation
 */
class BDT_SimpleDiff : public BulkDomainDescription
{
public:

    //! Constructor
    BDT_SimpleDiff(DomainLayout* dl_ptr);

    //! Constructor
    /*!
     * In the constructor, we have typically been laying out what the unknowns are
     * and what the equations are, that are solved within the domain.
     *
     */
    BDT_SimpleDiff(DomainLayout* dl_ptr, int id);

    //! Destructor
    virtual
    ~BDT_SimpleDiff();

    //! Copy Constructor
    /*!
     * @param r Object to be copied
     */
    BDT_SimpleDiff(const BDT_SimpleDiff& r);

    //! Assignment operator
    /*!
     * @param r    Object to be copied
     * @return     Returns a changeable reference to the current object
     */
    BDT_SimpleDiff&
    operator=(const BDT_SimpleDiff& r);

    //! Determine the list of Equations and Variables
    /*!
     *  This routine is responsible for setting the variables:
     *    - VariableNameList
     *    - EquationNameList
     */
    virtual void
    SetEquationsVariablesList();

    //! Malloc and Return the object that will calculate the residual efficiently
    /*!
     * @return  Returns a pointer to the object that will calculate the residual
     *          efficiently
     */
    virtual BulkDomain1D* mallocDomain1D();


    //! Equation type and var type to apply them
    std::vector<VarType> EquationID;


};
//=====================================================================================================================
//=====================================================================================================================
//=====================================================================================================================

//! This class consists of one species diffusing using time dependence
/*!
 *  This class is used to test the implementation
 */
class BDT_SimpleTDDiff : public BulkDomainDescription
{
public:

    //! Constructor
    BDT_SimpleTDDiff(DomainLayout* dl_ptr);

    //! Constructor
    /*!
     * In the constructor, we have typically been laying out what the unknowns are
     * and what the equations are, that are solved within the domain.
     *
     */
    BDT_SimpleTDDiff(DomainLayout* dl_ptr, int id);

    //! Destructor
    virtual
    ~BDT_SimpleTDDiff();

    //! Copy Constructor
    /*!
     * @param r Object to be copied
     */
    BDT_SimpleTDDiff(const BDT_SimpleTDDiff& r);

    //! Assignment operator
    /*!
     * @param r    Object to be copied
     * @return     Returns a changeable reference to the current object
     */
    BDT_SimpleTDDiff&
    operator=(const BDT_SimpleTDDiff& r);

    //! Determine the list of Equations and Variables
    /*!
     *  This routine is responsible for setting the variables:
     *    - VariableNameList
     *    - EquationNameList
     */
    virtual void
    SetEquationsVariablesList();

    //! Malloc and Return the object that will calculate the residual efficiently
    /*!
     * @return  Returns a pointer to the object that will calculate the residual
     *          efficiently
     */
    virtual BulkDomain1D* mallocDomain1D();


    //! Equation type and var type to apply them
    std::vector<VarType> EquationID;


};


}

#endif /* M1D_BULKDOMAINTYPES_H_ */
