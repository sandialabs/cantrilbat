/**
 * @file m1d_BDT_porousSeparator_LiIon.h
 */

/*
 * Copywrite 2013 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef M1D_BDT_PORSEPARATOR_LIION_H_
#define M1D_BDT_PORSEPARATOR_LIION_H_

#include "m1d_BulkDomainDescription.h"
#include "m1d_BDD_porousFlow.h"

namespace Cantera
{
class Transport;
class ThermoPhase;
}

namespace m1d
{

//! This class consists of multiple species diffusing in a time
//! dependent manner.  There is a net flow and a net electric current.
/*!
 *  This class is used to test the implementation
 */
class BDT_porSeparator_LiIon : public BDD_porousFlow
{
public:

    //! Constructor
    /*!
     * This constructor constructs the bulk domain from a MultiPhase object.
     *
     * In the constructor, we have typically been laying out what the unknowns are
     * and what the equations are, that are solved within the domain.
     *
     * @param dl_ptr   Pointer to the domain layout object
     */
    BDT_porSeparator_LiIon(DomainLayout* dl_ptr);

    //! Destructor
    virtual
    ~BDT_porSeparator_LiIon();

    //! Copy Constructor
    /*!
     * @param r Object to be copied
     */
    BDT_porSeparator_LiIon(const BDT_porSeparator_LiIon& r);

    //! Assignment operator
    /*!
     * @param r    Object to be copied
     * @return     Returns a changeable reference to the current object
     */
    BDT_porSeparator_LiIon&
    operator=(const BDT_porSeparator_LiIon& r);

    //! Malloc and Return the object that will calculate the residual efficiently
    /*!
     * @return  Returns a pointer to the object that will calculate the residual
     *          efficiently
     */
    virtual BulkDomain1D*
    mallocDomain1D();

    // --------------------------------------------------------------------------------------------
    //                          DATA
    // --------------------------------------------------------------------------------------------

};
//=====================================================================================================================
}
//=====================================================================================================================
#endif /* M1D_BULKDOMAINTYPES_H_ */
