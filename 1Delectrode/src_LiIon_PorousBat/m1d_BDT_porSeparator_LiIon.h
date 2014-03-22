/**
 * @file m1d_BDT_porousSeparator_LiIon.h
 */

/*
 *   $Id: m1d_BDT_porSeparator_LiIon.h 5 2012-02-23 21:34:18Z hkmoffa $
 */

#ifndef M1D_BDT_PORSEPARATOR_LIION_H_
#define M1D_BDT_PORSEPARATOR_LIION_H_

#include "m1d_BulkDomainDescription.h"

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
class BDT_porSeparator_LiIon : public BulkDomainDescription
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

    //! Pointer to the thermo object for the molten salt
    /*!
     *   We own this object
     */
    Cantera::ThermoPhase* ionicLiquid_;

    //! Pointer to the transport object for the molten salt
    /*!
     * We own this object
     */
    Cantera::Transport* trans_;
};
//=====================================================================================================================
}
//=====================================================================================================================
#endif /* M1D_BULKDOMAINTYPES_H_ */
