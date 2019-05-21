/**
 * @file m1d_SurfDomainTypes.h
 */

#ifndef M1D_SDD_FLATANODE_H_
#define M1D_SDD_FLATANODE_H_

#include "m1d_SDD_Mixed.h"

#include "zuzax/thermo/IonsFromNeutralVPSSTP.h"  // ion properties

namespace Zuzax
{
class Electrode;
}
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{
//==================================================================================================================================
//! This class specifies that all equations are handled
//! by a simple Dirichlet condition
/*!
 *  This is a substitute class for a full electrode as well
 */
class SDD_FlatAnode : public SDD_Mixed
{
public:

    //! Constructor
    /*!
     *   We construct the object but don't actually specify any Dirichlet conditions.
     *   Later we can add dirichlet conditions into the object.
     *
     * In the constructor, we have typically been laying out what the unknowns are
     * and what the equations are, that are solved within the domain.
     *
     * @param dl_ptr  Domain Layout object that owns this description.
     */
    SDD_FlatAnode(DomainLayout* dl_ptr, int position);

    //! Destructor
    virtual
    ~SDD_FlatAnode();

    //! Copy Constructor
    /*!
     * @param r Object to be copied
     */
    SDD_FlatAnode(const SDD_FlatAnode& r);

    //! Assignment operator
    /*!
     * @param r    Object to be copied
     * @return     Returns a changeable reference to the current object
     */
    SDD_FlatAnode&
    operator=(const SDD_FlatAnode& r);

    //! Set the equation and variables list
    /*!
     *  This routine is responsible for setting the variables:
     *    - VariableNameList
     *    - EquationNameList
     */
    virtual void SetEquationsVariablesList();

    //! Set the equation description
    /*!
     *  This routine is responsible for setting the variables:
     *    - NumEquationsPerNode
     *    - EquationIndexStart_EqName
     */
    virtual void SetEquationDescription();

    //! Malloc and Return the object that will calculate the residual efficiently
    /*!
     *
     * @return  Returns a pointer to the object that will calculate the residual
     *          efficiently
     */
    virtual SurDomain1D* mallocDomain1D();

    //! top or bottom of the domain
    /*!
     *   0 - top, right
     *   1 - bottom, left
     */
    int m_position;

    //! Pointer to the electrode object
    /*!
     *   We own the electrode object
     */
    Zuzax::Electrode* ElectrodeA_;

    //! Pointer to the thermo object for the molten salt
    /*!
     */
    Zuzax::IonsFromNeutralVPSSTP *ionicLiquidIFN_;

    //! Make the SurDomain1D class a friend so that it can access all of the stuff in this class
    friend class SurDomain_FlatLiSiAnode;
};
//==================================================================================================================================
}
//----------------------------------------------------------------------------------------------------------------------------------
#endif 
