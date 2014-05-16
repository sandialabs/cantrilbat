/**
 * @file ReactingVolDomain.h
 *
 * $Author: hkmoffa $
 * $Revision: 571 $
 * $Date: 2013-03-26 10:44:21 -0600 (Tue, 26 Mar 2013) $
 */


/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef EXTRAGLOBALRXN_H
#define EXTRAGLOBALRXN_H

#include "cantera/kinetics/InterfaceKinetics.h"
#include <string>
#include <vector>

namespace Cantera
{


//!  Class describing a global reaction, which is defined as
//!  a linear combination of elementary reactions creating a global stoichiometric result
/*!
 *   This is useful for defining thermodynamics of global processes that occur
 *   on a surface or in a homogeneous phase.
 */
class ExtraGlobalRxn
{

public:
    //!  Constructor takes a default kinetics pointer
    ExtraGlobalRxn(Kinetics* k_ptr);

    virtual ~ExtraGlobalRxn();

    void setupElemRxnVector(double* RxnVector,
                            int specialSpecies = -1);
    std::string reactionString();
    double deltaSpecValue(double* speciesVectorProperty);

    std::vector<int>& reactants();
    std::vector<int>& products();
    bool isReversible();

    double ROPValue(double* ROPKinVector);
    double FwdROPValue(double* FwdROPElemKinVector, double* RevROPElemKinVector);
    double RevROPValue(double* FwdROPElemKinVector, double* RevROPElemKinVector);

    double reactantStoichCoeff(int kKin);
    double productStoichCoeff(int kKin);
    bool m_ThisIsASurfaceRxn;
    double deltaRxnVecValue(double* rxnVectorProperty);

    //! This kinetics operator is associated with just one
    //! homogeneous phase, associated with tpList[0] phase
    /*!
     * This object owns the Kinetics object
     */
    Cantera::Kinetics* m_kinetics;

    //! This kinetics operator is associated with multiple
    //! homogeneous and surface phases.
    /*!
     * This object owns the Kinetics object
     */
    Cantera::InterfaceKinetics* m_InterfaceKinetics;

    int m_nKinSpecies;
    int m_nReactants;
    std::vector<int> m_Reactants;
    std::vector<doublereal> m_ReactantStoich;

    int m_nProducts;
    std::vector<int> m_Products;
    std::vector<doublereal> m_ProductStoich;

    int m_nNetSpecies;
    std::vector<int> m_NetSpecies;
    std::vector<doublereal> m_netStoich;

    int m_nRxns;
    std::vector<doublereal> m_ElemRxnVector;

    int m_SpecialSpecies;
    bool m_SpecialSpeciesProduct;
    int m_SS_index;

    int iphaseKin;
    bool m_ok;
    bool m_reversible;


};
}
#endif
