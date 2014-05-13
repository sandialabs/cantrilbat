/*
 *  @file EState.h
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ESTATE_RADIALDISTRIB_H
#define _ESTATE_RADIALDISTRIB_H

#include "EState.h"

namespace Cantera
{

class Electrode;
class Electrode_SimpleDiff;
class Electrode_DiffTALE;


//! Child Class for the EState class concept refering to SimpleDiff Electrode object.
/*!
 *    This class adds the inner_platNum and second_platNum records to the output file.
 *
 *   It's an example of how to write child classes for EState.
 */
class EState_RadialDistrib : public EState
{

public:

    //! Default constructor for the base object
    EState_RadialDistrib();

    //! Default constructor for the base object
    EState_RadialDistrib(std::string electrodeType);

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    EState_RadialDistrib(const EState_RadialDistrib& right);

    //! Destructor is virtual
    virtual ~EState_RadialDistrib();

    //! Assignment operator
    /*!
     *  @param right  object to be duplicated
     */
    EState_RadialDistrib& operator=(const EState_RadialDistrib& right);

    //! Duplicator function for this class
    /*!
     *  @return Returns a duplication of the current state as a pointer to the base class
     */
    virtual EState* duplMyselfAsEState() const;

    //! Initialize the object based on an electrode Base class
    /*!
     *   This call will initialize all of the arrays within this class.
     *   All of the species and phase identification information is created and the class is
     *   readied for use as a state maintainer.
     *
     *  @param e    Pointer to the electrode base class
     */
    virtual int initialize(const Cantera::Electrode* const e);

    //! Create an indentification XML_Node element for this Electrode EState object
    /*!
     *  @return Returns a malloced XML_Node tree.
     */
    XML_Node* writeIdentificationToXML() const;

    //! Write the ElectrodeState to an XML_Node tree
    /*!
     *  (virtual function)
     *
     *  virtual function, because the base class is called from general code, allowing
     *      the child classes to be invoked.
     *
     *  @return pointer to the XML_Node tree
     */
    virtual XML_Node* writeStateToXML() const;

    //! Read the state from the XML_Node  given by the argument
    /*!
     *  @param xmlRoot   Root of the xml tree to get the information from
     */
    virtual void readStateFromXML(const XML_Node& xmlEState);

    void copyElectrode_SimpleDiff_intoState(const Cantera::Electrode_SimpleDiff* const emp);
    void copyElectrode_DiffTALE_intoState(const Cantera::Electrode_DiffTALE* const emp);

    //! Set the State of this object from the state of the Electrode object
    /*!
     *  (virtual function)
     *
     *  virtual function, because the base class is called from general code, allowing
     *      the child classes to be invoked.
     *
     *  This function takes the electrode objects _final_ state and copies it into this object.
     *  This function must be carried out before the XML_Node tree is written.
     *
     *  @param e   Pointer to the Electrode object. Note, this class may use dynamic casting
     *             to choose a child object, and then may invoke an error if the match isn't
     *             correct.
     */
    virtual void copyElectrode_intoState(const Cantera::Electrode* const e);

    void setStateElectrode_SimpleDiff_fromEState(Cantera::Electrode_SimpleDiff* const e) const;
    void setStateElectrode_DiffTALE_fromEState(Cantera::Electrode_DiffTALE* const e) const;


    //! Set the state of the Electrode from the state of this object
    /*!
     *   This function will set the _final_ state of the Electrode object.
     *   This function will also set all of the internal degrees of freedom such that
     *   the next time step taken by the object will result in the same values and time step
     *   that would have occurred if there hadn't been a restart.
     *
     *  @param  e  changeable pointer to the electrode object.
     */
    virtual void setStateElectrode_fromEState(Cantera::Electrode*   const e) const;

    /* --------------------------------------------------------------------------------------  */
protected:

    int numRCells_;

    //! Define the number of species that are defined to have radially distributed distributions
    //! within the solid
    /*!
     *   Note, for this object there is only one radial distribution. This concept will be enhanced
     *   in later formulations.
     *
     *   There are a few arrays which have numKRSpecies_ as their inner loop dimension.
     */
    int numKRSpecies_;

    //! Number of phases which have radial distributions of their species
    int numSPhases_;

    //  std::vector<int> numSpeciesPerSPhase_;



    //! Node position of the mesh - final
    std::vector<doublereal> rnodePos_;  

    std::vector<doublereal> cellBoundR_;

 
    //! Total concentration of each of the solid phases that are distributed - global init state
    /*!
     *  concTot_SPhase_Cell_init_init_[numSPhases_ * iCell + iJRPh]
     */
    std::vector<double> concTot_SPhase_Cell_;

    //! total concentration of the solid phases that are distributed - init state
    /*!
     *     concKRSpecies_Cell_init_init_[numKRSpecies_ * iCell + iKRSpecies]
     */
    std::vector<double> concKRSpecies_Cell_;

    //! Vector of solid species defined on the grid of the spherical particle
    /*!
     *
     *   The inner loop is over the number of species that are defined to have radially dependent
     *   distributions, numKRSpecies;
     *   init_init state value
     */
    std::vector<double> spMoles_KRsolid_Cell_;

    //! Mole fraction of the solid phase species that are distributed = final state
    /*!
     *
     */
    //std::vector<double> spMf_KRSpecies_Cell_;

    //! This integer describes if the system is current on a Region boundary at the end of a subgrid
    //! integration step
    /*!
     *   We define a region boundary here iff all cells are on that boundary
     *  If it isn't, we set it to -1. If it is, we set it to the region boundary index
     */
    int onRegionBoundary_;

    //! Statement that the Electrode class can access any information in this class
    /*!
     *  NOTE, I'm not sure that this direction of access is needed ATM.
     */
    friend class Cantera::Electrode;
    friend class Cantera::Electrode_SimpleDiff;
};

}
#endif
