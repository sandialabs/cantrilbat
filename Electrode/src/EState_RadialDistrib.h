/*
 *  @file EState.h
 */

/*
 * $Id: EState_RadialDiffusion.h 571 2013-03-26 16:44:21Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _ESTATE_RADIALDIFFUSION_H
#define _ESTATE_RADIALDIFFUSION_H


#include "EState.h"

namespace Cantera
{

class Electrode;
class Electrode_SimpleDiff;


//! Child Class for the EState class concept refering to SimpleDiff Electrode object.
/*!
 *    This class adds the inner_platNum and second_platNum records to the output file.
 *
 *   It's an example of how to write child classes for EState.
 */
class EState_RadialDiffusion  : public EState
{

public:

    //! Default constructor for the base object
    EState_RadialDiffusion();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    EState_RadialDiffusion(const EState_RadialDiffusion& right);

    //! Destructor is virtual
    virtual ~EState_RadialDiffusion();

    //! Assignment operator
    /*!
     *  @param right  object to be duplicated
     */
    EState_RadialDiffusion& operator=(const EState_RadialDiffusion& right);

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
    int initialize(const Cantera::Electrode_SimpleDiff* const e);

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

    //! Mapping to the first inner radius
    /*!
     *  Note, this can be -1, which indicates that there isn't an inner core region.
     *  Note both inner_platNum_ and second_platNum_ can't both be -1; that would be an error condition
     */
    int inner_platNum_;

    //! Mapping to the second inner radius
    /*!
     *  Note, this can be -1, which indicates that there isn't an outer annular region
     */
    int second_platNum_;

    //! Statement that the Electrode class can access any information in this class
    /*!
     *  NOTE, I'm not sure that this direction of access is needed ATM.
     */
    friend class Cantera::Electrode;
    friend class Cantera::Electrode_SimpleDiff;
};

}
#endif
