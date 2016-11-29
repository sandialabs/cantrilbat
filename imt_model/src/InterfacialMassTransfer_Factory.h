/**
 *  @file ThermoFactory.h
 *     Headers for the factory class that can create known %ThermoPhase objects
 *     (see \ref thermoprops and class \link ZZCantera::ThermoFactory ThermoFactory\endlink).
 *
 */

/*
 * $Author: hkmoffa $
 * $Revision: 507 $
 * $Date: 2013-01-07 15:48:29 -0700 (Mon, 07 Jan 2013) $
 */

// Copyright 2001  California Institute of Technology
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tok_input_util.h"


#include "mdp_allo.h"
#include "cantera/equilibrium.h"

#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/thermo/FixedChemPotSSTP.h"

#include "cantera/numerics/solveProb.h"
#include "cantera/numerics/BEulerInt.h"

#include "cantera/solvers.h"

//#include "PhaseList.h"
#include "BlockEntryGlobal.h"
#include "InterfacialMassTransfer.h"
#include "InterfacialMassTransfer_input.h"
#include "ApplBase_print.h"


#ifndef INTERFACIALMASSTRANSFER_FACTORY_H
#define INTERFACIALMASSTRANSFER_FACTORY_H

//#include "ThermoPhase.h"
//#include "xml.h"

#if defined(THREAD_SAFE_CANTERA)
#include <boost/thread/mutex.hpp>
#endif

#include "cantera/base/FactoryBase.h"

using namespace std;

#ifdef useZuzaxNamespace
namespace Zuzax
#else
#ifdef useZuzaxNamespace
namespace Zuzax
#else
namespace Cantera
#endif
#endif 
{

  class SpeciesThermoFactory;
  class VPSSMgr;


  //! Enum to String routine for the enum IMT_Types_Enum 
  /*!
   *  @param etype The model of the electrode
   *
   *  @return Returns the characteristic string for that interfacial mass transport Model
   */
  std::string to_string(const IMT_Types_Enum & etype);

  //! String to Enum Routine for the enum IMT_Types_Enum 
  /*!
   *  Matches are first made using case. Then, they are made by ignoring case
   *
   *  @param       input_string
   *
   *  @return      Returns the Enum type for the string
   */
  IMT_Types_Enum string_to_IMT_Types_Enum(const std::string & input_string);


  //! Factory class for thermodynamic property managers.
  /*!
   * This class keeps a list of the known ThermoPhase classes, and is
   * used to create new instances of these classes.
   */
  class InterfacialMassTransfer_Factory : public FactoryBase {

  private:
    //! Private constructor
    InterfacialMassTransfer_Factory();
       
    //!  Private virtual destructor
    /*!
     * We do not delete statically created single instance of this
     * class here, because it would create an infinite loop if
     * destructor is called for that single instance.
     */
    virtual ~InterfacialMassTransfer_Factory();

  public:

    //! Static function that creates a static instance of the factory.
    static InterfacialMassTransfer_Factory* factory() {
#if defined(THREAD_SAFE_CANTERA)
      boost::mutex::scoped_lock lock(thermo_mutex);
#endif
      if (!s_factory) s_factory = new InterfacialMassTransfer_Factory;
      return s_factory;
    }

    //! delete the static instance of this factory
    virtual void deleteFactory() {
#if defined(THREAD_SAFE_CANTERA)
      boost::mutex::scoped_lock lock(thermo_mutex);
#endif
      if (s_factory) {
	delete s_factory;
	s_factory = 0;
      }
    }



    //! Create a new thermodynamic property manager.
    /*!
     * @param model  String to look up the model against
     *
     * @return 
     *   Returns a pointer to a new ThermoPhase instance matching the
     *   model string. Returns NULL if something went wrong.
     *   Throws an exception UnknownThermoPhaseModel if the string
     *   wasn't matched.
     */
    virtual InterfacialMassTransfer* newInterfacialMassTransferObject(std::string model);

    //!    Create a new IMT_KEY_INPUT Object given a model name
    /*
     * @param model   String to look up the model 
     * @param f       ThermoFactor instance to use in matching the string
     *
   * @return 
   *   Returns a pointer to a new IMT_KEY_INPUT instance matching the
   *   model string for the imt object. Returns NULL if something went wrong.
   *   Throws an exception  if the string
   *   wasn't matched.
   */
  IMT_KEY_INPUT * newIMTKeyInputObject(std::string model);
  
  public:
    //! static member of a single instance
    static InterfacialMassTransfer_Factory * s_factory;


#if defined(THREAD_SAFE_CANTERA)
    //! Decl for locking mutex for thermo factory singelton
    static boost::mutex electrode_mutex;
#endif

  };

  
  //!  Create a new thermo manager instance.
  /*!
   * @param model   String to look up the model against
   * @param f       ThermoFactor instance to use in matching the string
   *
   * @return 
   *   Returns a pointer to a new ThermoPhase instance matching the
   *   model string. Returns NULL if something went wrong.
   *   Throws an exception UnknownThermoPhaseModel if the string
   *   wasn't matched.
   */ 
  inline InterfacialMassTransfer * newInterfacialMassTransferObject(std::string model,  
					InterfacialMassTransfer_Factory* f=0) {
    if (f == 0) {
      f = InterfacialMassTransfer_Factory::factory();
    }
    return f->newInterfacialMassTransferObject(model);
  }
  
  //!  Create a new IMT_KEY_INPUT Object
  /*!
   * @param model   String to look up the model against
   * @param f        instance to use in matching the string
   *
   * @return 
   *   Returns a pointer to a new IMT_KEY_INPUT instance matching the
   *   model string for the IMT object. Returns NULL if something went wrong.
   *   Throws an exception  if the string
   *   wasn't matched.
   */ 
   IMT_KEY_INPUT * newIMTKeyInputObject(std::string model, InterfacialMassTransfer_Factory* f = 0);

}
#endif 
