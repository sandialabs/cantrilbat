/*
 * $Id: mpequil_input.h 502 2013-01-07 22:25:47Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _MPEQUIL_INPUT_H
#define _MPEQUIL_INPUT_H

#include "cantera/equilibrium.h"
#include "tok_input_util.h"
#include <string>
#include <vector>

#ifdef useZuzaxNamespace
#define ZZCantera Zuzax
#else
#define ZZCantera Cantera
#endif

/*
 *-----------------------------------------------------------------------------
 *
 * Include file containing constant declarations for inputs to 
 * mpequil
 *
 *-----------------------------------------------------------------------------
 */
#define MPEQUIL_MAX_NAME_LEN_P1 81
#define MPEQUIL_MAX_NAME_LEN    80

#define MPEQUIL_SUCCESS 0


//! storage for Command file input
/*!
 * This is the current command file specification
 *                       of the problem statement.
 */
class MPEQUIL_KEY_INPUT {
public:
    MPEQUIL_KEY_INPUT() ;
    ~MPEQUIL_KEY_INPUT();

    std::string CanteraFN1;
    int NumberCanteraFiles;
    char **CanteraFileNames;
    double Temperature;
    double Pressure;
    double Vol;
    double *MoleNumber;
    double * MoleNumberIG;
    int    *PhaseInclude;
    int    ProblemType;
    char **SpeciesNames;
    char **PhaseNames;
    char **ElementNames;
    double *ElementAbundances;
    bool  specifiedElementAbundances;
    std::string Title;
    void InitForInput(ZZCantera::MultiPhase *);
	
};
extern MPEQUIL_KEY_INPUT PO;

/************************************************************************
 * Complete problem statement
 */
class MPEQUIL_INPUT {
public:
  MPEQUIL_INPUT();
  ~MPEQUIL_INPUT();

  ZZCantera::ThermoPhase **tplist;

  ZZCantera::MultiPhase *m_mp;

  //* Integer representing the Problem type.
  /*!
   *  The identity of  what is held constant. Currently, 
   *   T and P are held constant, and this input is ignored 
   */
  int    prob_type;

  int     nspecies; /* nspecies   = Total number of species in the 
		     *              problems*/
 
  int     ne;       /* ne         = Number of elements in the problem  */
 
  int     nphase;   /* NPhase     = Number of phases in the problem    */
 
  double   T;       /* T          = Temperature (K)                    */
  double   Pres;    /* Pres       = Pressure (MKS)                     */
  double   Vol;     /* Vol        = Volume (cm^3)                      */
  double  *VolPM;   /* VolPM[k]   = Partial Molar Volumes of species   */
  double  *spMoles;
  double  *spMf;
  double  *spChemPot;

  //! Vector containing the number of moles of each element
  /*! 
   * It has a length equal to ne, the number of elements that are in the
   * multiphase object and the order is the same as the elements in the
   * multiphase object.
   *
   * If specifiedElementAbundances is true, then these numbers are not
   * calculated from the species moles, and actually become the 
   * controlling input into the equilibrium calculation. 
   */
  double  *elementMoles;

  double  *phaseMoles;
  char    *Title;

  //! Boolean indicating whether this routine should try to find its
  //! own initial estimate. 
  /*!
   *  If true, then we start from scratch in terms of the starting mole
   *  numbers. If false, then we use the mole numbers of the species as
   *  a starting point.
   */
  int     iest;
  bool  specifiedElementAbundances;

};

extern MPEQUIL_KEY_INPUT Key;
extern int mpequil_input(MPEQUIL_INPUT *, std::string commandFile);


#endif 
/*****************************************************************************/
