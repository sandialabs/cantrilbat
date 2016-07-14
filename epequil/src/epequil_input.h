/*
 * $Id: epequil_input.h 508 2013-01-07 22:54:04Z hkmoffa $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#ifndef _EPEQUIL_INPUT_H
#define _EPEQUIL_INPUT_H

#include "cantera/equilibrium.h"
#include "tok_input_util.h"
#include <string>
/*
 *-----------------------------------------------------------------------------
 *
 * Include file containing constant declarations for inputs to 
 * epequil
 *
 *-----------------------------------------------------------------------------
 */
#define EPEQUIL_MAX_NAME_LEN_P1 81
#define EPEQUIL_MAX_NAME_LEN    80

#define EPEQUIL_SUCCESS 0

#ifdef useZuzaxNamespace
#define ZZCantera Zuzax
#else
#define ZZCantera Cantera
#endif


/*
 * Command file input -> This is the current command file specification
 *                       of the problem statement.
 */
class EPEQUIL_KEY_INPUT {
public:
    EPEQUIL_KEY_INPUT() ;
    ~EPEQUIL_KEY_INPUT();

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
extern EPEQUIL_KEY_INPUT PO;

/************************************************************************
 * Complete problem statement
 */
class EPEQUIL_INPUT {
public:
    EPEQUIL_INPUT();
    ~EPEQUIL_INPUT();

    ZZCantera::ThermoPhase **tplist;

    ZZCantera::MultiPhase *m_mp;
    int    prob_type; /* prob_type = Problem type. I.e., the identity of 
		       *             what is held constant. Currently, 
		       *             T and P are held
		       *             constant, and this input is ignored */
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
    double  *elementMoles;
    double  *phaseMoles;
    char    *Title;
    int     iest;
    bool  specifiedElementAbundances;

    /*********************************************************************/
    /*
     * Functions associated with classes
     */
 
};

extern EPEQUIL_KEY_INPUT Key;
extern int epequil_input(EPEQUIL_INPUT *, std::string commandFile);


/*****************************************************************************/
/*****************************************************************************/
#endif 
/*****************************************************************************/
