/**
 *  @file example2.cpp
 *
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cttables_kin.h"
#include "cantera/kinetics/solveSP.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

using std::cout;
using std::endl;

//======================================================================================================================
/*
 * CONVENTION ON THE CURRENT USED IN CTTABLES
 *----------------------------------------------------
 *
 *  The current always refers to the current from the metal into the solution. i is positive when the overptential
 *  of the reaction is positive. When i is positive, electrons in the metal are created.
 *  When the overpotential is positive the potential of the metal is higher than the solution.
 *  Making the metal potential postivie means that the electrons have a lower chemical potential thus
 *  favoring their formation.
 *
 *  Some useful equations:
 *
 *   Voltage = Phi_Metal - Phi_Soln
 *
 *   Overpotential = Voltage - Voltage_eq
 *
 *   I = d[e-] / dt
 *
 */
//======================================================================================================================

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/*
 * Calculate the arrhenius parameters from a fit from two points.
 * Ea/R is returned. This has units of Kelvin.
 *
 */
void calcArrhParam(double k1, double T1,
                   double k2, double T2,
                   double& Afac, double& EadivR)
{
    if (T1 <= 0.0 || T2 <= 0.0) {
        Afac = k1;
        EadivR = 0.0;
        return;
    }
    if (T1 == T2) {
        Afac = k1;
        EadivR = 0.0;
        return;
    }
    if (k1 <= 0.0 || k2 <= 0.0) {
        Afac = k1;
        EadivR = 0.0;
        return;
    }
    EadivR = log(k1/k2) / (1.0/T2 - 1.0/T1);
    Afac = exp(log(k1) + EadivR/T1);
    return;
}

//!  Fit two current vs voltage points to an exponential relation
//!  based on the overpotential
/*!
 *   if iform = 1, we solve for the anodic branch
 *
 *     i = +- i0 exp (alpha_a F nu / RT)
 *
 *         where nu = V - Erxn
 *               F = faraday's constant
 *               R = Gas Constant
 *               T = temperature
 *     -> solve for i0 and alpha_a
 *
 *  if iform = 0, we solve for the cathodic branch
 *
 *     i = +- i0 exp (- alpha_c F nu / RT)
 *
 *         where nu = V - Erxn
 *               F = faraday's constant
 *               R = Gas Constant
 *               T = temperature
 *     -> solve for i0 and alpha_c
 *
 *  Currently, the signs are passed through unchanged
 *
 */
void calcBF(double i1, double V1, double i2,
            double V2, double Erxn, int iform, double FdRT,
            double& i0, double& alpha)
{

    double n1 = V1 - Erxn;
    double n2 = V2 - Erxn;
    double sgn = 1.0;
    AssertTrace(iform == 1 || iform == 0);
    if (iform == 1) {
        // anodic branch
        if (i1 < 0.0) {
            if (i2 < 0.0) {
                sgn = -1.0;
                i1 = -i1;
                i2 = -i2;
            } else {
                throw CanteraError("calcBF"," two currents have different signs");
            }
        }
        double alpha_prime = (log(i1) - log(i2)) / (n1 - n2);
        alpha = alpha_prime / FdRT;
        i0 = exp(log(i1) - n1 * alpha_prime);
    } else {
        // Cathodic branch
        if (i1 < 0.0) {
            if (i2 < 0.0) {
                sgn = -1.0;
                i1 = -i1;
                i2 = -i2;
            } else {
                throw CanteraError("calcBF", " two currents have different signs");
            }
        }
        double alpha_prime = - (log(i1) - log(i2)) / (n1 - n2);
        alpha = alpha_prime / FdRT;
        i0 = exp(log(i1)  + n1 * alpha_prime);
    }
    i0 = i0 * sgn;

    return;
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 * Setup the tables for the reactions.
 */
void
getKineticsTables(TemperatureTable& TT, PhaseList* pl,
                  ZZCantera::Kinetics& kin,
                  DenseMatrix& kfwd_Table,
                  DenseMatrix& krev_Table,
                  DenseMatrix& deltaG_Table,
                  DenseMatrix& deltaH_Table,
                  DenseMatrix& deltaS_Table,
                  DenseMatrix& Afwd_Table,
                  DenseMatrix& EafwddivR_Table,
                  DenseMatrix& Arev_Table,
                  DenseMatrix& EarevdivR_Table,
                  DenseMatrix& kfwdPrime_Table,
                  DenseMatrix& krevPrime_Table)
{
    int i, j, iph;
    int nReactions = kin.nReactions();
    int nPhases = kin.nPhases();
    double* Rarray = new double [nReactions];
    double deltaT = 0.01;
    double Afac, EadivR, csc;
    int kindex;

    double T, Tdelta;
    for (i = 0; i < TT.size(); i++) {
        T = TT[i];
        for (iph = 0; iph < nPhases; iph++) {
            ThermoPhase& gRef = kin.thermo(iph);
            size_t kstart = pl->globalSpeciesIndex(&gRef);
            gRef.setState_TPX(T, BG.Pressure, BG.XmolPLSpecVec +kstart);
        }

        kin.getFwdRateConstants(Rarray);
        for (j = 0; j < nReactions; j++) {
            kfwd_Table(i, j) = Rarray[j];
        }

        kin.getRevRateConstants(Rarray);
        for (j = 0; j < nReactions; j++) {
            krev_Table(i, j) = Rarray[j];
        }

        kin.getDeltaSSGibbs(Rarray);
        for (j = 0; j < nReactions; j++) {
            deltaG_Table(i, j) = Rarray[j];
        }

        kin.getDeltaSSEnthalpy(Rarray);
        for (j = 0; j < nReactions; j++) {
            deltaH_Table(i, j) = Rarray[j];
        }

        kin.getDeltaSSEntropy(Rarray);
        for (j = 0; j < nReactions; j++) {
            deltaS_Table(i, j) = Rarray[j];
        }

        Tdelta = T + deltaT;
        for (iph = 0; iph < nPhases; iph++) {
            ThermoPhase& gRef = kin.thermo(iph);
            int kstart = pl->globalSpeciesIndex(&gRef);

            gRef.setState_TPX(Tdelta, BG.Pressure, BG.XmolPLSpecVec +kstart);
        }


        kin.getFwdRateConstants(Rarray);
        for (j = 0; j < nReactions; j++) {
            calcArrhParam(kfwd_Table(i, j), T, Rarray[j], Tdelta, Afac, EadivR);
            Afwd_Table(i,j) = Afac;
            EafwddivR_Table(i,j) = EadivR;
        }

        kin.getRevRateConstants(Rarray);
        for (j = 0; j < nReactions; j++) {
            calcArrhParam(krev_Table(i, j), T, Rarray[j], Tdelta, Afac, EadivR);
            Arev_Table(i,j) = Afac;
            EarevdivR_Table(i,j) = EadivR;
        }

        for (j = 0; j < nReactions; j++) {
            /*
             * reactants
             */
            const vector<size_t>& reactants = kin.reactants(j);
            int nReac = reactants.size();
            /*
             * products
             */
            const vector<size_t>& products = kin.products(j);
            int nProd = products.size();

            kfwdPrime_Table(i, j) = kfwd_Table(i,j);
            for (int k = 0; k < nReac; k++) {
                kindex = reactants[k];
                int iphase = kin.speciesPhaseIndex(kindex);
                ThermoPhase& tpRef = kin.thermo(iphase);

                int kLocal = kindex - kin.kineticsSpeciesIndex(0, iphase);
                csc = tpRef.standardConcentration(kLocal);
                kfwdPrime_Table(i, j) *= csc;
            }

            krevPrime_Table(i, j) = krev_Table(i,j);
            for (int k = 0; k < nProd; k++) {
                kindex = products[k];
                int iphase = kin.speciesPhaseIndex(kindex);
                ThermoPhase& tpRef = kin.thermo(iphase);
                int kLocal = kindex - kin.kineticsSpeciesIndex(0, iphase);
                csc = tpRef.standardConcentration(kLocal);
                krevPrime_Table(i, j) *= csc;
            }
        }
    }

    for (iph = 0; iph < nPhases; iph++) {
        ThermoPhase& gRef = kin.thermo(iph);
        size_t kstart = pl->globalSpeciesIndex(&gRef);
        gRef.setState_TPX(BG.Temperature, BG.Pressure, BG.XmolPLSpecVec + kstart);
    }

    delete [] Rarray;
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

std::string iunit_string(int i)
{
    std::string val;
    switch (i) {
    case 0:
        val = "kmol";
        break;
    case 1:
        val = "m";
        break;
    case 2:
        val = "Kg";
        break;
    case 3:
        val = "Pa";
        break;
    case 4:
        val = "K";
        break;
    case 5:
        val = "s";
        break;
    default:
        val = "";
        break;
    }
    return val;
}

bool isInteger(double num)
{
    int inum = (int) num;
    double rnum = inum;
    if (rnum == num) {
        return true;
    }
    return false;
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

std::string formUnitsString(double unitsA[6])
{
    std::string num;
    std::string denom;
    std::string val;
    int numNum = 0;
    int numDenom = 0;
    for (int i = 0; i < 6; i++) {
        if (unitsA[i] > 0.0) {
            numNum++;
        }
        if (unitsA[i] < 0.0) {
            numDenom++;
        }
    }
    if ((numNum + numDenom) == 0) {
        return "";
    }
    if (numNum > 0) {
        int itmp = 0;
        for (int i = 0; i < 6; i++) {
            if (unitsA[i] > 0.0) {
                val = iunit_string(i);
                if (itmp > 0) {
                    num = num + " ";
                }
                num += val;
                if (isInteger(unitsA[i])) {
                    int np = (int) unitsA[i];
                    if (np != 1) {
                        num += "^" + int2str(np);
                    }
                } else {
                    if (unitsA[i] == 0.5) {
                        num += "^0.5";
                    } else {
                        num += "^" + fp2str(unitsA[i]);
                    }
                }
                itmp++;
            }
        }
    } else {
        num = "1";
    }
    if (numDenom > 0) {
        int idtmp = 0;
        for (int i = 0; i < 6; i++) {
            if (unitsA[i] < 0.0) {
                val = iunit_string(i);
                double dp = - unitsA[i];
                if (idtmp > 0) {
                    denom = denom + " ";
                }
                denom += val;
                if (isInteger(dp)) {
                    int np = (int) dp;
                    if (np != 1) {
                        denom += "^" + int2str(np);
                    }
                } else {
                    if (dp == 0.5) {
                        denom += "^0.5";
                    } else {
                        denom += "^" + fp2str(dp);
                    }
                }
                idtmp++;
            }
        }
    }
    if (numNum > 0 && numDenom > 0) {
        if (numNum > 1) {
            std::string tmp = "(" + num + ")";
            num = tmp;
        }
        if (numDenom > 1) {
            std::string tmp = "(" + denom + ")";
            denom = tmp;
        }
    }
    if (numDenom > 0) {
        val = num + "/" + denom;
    } else {
        val = num;
    }
    return val;
}

//==================================================================================================================================
/*
 *  This routine will print out a table of information a single
 *  reaction, j.
 *
 *
 */
void printKineticsTable(PhaseList* pl, int j,
                        TemperatureTable& TT,
                        ZZCantera::Kinetics& kin,
                        DenseMatrix& kfwd_Table,
                        DenseMatrix& krev_Table,
                        DenseMatrix& deltaG_Table,
                        DenseMatrix& deltaH_Table,
                        DenseMatrix& deltaS_Table,
                        DenseMatrix& Afwd_Table,
                        DenseMatrix& EafwddivR_Table,
                        DenseMatrix& Arev_Table,
                        DenseMatrix& EarevdivR_Table,
                        DenseMatrix& kfwdPrime_Table,
                        DenseMatrix& krevPrime_Table,
                        ZZCantera::RxnMolChange* rmc)
{
    /*
     *  Get the species data object from the Mixture object
     *  this is defined in the constituents.h file, and is
     *  inherited by Mixture through BaseMix
     */

    int tableWidth = 128;

    /*
     * Conversion factor between Joules/kmol, cantera's units,
     * and kcal / gmol.
     */
    double conv = 1.0 / (4.184E6);
    double convS = conv * 1.0E3;
    double Rkcal = 1.98721E-3;
    // Possible change the conv units to kJ / gmol
    if (UIO.unitDef == UNITS_KJOULE) {
        conv = 1.0 / 1.0E6;
        convS = conv * 1.0E3;
        Rkcal = 8.314472E-3;
    }


    ThermoPhase* gThermo = &kin.thermo(0);
    /*
     *  Dump out all of the information about the reaction
     */
    cout << endl;
    print_char('=', tableWidth);
    cout << endl;
    cout << endl;
    cout << "INFORMATION TABLE FOR Reaction \""  << j;
    cout << "\" IN PHASE \"";
    cout << gThermo->id();
    cout << "\"" << endl;
    cout << endl;

    /*
     * Reaction String
     */
    std::string rs = kin.reactionString(j);
    dnt(1);
    cout << rs << endl;
    cout << endl;

    /*
     * Print out reactants
     */
    const vector<size_t>& reactants = kin.reactants(j);
    int nReac = reactants.size();
    dnt(1);
    if (nReac == 1) {
        cout << "There is one reactant: ";
    } else {
        cout << "There are " << nReac << " reactants: ";
    }
    for (int k = 0; k < nReac; k++) {
        int kindex = reactants[k];
        cout << kin.kineticsSpeciesName(kindex) << " ";
    }
    cout << endl;


    /*
     * Print out products
     */
    const std::vector<size_t>& products = kin.products(j);
    int nProd = products.size();
    dnt(1);
    if (nProd == 1) {
        cout << "There is one product: ";
    } else {
        cout << "There are " << nProd << " products: ";
    }
    for (int k = 0; k < nProd; k++) {
        int kindex = products[k];
        cout << kin.kineticsSpeciesName(kindex) << " ";
    }
    cout << endl;

    /*
     *  reaction type
     */
    double pCurrent = gThermo->pressure();
    int rType = kin.reactionType(j);
    dnt(1);
    cout << "Reaction Type = ";
    switch (rType) {
    case ELEMENTARY_RXN:
        cout << "Elementary Rxn" << endl;
        break;
    case THREE_BODY_RXN:

        cout << "Three Body Rxn: "<< endl;
        dnt(1);
        cout <<"        Rate Constants below include third body collision partner "
             << endl;
        dnt(1);
        cout <<"        and are effected by current pressure, "
             << pCurrent << " Pa " << endl;
        break;
    case FALLOFF_RXN:
        cout << "Falloff Rxn" << endl;
        dnt(1);
        cout <<"        Rate Constants below include third body collision partner effects"
             << endl;
        dnt(1);
        cout <<"        and are effected by current pressure, "
             << pCurrent << " Pa " << endl;
        break;
    case CHEMACT_RXN:
        cout << "Chemically activated rxn" << endl;
        dnt(1);
        cout <<"        Rate Constants below include third body collision partner effects"
             << endl;
        dnt(1);
        cout <<"        and are effected by current pressure, "
             << pCurrent << " Pa " << endl;
        break;
    case SURFACE_RXN:
        cout << "Surface rxn" << endl;
        break;
    case BUTLERVOLMER_NOACTIVITYCOEFFS_RXN:
        cout << "Butler Volmer rxn without Activity Coefficients" << endl;
        break;
    case BUTLERVOLMER_RXN:
        cout << "Butler Volmer rxn" << endl;
        break;
    case SURFACEAFFINITY_RXN:
        cout << "Surface Affinity rxn" << endl;
        break;
    case EDGE_RXN:
        cout << "Edge Rxn" << endl;
        break;
    case GLOBAL_RXN:
        cout << "Global Rxn" << endl;
        break;
    default:
        cout << "Unknown_Type - (" << rType << ")" << endl;
        break;
    }


    /*
     * Reaction Reversible?
     */
    bool rReversible = kin.isReversible(j);
    dnt(1);
    if (rReversible) {
        cout << "Reaction is reversible" << endl;
    } else {
        cout << "Reaction is not reversible" << endl;
    }

    cout << endl;



    /*
     * Figure out the units for the forward reaction rate constant
     * First start out with the units for the rop
     */
    int nphase = kin.nPhases();
    int ndim = 3;
    for (int iph = 0; iph < nphase; iph++) {
        ThermoPhase& tpRef = kin.thermo(iph);
        int idim = tpRef.nDim();
        if (idim < ndim) {
            ndim = idim;
        }
    }
    double unitsROP[6] = { 1.0, double(-ndim), 0.0, 0.0, 0.0, -1.0 };
    double unitskfwd[6];
    double unitskrev[6] ;
    double unitsSpecies[6];
    for (int i = 0; i < 6; i++) {
        unitskfwd[i] = unitsROP[i];
        unitskrev[i] = unitsROP[i];
    }
    for (int k = 0; k < nReac; k++) {
        int kindex = reactants[k];
        int iphase = kin.speciesPhaseIndex(kindex);
        ThermoPhase& tpRef = kin.thermo(iphase);
        tpRef.getUnitsStandardConc(unitsSpecies, kindex);
        for (int i = 0; i < 6; i++) {
            unitskfwd[i] -= unitsSpecies[i];
        }
    }
    for (int k = 0; k < nProd; k++) {
        int kindex = products[k];
        int iphase = kin.speciesPhaseIndex(kindex);
        ThermoPhase& tpRef = kin.thermo(iphase);
        tpRef.getUnitsStandardConc(unitsSpecies, kindex);
        for (int i = 0; i < 6; i++) {
            unitskrev[i] -= unitsSpecies[i];
        }
    }

    /****************** phases change table ****************************/
    int pwidth = 52;
    if (rmc->m_ChargeTransferInRxn != 0.0) {
        dnt(1);
        printf("This reaction is charge transfer reaction, nStoich = %g\n",
               rmc->m_ChargeTransferInRxn);
        dnt(1);
        printf("Electrochemical Beta = %g\n", rmc->m_beta);
        pwidth = 86;
    }

    dnt(2);
    print_char('-', pwidth);
    printf("\n");

    dnt(2);
    printf("|              Phase | Dim Delta_Moles  Delta_Mass |");
    if (rmc->m_ChargeTransferInRxn != 0.0) {
        printf("Electron_Transfer Base_Potential |");
    }
    printf("\n");
    for (size_t iph = 0; iph < rmc->m_nPhases; iph++) {
        ThermoPhase& tpRef = kin.thermo(iph);
        std::string sname = tpRef.name();
        dnt(2);
        printf("|%20s|", sname.c_str());
        printf("%4d",    rmc->m_phaseDims[iph]);
        printf("%12.5g", rmc->m_phaseMoleChange[iph]);
        printf("%12.5g", rmc->m_phaseMassChange[iph]);
        if (rmc->m_ChargeTransferInRxn != 0.0) {
            printf(" |");
            printf("  %12.5g", rmc->m_phaseChargeChange[iph]);
            ThermoPhase& tpp = kin.thermo(iph);
            double volts = tpp.electricPotential();
            printf("      %12.5g", volts);
        }
        printf(" |\n");
    }

    dnt(2);
    print_char('-', pwidth);
    printf("\n\n");

    if (rType == SURFACEAFFINITY_RXN) {
        printAffinityHeader(rmc, pl, j, TT, kin, kfwd_Table, krev_Table,  deltaG_Table, deltaH_Table, deltaS_Table,
                            Afwd_Table, EafwddivR_Table,   Arev_Table, EarevdivR_Table,
                            kfwdPrime_Table, krevPrime_Table, &unitskfwd[0], &unitskrev[0]);
    }


    /*******************************************************************/


    std::string units_kfwd = formUnitsString(unitskfwd);
    std::string units_krev = formUnitsString(unitskrev);
    dnt(1);
    cout << "units_kfwd = Units for forward reaction rate constant = "
         << units_kfwd << endl;
    dnt(1);
    cout << "units_krev = Units for reverse reaction rate constant = "
         << units_krev << endl;

    /*
     * Dump out the kinetics table as a function of the
     * temperature
     */
    cout << endl;
    cout << "|";
    print_char('-', tableWidth);
    cout << "|" << endl;
    cout << "|        |                              "
         << " |                                 ";
    if (rReversible) {
        cout << "|                               ";
    } else {
        cout << "| -- REVERSE RXN NOT IN MECH -- ";
    }
    cout << "|                     ";
    cout << "|" << endl;
    cout << "|  Temp  |   kfwd      Afwd     Ea_fwd  "
         << " |  DeltaG0    DeltaH0    DeltaS0  "
         << "|  k_rev     A_rev     Ea_rev   "
         << "|  kstar   kstar_rev  ";

    cout << "|" << endl;


    if (UIO.unitDef == UNITS_KJOULE) {
        cout << "|   (K)  |   (units_kfwd)       (kJ/gmol)"
             << "|  (kJ/gmol) (kJ/gmol)  (J/gmolK) "
             << "|  (units_krev)       (kJ/gmol) |"
             << "    (kmol / m^" << ndim << " s)   |";
        cout << endl;
    } else {
        cout << "|   (K)  |   (units_kfwd)     (kcal/gmol)"
             << "|(kcal/gmol)(kcal/gmol)(cal/gmolK)"
             << "|  (units_krev)      (kcal/gmol)|"
             << "    (kmol / m^" << ndim << " s)   |";
        cout << endl;
    }

    cout << "|--------|----------"
         << "---------------------|-----------------"
         << "----------------|-----------------------------"
         << "--|---------------------|" << endl;

    for (int i = 0; i < TT.size(); i++) {
        cout << "|";
        pr_df(TT[i], 7, 1);
        cout << " |";
        pr_de(kfwd_Table(i,j), 10, 2);
        pr_de(Afwd_Table(i,j), 10, 2);
        pr_dg(EafwddivR_Table(i,j) * Rkcal, 10, 2);
        cout << " | ";
        pr_dg(deltaG_Table(i,j) * conv, 10, 2);
        pr_dg(deltaH_Table(i,j) * conv, 10, 2);
        pr_dg(deltaS_Table(i,j) * convS, 10, 2);
        cout << "  |";
        pr_de(krev_Table(i,j), 10, 2);
        pr_de(Arev_Table(i,j), 10, 2);
        pr_dg(EarevdivR_Table(i,j) * Rkcal, 10, 2);
        cout << " |";
        pr_de(kfwdPrime_Table(i,j), 10, 2);
        pr_de(krevPrime_Table(i,j), 10, 2);
        cout << " |" << endl;
    }
    cout << "|" ;
    print_char('-', tableWidth);
    cout << "|" << endl;
}
/*************************************************************************/
/**********************************************************************************************************************************/
/*************************************************************************/

void doKineticsTablesHetero(PhaseList* pl,  InterfaceKinetics* gKinetics, TemperatureTable& TT)
{

    size_t nReactions = gKinetics->nReactions();
    // int nPhases = gKinetics->nPhases();
    DenseMatrix kfwd_Table(TT.size(), nReactions);
    DenseMatrix krev_Table(TT.size(), nReactions);
    DenseMatrix deltaG_Table(TT.size(), nReactions);
    DenseMatrix deltaH_Table(TT.size(), nReactions);
    DenseMatrix deltaS_Table(TT.size(), nReactions);
    DenseMatrix Afwd_Table(TT.size(), nReactions);
    DenseMatrix EafwddivR_Table(TT.size(), nReactions);
    DenseMatrix Arev_Table(TT.size(), nReactions);
    DenseMatrix EarevdivR_Table(TT.size(), nReactions);
    DenseMatrix kfwdPrime_Table(TT.size(), nReactions);
    DenseMatrix krevPrime_Table(TT.size(), nReactions);


    getKineticsTables(TT, pl, *gKinetics, kfwd_Table, krev_Table, deltaG_Table, deltaH_Table, deltaS_Table,
                      Afwd_Table, EafwddivR_Table, Arev_Table, EarevdivR_Table, kfwdPrime_Table, krevPrime_Table);

    vector<RxnMolChange*> rmcVector;
    rmcVector.resize(nReactions,0);

    for (size_t i = 0; i < nReactions; i++) {
        rmcVector[i] = new RxnMolChange(gKinetics, i);

        printKineticsTable(pl, i, TT, *gKinetics, kfwd_Table, krev_Table,
                           deltaG_Table, deltaH_Table, deltaS_Table,
                           Afwd_Table, EafwddivR_Table,
                           Arev_Table, EarevdivR_Table,
                           kfwdPrime_Table, krevPrime_Table, rmcVector[i]);


        if ((rmcVector[i])->m_ChargeTransferInRxn != 0.0) {
            processCurrentVsPotTable(rmcVector[i],
                                     pl, i, TT,
                                     *gKinetics, kfwd_Table, krev_Table,
                                     deltaG_Table, deltaH_Table, deltaS_Table,
                                     Afwd_Table, EafwddivR_Table,
                                     Arev_Table, EarevdivR_Table,
                                     kfwdPrime_Table, krevPrime_Table);
        }

        int rType = gKinetics->reactionType(i);
        if (rType == SURFACEAFFINITY_RXN) {
            processAffinityTable(rmcVector[i], pl, i, TT, *gKinetics, kfwd_Table, krev_Table,  deltaG_Table, deltaH_Table,
                                 deltaS_Table,
                                 Afwd_Table, EafwddivR_Table,   Arev_Table, EarevdivR_Table,
                                 kfwdPrime_Table, krevPrime_Table);
        }


    }
    setAllBathSpeciesConditions(pl);

    for (int iextra = 0; iextra < IOO.numExtraGlobalRxns; iextra++) {
        struct EGRInput* egr_ptr = IOO.m_EGRList[iextra];
        ExtraGlobalRxn* egr = new ExtraGlobalRxn(*gKinetics);
        double* RxnVector  = new double[nReactions];
        for (size_t i = 0; i < nReactions; i++) {
            RxnVector[i] = 0.0;
        }
        for (int iErxn = 0; iErxn < egr_ptr->m_numElemReactions; iErxn++) {
            struct ERSSpec* ers_ptr = egr_ptr->m_ERSList[iErxn];
            RxnVector[ers_ptr->m_reactionIndex] = ers_ptr->m_reactionMultiplier;
        }
        egr->setupElemRxnVector(RxnVector,
                                egr_ptr->m_SS_KinSpeciesKindex);

        RxnMolChange* rmcEGR = new RxnMolChange(gKinetics, egr);
        RxnTempTableStuff* rts = new RxnTempTableStuff(-1,0);
        getGERKineticsTables(TT, pl, *gKinetics, *egr, *rts);
        setAllBathSpeciesConditions(pl);
        printAllBathSpeciesConditions(pl);
        printGERKineticsTable(pl, -1, TT, *gKinetics, *egr, rmcEGR, *rts);
        if ((rmcEGR)->m_ChargeTransferInRxn != 0.0) {
            processGERCurrentVsPotTable(rmcEGR, pl, 0, TT,
                                        *gKinetics, *egr, *rts);
        }

        delete egr;
        egr=0;
        delete [] RxnVector;
        delete rmcEGR;
        rmcEGR = 0;
        delete rts;
        rts=0;
    }


    for (size_t i = 0; i < nReactions; i++) {
        delete  rmcVector[i];
    }
}
//==================================================================================================================================
void doKineticsTablesHomog(PhaseList* pl, ZZCantera::Kinetics* gKinetics, TemperatureTable& TT)
{
    int nReactions = gKinetics->nReactions();
    DenseMatrix kfwd_Table(TT.size(), nReactions);
    DenseMatrix krev_Table(TT.size(), nReactions);
    DenseMatrix deltaG_Table(TT.size(), nReactions);
    DenseMatrix deltaH_Table(TT.size(), nReactions);
    DenseMatrix deltaS_Table(TT.size(), nReactions);
    DenseMatrix Afwd_Table(TT.size(), nReactions);
    DenseMatrix EafwddivR_Table(TT.size(), nReactions);
    DenseMatrix Arev_Table(TT.size(), nReactions);
    DenseMatrix EarevdivR_Table(TT.size(), nReactions);
    DenseMatrix kfwdPrime_Table(TT.size(), nReactions);
    DenseMatrix krevPrime_Table(TT.size(), nReactions);

    getKineticsTables(TT, pl,
                      *gKinetics, kfwd_Table, krev_Table,
                      deltaG_Table, deltaH_Table, deltaS_Table,
                      Afwd_Table, EafwddivR_Table,
                      Arev_Table, EarevdivR_Table,
                      kfwdPrime_Table, krevPrime_Table);
    vector<RxnMolChange*> rmcVector;
    rmcVector.resize(nReactions,0);
    for (int i = 0; i < nReactions; i++) {
        rmcVector[i] = new RxnMolChange(gKinetics, i);
        printKineticsTable(pl, i, TT,
                           *gKinetics, kfwd_Table, krev_Table,
                           deltaG_Table, deltaH_Table, deltaS_Table,
                           Afwd_Table, EafwddivR_Table,
                           Arev_Table, EarevdivR_Table,
                           kfwdPrime_Table, krevPrime_Table, rmcVector[i]);
    }

    for (int i = 0; i < nReactions; i++) {
        delete  rmcVector[i];
    }
}
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/


RxnMolChangeLocal::RxnMolChangeLocal(Kinetics* kinPtr, int irxn) :
    m_nPhases(0),
    m_kinBase(kinPtr),
    m_iRxn(irxn),
    m_ChargeTransferInRxn(0.0),
    m_beta(0.0),
    m_egr(0)
{
    int nReac = kinPtr->nReactions();
    int iph;
    AssertTrace(irxn >= 0);
    AssertTrace(irxn < nReac);

    m_nPhases = kinPtr->nPhases();

    m_phaseMoleChange.resize(m_nPhases, 0.0);
    m_phaseMassChange.resize(m_nPhases, 0.0);
    m_phaseChargeChange.resize(m_nPhases, 0.0);
    m_phasePotentials.resize(m_nPhases, 0.0);
    m_phaseTypes.resize(m_nPhases, 0);
    m_phaseDims.resize(m_nPhases, 0);

    int m_kk = kinPtr->nTotalSpecies();

    for (int kKin = 0; kKin < m_kk; kKin++) {
        iph =  m_kinBase->speciesPhaseIndex(kKin);
        ThermoPhase& tpRef = m_kinBase->thermo(iph);
        int kLoc = kKin - m_kinBase->kineticsSpeciesIndex(0, iph);
        double rsc = m_kinBase->reactantStoichCoeff(kKin, irxn);
        double psc = m_kinBase->productStoichCoeff(kKin, irxn);
        double nsc = psc - rsc;
        m_phaseMoleChange[iph] += (nsc);
        double mw = tpRef.molecularWeight(kLoc);
        m_phaseMassChange[iph] += (nsc) * mw;
        double chg = tpRef.charge(kLoc);
        m_phaseChargeChange[iph] += nsc * chg;
    }

    for (iph = 0; iph < m_nPhases; iph++) {
        ThermoPhase& tpRef = m_kinBase->thermo(iph);
        m_phasePotentials[iph] =  tpRef.electricPotential();
        m_phaseDims[iph] = tpRef.nDim();
        m_phaseTypes[iph] = tpRef.eosType();
        if (m_phaseChargeChange[iph] != 0.0) {
            double tmp = fabs(m_phaseChargeChange[iph]);
            if (tmp >  m_ChargeTransferInRxn) {
                m_ChargeTransferInRxn = tmp;
            }
        }
    }

    if (m_ChargeTransferInRxn) {
        InterfaceKinetics* iK = dynamic_cast<InterfaceKinetics*>(kinPtr);
        if (iK) {
            m_beta = iK->electrochem_beta(irxn);
        } else {
            throw CanteraError("RxnMolChange", "unknown condition on charge");
        }
    }

}
//==================================================================================================================================
RxnMolChangeLocal::RxnMolChangeLocal(Kinetics* kinPtr, ExtraGlobalRxn* egr) :
    m_nPhases(0),
    m_kinBase(kinPtr),
    m_iRxn(-1),
    m_ChargeTransferInRxn(0.0),
    m_beta(0.0),
    m_egr(egr)
{
    int iph;
    AssertTrace(egr != 0);

    m_nPhases = kinPtr->nPhases();

    m_phaseMoleChange.resize(m_nPhases, 0.0);
    m_phaseMassChange.resize(m_nPhases, 0.0);
    m_phaseChargeChange.resize(m_nPhases, 0.0);
    m_phasePotentials.resize(m_nPhases, 0.0);
    m_phaseTypes.resize(m_nPhases, 0);
    m_phaseDims.resize(m_nPhases, 0);

    int m_kk = kinPtr->nTotalSpecies();

    for (int kKin = 0; kKin < m_kk; kKin++) {
        iph =  m_kinBase->speciesPhaseIndex(kKin);
        ThermoPhase& tpRef = m_kinBase->thermo(iph);
        int kLoc = kKin - m_kinBase->kineticsSpeciesIndex(0, iph);
        double rsc = egr->reactantStoichCoeff(kKin);
        double psc = egr->productStoichCoeff(kKin);
        double nsc = psc - rsc;
        m_phaseMoleChange[iph] += (nsc);
        double mw = tpRef.molecularWeight(kLoc);
        m_phaseMassChange[iph] += (nsc) * mw;
        double chg = tpRef.charge(kLoc);
        m_phaseChargeChange[iph] += nsc * chg;
    }

    for (iph = 0; iph < m_nPhases; iph++) {
        ThermoPhase& tpRef = m_kinBase->thermo(iph);
        m_phasePotentials[iph] =  tpRef.electricPotential();
        m_phaseDims[iph] = tpRef.nDim();
        m_phaseTypes[iph] = tpRef.eosType();
        if (m_phaseChargeChange[iph] != 0.0) {
            double tmp = fabs(m_phaseChargeChange[iph]);
            if (tmp >  m_ChargeTransferInRxn) {
                m_ChargeTransferInRxn = tmp;
            }
        }
    }

    if (m_ChargeTransferInRxn) {
        InterfaceKinetics* iK = dynamic_cast<InterfaceKinetics*>(kinPtr);
        if (iK) {
            m_beta = 0.0;
        } else {
            throw CanteraError("RxnMolChange", "unknown condition on charge");
        }
    }

}


RxnMolChangeLocal::~RxnMolChangeLocal()
{

}
//==================================================================================================================================

void processCurrentVsPotTable(RxnMolChange* rmc,
                              PhaseList* pl, int irxn,
                              TemperatureTable& TT,
                              Kinetics& kin,
                              DenseMatrix& kfwd_Table,
                              DenseMatrix& krev_Table,
                              DenseMatrix& deltaG_Table,
                              DenseMatrix& deltaH_Table,
                              DenseMatrix& deltaS_Table,
                              DenseMatrix& Afwd_Table,
                              DenseMatrix& EafwddivR_Table,
                              DenseMatrix& Arev_Table,
                              DenseMatrix& EarevdivR_Table,
                              DenseMatrix& kfwdPrime_Table,
                              DenseMatrix& krevPrime_Table)
{
    int iph, k, e;
    ThermoPhase* tp = 0;
    ThermoPhase* tpMetal = 0;
    int nTotalSpecies = kin.nTotalSpecies();
    int nRxns         = kin.nReactions();
    vector<double> Rfwd(nTotalSpecies, 0.0);
    vector<double> Rrev(nTotalSpecies, 0.0);
    vector<double> Rnet(nTotalSpecies, 0.0);
    vector<double> deltaG(nRxns, 0.0);
    vector<double> deltaSSG(nRxns, 0.0);

    InterfaceKinetics* iK = dynamic_cast<InterfaceKinetics*>(&kin);
    if (!iK) {
        throw CanteraError("RxnMolChange", "unknown condition on charge");
    }
    /*
     * First thing is to figure out how to calculate the voltage
     *    V = phi_metal - phi_soln
     * Look for the electron species. The phase where this occurs will be
     * called the metal, iMetal;
     */
    int iMetal = -1;
    int iSoln = 0;
    int kElectron = -1;
    int nPhases = iK->nPhases();
    for (iph = 0; iph < nPhases; iph++) {
        tp = &(iK->thermo(iph));
        std::string pName = tp->id();
        int nSpecies = tp->nSpecies();
        int nElements = tp->nElements();
        int eElectron = tp->elementIndex("E");
        if (eElectron >= 0) {
            for (k = 0; k < nSpecies; k++) {
                if (fabs(tp->nAtoms(k,eElectron) - 1) < 1.0E-10) {
                    int ifound = 1;
                    for (e = 0; e < nElements; e++) {
                        if (fabs(tp->nAtoms(k,e)) > 1.0E-10) {
                            if (e != eElectron) {
                                ifound = 0;
                            }
                        }
                    }
                    if (ifound == 1) {
                        if (iMetal == -1) {
                            iMetal = iph;
                            kElectron = iK->kineticsSpeciesIndex(k, iph);
                        } else {
                            int pi = iK->phaseIndex(pName);
                            if (rmc->m_phaseChargeChange[pi] != 0.0) {
                                iMetal = iph;
                                kElectron = iK->kineticsSpeciesIndex(k, iph);
                            }
                        }
                    }
                }
            }
        }
        if (iph != iMetal) {
            if (rmc->m_phaseChargeChange[iph] != 0) {
                iSoln = iph;
            }
        }
    }
    //
    //  Set State back to reference conditions
    //
    for (iph = 0; iph < nPhases; iph++) {
        ThermoPhase& gRef = kin.thermo(iph);
        size_t kstart = pl->globalSpeciesIndex(&gRef);
        gRef.setState_TPX(BG.Temperature, BG.Pressure, BG.XmolPLSpecVec +kstart);
    }
    tp = &(iK->thermo(iMetal));

    double phi0Metal = (iK->thermo(iMetal)).electricPotential();
    //double phi0Metal = rmc->m_phasePotentials[iMetal];
    double phi0Soln = (iK->thermo(iSoln)).electricPotential();
    //double phi0Soln = rmc->m_phasePotentials[iSoln];
    double V0 = phi0Metal - phi0Soln;
    //double V0 = rmc->m_phasePotentials[iMetal] - rmc->m_phasePotentials[iSoln];
    tpMetal = &(iK->thermo(iMetal));
    ThermoPhase* tpSoln =  &(iK->thermo(iSoln));
    int iSurf = iK->reactionPhaseIndex();
    ThermoPhase* tps = &(iK->thermo(iSurf));
    SurfPhase* tpSurface = dynamic_cast<SurfPhase*>(tps);

    int nSurfSpecies = tpSurface->nSpecies();
    vector<double> thetaSurf(nSurfSpecies, 0.0);
    tpSurface->getCoverages(DATA_PTR(thetaSurf));
    double nStoichElectrons = - rmc->m_phaseChargeChange[iMetal];


    int pwidth = 100;
    int twidth = 93;
    dnt(1);
    print_char('-', pwidth);
    printf("\n");
    std::string rs = kin.reactionString(irxn);
    dnt(1);
    printf("Voltage Dependence for Rxn: %s", rs.c_str());
    printf("\n\n");
    dnt(2);
    printf("ThermoPhase index pertaining to the solution = %d\n", iSoln);
    dnt(2);
    printf("ThermoPhase index pertaining to the metal = %d\n", iMetal);
    dnt(2);
    printf("Kinetic Species Index pertaining to the electrons = %d\n", kElectron);
    dnt(2);
    printf("Input Electric Potential of the Metal = %g Volts\n", phi0Metal);
    dnt(2);
    printf("Input Electric Potential of the Soln  = %g Volts\n", phi0Soln);
    dnt(2);
    printf("Input Voltage = phi_Metal - phi_Soln  = %g Volts\n", V0);
    dnt(2);
    printf("nStoichElectrons (product electrons - reactant electrons) = %g \n", nStoichElectrons);
    /*
     * the 1.0E-4 is to change from m-2 to cm-2.
     */
    iK->getDeltaGibbs(DATA_PTR(deltaG));
    iK->getDeltaSSGibbs(DATA_PTR(deltaSSG));
    double deltaGrxn = deltaG[irxn];
    double deltaSSGrxn = deltaSSG[irxn];
    double Erxn =  deltaGrxn/Faraday/nStoichElectrons;
    double ESSrxn =  deltaSSGrxn/Faraday/nStoichElectrons;


    dnt(2);
    printf("Delta G (rxn)   = %13.5g J/kmol   --->   ", deltaGrxn);
    dnt(2);
    printf("E    = %g Volts", Erxn);

    printf("\n");
    dnt(2);
    printf("Delta G_SS (rxn)= %13.5g J/kmol   --->   ", deltaSSGrxn);
    dnt(2);
    printf("E_SS = %g Volts", ESSrxn);

    printf("\n");
    dnt(2);
    printf("Bath Gas conditions:\n");
    dnt(3);
    printf("T = %g K\n", tpSoln->temperature());
    dnt(3);
    printf("P = %g Pascal\n", tpSoln->pressure());
    dnt(3);
    printf("     SurfSpecName         Theta_k_Bath\n");
    for (k = 0; k < nSurfSpecies; k++) {
        std::string sName = tpSurface->speciesName(k);
        dnt(3);
        printf("%16s %12.5g\n", sName.c_str(), thetaSurf[k]);
    }
    printf("\n");

    double startingVoltage = -0.6;
    if (Erxn < startingVoltage) {
        startingVoltage = Erxn - 0.2;
        startingVoltage = (int(startingVoltage*10.0))/10.0;
    }

    /*
     *  We use this to determine the current. The current will be equal
     *  to the net rate of electron production in the metal. A positive
     *  current will have positive charge flow from the metal into
     *  the solution
     */
    double nF = - Faraday * rmc->m_phaseChargeChange[iMetal] * 1.0E-4;


    dnt(2);
    print_char('-', twidth);
    printf("\n");
    dnt(2);
    cout << "|  Voltage  |    Rfwd          Rrev         Rnet   "
         << "|        Ifwd         Irev         Inet  ";
    cout << "|" << endl;

    dnt(2);
    cout << "|   Volts   |             (kmol/m**2 s)            "
         << "|               (amps/cm2)               "
         << "|" << endl;
    dnt(2);
    print_char('-', twidth);
    printf("\n");
    double Voltage, phiMetal;

#ifdef DEBUG_JCH_ETABLE
    FILE* EFP = fopen("polarizationData.dat", "w");
    fprintf(EFP,"TITLE = \"EMF  = %g Volts\"\n", Erxn);
    fprintf(EFP,"VARIABLES =  \"Voltage\"  \"Rfwd\"  \"Rrev\"  \"Rnet\"  \"Ifwd\"   \"Irev\"  \"Inet\"\n");
#endif

    if (IOO.VVincEeq) {
        VV_ptr->AddEeq(Erxn);
    }
    if (IOO.VVincEzero) {
        VV_ptr->AddEzero(ESSrxn);
    }

#ifdef CENTERED_VOLTAGE_TABLE
    //Make a new Voltage Table around Erxn
    if (IOO.VVincEeq) {
        double dVtab[9] = { 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0 };
        for (int i = 0; i < 9; i++) {
            VV_ptr->addPoint(Erxn + dVtab[i]);
            VV_ptr->addPoint(Erxn - dVtab[i]);
        }
    }
#endif // CENTERED_VOLTAGE_TABLE

    int npts = VV_ptr->size();

    for (int iV = 0; iV < npts; iV++) {
        Voltage = (*VV_ptr)[iV];
        phiMetal = Voltage + (iK->thermo(iSoln)).electricPotential();
        //phiMetal = Voltage + rmc->m_phasePotentials[iSoln];
        tpMetal->setElectricPotential(phiMetal);

        iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
        iK->getRevRatesOfProgress(DATA_PTR(Rrev));
        iK->getNetRatesOfProgress(DATA_PTR(Rnet));
        dnt(2);
        printf("!%10.4g |%11.3g %12.3g %12.3g |%12.3g %12.3g  %12.3g |\n",
               Voltage, Rfwd[irxn], Rrev[irxn], Rnet[irxn], nF*Rfwd[irxn], - nF*Rrev[irxn], nF*Rnet[irxn]);
#ifdef DEBUG_JCH_ETABLE
        fprintf(EFP, "%10.4g  %11.5g %11.5g %11.5g %11.5g %11.5g %11.5g\n",
                Voltage, Rfwd[irxn], Rrev[irxn], Rnet[irxn], nF*Rfwd[irxn], - nF*Rrev[irxn], nF*Rnet[irxn]);
#endif

    }
#ifdef DEBUG_JCH_ETABLE
    fclose(EFP);
#endif

    dnt(2);
    print_char('-', twidth);
    printf("\n");

    dnt(2);
    printf(" TABLE that uses the Surface Solver to Solve for Surface Site Concentrations\n");


    int nSS = std::min(nSurfSpecies, 6);
    twidth = 54 + nSS * 11;
    dnt(2);
    print_char('-', twidth);
    printf("\n");
    dnt(2);
    cout << "|  Voltage  "
         << "|       Ifwd         Irev        Inet  ";
    cout << "| ";

    for (k = 0; k < nSS; k++) {
        std::string sname =  tpSurface->speciesName(k);
        printf("%10s ", sname.c_str());
    }
    cout <<"|" << endl;
    dnt(2);
    cout << "|   Volts   "
         << "|               (amps/cm2)             "
         << "| " ;
    for (k = 0; k < nSS; k++) {
        printf("           ");
    }
    cout <<"|" << endl;
    dnt(2);
    print_char('-', twidth);
    printf("\n");


    for (int iV = 0; iV < npts; iV++) {
        Voltage = (*VV_ptr)[iV];
        phiMetal = Voltage + (iK->thermo(iSoln)).electricPotential();
        //phiMetal = Voltage + rmc->m_phasePotentials[iSoln];
        tpMetal->setElectricPotential(phiMetal);
        iK->advanceCoverages(100.0);
        iK->solvePseudoSteadyStateProblem(SFLUX_RESIDUAL);
        iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
        iK->getRevRatesOfProgress(DATA_PTR(Rrev));
        iK->getNetRatesOfProgress(DATA_PTR(Rnet));

        dnt(2);
        printf("!%10.4g |%11.3g %12.3g %12.3g |",
               Voltage,
               nF*Rfwd[irxn], - nF*Rrev[irxn], nF*Rnet[irxn]);
        tpSurface->getCoverages(DATA_PTR(thetaSurf));
        for (k = 0; k < nSS; k++) {
            printf("%10.4g ", thetaSurf[k]);
        }
        printf(" |\n");
    }
    dnt(2);
    print_char('-', twidth);
    printf("\n");

    dnt(1);
    print_char('-', pwidth);
    printf("\n");
}


RxnTempTableStuff::RxnTempTableStuff(int irxn, int irxnGE) :
    m_irxn(irxn),
    m_irxnGlobalExtra(irxnGE)
{
}

RxnTempTableStuff::~RxnTempTableStuff()
{
}

/*
 *
 */
void
getGERKineticsTables(TemperatureTable& TT, PhaseList* pl,
                     ZZCantera::Kinetics& kin,
                     ExtraGlobalRxn& egr,
                     RxnTempTableStuff& rts)
{


    int i, iph;
    int nReactions = kin.nReactions();
    int nPhases = kin.nPhases();
    int nTotSpecies = kin.nTotalSpecies();
    double* Rarray = new double [nReactions];
    double* R2array = new double [nReactions];
    double* Sarray = new double [nTotSpecies];
    double deltaT = 0.001;
    double Afac, EadivR;

    vector<double>& kfwd_Table      = rts.kfwd_Table;
    vector<double>& krev_Table      = rts.krev_Table;
    vector<double>& deltaG_Table    = rts.deltaG_Table;
    vector<double>& deltaH_Table    = rts.deltaH_Table;
    vector<double>& deltaS_Table    = rts.deltaS_Table;
    vector<double>& deltaGss_Table    = rts.deltaGss_Table;
    vector<double>& deltaHss_Table    = rts.deltaHss_Table;
    vector<double>& deltaSss_Table    = rts.deltaSss_Table;
    vector<double>& Afwd_Table      = rts.Afwd_Table;
    vector<double>& EafwddivR_Table = rts.EafwddivR_Table;
    vector<double>& Arev_Table      = rts.Arev_Table;
    vector<double>& EarevdivR_Table = rts.EarevdivR_Table;
    vector<double>& kfwdPrime_Table = rts.kfwdPrime_Table;
    vector<double>& krevPrime_Table = rts.krevPrime_Table;
    vector<double>& NetROP_Table    = rts.NetROP_Table;
    vector<double>& FwdROP_Table    = rts.FwdROP_Table;
    vector<double>& RevROP_Table    = rts.RevROP_Table;
    vector<double>& Anet_Table      = rts.Anet_Table;
    vector<double>& EanetdivR_Table = rts.EanetdivR_Table;

    int m = TT.size();
    kfwd_Table.resize(m,0.0);
    krev_Table.resize(m,0.0);
    deltaG_Table.resize(m,0.0);
    deltaH_Table.resize(m,0.0);
    deltaS_Table.resize(m,0.0);
    deltaGss_Table.resize(m,0.0);
    deltaHss_Table.resize(m,0.0);
    deltaSss_Table.resize(m,0.0);
    Afwd_Table.resize(m,0.0);
    EafwddivR_Table.resize(m,0.0);
    Arev_Table.resize(m,0.0);
    EarevdivR_Table.resize(m,0.0);
    kfwdPrime_Table.resize(m,0.0);
    krevPrime_Table.resize(m,0.0);
    NetROP_Table.resize(m,0.0);
    FwdROP_Table.resize(m,0.0);
    RevROP_Table.resize(m,0.0);
    Anet_Table.resize(m,0.0);
    EanetdivR_Table.resize(m,0.0);

    InterfaceKinetics* iKin = dynamic_cast<InterfaceKinetics*>(&kin);

    setAllBathSpeciesConditions(pl);
    double T, Tdelta;
    for (i = 0; i < TT.size(); i++) {
        T = TT[i];
        for (iph = 0; iph < nPhases; iph++) {
            ThermoPhase& gRef = kin.thermo(iph);
            int kstart = pl->globalSpeciesIndex(&gRef);
            gRef.setState_TPX(T, BG.Pressure, BG.XmolPLSpecVec +kstart);
        }

        if (iKin) {
            iKin->advanceCoverages(100.0);
        }

        kin.getFwdRatesOfProgress(Rarray);
        kin.getRevRatesOfProgress(R2array);
        FwdROP_Table[i] = egr.FwdROPValue(Rarray, R2array);
        RevROP_Table[i] = egr.RevROPValue(Rarray, R2array);

        kin.getNetRatesOfProgress(Rarray);
        NetROP_Table[i] = egr.ROPValue(Rarray);

        kin.getDeltaSSGibbs(Rarray);
        deltaGss_Table[i] = egr.deltaRxnVecValue(Rarray);

        kin.getDeltaSSEnthalpy(Rarray);
        deltaHss_Table[i] = egr.deltaRxnVecValue(Rarray);

        kin.getDeltaSSEntropy(Rarray);
        deltaSss_Table[i] = egr.deltaRxnVecValue(Rarray);

        kin.getDeltaGibbs(Rarray);
        deltaG_Table[i] = egr.deltaRxnVecValue(Rarray);

        kin.getDeltaEnthalpy(Rarray);
        deltaH_Table[i] = egr.deltaRxnVecValue(Rarray);

        kin.getDeltaEntropy(Rarray);
        deltaS_Table[i] = egr.deltaRxnVecValue(Rarray);

        Tdelta = T + deltaT;
        for (iph = 0; iph < nPhases; iph++) {
            ThermoPhase& gRef = kin.thermo(iph);
            int kstart = pl->globalSpeciesIndex(&gRef);
            gRef.setState_TPX(Tdelta, BG.Pressure, BG.XmolPLSpecVec +kstart);
        }
        if (iKin) {
            iKin->advanceCoverages(100.0);
        }
        kin.getFwdRatesOfProgress(Rarray);
        kin.getRevRatesOfProgress(R2array);

        double FwdROP_delta = egr.FwdROPValue(Rarray, R2array);
        double RevROP_delta = egr.RevROPValue(Rarray, R2array);

        kin.getNetRatesOfProgress(Rarray);
        double NetROP_delta = egr.ROPValue(Rarray);

        calcArrhParam(FwdROP_Table[i], T, FwdROP_delta, Tdelta, Afac, EadivR);
        Afwd_Table[i] = Afac;
        EafwddivR_Table[i] = EadivR;

        calcArrhParam(RevROP_Table[i], T, RevROP_delta,  Tdelta, Afac, EadivR);
        Arev_Table[i] = Afac;
        EarevdivR_Table[i] = EadivR;

        if (NetROP_Table[i] >= 0.0) {
            calcArrhParam(NetROP_Table[i], T, NetROP_delta, Tdelta, Afac, EadivR);
            Anet_Table[i] = Afac;
        } else {
            calcArrhParam(-NetROP_Table[i], T, -NetROP_delta, Tdelta, Afac, EadivR);
            Anet_Table[i] = -Afac;
        }
        EanetdivR_Table[i] = EadivR;


    }


    setAllBathSpeciesConditions(pl);




    delete [] Sarray;
    delete [] Rarray;
    delete [] R2array;
}

//===============================================================================================
struct StoichVectors {
    StoichVectors() :
        kKinV(0),
        sCoeffV(0)
    {
    }
    size_t num() const
    {
        return kKinV.size();
    }
    void add(size_t k, double p, bool force=false)
    {
        if (p != 0.0 || force) {
            kKinV.push_back(k);
            sCoeffV.push_back(p);
        }
    }
    double coeff(size_t kFind) const
    {
        for (size_t k = 0; k < kKinV.size(); ++k) {
            if (kKinV[k] == kFind) {
                return sCoeffV[k];
            }
        }
        return 0.0;
    }
    size_t index(size_t kFind) const
    {
        for (size_t k = 0; k < kKinV.size(); ++k) {
            if (kKinV[k] == kFind) {
                return k;
            }
        }
        return npos;
    }
    std::vector<size_t> kKinV;
    std::vector<double> sCoeffV;
};
//==================================================================================================================================
void printAffinityHeader(ZZCantera::RxnMolChange* rmc, PhaseList* pl, int iRxn, TemperatureTable& TT,
                         Kinetics& kin,  DenseMatrix& kfwd_Table,   DenseMatrix& krev_Table,
                         DenseMatrix& deltaG_Table, DenseMatrix& deltaH_Table,  DenseMatrix& deltaS_Table,
                         DenseMatrix& Afwd_Table, DenseMatrix& EafwddivR_Table,  DenseMatrix& Arev_Table,
                         DenseMatrix& EarevdivR_Table, DenseMatrix& kfwdPrime_Table, DenseMatrix& krevPrime_Table,
                         double* unitskfwd, double* unitskrev)
{
    size_t nTotSpecies = kin.nTotalSpecies();

    double affinity;
    affinityRxnData aj;
    bool hasElectricalTerm;
    double betaf;
    double kc;
    //double phiSoln;
    double perturb;
    //double deltaM[10];
    ZZCantera::InterfaceKinetics* iKin = dynamic_cast<InterfaceKinetics*>(&kin);
    if (!iKin) {
        throw CanteraError("processAffinityTable", "failure");
    }


    double* actConc = new double [nTotSpecies];
    iKin->getActivityConcentrations(actConc);
    bool ok = iKin->getAffinityRxnFormulation(iRxn, affinity, kc, perturb, aj, hasElectricalTerm, betaf);
    if (!ok) {
        throw CanteraError("printAffinityHeader()",  "error calling getAffinityRxnFormulation");
    }

    const double affinityPowerFwd = aj.affinityPowerFwd;
    const double equilibPowerFwd = aj.equilibPowerFwd;
    const size_t numForwardRxnSpecies = aj.numForwardRxnSpecies;
    const size_t numRevFwdRxnSpecies = aj.numRevFwdRxnSpecies;
    const std::vector<size_t>& affinSpec_FRC = aj.affinSpec_FRC;
    const std::vector<double>&  affinOrder_FRC = aj.affinOrder_FRC;
    const std::vector<size_t>& affinSpec_RFRC = aj.affinSpec_RFRC;
    const std::vector<double>&  affinOrder_RFRC = aj.affinOrder_RFRC;

    dnt(1);
    printf("Forward Affinity Power = %g \n", affinityPowerFwd);
    dnt(1);
    printf("Equlibrium Power Value = %g \n", equilibPowerFwd);
    dnt(1);
    printf("Number of Forward Rxn Species  = %d \n", (int) numForwardRxnSpecies);


    size_t jRevIndex = aj.index_RevRateConstant;
    if (jRevIndex != npos) {
        dnt(1);
        printf("Has a different reverse reaction rate constant, index  = %d \n", (int) jRevIndex);
    }

    StoichVectors rstoichVec;
    StoichVectors pstoichVec;
    for (size_t k = 0; k < nTotSpecies; ++k) {
        double val = kin.reactantStoichCoeff(k, iRxn);
        if (val != 0.0) {
            rstoichVec.add(k, val);
        }
        val = kin.productStoichCoeff(k, iRxn);
        if (val != 0.0) {
            pstoichVec.add(k, val);
        }
    }
    dnt(2);
    printf("        Rate form in forward direction: ");
    iKin->write_one(iRxn, cout, false);

    dnt(2);
    printf("        Reaction Orders in the Forwards Direction Far from Equilibrium:\n");
    int pwidth = 111;
    dnt(2);
    print_char('-', pwidth);
    printf("\n");
    dnt(2);
    printf("|      Species  Name | FStoichCoeffient | ForwardRxnPower|StandardConc|  moleFrac  |  ActCoeff  | RateContrib|\n");
    double uA[6];

    std::vector<bool> kStoichHandled(rstoichVec.num(), false);
    for (size_t kk = 0; kk <  numForwardRxnSpecies; ++kk) {
        size_t kKin = affinSpec_FRC[kk];
        std::string kname = kin.kineticsSpeciesName(kKin);
        double stoichV = kin.reactantStoichCoeff(kKin, iRxn);
        double forwP = affinOrder_FRC[kk];

        size_t kThermo;
        size_t n = kin.speciesPhaseSpeciesIndex(kKin, kThermo);
        double sConc = kin.thermo(n).standardConcentration(kThermo);
        ThermoPhase& tp = kin.thermo(n);

        tp.getUnitsStandardConc(uA, kThermo);
        for (size_t jj = 0; jj < 6; ++jj) {
            unitskfwd[jj] += forwP * uA[jj];
        }

        double xMol = tp.moleFraction(kThermo);
        vector<double> ac(tp.nSpecies());
        tp.getActivityCoefficients(DATA_PTR(ac));
        double a = ac[kThermo] * xMol * sConc;
        double value = pow(a, forwP);

        dnt(2);
        printf("| %17.17s  |     %9.4f    |     %9.4f   |%11.4E |%11.4E |%11.4E |%11.4E |\n",
               kname.c_str(), stoichV, forwP, sConc, xMol, ac[kThermo], value);
        size_t index = rstoichVec.index(kKin);
        if (index != npos) {
            kStoichHandled[index] = true;
        }
    }
    for (size_t  kk = 0; kk < rstoichVec.num(); ++kk) {
        if (! kStoichHandled[kk]) {
            size_t kKin = rstoichVec.kKinV[kk];
            std::string kname = kin.kineticsSpeciesName(kKin);
            double stoichV = rstoichVec.sCoeffV[kk];
            double forwP = 0.0;
            size_t kThermo;
            size_t n = kin.speciesPhaseSpeciesIndex(kKin, kThermo);
            ThermoPhase& tp = kin.thermo(n);
            double sConc = tp.standardConcentration(kThermo);
            double xMol = tp.moleFraction(kThermo);
            vector<double> ac(tp.nSpecies());
            tp.getActivityCoefficients(DATA_PTR(ac));
            dnt(1);
            printf("| %17.17s  |     %9.4f    |     %9.4f   |%11.4E |%11.4E |%11.4E |     NA     |\n",
                   kname.c_str() , stoichV, forwP,  sConc, xMol, ac[kThermo]);
        }
    }
    dnt(2);
    print_char('-', pwidth);
    printf("\n");
    for (size_t  kk = 0; kk < rstoichVec.num(); ++kk) {
        kStoichHandled[kk] = false;
    }
    //
    // Here we will assume that 1-exp(-A/RT) < 0 and that the reverse direction is larger
    //
    dnt(2);
    printf("        Rate form in reverse direction: ");
    iKin->write_one(iRxn, cout, false);
    dnt(2);
    printf("        Reaction Orders in the Reverse Direction Far from Equilibrium:\n");
    dnt(2);
    print_char('-', pwidth);
    printf("\n");
    dnt(2);
    printf("|      Species  Name |NetStoichCoeffient| ForwardRxnPower | EffReversePower|StandardConc|  moleFrac  |  ActCoeff  |RevRateContr|\n");
    double multiplier = 1.0 * aj.affinityPowerFwd *  aj.equilibPowerFwd;
    std::vector<bool> kpStoichHandled(pstoichVec.num(), false);
    for (size_t kk = 0; kk < numRevFwdRxnSpecies; ++kk) {
        size_t kKin = affinSpec_RFRC[kk];
        std::string kname = kin.kineticsSpeciesName(kKin);
        //double stoichV = kin.productStoichCoeff(kKin, iRxn) - kin.reactantStoichCoeff(kKin, iRxn);
        double forwP = affinOrder_RFRC[kk];
        double netStoich = kin.productStoichCoeff(kKin, iRxn) - kin.reactantStoichCoeff(kKin, iRxn);
        double effReversePower = forwP + netStoich * multiplier;

        size_t kThermo;
        size_t n = kin.speciesPhaseSpeciesIndex(kKin, kThermo);
        ThermoPhase& tp = kin.thermo(n);
        double sConc = tp.standardConcentration(kThermo);
        double xMol = tp.moleFraction(kThermo);
        vector<double> ac(tp.nSpecies());
        tp.getActivityCoefficients(DATA_PTR(ac));
        double act = ac[kThermo] * xMol * sConc;
        double actConc = ac[kThermo] * xMol * sConc;
        double value = pow(act, netStoich * multiplier) * pow(actConc, forwP);

        tp.getUnitsStandardConc(uA, kThermo);
        for (size_t jj = 0; jj < 6; ++jj) {
            unitskrev[jj] += forwP * uA[jj];
        }

        dnt(2);
        printf("| %17.17s  |     %9.4f    |     %9.4f   |    %9.4f   |%11.4E |%11.4E |%11.4E ",
               kname.c_str(), netStoich, forwP, effReversePower, sConc, xMol, ac[kThermo]);
        if (effReversePower == 0.0 && value == 1.0) {
            printf("|     NA     |\n");
        } else {
            printf("|%11.4E |\n", value);
        }
        size_t index = pstoichVec.index(kKin);
        if (index != npos) {
            kpStoichHandled[index] = true;
        }
        index = rstoichVec.index(kKin);
        if (index != npos) {
            kStoichHandled[index] = true;
        }
    }
    for (size_t  kk = 0; kk < pstoichVec.num(); ++kk) {
        if (! kpStoichHandled[kk]) {
            size_t kKin = pstoichVec.kKinV[kk];
            std::string kname = kin.kineticsSpeciesName(kKin);
            // double stoichV = rstoichVec.sCoeffV[kk];
            double forwP = 0.0;
            double netStoich = kin.productStoichCoeff(kKin, iRxn) - kin.reactantStoichCoeff(kKin, iRxn);
            double effReversePower = forwP + netStoich * multiplier;
            size_t kThermo;
            size_t n = kin.speciesPhaseSpeciesIndex(kKin, kThermo);
            ThermoPhase& tp = kin.thermo(n);
            double sConc = tp.standardConcentration(kThermo);
            double xMol = tp.moleFraction(kThermo);
            vector<double> ac(tp.nSpecies());
            tp.getActivityCoefficients(DATA_PTR(ac));

            double act = ac[kThermo] * xMol;
            double value = pow(act, effReversePower);
            dnt(2);
            printf("| %17.17s  |     %9.4f    |     %9.4f   |    %9.4f   |%11.4E |%11.4E |%11.4E ",
                   kname.c_str(), netStoich, forwP, effReversePower, sConc, xMol, ac[kThermo]);
            if (effReversePower == 0.0 && value == 1.0) {
                printf("|   NA     |\n");
            } else {
                printf("|%11.4E |\n", value);
            }
        }
    }
    for (size_t  kk = 0; kk < rstoichVec.num(); ++kk) {
        if (! kStoichHandled[kk]) {
            size_t kKin = rstoichVec.kKinV[kk];
            std::string kname = kin.kineticsSpeciesName(kKin);
            // double stoichV = rstoichVec.sCoeffV[kk];
            double forwP = 0.0;
            double netStoich = kin.productStoichCoeff(kKin, iRxn) - kin.reactantStoichCoeff(kKin, iRxn);
            double effReversePower = forwP + netStoich * multiplier;
            size_t kThermo;
            size_t n = kin.speciesPhaseSpeciesIndex(kKin, kThermo);
            ThermoPhase& tp = kin.thermo(n);
            double sConc = tp.standardConcentration(kThermo);
            double xMol = tp.moleFraction(kThermo);
            vector<double> ac(tp.nSpecies());
            tp.getActivityCoefficients(DATA_PTR(ac));

            double act = ac[kThermo] * xMol;
            double value = pow(act, effReversePower);
            dnt(2);
            printf("| %17.17s  |     %9.4f    |     %9.4f   |    %9.4f   |%11.4E |%11.4E |%11.4E ",
                   kname.c_str(), netStoich, forwP, effReversePower, sConc, xMol, ac[kThermo]);
            if (effReversePower == 0.0) {
                printf("|   NA     |\n");
            } else {
                printf("|%11.4E |\n", value);
            }
        }
    }


    dnt(2);
    print_char('-', pwidth);
    printf("\n");



    delete [] actConc;

}
//==================================================================================================================================
void processAffinityTable(ZZCantera::RxnMolChange* rmc, PhaseList* pl, int irxn, TemperatureTable& TT,
                          Kinetics& kin,  DenseMatrix& kfwd_Table,   DenseMatrix& krev_Table,
                          DenseMatrix& deltaG_Table, DenseMatrix& deltaH_Table,  DenseMatrix& deltaS_Table,
                          DenseMatrix& Afwd_Table, DenseMatrix& EafwddivR_Table,  DenseMatrix& Arev_Table,
                          DenseMatrix& EarevdivR_Table, DenseMatrix& kfwdPrime_Table,   DenseMatrix& krevPrime_Table)
{
    int nTotSpecies = kin.nTotalSpecies();

    double affinity;
    affinityRxnData aj;
    bool hasElectricalTerm;
    double betaf;
    double kc;
    //double phiSoln;
    double perturb;
    //double deltaM[10];
    ZZCantera::InterfaceKinetics* iKin = dynamic_cast<InterfaceKinetics*>(&kin);
    if (!iKin) {
        throw CanteraError("processAffinityTable", "failure");
    }


    double* actConc = new double [nTotSpecies];
    iKin->getActivityConcentrations(actConc);
    bool ok = iKin->getAffinityRxnFormulation(irxn, affinity, kc, perturb, aj, hasElectricalTerm, betaf);
    if (!ok) {
        throw CanteraError(" processAffinityTable()",
                           " >getAffinityRxnFormulation() failed ");
    }
    /*
     const double affinityPowerFwd = aj.affinityPowerFwd;
     const double equilibPowerFwd = aj.equilibPowerFwd;
     const size_t numForwardRxnSpecies = aj.numForwardRxnSpecies;
     const std::vector<size_t>& affinSpec_FRC = aj.affinSpec_FRC;
     const std::vector<double>&  affinOrder_FRC = aj.affinOrder_FRC;
    */

    delete [] actConc;
}
//==================================================================================================================================
/*
 *  This routine will print out a table of information a single
 *  reaction, j.
 *
 *
 */
void printGERKineticsTable(PhaseList* pl, int iGER,
                           TemperatureTable& TT,
                           ZZCantera::Kinetics& kin,
                           ExtraGlobalRxn& egr,
                           RxnMolChange* rmc,
                           RxnTempTableStuff& rts)
{

    // vector<double> & kfwd_Table      = rts.kfwd_Table;
    // vector<double> & krev_Table      = rts.krev_Table;
    // vector<double> & deltaG_Table    = rts.deltaG_Table;
    // vector<double> & deltaH_Table    = rts.deltaH_Table;
    // vector<double> & deltaS_Table    = rts.deltaS_Table;
    vector<double>& deltaGss_Table  = rts.deltaGss_Table;
    vector<double>& deltaHss_Table  = rts.deltaHss_Table;
    vector<double>& deltaSss_Table  = rts.deltaSss_Table;
    vector<double>& Afwd_Table      = rts.Afwd_Table;
    vector<double>& EafwddivR_Table = rts.EafwddivR_Table;
    vector<double>& Arev_Table      = rts.Arev_Table;
    vector<double>& EarevdivR_Table = rts.EarevdivR_Table;
    // vector<double> & kfwdPrime_Table = rts.kfwdPrime_Table;
    // vector<double> & krevPrime_Table = rts.krevPrime_Table;
    vector<double>& NetROP_Table    = rts.NetROP_Table;
    vector<double>& FwdROP_Table    = rts.FwdROP_Table;
    vector<double>& RevROP_Table    = rts.RevROP_Table;

    vector<double>& Anet_Table      = rts.Anet_Table;
    vector<double>& EanetdivR_Table = rts.EanetdivR_Table;

    /*
     *  Get the species data object from the Mixture object
     *  this is defined in the constituents.h file, and is
     *  inherited by Mixture through BaseMix
     */

    int tableWidth = 138;

    /*
     * Conversion factor between Joules/kmol, cantera's units,
     * and kcal / gmol.
     */
    double conv = 1.0 / (4.184E6);
    double convS = conv * 1.0E3;
    double Rkcal = 1.98721E-3;

    // Possible change the conv units to kJ / gmol
    if (UIO.unitDef == UNITS_KJOULE) {
        conv = 1.0 / 1.0E6;
        convS = conv * 1.0E3;
        Rkcal = 8.314472E-3;
    }

    ThermoPhase* gThermo = &kin.thermo(0);
    /*
     *  Dump out all of the information about the reaction
     */
    cout << endl;
    print_char('=', tableWidth);
    cout << endl;
    cout << endl;
    cout << "INFORMATION TABLE FOR Global Extra Reaction Reaction \""  << iGER;
    cout << "\" IN PHASE \"";
    cout << gThermo->id();
    cout << "\"" << endl;
    cout << endl;

    /*
     * Reaction String
     */
    std::string rs = egr.reactionString();
    dnt(1);
    cout << rs << endl;
    cout << endl;

    /*
     * Print out reactants
     */
    const std::vector<size_t>& reactants = egr.reactants();
    int nReac = reactants.size();
    dnt(1);
    if (nReac == 1) {
        cout << "There is one reactant: ";
    } else {
        cout << "There are " << nReac << " reactants: ";
    }
    for (int k = 0; k < nReac; k++) {
        int kindex = reactants[k];
        cout << kin.kineticsSpeciesName(kindex) << " (" << egr.m_ReactantStoich[k] << ") ";
    }
    cout << endl;


    /*
     * Print out products
     */
    const std::vector<size_t>& products = egr.products();
    int nProd = products.size();
    dnt(1);
    if (nProd == 1) {
        cout << "There is one product: ";
    } else {
        cout << "There are " << nProd << " products: ";
    }
    for (int k = 0; k < nProd; k++) {
        int kindex = products[k];
        cout << kin.kineticsSpeciesName(kindex) << " (" << egr.m_ProductStoich[k] << ") ";
    }
    cout << endl;

    /*
     *  reaction type
     */
    // double pCurrent = gThermo->pressure();

    dnt(1);
    cout << "Reaction Type = ";
    cout << "Global Extra Rxn - not in mechanism" << endl;

    /*
     * Reaction Reversible?
     */
    bool rReversible = egr.isReversible();
    dnt(1);
    if (rReversible) {
        cout << "Reaction is reversible" << endl;
    } else {
        cout << "Reaction is not reversible" << endl;
    }
    cout << endl;

    /*
     * Figure out the units for the forward reaction rate constant
     * First start out with the units for the ROP
     */
    int nphase = kin.nPhases();
    int ndim = 3;
    for (int iph = 0; iph < nphase; iph++) {
        ThermoPhase& tpRef = kin.thermo(iph);
        int idim = tpRef.nDim();
        if (idim < ndim) {
            ndim = idim;
        }
    }
    double unitsROP[6] = { 1.0, double(-ndim), 0.0, 0.0, 0.0, -1.0 };
    double unitskfwd[6];
    double unitskrev[6] ;
    double unitsSpecies[6];
    for (int i = 0; i < 6; i++) {
        unitskfwd[i] = unitsROP[i];
        unitskrev[i] = unitsROP[i];
    }
    for (int k = 0; k < nReac; k++) {
        int kindex = reactants[k];
        int iphase = kin.speciesPhaseIndex(kindex);
        ThermoPhase& tpRef = kin.thermo(iphase);
        tpRef.getUnitsStandardConc(unitsSpecies, kindex);
        for (int i = 0; i < 6; i++) {
            unitskfwd[i] -= unitsSpecies[i];
        }
    }
    for (int k = 0; k < nProd; k++) {
        int kindex = products[k];
        int iphase = kin.speciesPhaseIndex(kindex);
        ThermoPhase& tpRef = kin.thermo(iphase);
        tpRef.getUnitsStandardConc(unitsSpecies, kindex);
        for (int i = 0; i < 6; i++) {
            unitskrev[i] -= unitsSpecies[i];
        }
    }

    /****************** phases change table ****************************/
    int pwidth = 52;
    if (rmc->m_ChargeTransferInRxn != 0.0) {
        dnt(1);
        printf("This reaction is a charge transfer reaction, n = %g\n",
               rmc->m_ChargeTransferInRxn);
        pwidth = 86;
    }

    dnt(2);
    print_char('-', pwidth);
    printf("\n");

    dnt(2);
    printf("|              Phase | Dim Delta_Moles  Delta_Mass |");
    if (rmc->m_ChargeTransferInRxn != 0.0) {
        printf("Electron_Transfer Base_Potential |");
    }
    printf("\n");
    for (size_t iph = 0; iph < rmc->m_nPhases; iph++) {
        ThermoPhase& tpRef = kin.thermo(iph);
        std::string sname = tpRef.name();
        dnt(2);
        printf("|%20s|", sname.c_str());
        printf("%4d",    rmc->m_phaseDims[iph]);
        printf("%12.5g", rmc->m_phaseMoleChange[iph]);
        printf("%12.5g", rmc->m_phaseMassChange[iph]);
        if (rmc->m_ChargeTransferInRxn != 0.0) {
            printf(" |");
            printf("  %12.5g", rmc->m_phaseChargeChange[iph]);
            ThermoPhase& tpp = (kin.thermo(iph));
            double volts = tpp.electricPotential();
            printf("      %12.5g", volts);
        }
        printf(" |\n");
    }

    dnt(2);
    print_char('-', pwidth);
    printf("\n\n");


    /*******************************************************************/

    std::string units_kfwd = formUnitsString(unitskfwd);
    std::string units_krev = formUnitsString(unitskrev);
    dnt(1);
    std::string units_ROP =  "kmol / m^" + int2str(ndim) + " s";
    cout << "units_ROP = Units for Global Reaction Rate of progress = "
         << units_ROP << endl;


    /*
     * Dump out the kinetics table as a function of the
     * temperature
     */
    cout << endl;
    cout << "|";
    print_char('-', tableWidth);
    cout << "|" << endl;
    cout << "|        |                              "
         << " |                                 ";
    if (rReversible) {
        cout << "|                               ";
    } else {
        cout << "| -- REVERSE RXN NOT IN MECH -- ";
    }
    cout << "|                               ";
    cout << "|" << endl;
    cout << "|  Temp  |  FwdROP     Afwd     Ea_fwd  "
         << " |  DeltaG0    DeltaH0    DeltaS0  "
         << "|  RevROP    A_rev     Ea_rev   "
         << "|  NetROP       Anet     Ea_net ";
    cout << "|" << endl;
    if (UIO.unitDef == UNITS_KJOULE) {
        cout << "|   (K)  |    (units_ROP)       (kJ/gmol)"
             << "|  (kJ/gmol) (kJ/gmol)  (J/gmolK) "
             << "|  (units_ROP)        (kJ/gmol) |"
             << "   (units_ROP)       (kJ/gmol) |";
        cout << endl;
    } else {
        cout << "|   (K)  |    (units_ROP)     (kcal/gmol)"
             << "|(kcal/gmol)(kcal/gmol)(cal/gmolK)"
             << "|  (units_ROP)      (kcal/gmol) |"
             << "   (units_ROP)     (kcal/gmol) |";
        cout << endl;
    }
    cout << "|--------|----------"
         << "---------------------|-----------------"
         << "----------------|-----------------------------"
         << "--|-------------------------------|" << endl;

    for (int i = 0; i < TT.size(); i++) {
        cout << "|";
        pr_df(TT[i], 7, 1);
        cout << " |";
        pr_de(FwdROP_Table[i], 10, 2);
        pr_de(Afwd_Table[i], 10, 2);
        pr_dg(EafwddivR_Table[i] * Rkcal, 10, 2);
        cout << " | ";
        pr_dg(deltaGss_Table[i] * conv, 10, 2);
        pr_dg(deltaHss_Table[i] * conv, 10, 2);
        pr_dg(deltaSss_Table[i] * convS, 10, 2);
        cout << "  |";
        pr_de(RevROP_Table[i], 10, 2);
        pr_de(Arev_Table[i], 10, 2);
        pr_dg(EarevdivR_Table[i] * Rkcal, 10, 2);
        cout << " |";
        pr_de(NetROP_Table[i], 10, 2);
        pr_de(Anet_Table[i], 10, 2);
        pr_dg(EanetdivR_Table[i] * Rkcal, 10, 2);
        cout << " |" << endl;
    }
    cout << "|" ;
    print_char('-', tableWidth);
    cout << "|" << endl;
}

/*************************************************************************/
void processGERCurrentVsPotTable(RxnMolChange* rmc,
                                 PhaseList* pl, int iGERrxn,
                                 TemperatureTable& TT,
                                 Kinetics& kin,
                                 ExtraGlobalRxn& egr,
                                 RxnTempTableStuff& rts)
{
    int iph, k, e;
    ThermoPhase* tp = 0;
    ThermoPhase* tpMetal = 0;
    //int nTotalSpecies = kin.nTotalSpecies();
    int nRxns         = kin.nReactions();
    vector<double> Rfwd(nRxns, 0.0);
    vector<double> Rrev(nRxns, 0.0);
    vector<double> Rnet(nRxns, 0.0);
    vector<double> deltaG(nRxns, 0.0);
    vector<double> deltaSSG(nRxns, 0.0);

    InterfaceKinetics* iK = dynamic_cast<InterfaceKinetics*>(&kin);
    if (!iK) {
        throw CanteraError("RxnMolChange", "unknown condition on charge");
    }
    /*
     * First thing is to figure out how to calculate the voltage
     *    V = phi_metal - phi_soln
     *    Look for the electron species. The phase where this occurs will be
     *    called the metal, iMetal;
     *
     *    The species index within the kinetics object for the electron will
     *    be called kElectron
     */
    // Phase index for the metal
    int iMetal = -1;
    // Phase index for the solution
    int iSoln = 0;
    // Species index for the electron
    int kElectron = -1;
    int nPhases = iK->nPhases();

    for (iph = 0; iph < nPhases; iph++) {
        tp = &(iK->thermo(iph));
        std::string pName = tp->id();
        int nSpecies = tp->nSpecies();
        int nElements = tp->nElements();
        int eElectron = tp->elementIndex("E");
        if (eElectron >= 0) {
            for (k = 0; k < nSpecies; k++) {
                if (fabs(tp->nAtoms(k,eElectron) - 1) < 1.0E-10) {
                    int ifound = 1;
                    for (e = 0; e < nElements; e++) {
                        if (fabs(tp->nAtoms(k,e)) > 1.0E-10) {
                            if (e != eElectron) {
                                ifound = 0;
                            }
                        }
                    }
                    if (ifound == 1) {
                        if (iMetal == -1) {
                            iMetal = iph;
                            kElectron = iK->kineticsSpeciesIndex(k, iph);
                        } else {
                            int pi = iK->phaseIndex(pName);
                            if (rmc->m_phaseChargeChange[pi] != 0.0) {
                                iMetal = iph;
                                kElectron = iK->kineticsSpeciesIndex(k, iph);
                            }
                        }
                    }
                }
            }
        }
        if (iph != iMetal) {
            if (rmc->m_phaseChargeChange[iph] != 0.0) {
                iSoln = iph;
            }
        }
    }
    tpMetal = &(iK->thermo(iMetal));
    ThermoPhase* tpSoln =  &(iK->thermo(iSoln));
    double nStoichElectrons = - rmc->m_phaseChargeChange[iMetal];
    double phi0Metal = tpMetal->electricPotential();
    //double phi0Metal = rmc->m_phasePotentials[iMetal];
    double phi0Soln = tpSoln->electricPotential();
    //double phi0Soln = rmc->m_phasePotentials[iSoln];
    double V0 = phi0Metal -  phi0Soln;
    int iSurf = iK->reactionPhaseIndex();
    ThermoPhase* tps = &(iK->thermo(iSurf));
    SurfPhase* tpSurface = dynamic_cast<SurfPhase*>(tps);

    int nSurfSpecies = tpSurface->nSpecies();
    vector<double> thetaSurf(nSurfSpecies, 0.0);
    tpSurface->getCoverages(DATA_PTR(thetaSurf));

    int pwidth = 100;
    int twidth = 93;
    dnt(1);
    print_char('-', pwidth);
    printf("\n");
    std::string rs = egr.reactionString();
    dnt(1);
    printf("Voltage Dependence for Rxn: %s", rs.c_str());
    printf("\n\n");
    dnt(2);
    printf("ThermoPhase index pertaining to the solution = %d\n", iSoln);
    dnt(2);
    printf("ThermoPhase index pertaining to the metal = %d\n", iMetal);
    dnt(2);
    printf("Kinetic Species Index for the electrons = %d\n", kElectron);
    dnt(2);
    printf("Input Electric Potential of the Metal = %g Volts\n", phi0Metal);
    dnt(2);
    printf("Input Electric Potential of the Soln  = %g Volts\n", phi0Soln);
    dnt(2);
    printf("Input Voltage                         = %g Volts\n", V0);
    dnt(2);
    printf("nStoicElectrons (product electrons - reactant electrons) = %g \n", nStoichElectrons);

    // dnt(2); printf("Number of reactant Electrons          = %g \n", nElectrons);
    /*
     * the 1.0E-4 is to change from m-2 to cm-2.
     */
    iK->getDeltaGibbs(DATA_PTR(deltaG));
    double deltaGrxn = egr.deltaRxnVecValue(DATA_PTR(deltaG));

    iK->getDeltaSSGibbs(DATA_PTR(deltaSSG));
    double deltaSSGrxn = egr.deltaRxnVecValue(DATA_PTR(deltaSSG));


    double Erxn = deltaGrxn/Faraday/nStoichElectrons;
    double ESSrxn = deltaSSGrxn/Faraday/nStoichElectrons;
    dnt(2);
    printf("Delta G (rxn)   = %13.5g J/kmol   --->   ", deltaGrxn);
    dnt(2);
    printf("E    = %g Volts\n", Erxn);
    dnt(2);
    printf("Delta G_SS (rxn)= %13.5g J/kmol   --->   ", deltaSSGrxn);
    dnt(2);
    printf("E_SS = %g Volts\n", ESSrxn);

    dnt(2);
    printf("Bath Conditions:\n");
    dnt(3);
    printf("T = %g K\n", tpSoln->temperature());
    dnt(3);
    printf("P = %g Pascal\n", tpSoln->pressure());
    dnt(3);
    printf("     SurfSpecName         Theta_k_Bath\n");
    for (k = 0; k < nSurfSpecies; k++) {
        std::string sName = tpSurface->speciesName(k);
        dnt(3);
        printf("%16s %12.5g\n", sName.c_str(), thetaSurf[k]);
    }
    printf("\n");



    /*
     *  n here is defined as the stoichiometric coefficient for the electron in the global
     *  reaction. Therefore n = -1 means the reaction is written in the cathodic direction,
     *  where electrons are consumed.
     *  n = 1 means that the reaction is written in the anodic direction, where electrons
     *  are created.
     */
    double nF = - Faraday * rmc->m_phaseChargeChange[iMetal] * 1.0E-4;

    double Voltage, phiMetal;

#ifdef DEBUG_HKM_ETABLE
    FILE* EFP = fopen("polarizationData.csv", "w");
#endif
#ifdef DEBUG_JCH_ETABLE
    FILE* EFP = fopen("polarizationData_GER.dat", "w");
    fprintf(EFP,"TITLE = \"EMF  = %g Volts\"\n", Erxn);
#endif
    dnt(2);
    printf(" TABLE that uses the Surface Solver to Solve for Surface Site Concentrations\n");


    int nSS = std::min(nSurfSpecies, 6);
    twidth = 104 + nSS * 11;
    dnt(2);
    print_char('-', twidth);
    printf("\n");
    dnt(2);
    cout << "|  Voltage  "
         << "|       Ifwd   Iexchange_f   Beta_f   "
         << "|       Irev   Iexchange_r   Beta_r   "
         << "|    Inet     "
         << "| ";

    for (k = 0; k < nSS; k++) {
        std::string sname =  tpSurface->speciesName(k);
        printf("%10s ", sname.c_str());
    }
    cout <<"|" << endl;

#ifdef DEBUG_HKM_ETABLE
    fprintf(EFP,"  Voltage,    Ifwd,   Iexchange_f,   Beta_f,   "
            "   Irev,   Iexchange_r,   Beta_r,      Inet,  ");
    for (k = 0; k < nSS; k++) {
        std::string sname =  tpSurface->speciesName(k);
        fprintf(EFP, "%10s ", sname.c_str());
        if (k < nSS - 1) {
            fprintf(EFP,", ");
        }
    }
    fprintf(EFP,"\n");
#endif
#ifdef DEBUG_JCH_ETABLE
    fprintf(EFP,"VARIABLES =  \"Voltage\"    \"Ifwd\"   \"Iexchange_f\"   \"Beta_f\"   "
            "   \"Irev\"   \"Iexchange_r\"  \"Beta_r\"      \"Inet\"  ");
    for (k = 0; k < nSS; k++) {
        std::string sname =  tpSurface->speciesName(k);
        fprintf(EFP, "\"%10s\" ", sname.c_str());
    }
    fprintf(EFP,"\n");
#endif

    dnt(2);
    std::cout << "|   Volts   "
         << "|         (amps/cm2)                  "
         << "|         (amps/cm2)                  "
         << "|             "
         << "| ";
    for (k = 0; k < nSS; k++) {
        printf("           ");
    }
    cout <<"|" << endl;
    dnt(2);
    print_char('-', twidth);
    printf("\n");

    if (IOO.VVincEeq) {
        VV_ptr->AddEeq(Erxn);
    }
    if (IOO.VVincEzero) {
        VV_ptr->AddEzero(ESSrxn);
    }

#ifdef CENTERED_VOLTAGE_TABLE
    //Make a new Voltage Table around Erxn
    if (IOO.VVincEeq) {
        double dVtab[7] = { 0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5 };
        for (int i = 0; i < 7; i++) {
            VV_ptr->addPoint(Erxn + dVtab[i]);
            VV_ptr->addPoint(Erxn - dVtab[i]);
        }
    }
#endif // CENTERED_VOLTAGE_TABLE

    int npts = VV_ptr->size();

    for (int iV = 0; iV < npts; iV++) {
        Voltage = (*VV_ptr)[iV];
        phiMetal = (*VV_ptr)[iV] + (iK->thermo(iSoln)).electricPotential();
        //phiMetal = (*VV_ptr)[iV] + rmc->m_phasePotentials[iSoln];
        if (fabs(Voltage) < 1.0E-10) {
            Voltage = 0.0;
        }
        tpMetal->setElectricPotential(phiMetal);
        iK->advanceCoverages(1000.0);

        iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
        iK->getRevRatesOfProgress(DATA_PTR(Rrev));
        double FwdROP = egr.FwdROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
        double RevROP = egr.RevROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));

        iK->getNetRatesOfProgress(DATA_PTR(Rnet));
        double NetROP = egr.ROPValue(DATA_PTR(Rnet));

        double deltaPhi = 1.0E-4;
        double phiMetalDelta = phiMetal + deltaPhi;
        double  VoltageDelta = Voltage + deltaPhi;
        tpMetal->setElectricPotential(phiMetalDelta);
        iK->advanceCoverages(1000.0);
        iK->getFwdRatesOfProgress(DATA_PTR(Rfwd));
        iK->getRevRatesOfProgress(DATA_PTR(Rrev));

        double FwdROP_delta = egr.FwdROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
        double RevROP_delta = egr.RevROPValue(DATA_PTR(Rfwd), DATA_PTR(Rrev));
        double i0_f;
        double alpha_f;
        double i0_r;
        double alpha_r;
        double FdRT = Faraday/(GasConstant * tpSoln->temperature());
        if (nF > 0.0) {
            calcBF(nF*FwdROP, Voltage, nF*FwdROP_delta, VoltageDelta,  Erxn, 1, FdRT,
                   i0_f, alpha_f);
        } else {
            calcBF(nF*FwdROP, Voltage, nF*FwdROP_delta, VoltageDelta,  Erxn, 0, FdRT,
                   i0_f, alpha_f);
        }
        if (nF > 0.0) {
            calcBF(nF*RevROP, Voltage, nF*RevROP_delta, VoltageDelta,  Erxn, 0, FdRT,
                   i0_r, alpha_r);
        } else {
            calcBF(nF*RevROP, Voltage, nF*RevROP_delta, VoltageDelta,  Erxn, 1, FdRT,
                   i0_r, alpha_r);
        }

        dnt(2);
        printf("!%10.4g | %11.5g %11.5g %11.5g | %11.5g %11.5g %11.5g | %11.5g |",
               Voltage, nF*FwdROP, i0_f, alpha_f, -nF*RevROP, i0_r, alpha_r, nF*NetROP);

#ifdef DEBUG_HKM_ETABLE
        fprintf(EFP, "%10.4g , %11.5g, %11.5g, %11.5g, %11.5g, %11.5g, %11.5g, %11.5g,",
                Voltage, nF*FwdROP, i0_f, alpha_f, -nF*RevROP, i0_r, alpha_r, nF*NetROP);
#endif
#ifdef DEBUG_JCH_ETABLE
        fprintf(EFP, "%10.4g  %11.5g %11.5g %11.5g %11.5g %11.5g %11.5g %11.5g",
                Voltage, nF*FwdROP, i0_f, alpha_f, -nF*RevROP, i0_r, alpha_r, nF*NetROP);
#endif

        tpSurface->getCoverages(DATA_PTR(thetaSurf));
        for (k = 0; k < nSS; k++) {
            printf("%10.5g ", thetaSurf[k]);
        }
        printf(" |\n");

#ifdef DEBUG_HKM_ETABLE
        for (k = 0; k < nSS; k++) {
            fprintf(EFP,"%10.5g ", thetaSurf[k]);
            if (k < nSS - 1) {
                fprintf(EFP,", ");
            }
        }
        fprintf(EFP,"\n");
#endif
#ifdef DEBUG_JCH_ETABLE
        for (k = 0; k < nSS; k++) {
            fprintf(EFP,"%10.5g ", thetaSurf[k]);
        }
        fprintf(EFP,"\n");
#endif
    }
    dnt(2);
    print_char('-', twidth);
    printf("\n");



    dnt(1);
    print_char('-', pwidth);
    printf("\n");
#ifdef DEBUG_HKM_ETABLE
    fclose(EFP);
#endif
#ifdef DEBUG_JCH_ETABLE
    fclose(EFP);
#endif

}

