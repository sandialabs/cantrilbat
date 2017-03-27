/**
 *  @file cttables.cpp
 *
 *  $Id: cttables.cpp 554 2013-03-01 23:28:07Z hkmoffa $
 */


/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


#include "cttables.h"

#include "cttInput.h"
#include "cttables_kin.h"
#include "BlockEntryGlobal.h"
#include "cantera/base/clockWC.h"

// Kinetics includes
#include "importAllCTML.h"
#include "importPL.h"
#include "ReactingVolDomain.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

using namespace std;
using namespace BEInput;
using namespace mdpUtil;

std::string InputFile;
std::string TransportFile;
std::string LogFile;
Transport* GTran         = 0;
TemperatureTable* TT_ptr = 0;
VoltageTable* VV_ptr = 0;

Bath BG;
UnitsIO UIO;
int DebugPrinting(0);

/**!
 * Driving object is a cantera Kinetics object.
 * If false, driving object is a Cantera ThermoPhase object
 */
bool TopIsKineticsObject = false;

std::vector<double> S298;

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void print_char(const char c, const int nTimes)
{
    for (int i = 0; i < nTimes; i++) {
        cout << c;
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void print_bool(const bool b)
{
    if (b) {
        cout << "yes";
    } else {
        cout << "no ";
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void pr_if(const int i, const int w)
{
    cout.width(w);
    cout << i;
}
//===================================================================================================================================
// Print a string with a minimum string width
//    -> note the maximum string width is not controlled
void pr_sf(const std::string s, const int w)
{
    int sz = s.size();
    if (sz < w) {
        int num = w - sz;
        for (int i = 0; i < num; i++) {
            cout << " ";
        }
    }
    cout << s;
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
void pr_sf_lj(const string s, const int w, const int crop)
{
    int sz = s.size();
    if (crop && sz > w) {
        const char* pos = s.c_str();
        for (int i = 0; i < w; i++, pos++) {
            cout << *pos;
        }
    } else {
        cout << s;
    }
    if (sz < w) {
        int num = w - sz;
        for (int i = 0; i < num; i++) {
            cout << " ";
        }
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/*
 * pr_df():
 *   print with a fixed precision and width in fixed format.
 *      d = value
 *      w = width
 *      p = precision
 */
void pr_df(const double d, const int w, const int p)
{
    int wmax = w;
    if (d < 0.0) {
        wmax--;
    }
    double dlnd10 = 0.0;
    double dabs = fabs(d);
    if (dabs >=10.0) {
        dlnd10 = log10(dabs);
    }
    int idlnd10 = (int) dlnd10;
    int pmax = wmax - 2 - idlnd10;
    int puse = p;
    if (puse > pmax) {
        puse = pmax;
    }
    printf("%*.*f", w, puse, d);
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/**
 * pr_dfp()
 *   Print with a fixed precision and a variable width in a
 *   fixed format (i.e., non scientific notation)
 */

void pr_dfp(const double d, const int p)
{
    //int pp = cout.precision(p);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout << d;
    cout.precision(6);
    cout.setf(ios_base::fmtflags(0), ios_base::floatfield);
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void pr_de(const double d, const int w, const int p)
{
    cout.setf(ios_base::scientific | ios_base::uppercase);
    cout.width(w);
    int pp = cout.precision(p);
    cout <<  d;
    pp = cout.precision(pp);
    cout.unsetf(ios_base::scientific);
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void pr_dg(const double d, const int w, const int p)
{
    if (d == 0.0) {
        pr_df(d, w, p);
    } else if (fabs(d) < 0.1) {
        pr_de(d, w, p);
    } else if (fabs(d) > 1000.) {
        pr_de(d, w, p);
    } else {
        pr_df(d, w, p);
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void dnt(const int i)
{
    if (i == 0) {
        /* Section headings */
        print_char(' ', 5);
    } else if (i == 1) {
        /* informative text below section headings */
        print_char(' ', 8);
    } else if (i == 2) {
        /* indentation of tables themselves */
        print_char(' ', 15);
    } else if (i == 3) {
        /* indentation for txt pertaining to 1 */
        print_char(' ', 18);
    } else if (i == 4) {
        /* indentation for small tables */
        print_char(' ', 10);
    } else {
        print_char(' ', 8);
    }
}
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

void print_map(const map<string,double>& m, const string& prefix)
{
    if (prefix.size() > 0) {
        cout << prefix;
    }
    map<string,double>::const_iterator it;
    for (it = m.begin(); it != m.end(); it++) {
        if (it != m.begin()) {
            cout << " ";
        }
        cout << "(";
        cout << it->first << ", " << it->second << ")";
    }
}
/**********************************************************************/
/******* bltest_blessed.out
***************************************************************/
/**********************************************************************/
/*
 *  print_phase:
 *        Phase is an enumerated int, as currently defined in
 *        Cantera. This program interprets the int, and prints
 *        a string
 */
//void print_phase(const int p) {
// if      (p == GAS)         cout << "GAS";
// else if (p == LIQUID)      cout << "LIQUID";
// else if (p == SOLID)       cout << "SOLID";
// else if (p == PURE_FLUID) cout << "PURE_FLUID";
// else {
//   cout << "internal error" << endl;
//   exit(-1);
// }
//}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 *  This routine will print out a table of information about
 *  a species in an ideal gas thermo phse. It explicitly
 *  assumes that a multitransport object has been created for
 *  the phase, and it presumes a NASA polynomial form for the
 *  species thermodynamics.
 */
void printIdealGasSpeciesTable(ThermoPhase& g,
                               int k, TemperatureTable& TT,
                               DenseMatrix& Cp_Table,
                               DenseMatrix& Hrel_Table,
                               DenseMatrix& Grel_Table,
                               DenseMatrix& S_Table,
                               bool haveSpeciesTransportProps,
                               DenseMatrix& Visc_Table,
                               DenseMatrix& Cond_Table,
                               DenseMatrix& Diff_Table,
                               vector<double> H298,
                               DenseMatrix& Urel_Table,
                               vector<double> U298)
{
    int tModel = 0;
    /*
     *  Get the species data object from the Mixture object
     *  this is defined in the constituents.h file, and is
     *  inherited by Mixture through BaseMix
     */
    string sName = g.speciesName(k);
    /*
     *  Dump out all of the information about the species
     */
    cout << endl;
    print_char('=', 120);
    cout << endl;
    cout << "INFORMATION TABLE FOR IDEAL GAS SPECIES \""  << sName;
    cout << "\" IN PHASE \"";
    cout << g.id();
    cout << "\"" << endl;
    dnt(1);
    cout << "Overall, this is the " << k+1
         <<"th species in the mechanism" << endl;
    dnt(1);
    cout << "It is the " << k+1
         <<"th species in the phase" << endl;

    double elementEntropyTotal = entropyElem298(&g, (size_t) k);
    double ch = g.charge(k);
    if (ch != 0.0) {
        double entCh = g.entropyCharge(ch);
        elementEntropyTotal += entCh;
    }

    dnt(1);
    cout << "Elemental Composition:" << endl;
    for (int m = 0; m < (int) g.nElements(); m++) {
        double na = g.nAtoms(k, m);
        if (na != 0.0) {
            dnt(3);
            cout << g.elementName(m) << ": " << na <<endl;
        }
    }
    dnt(1);
    cout << "Electronic Charge = ";
    pr_dfp(g.charge(k), 2);
    cout << endl;
    double mw = g.molecularWeight(k);
    if (mw == Tiny) {
        dnt(1);
        cout << "Molecular Weight = 0 (set to " << Tiny << " for massF conversion)"
             << " gm/mol" << endl;
    } else {
        dnt(1);
        cout << "Molecular Weight = " << mw << " gm/mol" << endl;
    }
    if (GTran) {
        tModel = GTran->model();
    }
    if (tModel == cMulticomponent || tModel == CK_Multicomponent) {
        if (haveSpeciesTransportProps) {

//TODO: figure out how to do this w/2.0
            MultiTransport* mt = dynamic_cast<MultiTransport*>(GTran);
            //GasTransportData td = mt->getGasTransportData(k);
            dnt(1);
            double wellDepth = mt->m_eps[k] / Boltzmann;
            cout << "L-J Potential Well Depth = " << wellDepth << " K" << endl;
            dnt(1);
            double diameter = mt->m_sigma[k] * 1.0E10;
            cout << "L-J collision diameter = "
                 << diameter << " Angstroms" << endl;
            dnt(1);
            double dipoleMoment =  1.0E25 / SqrtTen * mt->m_dipole(k,k);
            cout << "Dipole Moment = " << dipoleMoment << " Debeye" << endl;
            dnt(1);
            double polarizability = mt->m_alpha[k] * 1.0E30;
            cout << "Polarizability = ";
            pr_dfp(polarizability, 2);
            cout << " Angstroms**3" << endl;
            dnt(1);
            double rotRelaxNumber = mt->m_zrot[k];
            cout << "Rotational Collision Number at 298 K = "
                 << rotRelaxNumber << endl;
            dnt(1);
            double crot = mt->m_crot[k];
            if (crot == 0.0) {
                cout <<"This molecule is monatomic" << endl;
            } else if (crot == 1.0) {
                cout <<"This molecule is linear" << endl;
            } else if (crot == 1.5) {
                cout <<"This molecule is nonlinear" << endl;
            } else {
                throw CanteraError("", "Unknown crot = " + ZZCantera::fp2str(crot));
            }

        }
    }
    if (tModel == cMixtureAveraged || tModel == CK_MixtureAveraged) {
        if (haveSpeciesTransportProps) {
            MixTransport* mt = dynamic_cast<MixTransport*>(GTran);

            //	GasTransportData td = mt->getGasTransportData(k);
            dnt(1);
            double wellDepth = mt->m_eps[k] / Boltzmann;
            cout << "L-J Potential Well Depth = "
                 << wellDepth << " K" << endl;
            dnt(1);
            double diameter = mt->m_sigma[k] * 1.0E10;
            cout << "L-J collision diameter = " << diameter
                 << " Angstroms" << endl;
            dnt(1);
            double dipoleMoment = 1.0E25 / SqrtTen * mt->m_dipole(k,k);
            cout << "Dipole Moment = " << dipoleMoment << " Debeye" << endl;
            dnt(1);
            double polarizability = mt->m_alpha[k] * 1.0E30;
            cout << "Polarizability = ";
            pr_dfp(polarizability, 2);
            cout << " Angstroms**3" << endl;
            dnt(1);
            double rotRelaxNumber = mt->m_zrot[k];
            cout << "Rotational Collision Number at 298 K = "
                 <<  rotRelaxNumber << endl;
            dnt(1);
            double crot = mt->m_crot[k];
            if (crot == 0.0) {
                cout <<"This molecule is monatomic" << endl;
            } else if (crot == 1.0) {
                cout <<"This molecule is linear" << endl;
            } else if (crot == 1.5) {
                cout <<"This molecule is nonlinear" << endl;
            } else {
                throw CanteraError("", "Unknown crot = " + ZZCantera::fp2str(crot));
            }

        }
    }
    /*
     * Print out the Heat of Formation at 298.15 K
     */
    dnt(1);
    cout << "Heat of formation (298.15K) = ";
    pr_dfp(H298[k], 4);
    printf(" %8s\n", UIO.sGibbs.c_str());
    dnt(1);
    double mu_o = H298[k];
    cout << "mu_0(298.15K) = H298[k]  " << mu_o << endl;

    double deltaGf = 1.0E6 * H298[k] - 298.15 * S298[k] + 298.15 * elementEntropyTotal;
    deltaGf /= 1.0E6;
    dnt(1);
    cout << "DeltaGf (298.15K) = ";
    pr_dfp(deltaGf, 4);
    printf(" %8s\n", UIO.sGibbs.c_str());
    dnt(1);

    /*
     * Optionally print out the internal energy at 298.15 K
     */
    if (IOO.IntEnergyColumn) {
        dnt(1);
        cout << "Int Energy (298.15K, refPres) = ";
        pr_dfp(U298[k], 4);
        printf(" %8s\n", UIO.sGibbs.c_str());
    }

    /*
     * Print out the reference pressure
     */
    SpeciesThermo& sThermo = g.speciesThermo();
    double presRef = sThermo.refPressure(k);
    double presRefPhase = g.refPressure();
    if (IOO.OutputUnits == UNITS_KCAL_CGS) {
        if (doubleEqual(presRef, presRefPhase)) {
            dnt(1);
            cout << "Phase and Species Reference Pressure = "
                 << presRefPhase *10. << " erg cm-2" << endl;
        } else {
            dnt(1);
            cout << "Phase Reference Pressure = "
                 << presRefPhase *10. << " erg cm-2" << endl;
            dnt(1);
            cout << "Species Reference Pressure = "
                 << presRef  *10. << " erg cm-2" << endl;
        }
    } else  if (IOO.OutputUnits == UNITS_KJOULE) {
        if (doubleEqual(presRef, presRefPhase)) {
            dnt(1);
            cout << "Phase and Species Reference Pressure = "
                 << presRefPhase << " Pa" << endl;
        } else {
            dnt(1);
            cout << "Phase Reference Pressure = "
                 << presRefPhase << " Pa" << endl;
            dnt(1);
            cout << "Species Reference Pressure = "
                 << presRef  << " Pa" << endl;
        }
    }
    double minTemp = g.minTemp(k);
    double maxTemp = g.maxTemp(k);
    dnt(1);
    printf("Minimum Temperature = %g K\n", minTemp);
    dnt(1);
    printf("Maximum Temperature = %g K\n", maxTemp);

    printThermoCoeffSpecies(&g, k);

    /*
     *  Dump out a species thermo and transport table
     */
    int tableWidth = 120;
    print_char('-', 120);
    if (IOO.IntEnergyColumn) {
        tableWidth += 15;
        print_char('-', 15);
    }
    if (IOO.ChemPotColumn) {
        tableWidth += 15;
        print_char('-', 15);
    }
    printf("\n");

    double ptable = presRef;
    if (IOO.UseRefPressureInThermoTables) {
        printf("|------------ Thermo Functions for Reference Pressure");
    } else {
        printf("|----------- Thermo Functions for BathSpecies Pressure");
        ptable = BG.Pressure;
    }
    if (IOO.OutputUnits == UNITS_KCAL_CGS) {
        printf(" %11.3g erg cm-2", ptable * 10.0);
    } else {
        printf(" %11.3g Pa -----", ptable);
    }
    for (int i = 0; i < tableWidth-79; i++) {
        printf("-");
    }
    printf("----|\n");

    if (haveSpeciesTransportProps) {
        printf("|     Temp  |   (H-H298)        (G-H298)            Cp      "
               "   S       ");
        if (IOO.IntEnergyColumn) {
            printf("|  (U - U298)  ");
        }
        if (IOO.ChemPotColumn) {
            printf("|     G_abs    ");
        }
        printf("|  Viscosity  Therm_Cond   Dif_Co_with_");;
        string tmp = g.speciesName(BG.BathSpeciesID);
        pr_sf_lj(tmp, 9, 1);
        cout << "|" << endl;
        if (IOO.OutputUnits == UNITS_KCAL_CGS) {
            printf("|      (K)  |  (kcal/mol)     (kcal/mol)     (cal/mol*K)    "
                   "(cal/mol*K)");
            if (IOO.IntEnergyColumn) {
                printf("|  (kcal/mol)  ");
            }
            if (IOO.ChemPotColumn) {
                printf("|  (kcal/mol)  ");
            }
            printf("|   (gm/cm*sec) (erg/cm*sec*K)   (cm**2/sec)    ");
            printf("|\n");
        } else  if (IOO.OutputUnits == UNITS_KJOULE) {
            printf("|      (K)  |  (kJ/gmol)      (kJ/gmol)      (J/gmol*K)     "
                   "(J/gmol*K) ");
            if (IOO.IntEnergyColumn) {
                printf("|  (kJ/gmol)  ");
            }
            if (IOO.ChemPotColumn) {
                printf("|  (kJ/gmol)  ");
            }
            printf("|   (kg/m*sec)    (J/m*sec*K)     (m**2/sec)    ");
            printf("|\n");
        }
        printf("|-----------|-----------------------------------------------"
               "-----------");
        if (IOO.IntEnergyColumn) {
            printf("|--------------");
        }
        if (IOO.ChemPotColumn) {
            printf("|--------------");
        }
        printf("|-----------------------------------------------");
        printf("|\n");
    } else {
        printf("|     Temp  |   (H-H298)        (G-H298)            Cp      "
               "   S       ");
        if (IOO.IntEnergyColumn) {
            printf("|  (U - U298)    ");
        }
        if (IOO.ChemPotColumn) {
            printf("|      G_abs     ");
        }
        printf("|\n");
        if (IOO.OutputUnits == UNITS_KCAL_CGS) {
            printf("|      (K)  |  (kcal/mol)     (kcal/mol)     (cal/mol*K)    "
                   "(cal/mol*K)");
            if (IOO.IntEnergyColumn) {
                printf("|   (kcal/mol)  ");
            }
            if (IOO.ChemPotColumn) {
                printf("|   (kcal/mol)  ");
            }
            printf("|\n");
        } else  if (IOO.OutputUnits == UNITS_KJOULE) {
            printf("|      (K)  |   (kJ/gmol)      (kJ/gmol)      (J/gmol*K)    "
                   "(J/gmol*K) ");
            if (IOO.IntEnergyColumn) {
                printf("|   (kJ/gmol)  ");
            }
            if (IOO.ChemPotColumn) {
                printf("|   (kJ/gmol)  ");
            }
            printf("|\n");
        }
        printf("|-----------|-----------------------------------------------"
               "-----------");
        if (IOO.IntEnergyColumn) {
            printf("|--------------");
        }
        if (IOO.ChemPotColumn) {
            printf("|--------------");
        }
        printf("|\n");
    }

    bool outOfRange = false;
    for (int i = 0; i < TT.size(); i++) {
        outOfRange = false;
        if (TT[i] + 1.0E-3 < minTemp) {
            outOfRange = true;
        }
        if (TT[i] - 1.0E-3 > maxTemp) {
            outOfRange = true;
        }
        if (outOfRange) {
            cout << "|**";
            pr_df(TT[i], 8, 2);
            cout << "*|";
        } else {
            cout << "|  ";
            pr_df(TT[i], 8, 2);
            cout << " |";
        }
        pr_df(Hrel_Table(i,k), 11, 4);
        pr_df(Grel_Table(i,k), 16, 4);
        pr_df(Cp_Table(i,k), 15, 4);
        pr_df(S_Table(i,k), 13, 4);
        printf("   ");
        if (IOO.IntEnergyColumn) {
            printf("|");
            double valu = Urel_Table(i,k);
            pr_df(valu, 13, 4);
            printf(" ");
        }
        if (IOO.ChemPotColumn) {
            printf("|");
            double valmu = Grel_Table(i,k) + H298[k];
            pr_df(valmu, 13, 4);
            printf(" ");
        }
        printf("| ");
        if (haveSpeciesTransportProps) {
            pr_de(Visc_Table(i,k), 12, 3);
            pr_de(Cond_Table(i,k), 14, 3);
            pr_dg(Diff_Table(i,k), 14, 3);
            cout << "      |";
        }
        printf("\n");
    }
    print_char('-', 120);
    if (IOO.IntEnergyColumn) {
        print_char('-', 15);
    }
    if (IOO.ChemPotColumn) {
        print_char('-', 15);
    }
    printf("\n");
}

//=====================================================================================================================================
/**
 *  This routine will print out a table of information about
 *  a species in an ideal gas thermo phse. It explicitly
 *  assumes that a multitransport object has been created for
 *  the phase, and it presumes a NASA polynomial form for the
 *  species thermodynamics.
 */
void printThermoPhaseSpeciesTable(ThermoPhase* g_ptr,
                                  int k, TemperatureTable& TT,
                                  DenseMatrix& Cp_Table,
                                  DenseMatrix& Hrel_Table,
                                  DenseMatrix& Grel_Table,
                                  DenseMatrix& S_Table,
                                  bool haveSpeciesTransportProps,
                                  DenseMatrix& Visc_Table,
                                  DenseMatrix& Cond_Table,
                                  DenseMatrix& Diff_Table,
                                  vector<double> H298,
                                  DenseMatrix& Urel_Table,
                                  vector<double> U298)
{
    /*
     *  Get the species data object from the Mixture object
     *  this is defined in the constituents.h file, and is
     *  inherited by Mixture through BaseMix
     */
    string sName = g_ptr->speciesName(k);
    int tableWidth = 72;
    if (haveSpeciesTransportProps) {
        tableWidth += 30;
    }
    if (IOO.IntEnergyColumn) {
        tableWidth += 15;
    }
    if (IOO.ChemPotColumn) {
        tableWidth += 15;
    }

    /*
     *  Dump out all of the information about the species
     */
    cout << endl;
    print_char('=', tableWidth);
    cout << endl;
    cout << "INFORMATION TABLE FOR SPECIES \""  << sName;
    cout << "\" IN PHASE \"";
    cout << g_ptr->id();
    cout << "\"" << endl;
    dnt(1);
    cout << "Overall, this is the " << k+1
         <<"th species in mechanism" << endl;
    dnt(1);
    cout << "It is the " << k+1
         <<"th species in the phase" << endl;

    double elementEntropyTotal = entropyElem298(g_ptr, (size_t) k);
    dnt(1);
    cout << "Elemental Composition:               ElementEntropy298" << endl;
    for (size_t m = 0; m < g_ptr->nElements(); m++) {
        double na = g_ptr->nAtoms(k, m);
        if (na != 0.0) {
            double se = g_ptr->entropyElement298(m, true);
            dnt(3);
            pr_sf(g_ptr->elementName(m), 3);
            cout << ": " << na << "      |              ";
            if (se == ENTROPY298_UNKNOWN) {
                printf("[UNAVAILABLE] ");
            } else {
                pr_df(se * na / 1.0E3, 13, 3);
            }
            printf("\n");
        }
    }
    double ch = g_ptr->charge(k);
    if (ch != 0.0 &&  elementEntropyTotal != ENTROPY298_UNKNOWN) {
        double entCh = g_ptr->entropyCharge(ch);
        dnt(3);
        pr_sf("E-", 3);
        cout << ": " << ch << "      |              ";
        pr_df(-entCh / 1.0E3, 13, 3);
        printf("\n");
        elementEntropyTotal += entCh;
    }
    dnt(3);
    printf("                           ");
    if (elementEntropyTotal == ENTROPY298_UNKNOWN) {
        printf("    [UNAVAILABLE] ");
    } else {
        pr_df(elementEntropyTotal / 1.0E3, 13, 3);
    }
    printf(" J/gmol/K \n");

    dnt(1);
    cout << "Electronic Charge = ";
    pr_dfp(g_ptr->charge(k), 2);
    cout << endl;
    double mw = g_ptr->molecularWeight(k);
    if (mw == Tiny) {
        dnt(1);
        cout << "Molecular Weight = 0 (set to " << Tiny << " for massF conversion)"
             << " gm/mol" << endl;
    } else {
        dnt(1);
        cout << "Molecular Weight = " << mw << " gm/mol" << endl;
    }


    /*
     * Print out the Heat of Formation at 298.15 K
     */
    dnt(1);
    cout << "Heat of formation (298.15K) = ";
    pr_dfp(H298[k], 4);
    printf(" %8s\n", UIO.sGibbs.c_str());
    /*
     * Optionally print out the internal energy at 298.15 K
     */
    if (IOO.IntEnergyColumn) {
        dnt(1);
        cout << "Int Energy (298.15K, refPres) = ";
        pr_dfp(U298[k], 4);
        printf(" %8s\n", UIO.sGibbs.c_str());
    }

    double deltaGf = 1.0E6 * H298[k] - 298.15 * S298[k] + 298.15 * elementEntropyTotal;
    deltaGf /= 1.0E6;
    dnt(1);
    cout << "DeltaGf (298.15K) = ";
    pr_dfp(deltaGf, 4);
    printf(" %8s\n", UIO.sGibbs.c_str());
    dnt(1);


    /*
     * Print out the reference pressure
     */
    g_ptr->setTemperature(298.15);  // temp fix for variable reference pressures
    SpeciesThermo& sThermo = g_ptr->speciesThermo();
    double presRef = sThermo.refPressure(k);
    double presRefPhase = g_ptr->refPressure();

    if (IOO.OutputUnits == UNITS_KCAL_CGS) {
        if (doubleEqual(presRef, presRefPhase)) {
            dnt(1);
            cout << "Phase and Species Reference Pressure = "
                 << presRefPhase *10. << " erg cm-2" << endl;
        } else {
            dnt(1);
            cout << "Phase Reference Pressure = "
                 << presRefPhase *10. << " erg cm-2" << endl;
            dnt(1);
            cout << "Species Reference Pressure = "
                 << presRef  *10. << " erg cm-2" << endl;
        }
    } else  if (IOO.OutputUnits == UNITS_KJOULE) {
        if (doubleEqual(presRef, presRefPhase)) {
            dnt(1);
            cout << "Phase and Species Reference Pressure = "
                 << presRefPhase << " Pa" << endl;
        } else {
            dnt(1);
            cout << "Phase Reference Pressure = "
                 << presRefPhase << " Pa" << endl;
            dnt(1);
            cout << "Species Reference Pressure = "
                 << presRef  << " Pa" << endl;
        }
    }
    double minTemp = g_ptr->minTemp(k);
    double maxTemp = g_ptr->maxTemp(k);
    dnt(1);
    printf("Minimum Temperature = %g K\n", minTemp);
    dnt(1);
    printf("Maximum Temperature = %g K\n", maxTemp);


    /*
     * Add to the species table all of the parameter
     * information about
     * the calculation of the partial molar volume and other
     * related volumetric information for this species in the
     * current phase.
     */
    g_ptr->setState_TP(BG.Temperature, BG.Pressure);
    printVolSpecies(g_ptr, k);

    /*
     * Dump the coefficients out
     */
    printThermoCoeffSpecies(g_ptr, k);

    /*
     *  Dump out a species thermo and transport table
     */
    print_char('-', tableWidth);
    printf("\n");
    double ptable;
    if (IOO.UseRefPressureInThermoTables) {
        printf("|------------ Thermo Functions for Reference Pressure");
        ptable = presRef;
    } else {
        printf("|----------- Thermo Functions for BathSpecies Pressure");
        ptable = BG.Pressure;
    }
    if (IOO.OutputUnits == UNITS_KCAL_CGS) {
        printf(" %11.3g erg cm-2", ptable * 10.0);
    } else {
        printf(" %11.3g Pa -----", ptable);
    }

    for (int i = 0; i < tableWidth-79; i++) {
        printf("-");
    }
    printf("----|\n");
    print_char('-', tableWidth);
    printf("\n");
    if (haveSpeciesTransportProps) {
        printf("|     Temp  |   (H-H298)        (G-H298)            Cp      "
               "   S       ");
        if (IOO.IntEnergyColumn) {
            printf("|  (U - U298)   ");
        }
        if (IOO.ChemPotColumn) {
            printf("|     G_abs     ");
        }
        printf("|  Viscosity  Therm_Cond   Dif_Co_with_");
        string tmp = g_ptr->speciesName(BG.BathSpeciesID);
        pr_sf_lj(tmp, 9, 1);
        printf("|\n");
        if (IOO.OutputUnits == UNITS_KCAL_CGS) {
            printf("|      (K)  |  (kcal/mol)     (kcal/mol)     (cal/mol*K)    "
                   "(cal/mol*K)");
            if (IOO.IntEnergyColumn) {
                printf("|  (kcal/mole)  ");
            }
            if (IOO.ChemPotColumn) {
                printf("|  (kcal/mole)  ");
            }
            printf("|  (gm/cm*sec) (erg/cm*sec*K)   (cm**2/sec)    ");
            printf("|\n");
        } else  if (IOO.OutputUnits == UNITS_KJOULE) {
            printf("|      (K)  |  (kJ/gmol)      (kJ/gmol)      (J/gmol*K)     "
                   "(J/gmol*K) ");
            if (IOO.IntEnergyColumn) {
                printf("|  (kJ/gmol)    ");
            }
            if (IOO.ChemPotColumn) {
                printf("|  (kJ/gmol)    ");
            }
            printf("|  (kg/m*sec)    (J/m*sec*K)     (m**2/sec)    ");
            printf("|\n");
        }
    } else {
        printf("|     Temp  |   (H-H298)        (G-H298)            Cp      "
               "   S       ");
        if (IOO.IntEnergyColumn) {
            printf("|  (U - U298)  ");
        }
        if (IOO.ChemPotColumn) {
            printf("|    G_abs     ");
        }
        printf("|\n");
        if (IOO.OutputUnits == UNITS_KCAL_CGS) {
            printf("|      (K)  |  (kcal/mol)     (kcal/mol)     (cal/mol*K)    "
                   "(cal/mol*K)");
            if (IOO.IntEnergyColumn) {
                printf("|  (kcal/mole) ");
            }
            if (IOO.ChemPotColumn) {
                printf("|  (kcal/mole) ");
            }
            printf("|\n");
        } else  if (IOO.OutputUnits == UNITS_KJOULE) {
            printf("|      (K)  |   (kJ/gmol)      (kJ/gmol)      (J/gmol*K)    "
                   "(J/gmol*K) ");
            if (IOO.IntEnergyColumn) {
                printf("| (kJ/gmol)    ");
            }
            if (IOO.ChemPotColumn) {
                printf("| (kJ/gmol)    ");
            }
            printf("|\n");
        }
        cout << "|-----------|-----------------------------------------------"
             << "-----------";
        if (IOO.IntEnergyColumn) {
            printf("|--------------");
        }
        if (IOO.ChemPotColumn) {
            printf("|--------------");
        }
        printf("|\n");
    }
    bool outOfRange = false;
    for (int i = 0; i < TT.size(); i++) {
        outOfRange = false;
        if (TT[i] + 1.0E-3 < minTemp) {
            outOfRange = true;
        }
        if (TT[i] - 1.0E-3 > maxTemp) {
            outOfRange = true;
        }
        if (outOfRange) {
            cout << "|**";
            pr_df(TT[i], 8, 2);
            cout << "*|";
        } else {
            cout << "|  ";
            pr_df(TT[i], 8, 2);
            cout << " |";
        }
        pr_df(Hrel_Table(i,k), 11, 4);
        pr_df(Grel_Table(i,k), 16, 4);
        pr_df(Cp_Table(i,k), 15, 4);
        pr_df(S_Table(i,k), 13, 4);
        printf("   ");
        if (IOO.IntEnergyColumn) {
            printf("|");
            double valu = Urel_Table(i,k);
            pr_df(valu, 13, 4);
            printf(" ");
        }
        if (IOO.ChemPotColumn) {
            printf("|");
            double valmu = Grel_Table(i,k) + H298[k];
            pr_df(valmu, 13, 4);
            printf(" ");
        }
        printf("| ");
        if (haveSpeciesTransportProps) {
            pr_de(Visc_Table(i,k), 12, 3);
            pr_de(Cond_Table(i,k), 14, 3);
            pr_dg(Diff_Table(i,k), 14, 3);
            printf("      |");
        }
        cout << endl;
    }
    print_char('-', tableWidth);
    printf("\n");
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/*
 * Get the Thermodynamic tables for all species.
 * Note, I have switched over to using functions from the
 * ThermoPhase base class exclusively. It works! and therefore,
 * I can use this for other ThermoPhase derived classes.
 * This works at either the reference pressure for all species
 * or a single pressure specified by the input file.
 */
void
getThermoTables(TemperatureTable& TT, DenseMatrix& Cp_Table,
                DenseMatrix& Hrel_Table, DenseMatrix& Grel_Table,
                DenseMatrix& S_Table, vector<double>& H298,
                ThermoPhase& g, DenseMatrix& Urel_Table,
                vector<double>& U298)
{

    int nsp = g.nSpecies();
    double T298 = 298.15;
    double* h_RT = mdp_alloc_dbl_1(nsp, 0.0);
    double* u_RT = mdp_alloc_dbl_1(nsp, 0.0);
    double* S_R  = mdp_alloc_dbl_1(nsp, 0.0);
    double* g_RT = mdp_alloc_dbl_1(nsp, 0.0);
    double* cp_R = mdp_alloc_dbl_1(nsp, 0.0);

    /*
     * Get units converters
     *     Gibbs Functions:
     */
    double RGibbs = UIO.mGibbs;

    /*
     * First we load a vector of standard state properties at 298.15.
     * We are dropping down here to the ThermoPhase class to call
     * their thermodyanmics routine. We don't need the concept
     * of a mixture thermo here, except for the fact that for every
     * phase, there is a unique definition of what the standard
     * state is.
     */
    double P = BG.Pressure;
    if (IOO.UseRefPressureInThermoTables) {
        P = g.refPressure();
    }
    double* y =  mdp_alloc_dbl_1(nsp, 0.0);
    g.getReferenceComposition(y);

    /*
     * We need to set the state here. Note, mass fractions
     * don't matter, since we are only querying the
     * standard state of the species.
     */
    g.setState_TPX(T298, P, y);
    g.getEnthalpy_RT(h_RT);
    g.getEnthalpy_RT_ref(h_RT);
    g.getEntropy_R_ref(DATA_PTR(S298));
    g.getGibbs_RT_ref(g_RT);


    if (IOO.IntEnergyColumn) {
        g.getIntEnergy_RT_ref(u_RT);
    }
    for (int i = 0; i < nsp; i++) {
        H298[i] = h_RT[i] * T298 * RGibbs;
        U298[i] = u_RT[i] * T298 * RGibbs;
        S298[i] *= GasConstant;
    }

    g.getEntropy_R_ref(S_R);
    g.getGibbs_RT(g_RT);
    double T;
    for (int i = 0; i < TT.size(); i++) {
        T = TT[i];
        //T = 440;
        g.setState_TPX(T, P, y);
        //g.getEnthalpy_RT(h_RT);
        g.getEntropy_R_ref(S_R);
        g.getGibbs_RT_ref(g_RT);
        g.getEnthalpy_RT_ref(h_RT);
        if (IOO.IntEnergyColumn) {
            g.getIntEnergy_RT_ref(u_RT);
        }
        g.getCp_R_ref(cp_R);
        /*
         * Now, copy the results into the matrices, and
         * zero off the 298 K heat of formation. Also,
         * translate to real units.
         */
        for (int j = 0; j < nsp; j++) {
            Cp_Table(i, j)   = cp_R[j] * UIO.mEntropy;
            Hrel_Table(i, j) = h_RT[j] * RGibbs * T;
            Urel_Table(i, j) = u_RT[j] * RGibbs * T;
            Grel_Table(i, j) = g_RT[j] * RGibbs * T;
            S_Table(i, j)    = S_R[j]  * UIO.mEntropy;
            Hrel_Table(i, j) -= H298[j];
            Urel_Table(i, j) -= U298[j];
            Grel_Table(i, j) -= H298[j];
        }
    }

    if (IOO.UseRefPressureInThermoTables) {
        SpeciesThermo& sThermo = g.speciesThermo();
        for (int k = 0; k < nsp; k++) {
            double presRef = sThermo.refPressure(k);
            if (! doubleEqual(P, presRef)) {
                g.setState_TPX(T298, presRef, y);
                g.getEnthalpy_RT_ref(h_RT);
                H298[k] = h_RT[k] * T298 * RGibbs;
                g.getEntropy_R_ref(S_R);
                g.getGibbs_RT_ref(g_RT);
                double T;
                for (int i = 0; i < TT.size(); i++) {
                    T = TT[i];
                    g.setState_TPX(T, presRef, y);
                    g.getEntropy_R_ref(S_R);
                    g.getGibbs_RT_ref(g_RT);
                    g.getEnthalpy_RT_ref(h_RT);
                    if (IOO.IntEnergyColumn) {
                        g.getIntEnergy_RT_ref(u_RT);
                    }
                    g.getCp_R_ref(cp_R);
                    /*
                     * Now, copy the results into the matrices, and
                     * zero off the 298 K heat of formation. Also,
                     * translate to real units.
                     */
                    Cp_Table(i, k)   = cp_R[k] * UIO.mEntropy;
                    Hrel_Table(i, k) = h_RT[k] * RGibbs * T;
                    Urel_Table(i, k) = u_RT[k] * RGibbs * T;
                    Grel_Table(i, k) = g_RT[k] * RGibbs * T;
                    S_Table(i, k)    = S_R[k]  * UIO.mEntropy;
                    Hrel_Table(i, k) -= H298[k];
                    Urel_Table(i, k) -= U298[k];
                    Grel_Table(i, k) -= H298[k];
                }
            }
        }
    }

    mdp_safe_free((void**) &y);
    mdp_safe_free((void**) &h_RT);
    mdp_safe_free((void**) &u_RT);
    mdp_safe_free((void**) &S_R);
    mdp_safe_free((void**) &g_RT);
    mdp_safe_free((void**) &cp_R);
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 * Setup the transport tables for gases. This function assumes
 * that the phase is a gas, and that the Multitransport class
 * handles the calculation of transport properties.
 */
void
getGasTransportTables(TemperatureTable& TT, ThermoPhase& g,
                      DenseMatrix& Visc_Table,
                      DenseMatrix& Cond_Table,
                      DenseMatrix& Diff_Table)
{

    Transport* t = GTran;
    MultiTransport* mt = dynamic_cast<MultiTransport*>(t);
    int nsp = g.nSpecies();
    int nsp2 = nsp * nsp;
    int bgs = BG.BathSpeciesID;
    double* darray = new double [nsp2];
    double* svisc = new double [nsp];
    double* y = new double[nsp];
    for (int i = 0; i < nsp; i++) {
        y[i] = 0.0;
    }

    double T, kspec;
    for (int i = 0; i < TT.size(); i++) {
        T = TT[i];
        g.setState_TPX(T, BG.Pressure, BG.Xmol);
        /*
         * Need to force an update of the thermo functions
         * manually, here, for it to work.
         */
        /* g.updateTransport_T(); */
        mt->getSpeciesViscosities(svisc);
        mt->getBinaryDiffCoeffs(nsp, darray);
        /*
         * Now store the results in the transport tables
         * Have to calculate the thermal conductivity by
         * changing the bath gas to pure species, k. Takes
         * a lot of work.
         *  Note change to cgs units.
         */
        for (int j = 0; j < nsp; j++) {
            Diff_Table(i, j)   = darray[bgs * nsp + j] * 1.0E4;
            Visc_Table(i, j) = svisc[j] * 10.;
            y[j] = 1.0;
            g.setState_TPY(T, BG.Pressure, y);
            kspec = mt->thermalConductivity();
            Cond_Table(i, j) = kspec * 1.0E5;
            y[j] = 0.0;
        }
    }
    g.setState_TPX(BG.Temperature, BG.Pressure, BG.Xmol);
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/**
 * Setup the transport tables for gases. This function assumes
 * that the phase is a gas, and that the Multitransport class
 * handles the calculation of transport properties.
 */
void
getGenericTransportTables(TemperatureTable& TT, ThermoPhase& g,
                          DenseMatrix& Visc_Table,
                          DenseMatrix& Cond_Table,
                          DenseMatrix& Diff_Table)
{
    int i, k;
    Transport* t = GTran;
    if (!t) {
        cout <<"getGenericTransportTables: Error no transport model initialized"
             << endl;
        exit(-1);
    }
    int tModel = t->model();
    int nsp = g.nSpecies();
    int nsp2 = nsp * nsp;
    int bgs = BG.BathSpeciesID;
    double* darray = new double [nsp2];
    double* svisc = new double [nsp];
    double* y = new double[nsp];
    for (int i = 0; i < nsp; i++) {
        y[i] = 0.0;
    }

    double T, kspec;
    for (i = 0; i < TT.size(); i++) {
        T = TT[i];
        g.setState_TPX(T, BG.Pressure, BG.Xmol);
        /*
         * Need to force an update of the thermo functions
         * manually, here, for it to work.
         */
        /*
         * To get the "Species Viscosities", we set the mole fraction
         * of the phase to each of the pure species individually.
         */
        for (k = 0; k < nsp; k++) {
            y[k] = 1.0;
            g.setState_TPX(T, BG.Pressure, y);
            svisc[k] = t->viscosity();
            Visc_Table(i, k) = svisc[k] * 10.;
            y[k] = 0.0;
        }


        switch (tModel) {
        case cMulticomponent:
        case CK_Multicomponent:
        case cMixtureAveraged:
        case CK_MixtureAveraged:
            t->getBinaryDiffCoeffs(nsp, darray);
            for (k = 0; k < nsp; k++) {
                Diff_Table(i, k) = darray[bgs * nsp + k] * 1.0E4;
            }
            break;
        case cSolidTransport:
        default:
            y[BG.BathSpeciesID] = 1.0;
            g.setState_TPX(T, BG.Pressure, y);
            y[k] = 0.0;
            t->getMixDiffCoeffs(darray);
            for (k = 0; k < nsp; k++) {
                Diff_Table(i, k) = darray[k] * 1.0E4;
            }
            break;
        }

        /*
         * Now store the results in the transport tables
         * Have to calculate the thermal conductivity by
         * changing the bath gas to pure species, k. Takes
         * a lot of work.
         *  Note change to cgs units.
         */
        for (int j = 0; j < nsp; j++) {
            y[j] = 1.0;
            g.setState_TPY(T, BG.Pressure, y);
            kspec = t->thermalConductivity();
            Cond_Table(i, j) = kspec * 1.0E5;
            y[j] = 0.0;
        }
    }
    g.setState_TPX(BG.Temperature, BG.Pressure, BG.Xmol);
}
//======================================================================================================================
void
setAllBathSpeciesConditions(PhaseList* pl)
{
    for (size_t iphase = 0; iphase < pl->nPhases(); iphase++) {
        ThermoPhase* gThermoMainPhase =  &(pl->thermo(iphase));
        setBathSpeciesConditions(*gThermoMainPhase, pl, 0);
    }
}
//======================================================================================================================
void
printAllBathSpeciesConditions(PhaseList* pl)
{
    for (size_t iphase = 0; iphase < pl->nPhases(); iphase++) {
        ThermoPhase* gThermoMainPhase = &(pl->thermo(iphase));
        printBathSpeciesConditions(*gThermoMainPhase, pl, 1);
    }
}
//======================================================================================================================
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/*
 *  Set the bath gas for the one ThermoPhase. Note, we hardcode
 *  what we want here. An input deck might be the more appropriate
 *  place to set these properties.
 */
void
setBathSpeciesConditions(ThermoPhase& g, PhaseList* pl, int printLvl)
{

    //int nsp = g.nSpecies();

    int iph = pl->globalPhaseIndex(&g);
    if (iph < 0) {
        throw CanteraError(" ", "phase not found");
    }

    BG.BathSpeciesName = g.speciesName(BG.BathSpeciesIDVec[iph]);
    if (DebugPrinting) {
        cout << "Bath Species set to " << BG.BathSpeciesName << endl;
    }

    BG.MajorGasName = "O2";
    /*
     * Upcast to a constituents object, since speciesIndex
     * is defined multiple times
     */
    BG.MajorGasID = g.speciesIndex(BG.MajorGasName);

    /*
     * We need to set the state here. Note, mass fractions
     * don't matter, since we are only querying the
     * standard state of the species.
     */
    g.setState_TPX(BG.Temperature, BG.Pressure, BG.XmolPLPhases[iph]);

    /*
     * Set the electric potential
     */
    g.setElectricPotential(BG.PotentialPLPhases[iph]);

    /*
     * Print out table summarizing bath gas conditions
     */
    if (printLvl) {
        printBathSpeciesConditions(g, pl, printLvl);
    }
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
/*
 *  prints the bath gas for the one thermophase.
 */
void
printBathSpeciesConditions(ThermoPhase& g, PhaseList* pl, int printLvl)
{

    int nsp = g.nSpecies();

    int iph = pl->globalPhaseIndex(&g);
    if (iph < 0) {
        throw CanteraError(" ", "phase not found");
    }

    BG.BathSpeciesName = g.speciesName(BG.BathSpeciesIDVec[iph]);
    if (DebugPrinting) {
        cout << "Bath Species set to " << BG.BathSpeciesName << endl;
    }

    //BG.MajorGasName = "O2";
    /*
     * Upcast to a constituents object, since speciesIndex
     * is defined multiple times
     */
    BG.MajorGasID = g.speciesIndex(BG.BathSpeciesName);



    /*
     * Print out table summarizing bath gas conditions
     */
    if (printLvl) {
        double* C = new double [nsp];
        double* act = new double [nsp];
        g.getConcentrations(C);
        g.getActivities(act);
        size_t dim = g.nDim();
        dnt(0);
        print_char('=', 100);
        cout << endl;
        dnt(0);
        cout << " Bath Composition for Phase id, \"" << g.id()
             << "\", with name, \"" << g.name() << "\", with eosType, \""
             << eosTypeString(g.eosType()) << "\"" << endl;
        dnt(0);
        cout << endl;
        dnt(1);
        cout << "Total pressure = " << BG.Pressure * 760. / OneAtm;
        cout << " torr " << endl;
        dnt(1);
        cout << "Temperature (where needed) = " << BG.Temperature
             << " Kelvin" << endl;
        double volts = g.electricPotential();
        dnt(1);
        cout << "Voltage (where needed) = " << volts
             << " Volts" << endl;
        dnt(1);
        cout << "Carrier Species (used in diff. calcs) = "
             << BG.BathSpeciesName << endl;
        dnt(4);
        double cFac = 1.0E-3;
        if (dim == 3) {
            cout << "Number       Name      Mole_fraction  Concentration (gmol/cm**3)   Activities";
        } else if (dim == 2) {
            cFac = 1.0E-1;
            cout << "Number       Name      Mole_fraction  Concentration (gmol/cm**2)   Activities";
        } else if (dim == 1) {
            cFac = 1.0E1;
            cout << "Number       Name      Mole_fraction  Concentration (gmol/cm**1)   Activities";
        }
        cout << endl;
        dnt(4);
        cout << "-----------------------------------------------------------------------------";
        cout << endl;
        string spN;
        for (int k = 0; k < (int) g.nSpecies(); k++) {
            dnt(4);
            pr_if(k+1, 5);
            spN = g.speciesName(k);
            pr_sf(spN, 16);
            pr_df(BG.XmolPLPhases[iph][k], 16, 4);
            pr_de(C[k] * cFac, 16, 3);
            pr_de(act[k], 16, 3);
            cout << endl;
        }
        dnt(4);
        cout << "------------------------------------------------------------------------------";
        cout << endl;
        delete [] C;
        delete [] act;
    }
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

bool setupGasTransport(XML_Node* xmlPhase, ThermoPhase& g)

/*************************************************************
 *
 * setupGasTransport
 *
 *    This routine sets up the transport for the entire
 *    application
 *************************************************************/
{
    bool retn = true;
    double T = 300.;
    double P = 1.0E5;
    int numSpecies = g.nSpecies();
    double* y = new double[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        y[i] = 0.0;
    }
    y[0] = 1.0;
    g.setState_TPY(T, P, y);
    try {
        GTran = processExpandedTransport(xmlPhase, &g);
    } catch (CanteraError& ctex) {
        cout << "setupGasTransport ERROR!" << endl;
        showErrors(cout);
        retn = false;
    }
    delete [] y;
    if (!GTran) {
        retn = false;
    }
    return retn;
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

void print_Mixture_members(ThermoPhase& g)

/*************************************************************
 *
 * print_Mixture_members():
 *
 *    This routine dumps out information contained in the
 *    ThermoPhase class itself.
 *
 *************************************************************/
{
    print_char('-', 64);
    cout << endl;
    cout << " Dump of Mixture properties: " << endl;
    int numSpecies = g.nSpecies();
    double T = 300.;
    double P = 1.0E5;
    double* y = new double[numSpecies];
    for (int i = 0; i < numSpecies; i++) {
        y[i] = 0.0;
    }
    y[0] = 1.0;
    g.setState_TPY(T, P, y);
    if (GTran) {
        double visc = GTran->viscosity();
        cout << "\t viscosity = " << visc << endl;

        for (int i = 0; i < 10; i++) {
            T = 300. + 100. * i;
            g.setState_TPY(T, P, y);
            visc = GTran->viscosity();
            cout << "\t viscosity(" << T <<" K) = " << visc << endl;
        }
        print_char('-', 64);
        cout << endl;
    }
    delete []y;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

void printUsage()
{
    cout << "usage: cttables [-h] [-d DebugLvl] [cantera.xml] < command_file.txt "
         <<  endl;
    cout << "                   DebugLvl is the level of debug printing\n";
    cout << "      " << endl;
    cout << "         The following command prints out the command file text:\n\n";
    cout << "       cttables -help_cmdfile  <command_file.txt"
         << endl;
}
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/

int main(int argc, char** argv)
{
    FILE* inputFP = stdin;
    string xmlfile;
    int i;
    bool haveSpeciesTransportProps = false;
    bool skipTransport = false;
    bool printInputFormat = false;
    clockWC tt;
    string commandFile = "cttables.inp";
    ZZCantera::Kinetics* gKinetics = 0;
    // look for command-line options
    if (argc > 1) {
        string tok;
        for (int j = 1; j < argc; j++) {
            tok = string(argv[j]);
            if (tok[0] == '-') {
                int nopt = static_cast<int>(tok.size());
                for (int n = 1; n < nopt; n++) {
                    if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
                        printInputFormat = true;
                        break;
                    } else if (tok[n] == 'h') {
                        printUsage();
                        exit(1);
                    } else if (tok[n] == 'd') {
                        int lvl = 0;
                        if (j < (argc - 1)) {
                            string tokla = string(argv[j+1]);
                            if (strlen(tokla.c_str()) > 0) {
                                lvl = atoi(tokla.c_str());
                                n = nopt - 1;
                                j += 1;

                                DebugPrinting = lvl;

                            }
                        }
                    } else {
                        printUsage();
                        exit(1);
                    }
                }
            } else if (commandFile == "cttables.inp") {
                commandFile = tok;
            } else {
                printUsage();
                exit(1);
            }
        }
    }

    if (commandFile != "") {
        inputFP = fopen(commandFile.c_str(), "r");
        if (!inputFP) {
            printf("Can't open file, %s, bailing\n", commandFile.c_str());
            return -1;
        }
    }

    /*
     * General Catch block to trap Cantera Errors and print them
     */
    try {
        ReactingVolDomain* rVolDomain = new ReactingVolDomain();
        ThermoPhase* gThermoMainPhase = 0;
        XML_Node* xmlPhase = 0;


        /*
         * Print out the copywrite notice.
         */
        printf("\n");
        printf("Copywrite (2005) Sandia Corporation. Under the terms of \n");
        printf("Contract DE-Ac04-94AL85000, there is a non-exclusive license\n");
        printf("for use of this work by or on behalf of the U.S. Government.\n");
        printf("Export of this program may require a license from the United\n");
        printf("States Government.\n");

        BEInput::BlockEntry* cf = new BEInput::BlockEntry("command_file");
        PhaseList* pl = new PhaseList();
        /**
         * Setup and process the input deck from standard input.
         * -> note we attempt to make this work with gThermo = 0
         *    for documentation
         */
        setup_input_pass1(cf, gKinetics, gThermoMainPhase);

        if (printInputFormat) {
            printf("\n");
            printf("Command Line Format for Pass 1:\n\n");
            cf->print_usage(1);
        }

        bool ok = process_input(cf, commandFile, gKinetics, gThermoMainPhase, pl);
        if (!ok) {
            return -1;
        }

        /**
         * Setup and process the input deck for second time.
         * -> Might have to print out the input and quit as well.
         */
        setup_input_pass2(cf, gKinetics, gThermoMainPhase);

        if (printInputFormat) {
            printf("\n");
            printf("Command Line Format for Pass 2:\n\n");
            cf->print_usage(1);
        }

        ok = process_input(cf, commandFile, gKinetics, gThermoMainPhase, pl);
        if (!ok) {
            return -1;
        }

        int ifiles = 0;
        for (; IOO.CanteraFileNames[ifiles] != 0; ifiles++) {
        }
        if (ifiles != IOO.NumberCanteraFiles) {
            printf("Number of requested files differ\n");
            exit(-1);
        }


        /**
         * Read in all of the phase specifications from the cantera
         * input files into the PhaseList structure.
         */

        for (i = 0; i < IOO.NumberCanteraFiles; i++) {
            string jjj = IOO.CanteraFileNames[i];
            //if (!strcmp(PO.CanteraFileNames[i], "gas.xml")) {
            //printf("we are here\n");
            //}
            importAllCTMLIntoPhaseList(pl, jjj);
        }
        /*
         * Setup internally for next pass through the input file.
         */
        IOO.InitForInput(pl);

        /**
         * Setup and process the input deck for second time
         * -> Might have to print out the input and quit as well.
         */
        setup_input_pass3(cf, gKinetics, gThermoMainPhase, pl);

        if (printInputFormat) {
            printf("\n");
            printf("Command Line Format for Pass 3:\n\n");
            cf->print_usage(1);
            exit(0);
        }

        /*
         * Process the first pass of the input file ->
         *   We are just after the information needed to initialize
         *   the Cantera structures and size the problem
         */
        //cf->print_usage(0);
        ok = process_input(cf, commandFile, gKinetics, gThermoMainPhase, pl);
        //cf->print_usage(0);

        /*
         * OK, Find the first kinetics object
         */
        int numSurPhases = pl->nSurPhases();
        int numVolPhases = pl->nVolPhases();
        int isfound = -1;
        size_t ivfound = npos;
        for (i = 0; i < numSurPhases; i++) {
            if (pl->surPhaseHasKinetics(i)) {
                isfound = i;
                break;
            }
        }
        if (isfound == -1) {
            for (i = 0; i < numVolPhases; i++) {
                if (pl->volPhaseHasKinetics(i)) {
                    ivfound = i;
                    break;
                }
            }
        }


        /*
         * Import thermo and kinetics for multiple phases
         */
        if (isfound >= 0) {
            ok = rVolDomain->importSurKinFromPL(pl, isfound);
        } else {
            ok = rVolDomain->importVolKinFromPL(pl, ivfound);
        }
        if (!ok) {
            throw CanteraError("cttables main:", "rVolDomain returned an error");
        }


        gThermoMainPhase = rVolDomain->tpList[0];
        gKinetics = rVolDomain->m_kinetics;

        for (size_t iphase = 0; iphase < rVolDomain->m_NumPLPhases; iphase++) {
            gThermoMainPhase = rVolDomain->tpList[iphase];
            //
            // Change the default behavior of phases to print out INF's and not to change the numbers.
            //
            // Eventually we'll make a keyline command for this.
            //
            gThermoMainPhase->realNumberRangeBehavior_ = DONOTHING_CTRB;
            // gThermoMainPhase->realNumberRangeBehavior_ = CHANGE_OVERFLOW_CTRB;
            // gThermoMainPhase->realNumberRangeBehavior_ = THROWON_OVERFLOW_CTRB;
            size_t iph = pl->globalPhaseIndex(gThermoMainPhase);

            size_t nSpecies = gThermoMainPhase->nSpecies();
            std::string phaseBath = "Bath Specification for Phase ";
            std::string phaseNm = gThermoMainPhase->name();
            phaseBath += phaseNm;

            /*
             * Set up default bath gas conditions
             */
            mdp_safe_free((void**) &BG.Xmol);
            BG.Xmol = mdp_alloc_dbl_1(nSpecies, 0.0);

            mdp_safe_free((void**) &BG.Molalities);
            BG.Molalities = mdp_alloc_dbl_1(nSpecies,0.0);

            gThermoMainPhase->getMoleFractions(BG.Xmol);

            bool molVecSpecified = false;
            BlockEntry* pblock = cf->searchBlockEntry(phaseBath.c_str());
            if (pblock) {
                BlockEntry* pbsmf = pblock->searchBlockEntry("Bath Species Mole Fraction");
                if (pbsmf) {
                    if (pbsmf->get_NumTimesProcessed() > 0) {
                        molVecSpecified = true;
                    }
                }
            }

            if (molVecSpecified) {
                for (size_t k = 0; k < nSpecies; k++) {
                    BG.Xmol[k] = BG.XmolPLPhases[iph][k];
                }
            } else {
                for (size_t k = 0; k < nSpecies; k++) {
                    BG.XmolPLPhases[iph][k] =  BG.Xmol[k];
                }
            }

            bool molalVecSpecified = false;
            if (pblock) {
                BlockEntry* pbsmm =
                    pblock->searchBlockEntry("Bath Species Molalities");
                if (pbsmm) {
                    if (pbsmm->get_NumTimesProcessed() > 0) {
                        molalVecSpecified = true;
                    }
                }
            }
            if (molalVecSpecified) {
                MolalityVPSSTP* m_ptr = dynamic_cast<MolalityVPSSTP*>(gThermoMainPhase);
                if (m_ptr == 0) {
                    printf("Dynamic cast failed for some reason\n");
                    exit(-1);
                }
                m_ptr->setState_TPM(BG.Temperature, BG.Pressure, BG.MolalitiesPLPhases[iph]);
                m_ptr->getMoleFractions(BG.XmolPLPhases[iph]);

            }


            setBathSpeciesConditions(*gThermoMainPhase, pl, 0);
        }

        delete cf;
        cf = 0;
        /*
         *  Formulate a vector of temperatures to be used in
         *  the species and reaction tables
         */
        TT_ptr = new TemperatureTable(IOO.m_TTnpts, true, IOO.m_TTTlow,
                                      IOO.m_TTDeltaT, IOO.NumAddedTemps, IOO.AddedTemperatures);

        VV_ptr = new VoltageTable(IOO.m_VVnpts,  IOO.m_VVVlow, IOO.m_VVDeltaV,
                                  IOO.VVincZero,  IOO.NumAddedVoltages, IOO.AddedVoltages);

        for (size_t iphaseRVD = 0; iphaseRVD < rVolDomain->m_NumPLPhases; iphaseRVD++) {
            gThermoMainPhase = rVolDomain->tpList[iphaseRVD];

            size_t nSpecies = gThermoMainPhase->nSpecies();

            /*
             *  Import and set up the transport properties
             */
            if (IOO.SkipTransport) {
                skipTransport = true;
            }
            if (! skipTransport) {
                xmlPhase = & rVolDomain->thermo(0).xml();
                haveSpeciesTransportProps = setupGasTransport(xmlPhase, *gThermoMainPhase);
                if (!haveSpeciesTransportProps) {
                    skipTransport = true;
                }
            }

            /*
             * Print out the mixture members
             */
            if (DebugPrinting) {
                print_Mixture_members(*gThermoMainPhase);
            }
            /*
             * Set the bath gas conditions for the rest of the problem
             * Note, all results will be based on these conditions. Even
             * generated tables will be based upon these conditions,
             * except for the varying parameter that constitutes the
             * variation in the table.
             */
            setBathSpeciesConditions(*gThermoMainPhase, pl, 1);


            /*
             *  Go get all of the thermo quantities for the species tables.
             */
            DenseMatrix Cp_Table(TT_ptr->size(), nSpecies);
            DenseMatrix Hrel_Table(TT_ptr->size(), nSpecies);
            DenseMatrix Urel_Table(TT_ptr->size(), nSpecies);
            DenseMatrix Grel_Table(TT_ptr->size(), nSpecies);
            DenseMatrix S_Table(TT_ptr->size(), nSpecies);
            std::vector<double> H298(nSpecies);
            std::vector<double> U298(nSpecies);
            if (nSpecies > S298.size()) {
                S298.resize(nSpecies);
            }
            getThermoTables(*TT_ptr, Cp_Table, Hrel_Table, Grel_Table, S_Table, H298, *gThermoMainPhase, Urel_Table, U298);

            /*
             * Go get all of the transport quantities for the species
             * tables
             */
            DenseMatrix Visc_Table(TT_ptr->size(), nSpecies);
            DenseMatrix Cond_Table(TT_ptr->size(), nSpecies);
            DenseMatrix Diff_Table(TT_ptr->size(), nSpecies);
            if (! skipTransport) {
                getGenericTransportTables(*TT_ptr, *gThermoMainPhase, Visc_Table, Cond_Table, Diff_Table);
            }

            /*
             * Look up which routine to print out the species.
             *       cIdealGas
             */
            /*
             * Loop over the species printing out a table
             */
            if (DebugPrinting) {
                cout << "EOS type = " << gThermoMainPhase->eosType() << endl;
            }
            if (gThermoMainPhase->eosType() == cIdealGas) {
                if (DebugPrinting) {
                    cout << "PrintOut of Phase, \"" << gThermoMainPhase->id()
                         << "\", using the Ideal Gas Thermodynamics Functions:" << std::endl;
                }
                for (size_t k = 0; k < nSpecies; k++) {
                    printIdealGasSpeciesTable(*gThermoMainPhase, k, *TT_ptr, Cp_Table,
                                              Hrel_Table, Grel_Table, S_Table,haveSpeciesTransportProps,
                                              Visc_Table, Cond_Table, Diff_Table, H298, Urel_Table, U298);
                }
            } else {
                if (DebugPrinting) {
                    cout << "PrintOut of Phase, \"" << gThermoMainPhase->id()
                         << "\", using Generic Thermodynamics Functions:" << endl;
                }
                /*
                 * Loop over the species printing out a table
                 */
                for (size_t k = 0; k < nSpecies; k++) {
                    printThermoPhaseSpeciesTable(gThermoMainPhase, k, *TT_ptr, Cp_Table,
                                                 Hrel_Table, Grel_Table, S_Table, haveSpeciesTransportProps,
                                                 Visc_Table, Cond_Table, Diff_Table, H298, Urel_Table, U298);
                }
            }
        }
        if (rVolDomain->m_kinetics) {
            doKineticsTablesHomog(pl, gKinetics, *TT_ptr);
        }
        if (rVolDomain->m_InterfaceKinetics) {
            doKineticsTablesHetero(pl, rVolDomain->m_InterfaceKinetics, *TT_ptr);
        }

        fclose(inputFP);
        delete rVolDomain;
        rVolDomain = 0;
        delete pl;
        pl = 0;
        delete TT_ptr;
        TT_ptr = 0;
        appdelete();

    } catch (CanteraError) {
        showErrors(cout);
        return 0;
    }
    return 0;
}
//==================================================================================================================================

