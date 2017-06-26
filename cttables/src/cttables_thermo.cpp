/**
 *  @file cttables_thermo.cpp
 *
 *
 */

/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "cttables.h"

#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/base/PrintCtrl.h"

#ifdef useZuzaxNamespace
using namespace Zuzax;
#else
using namespace Cantera;
#endif

#include <cstdio>

using std::endl;
using std::cout;

//===================================================================================================================================
// Gather the entropy of the elements of a species at 298 K. This is useful for going back and forth from the
// gibbs free energy of formation and the absolute gibbs free energy in NIST format.
//
double entropyElem298(ZZCantera::ThermoPhase* g_ptr, size_t k)
{
    double se;
    double stotal = 0.0;
    for (size_t m = 0; m < g_ptr->nElements(); m++) {
        double na = g_ptr->nAtoms(k, m);
        if (na != 0.0) {
            se = g_ptr->entropyElement298(m, true);
            if (se == ENTROPY298_UNKNOWN) {
                return ENTROPY298_UNKNOWN;
            }
            stotal += se * na;
        }
    }
    return stotal;
}
//===================================================================================================================================
void printThermoCoeffSpecies(ThermoPhase* g_ptr, int k)
{

    int j, nint, nzones = 1, i;
    std::string sName = g_ptr->speciesName(k);
    double* cptr;
    int index;
    size_t kindex;
    int type;
    PDSS* pdss_ptr = nullptr;
    // Get the species thermo manager
    SpeciesThermo& sp = g_ptr->speciesThermo();
    VPStandardStateTP* vpss = nullptr;
    PDSS_enumType ptype;
    double DH0_tr_pr;
    double S0_tr_pr;
    double Mu0_tr_pr;
    double DG0_tr_pr;
    double dg_consistent, as, bs, cs;

    // Get the report Type
    int rt = sp.reportType(k);
    double c[200];

    // Get the coefficients of the parameterization
    int rt2;
    double minTemp, maxTemp, refPressure;
    sp.reportParams(k, rt2, c, minTemp, maxTemp, refPressure);
    if (rt != rt2) {
        printf("report types don't match ERROR\n");
        exit(-1);
    }
    double stotal = entropyElem298(g_ptr, k);

    // switch on the report type, printing out a formatted representation of the polynomials
    switch (rt) {
    case NASA2:
        dnt(1); printf("NASA Polynomial format: 2 zones\n");
        dnt(2); printf("Low  temperature polynomials: %g < T < %g: \n", minTemp, c[0]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,\n", c[3], c[4], c[5], c[6]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E\n", c[7], c[1], c[2]);
        dnt(2); printf("High temperature polynomials: %g < T < %g: \n", c[0], maxTemp);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,\n", c[10], c[11], c[12],  c[13]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E\n", c[14], c[8], c[9]);
        break;

    case SHOMATE2:
        dnt(1);
        printf("SHOMATE Polynomial format: 2 zones\n");
        dnt(2); printf("Low  temperature polynomials: %g < T < %g: \n", minTemp, c[0]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,n", c[1], c[2], c[3], c[4]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E\n", c[5],  c[6], c[7]);
        dnt(2); printf("High temperature polynomials: %g < T < %g: \n", c[0], maxTemp);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,n", c[8], c[9], c[10], c[11]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E\n", c[12], c[13], c[14]);
        break;

    case NASA1:
        dnt(1); printf("NASA Polynomial format: 1 zones\n");
        dnt(2); printf("temperature polynomials: %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,n", c[1], c[2],  c[3],  c[4]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E\n", c[5], c[6], c[7]);
        break;

    case SHOMATE1:
        dnt(1);
        printf("SHOMATE Polynomial format: 1 zones\n");
        dnt(2);
        printf("temperature polynomials: %g < T < %g: \n", minTemp, maxTemp);
        dnt(2);
        printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,/n", c[0], c[1], c[2], c[3]);
        dnt(2);
        printf("% 23.15E,% 23.15E,% 23.15E\n", c[4], c[5], c[6]);
        break;

    case CONSTANT_CP:
    case SIMPLE:
        dnt(1); printf("CONSTANT_CP format:\n");
        dnt(2); printf(" Valid Range: %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf(" at T = %5g K, H  = % 23.15E J kmol-1\n", c[0], c[1]);
        dnt(2); printf("                 S  = % 23.15E J kmol-1 K-1\n", c[2]);
        dnt(2); printf("                 Cp = % 23.15E J kmol-1 K-1\n", c[3]);
        break;

    case MU0_INTERP:
        dnt(1); printf("MU0_POLY format:\n");
        nint = (int) c[0];
        dnt(2); printf(" Valid Range: %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf(" T        mu_0\n");
        j = 2;
        for (i =0 ; i < nint; i++) {
            dnt(2); printf(" %g   %17.11g\n", c[j], c[j+1]);
            j += 2;
        }
        break;

    case NASA9:
        nzones = 1;
        dnt(1); printf("NASA9 Polynomial format: %d zones\n", nzones);
        cptr = c;
        dnt(2); printf("%2d:  Temperature polynomials: %g < T < %g: \n", 0, minTemp, maxTemp);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,\n", cptr[3], cptr[4],  cptr[5],  cptr[6]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,\n", cptr[7], cptr[8], cptr[9], cptr[10]);
        dnt(2); printf("% 23.15E\n", cptr[11]);
        break;

    case NASA9MULTITEMP:
        nzones = (int) c[0];
        dnt(1);
        printf("NASA9 Polynomial format: %d zones\n", nzones);
        index = 1;
        for (i = 0; i < nzones; i++) {
            minTemp = c[index];
            maxTemp = c[index+1];
            cptr = c+index+2;
            dnt(2); printf("%2d:  Temperature polynomials: %g < T < %g: \n", i, minTemp, maxTemp);
            dnt(2); printf("% 23.15E % 23.15E % 23.15E % 23.15E\n", cptr[0], cptr[1],  cptr[2],  cptr[3]);
            dnt(2); printf("% 23.15E % 23.15E % 23.15E % 23.15E\n", cptr[4], cptr[5], cptr[6], cptr[7]);
            dnt(2); printf("% 23.15E\n", cptr[8]);
            index += 11;
        }
        break;

    case MINEQ3:
        nzones = 1;
        DH0_tr_pr = c[7];
        S0_tr_pr = c[8];
        Mu0_tr_pr = c[9];
        DG0_tr_pr = c[10];
        as = c[0] / 4.184;
        bs = c[1] / (4.184 * 1000.);
        cs = c[4] * 1.0E6 / 4.184;

        if (stotal == ENTROPY298_UNKNOWN) {
            dg_consistent = ENTROPY298_UNKNOWN;
        } else {
            dg_consistent =  Mu0_tr_pr + 298.15 * (stotal);
        }
        dnt(1); printf("MinEQ3 format:  (a varient of Shomate1 format) \n");
        dnt(2); printf("temperature polynomials (Shomate Form): %g < T < %g: \n", minTemp, maxTemp);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E,% 23.15E,\n", c[0],  c[1], c[2], c[3]);
        dnt(2); printf("% 23.15E,% 23.15E,% 23.15E\n",  c[4],  c[5], c[6]);
        dnt(2); printf("  Delta G0_Tr_Pr = %16.9E cal/gmol\n", DG0_tr_pr / (4.184 * 1.0E3));
        dnt(2); printf("                 = %16.6g kJ /gmol\n", DG0_tr_pr/1.0E6);
        dnt(2); printf("  Delta H0_Tr_Pr = %16.9E cal/gmol\n", DH0_tr_pr / (4.184 * 1.0E3));
        dnt(2); printf("        S0_Tr_Pr = %16.9g cal/gmol/K\n", S0_tr_pr / (4.184 * 1.0E3));
        dnt(2); printf("                 = %16.6g  J /gmol/K\n", S0_tr_pr / 1.0E3);
        dnt(2); printf("       mu0_Tr_Pr = %16.6g kJ /gmol\n",   Mu0_tr_pr / 1.0E6);
        if (dg_consistent == ENTROPY298_UNKNOWN) {
            dnt(2); printf(" Delta G0_consis = [UNAVAILABE BECAUSE ENTROPY298 NOT INPUT]\n");
        } else {
            dnt(2); printf(" Delta G0_consis = %16.6g kJ /gmol\n", dg_consistent / 1.0E6);
        }
        dnt(2); printf("               a = %16.9g cal/gkmol/K\n", as);
        dnt(2); printf("               b = %16.9g cal/kmol/K2\n", bs);
        dnt(2); printf("               c = %16.9g cal-K/gmol\n", cs);
        break;

    case PDSS_TYPE:
        dnt(1); printf("Presure Dependent Standard State form\n");
        vpss = dynamic_cast<VPStandardStateTP*>(g_ptr);
        pdss_ptr = vpss->providePDSS(k);
        ptype = pdss_ptr->reportPDSSType();
        double minTemp, maxTemp, refPressure;
        pdss_ptr->reportParams(kindex, type, &c[0], minTemp, maxTemp, refPressure);
        dnt(2); printf("PDSS Species index    = %d\n", (int) kindex);
        dnt(2); printf("PDSS min Temperature  = %g\n", minTemp);
        dnt(2); printf("PDSS max Temperature  = %g\n", maxTemp);
        dnt(2); printf("PDSS reference pres   = %g\n", refPressure);

        switch (ptype) {
        case cPDSS_IDEALGAS:
            dnt(2); printf("Ideal Gas PDSS thermo type\n");
            break;
        case  cPDSS_CONSTVOL:
            dnt(2); printf("ConstVol PDSS thermo type\n");
            break;
        case cPDSS_MOLAL_CONSTVOL:
            dnt(2); printf("Molal ConstVol PDSS thermo type\n");
            break;
        case cPDSS_WATER:
            dnt(2); printf("Water PDSS thermo type\n");
            break;
        case cPDSS_MOLAL_HKFT:
            dnt(2); printf("PDSS_MOLAL_HKFT PDSS thermo type\n");
            dnt(3); printf("deltaG_formation_tr_pr = %g\n", c[0]);
            dnt(3); printf("deltaH_formation_tr_pr = %g\n", c[1]);
            dnt(3); printf("Mu0_tr_pr              = %g\n", c[2]);
            dnt(3); printf("Entrop_tr_pr           = %g\n", c[3]);
            dnt(3); printf("a1                     = %g\n", c[4]);
            dnt(3); printf("a2                     = %g\n", c[5]);
            dnt(3); printf("a3                     = %g\n", c[6]);
            dnt(3); printf("a4                     = %g\n", c[7]);
            dnt(3); printf("c1                     = %g\n", c[8]);
            dnt(3); printf("c2                     = %g\n", c[9]);
            dnt(3); printf("omega_pr_tr            = %g\n", c[10]);

            break;
        case cPDSS_UNDEF:
            dnt(2); printf("Undefined PDSS thermo type\n");
            break;
        default:
            dnt(2); printf("unknown species PDSS thermo type %d\n", (int) ptype);
        }

        break;
    default:
        printf("unknown species reference thermo type %d\n", rt);
    }
}
//==================================================================================================================================
/*
 *  This routine adds the molar volume information onto a species table.
 *  This is not part of the reference state.  This is part of the description of the phase.
 */
void printVolSpecies(ThermoPhase* const g_ptr, size_t k)
{
    std::string sName = g_ptr->speciesName(k);
    int eos = g_ptr->eosType();
    //int dim = g_ptr->nDim();
    size_t nS = g_ptr->nSpecies();
    std::vector<double> ssVol(nS, 0.0);
    g_ptr->getStandardVolumes(DATA_PTR(ssVol));
    std::vector<double> pmVol(nS, 0.0);
    g_ptr->getPartialMolarVolumes(DATA_PTR(pmVol));
    double currTemp = g_ptr->temperature();
    double currPres = g_ptr->pressure();

    dnt(1);
    printf("Partial    Molar Volume at (T,P) of (%g,%g) = %g m%1d/kmol\n", currTemp, currPres, pmVol[k], 3);
    dnt(1);
    printf("StandState Molar Volume at (T,P) of (%g,%g) = %g m%1d/kmol\n", currTemp, currPres, ssVol[k], 3);

    /*
     * Loop over the known types of eos's. Print out what's available.
     */
    if (eos == cIdealSolidSolnPhase0  || eos == cIdealSolidSolnPhase1 || eos == cIdealSolidSolnPhase2) {

        IdealSolidSolnPhase* issp = static_cast<IdealSolidSolnPhase*>(g_ptr);
        double mV = issp->speciesMolarVolume(k);
        dnt(1);
        cout << "species Molar Volume = " << mV << " m**3 / kmol\n";

        switch (eos) {
        case cIdealSolidSolnPhase0:
            dnt(1);
            cout << "standard Concentration value =  unity\n";
            break;
        case cIdealSolidSolnPhase1:
            dnt(1);
            cout << "standard Concentration value =  1/(molar volume) units of kmol/m3\n";
            break;
        case cIdealSolidSolnPhase2:
            dnt(1);
            cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3\n";
            break;
        }

        dnt(1);
        cout << "Pres dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1);
        cout << "                -> (H and G dependent): h(T) = u(T) + pv" << endl;

    } else if (eos == cStoichSubstance) {

        //StoichSubstance* issp = static_cast<StoichSubstance*>(g_ptr);
        double dens = g_ptr->density();
        double mw = g_ptr->molecularWeight(0);
        double mV = mw / (dens);
        dnt(1);
        cout << "StoichSubstance constant density = " << dens << " kg m-3 \n";
        dnt(1);
        cout << "species Molar Volume = " << mV << " m**3 / kmol <- constant\n";
        dnt(1);
        cout << "Press dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1);
        cout << "                 -> (H and G dependent): h(T) = u(T) + pv\n";
        dnt(1);
        cout << "standard Concentration value =  unity" << endl;

    } else if (eos == cIdealGas) {

        dnt(1);
        cout << "Volume Properties are from ideal gas law\n";
        dnt(1);
        cout << "Standard concentration value =  P / (R T) with units of kmol/m3" << endl;

    } else if (eos == cDebyeHuckel0  || eos == cDebyeHuckel1  || eos == cDebyeHuckel2) {

        double mV = pmVol[k];
        dnt(1);
        cout << "species partial Molar Volume = " << mV << " m**3 / kmol\n";
        switch (eos) {
        case cDebyeHuckel0:
            dnt(1);
            cout << "standard Concentration value =  unity\n";
            break;
        case cDebyeHuckel1:
            dnt(1);
            cout << "standard Concentration value =  1/(molar volume) units of kmol/m3\n";
            break;
        case cDebyeHuckel2:
            dnt(1);
            cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3\n";
            break;
        }
        dnt(1);
        cout << "Pres dependence -> (U, S, V, and Cp are independent)\n";
        dnt(1);
        cout << "                -> (H and G dependent): h(T) = u(T) + pv" << endl;

    } else if (eos == cHMWSoln0  || eos == cHMWSoln1  || eos == cHMWSoln2) {

        double mV = pmVol[k];
        dnt(1);
        cout << "species partial Molar Volume = " << mV << " m**3 / kmol\n";
        switch (eos) {
        case cHMWSoln0:
            dnt(1);
            cout << "standard Concentration value =  unity" << endl;
            break;
        case  cHMWSoln1 :
            dnt(1);
            cout << "standard Concentration value =  1/(molar volume) units of kmol/m3" << endl;
            break;
        case cHMWSoln2:
            dnt(1);
            cout << "standard Concentration value =  1/(solvent volume) units of kmol/m3" << endl;
            break;
        }

    } else if (eos == cSurf) {

        dnt(1);
        cout << "Volume Properties are from fixed surface site equation\n";
        SurfPhase* isurf = static_cast<SurfPhase*>(g_ptr);
        double sited = isurf->siteDensity();
        sited  *= Avogadro / 1.0E4;
        dnt(1);
        cout << "Surface Site Density = " << sited << " # per cm2 " << endl;

    } else {

        dnt(1);
        cout << "Molar Volume formulation for eos type " << eos << " is not known to cttables" << endl;

    }
}
//=================================================================================================================================
/**
 *  This routine will print out a table of information about a species in an ideal gas thermo phse. It explicitly
 *  assumes that a multitransport object has been created for  the phase, and it presumes a NASA polynomial form for the
 *  species thermodynamics.
 */
void printIdealGasSpeciesTable(Zuzax::ThermoPhase& g, int k, TemperatureTable& TT, Zuzax::DenseMatrix& Cp_Table,
                               Zuzax::DenseMatrix& Hrel_Table, Zuzax::DenseMatrix& Grel_Table, Zuzax::DenseMatrix& S_Table,
                               bool haveSpeciesTransportProps, Zuzax::DenseMatrix& Visc_Table, Zuzax::DenseMatrix& Cond_Table,
                               Zuzax::DenseMatrix& Diff_Table, std::vector<double> H298, Zuzax::DenseMatrix& Urel_Table,
                               std::vector<double> U298)
{
    PrintCtrl pc;
    int tModel = 0;
    /*
     *  Get the species data object from the Mixture object
     *  this is defined in the constituents.h file, and is
     *  inherited by Mixture through BaseMix
     */
    std::string sName = g.speciesName(k);
    /*
     *  Dump out all of the information about the species
     */
    printf("\n");
    print_char('=', 120);
    printf("\n");
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
    printf("Electronic Charge = %2f\n", g.charge(k));
    double mw = g.molecularWeight(k);
    if (mw == Tiny) {
        dnt(1);
        printf("Molecular Weight = 0 (set to %g for massF conversion) gm/mol\n", Tiny);
    } else {
        dnt(1);
        printf("Molecular Weight = %g gm/mol\n", mw);
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
            const MixTransport* mt = dynamic_cast<MixTransport*>(GTran);
            double wellDepth = mt->m_eps[k] / Boltzmann;
            dnt(1); printf("L-J Potential Well Depth = %g K\n", wellDepth);
            double diameter = mt->m_sigma[k] * 1.0E10;
            dnt(1); printf("L-J collision diameter = %g Angstroms\n", diameter);
            double dipoleMoment = 1.0E25 / SqrtTen * mt->m_dipole(k,k);
            dnt(1); printf("Dipole Moment = %g Debeye\n", dipoleMoment);
            double polarizability = mt->m_alpha[k] * 1.0E30;
            dnt(1); printf("Polarizability = %2f Angstroms**3\n", polarizability);
            double rotRelaxNumber = mt->m_zrot[k];
            dnt(1); printf("Rotational Collision Number at 298 K = %g\n", rotRelaxNumber);
            double crot = mt->m_crot[k];
            if (crot == 0.0) {
                dnt(1); printf("This molecule is monatomic\n");
            } else if (crot == 1.0) {
                dnt(1); printf("This molecule is linear\n");
            } else if (crot == 1.5) {
                dnt(1); printf("This molecule is nonlinear\n");
            } else {
                throw ZuzaxError("printTable()", "Unknown crot = " + ZZCantera::fp2str(crot));
            }
        }
    }
    /*
     * Print out the Heat of Formation at 298.15 K
     */
    dnt(1); printf("Heat of formation (298.15K) = %.4f %8s\n", H298[k], UIO.sGibbs.c_str());
    double mu_o = H298[k];

    double deltaGf = 1.0E6 * H298[k] - 298.15 * g_S298_refThermo[k] + 298.15 * elementEntropyTotal;
    deltaGf /= 1.0E6;
    dnt(1); printf("DeltaGf (298.15K) = %.4f %8s\n", deltaGf, UIO.sGibbs.c_str());
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
        std::string tmp = g.speciesName(BG.BathSpeciesID);
        pc.pr_str_MinMax(tmp, 9, 9);
        //pr_sf_lj(tmp, 9, 1);

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
    print_char('-', tableWidth);
    printf("\n");
}
//=====================================================================================================================================
