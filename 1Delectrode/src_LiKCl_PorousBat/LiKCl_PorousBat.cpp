/**
 * @file tddiff.cpp
 *
 */

/*
 *  $Id: LiKCl_PorousBat.cpp 593 2013-05-13 21:25:47Z hkmoffa $
 */


#include "Ifpack.h"
#include "AztecOO.h"
#include "m1d_defs.h"

#include "LiKCl_PorousBat.h"
#include "m1d_SolNonlinear_CurrentSolve.h"
#include "m1d_DomainLayout_LiKCl_PorousBat.h"
#include "m1d_ProblemStatementCell.h"
#include "BEulerInt_Battery.h"
#include "m1d_CanteraElectrodeGlobals.h"
#include "m1d_GlobalIndices.h"

#include "zuzax/base/Array.h"

//=====================================================================================
/*
 * HKM
 *
 *     Formulation of the structure of the problem.
 *
 *
 *    We will create a system for solving one dimensional problems based
 *    on using Trilinos and Cantera.
 *
 *    The matrix stencil for the calculations will be based on  a set of
 *    of domains that are connected with each other by tie-regions.
 *    The tie-regions may have a set of equations associated with
 *    the tie region. This set of equations will have dependencies
 *    on the neighboring nodes of the bounding domains.
 *
 *
 *
 *    |  Tie        |  Domain  | Tie      |  Domain  |  Tie     |
 *    |  Region     |  Zero    | Region   |  One     |  Region  |
 *    |  0          |          | 1        |          |  2       |
 *    |             |          |          |          |          |
 *    |  Ne^T_0     |  Ne^D_0  | Ne^T_1   |  Ne^D_1  |  Ne_T_2  |
 *
 *
 *
 *
 *   The system will be able to do solved using a MP formulation. However,
 *   we will not build into the formulation any kind of extensibility. All
 *   nodes will know everything about the problem. All nodes will have all
 *   data. MP will be achieved by breaking the nodes of the problem up
 *   into sections. The residual for each section will belong to one processor
 *   The matrix row for that section will belong to one processor.
 *   Matrix solves will occur via a MP method.
 *
 *
 *
 *  For our first step we will assume one domain, with no tie regions !
 */

using namespace std;
using namespace m1d;
using namespace beuler;

m1d::BatteryResidEval *PS_ptr = 0;

m1d::ProblemStatementCell PSinput;

int flagPrecipitation = -1;

void printUsage()
{
    cout << "usage: LiKCl_PorousBat [-h] [-help_cmdfile] [-d #] LiKCl_PorousBat.inp" << endl;
    cout << "    -h           help" << endl;
    cout << "    -d           #   : level of debug printing" << endl;
    cout << "  separator.inp  #   : command file" << endl;
    cout << endl;
}

//=====================================================================================

int main(int argc, char **argv)
{

    ofstream oftmp;
    m1d::stream0 w0;
#ifdef HAVE_MPI
    // Initialize MPI
    MPI_Init(&argc, &argv);
#endif
    try {
#ifdef HAVE_MPI
        // Create the One New communications object for global communications
        Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
        Epetra_SerialComm Comm;
#endif
        Comm_ptr = &Comm;
        /*
         * Process the command line arguments
         */
        int retn = 0;
        std::string commandFile = "LiKCl_PorousBat.inp";
        if (argc > 1) {
            std::string tok;
            for (int j = 1; j < argc; j++) {
                tok = string(argv[j]);
                if (tok[0] == '-') {
                    int nopt = static_cast<int>(tok.size());
                    for (int n = 1; n < nopt; n++) {
                        if (!strcmp(tok.c_str() + 1, "help_cmdfile")) {
                            m1d::PrintInputFormat = true;
                            break;
                        } else if (tok[n] == 'h') {
                            printUsage();
                            exit(1);
                        } else if (tok[n] == 'd') {
                            if (j < (argc - 1)) {
                                string tokla = string(argv[j + 1]);
                                if (strlen(tokla.c_str()) > 0) {
                                    n = nopt - 1;
                                    j += 1;
                                }
                            }
                        } else {
                            printUsage();
                            exit(1);
                        }
                    }
                } else if (commandFile == "" || commandFile == "LiKCl_PorousBat.inp") {
                    commandFile = tok;
                } else {
                    printUsage();
                    exit(1);
                }
            }
        }

        /*
         * Go get the problem description from the input file
         */
        PSCinput_ptr = &PSinput;
        retn = PSinput.parse_input_1(commandFile);
        if (retn == -1) {
            printf("exiting with error\n");
            exit(-1);
        }
        retn = PSinput.parse_input_2();
        if (retn == -1) {
            printf("exiting with error\n");
            exit(-1);
        }
        retn = PSinput.parse_input_3();
        if (retn == -1) {
            printf("exiting with error\n");
            exit(-1);
        }

        PSinput.readAnodeInputFile();
        PSinput.readCathodeInputFile();

        m1d::BatteryResidEval *ps = new BatteryResidEval(1.0E-13);
        PS_ptr = ps;

        /*
         *  Initialize the domain structure for the problem
         * Some hard coded things here until we get the parser
         * integrated.
         */

        //not used?
        /* map boundary conditions to their string label */
        /* since I haven't yet invented boundary conditions,
         * we'll just let them be ints
         */
        map<string, int> boundaryMap;
        boundaryMap["AnodeCollectorPlate"] = 0;
        boundaryMap["CathodeCollectorPlate"] = 1;

        /*
         *  In this step we assign all of the domains in the problem.
         *  We assign the number of nodes per domain.
         *  In setting up the domains, we assign the number of unknowns
         *  per node within the domain.
         */
        DomainLayout *dl = new DomainLayout_LiKCl_infPorousBat(&PSinput);

        /*
         *  Here, we assign the domainlayout to the problem solver
         *  object.
         */
        ps->specifyProblem(dl, &PSinput);

        /*
         *  Because we know the number of nodes and the number of
         *  unknowns per node, we can now lay out the global indexing
         *  for the problem.
         */
        ps->generateGlobalIndices();

        // (Now that MPI is initialized, we can access/use MPI_COMM_WORLD.)
        m1d::GlobalIndices &GI = *ps->GI_ptr_;
#ifdef HAVE_MPI
        MPI_Comm_size(MPI_COMM_WORLD, &(GI.NumProc));
        MPI_Comm_rank(MPI_COMM_WORLD, &(GI.MyProcID));
#else
        GI.NumProc = 1;
        GI.MyProcID = 0;
#endif

        /*
         *  Generate and fill up the local node vectors on this processor
         *    These indices have to do with accessing the local nodes on the processor
         *    At the end of this routine, the object LocalNodeIncides is fully formed
         *    on this processor.
         */
        ps->generateLocalIndices();

        /*
         * Malloc and start to prep the Domains that will be used for the efficient
         * calculation of the residual
         */
        generateDomain1D(ps);

        if (m1d::PrintInputFormat) {
            exit(0);
        }

        /*
         *  Set the reference temperature for all domains
         *  Do this after we have malloced the Domain1D structures
         */
        dl->set_TP_Reference(PSinput.TemperatureReference_, PSinput.PressureReference_);

        // bool copyMode = false;
        // ps->LRN_VBR_ptr_ = new m1d::LocalRowNodeVBRIndices(Comm_ptr, copyMode, ps->GI_ptr_, ps->LI_ptr_);

        // Declare a matrix and two vectors, to form a linear system: A * v = b.

        // The following will cause an error in the importer, because of an internal error in epetra
        //   -> const ElementSize is not handled correctly.
        // this will cause an error in the importer -> Epetra_Vector *v = new Epetra_Vector(*GI.GbEqnstoOwnedLcEqnsMap);
        Epetra_Vector *v = new Epetra_Vector(* ( (ps->LI_ptr_)->GbBlockNodeEqnstoLcBlockNodeEqnsColMap));
        Epetra_Vector *b = new Epetra_Vector(* ( (ps->LI_ptr_)->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap));
        Epetra_Vector *res = new Epetra_Vector(* ( (ps->LI_ptr_)->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap));
        Epetra_Vector *soln = new Epetra_Vector(* ( (ps->LI_ptr_)->GbBlockNodeEqnstoLcBlockNodeEqnsColMap));

        v->PutScalar(0.0);

        b->PutScalar(0.0);
        soln->PutScalar(0.0);

        ps->domain_prep();

        EpetraJac *jac = new EpetraJac(*ps);
        jac->allocateMatrix();
#ifdef DEBUG_MATRIX_STRUCTURE
        // this block always produces different results on every call. Therefore, it's not appropriate
        // for a test suite.
        print0_sync_start(true, w0, Comm);
        ostringstream ssSave;
        jac->queryMatrixStructure(ssSave);
        w0 << ssSave.str();
        w0 << endl;
        ssprint0(w0);
        print0_sync_end(true, w0, Comm);
#endif

        /*
         *  Set the initial conditions in the solution vector by calling the ProblemResidEval object
         */
        double t_init = 0.0;
        double delta_t = 1.0E-8;
        double delta_t_np1;
        ps->residSetupTmps();
        ps->initialConditions(false, soln, 0, t_init, delta_t, delta_t_np1);
        ps->residEval(res, false, soln, 0, 0.0, 0.0);
        //    We need to get rid of state information about time step in the electrode objects, because
        //    we will be manipulating the initial conditions in the Electrode objects.
        ps->revertToInitialGlobalTime();

#ifdef DEBUG_MATRIX_STRUCTURE
        print0_epMultiVector(*res, "Residual Value");
        print0_epMultiVector(*soln, "Solution Values");
#endif

        BEulerInt_Battery t1;

#ifdef DO_INIT_CALC
        /*
         * Section to turn on DAEINIT capability - experimental
         */
        char *resp_str = getenv("Bat_DODAEINIT");
        t1.m_doSpecialStartCalc = 1;
        if (resp_str) {
            if (resp_str[0] == 'n') {
                t1.m_doSpecialStartCalc = 0;
            }
        }
#else 
        t1.m_doSpecialStartCalc = 1;
#endif
        /*
         *  Initialize the time stepper by providing it with the ProblemResidEval object
         *  so that internal structures can be malloced.
         */
        t1.initializePRE(*ps);
        /*
         *  Set the nonlinear solver options used within the time stepper, and the
         *  initial DAE solve algorithm
         */
        int min_newt_its = 2;
        int matrixConditioning = 0;
        bool colScaling = 1;
        bool rowScaling = 1;
        t1.setNonLinOptions(min_newt_its, matrixConditioning, colScaling, rowScaling);

#ifdef DEBUG_HKM
        /*
        int printFlag = 10;
        t1.setPrintFlag(printFlag);

        SolNonlinear::s_print_NumJac = true;
        */
#endif

        t1.determineInitialConditions(PSinput.startTime_, delta_t);
        ps->revertToInitialGlobalTime();

        Epetra_Vector_Ghosted &solnInt = t1.solnVector();
        Epetra_Vector_Ghosted &solnDot = t1.solnDotVector();

        // Ok after the initialization of the solution vector, redo the calculation
        ps->residEval(res, true, &solnInt, &solnDot, 0.0, 0.0);
        ps->revertToInitialGlobalTime();
#ifdef DEBUG_HKM
        /*
        print0_epMultiVector(*res, "Residual Value");
        print0_epMultiVector(solnInt, "Solution Values");
        print0_epMultiVector(solnDot, "Solution Values");
        */
#endif

#ifdef DEBUG_MATRIX_STRUCTURE
        print0_epMultiVector(*res, "Residual Value");
        print0_epMultiVector(solnInt, "Solution Values");
        print0_epMultiVector(solnDot, "Solution Values");
#endif

        if (PSinput.rootFinderForConstantCurrent_ > 0) {
            SolGlobalNonlinear *sn = new SolNonlinear_CurrentSolve();
            t1.specifyNonLinearSolver(sn);
        }

        t1.setProblemType(BEULER_JAC_ANAL);
        // Epetra_Vector *av = new Epetra_Vector(*(t1.m_ownedMap));
        // av->PutScalar(PSinput.absTol_);
        const Epetra_Vector_Owned & abstol = ps->atolVector();
        t1.setTolerancesEpetra(PSinput.relTol_, abstol);
        //delete av;
        //av = 0;
        /***************************************************************************/
        // debugging section to do a matrix formation and a solve
#ifdef DEBUG_MATRIX_STRUCTURE
        ps->residEval(res, false, &solnInt, 0, 0.0, 0.0);
        jac->matrixEval(false, &solnInt, 0, 0.0, 0.0);

        print0_epVbrMatrix(*(jac->A_));

        for (int i = 0; i < ps->LI_ptr_->NumLcOwnedEqns; i++) {
            (*b)[i] = -(*res)[i];
        }

        // Call a function to solve the linear system using an iterative method from AztecOO.
        double linNorm;
        int linearIts;
        jac->solve(b, v, linearIts, linNorm, true);

        print0_epBlockMap((*soln).Map());
        print0_epBlockMap((*v).Map());
#endif
        /***************************************************************************/

        int numInitialConstantDeltaTSteps = 3;
        t1.setNumInitialConstantDeltaTSteps(numInitialConstantDeltaTSteps);
        int printNumToTout = 30;
        int printSolnInterval = 1;
        int printSolnFirstSteps = 30;
        int dumpJacobians = 0;
        t1.setPrintSolnOptions(printSolnInterval, printNumToTout, printSolnFirstSteps, dumpJacobians);
        t1.setPrintFlag(PSinput.TimeStepper_printLvl_, PSinput.NonlinSolver_printLvl_);
        t1.setMaxNumTimeSteps(PSinput.maxNumTimeSteps_);
        t1.setMaxStep(PSinput.MaxTimeStep_);
        t1.setMinStep(PSinput.MinTimeStep_);

        /*
         * Create a set of all of times that abrupt changes in the solution occur
         */
        std::set<double> stepTimes;
        stepTimes.insert(PSinput.startTime_);
        stepTimes.insert(PSinput.endTime_);
        std::set<double>::iterator step;

        /*
         *   Look at the boundary condition operator, and figure out when step changes in the
         *   boundary conditions occur.
         *   Cycle through the boundary condition times and add the times to the set, stepTimes
         */
        if (PSinput.BC_TimeDep_) {
            BoundaryCondition *BC = PSinput.BC_TimeDep_;
            double nextTime = BC->nextStep();
            //cycle through step times and add to the set of stepTimes
            do {
                if (nextTime < PSinput.endTime_) {
                    stepTimes.insert(nextTime);
                }
            } while ( (nextTime = BC->nextStep()) > 0);
            BC->resetSteps();
        }

        /*
         * Inform the time stepper of the changes in the times
         */
        std::vector<double> stepTimesV;
        for (step = stepTimes.begin(); step != stepTimes.end(); step++) {
            stepTimesV.push_back(*step);
        }
        t1.setTimeRegionBoundaries(stepTimesV);

        /*
         * Call the Time Integrator for each section of time
         */
        double TInit =  *step;
        double TFinal = *step;
        double TStop;
        step = stepTimes.begin();
        TInit = *step;
        step++;
        for (; step != stepTimes.end(); step++) {
            t1.setInitialTimeStep(PSinput.initialTimeStep_);
            TFinal = *step;
            fprintf(stderr, "BOUNDARY CONDITION time step from %f until %f\n", TInit, TFinal);
            TStop = t1.integratePRE(TFinal);
            if (TStop != TFinal) {
                printf(" Abnormal: Tstop %g not equal to Tfinal %g. Something happened\n", TStop, TFinal);
                break;
            }
            TInit = TFinal;
        }


        /*
         * Cleanup
         */
        Zuzax::appdelete();
        safeDelete(jac);
        safeDelete(ps);

        safeDelete(b);
        safeDelete(v);
        safeDelete(soln);
        safeDelete(res);

#undef PRECIPITATE
#ifdef  PRECIPITATE
        if ( PSinput.useDakota_ ) {
            std::ofstream dfp;
            dfp.open( PSinput.toDakotaFileName_.c_str() , std::ios_base::out );
            dfp << flagPrecipitation << std::endl;
            dfp.close();
        }
#endif // PRECIPITATE
#undef PRECIPITATE

    }
    /*
     *  If we have thrown a ZuzaxError of some sort, we catch it here and print out
     *  an informative error message.
     */
    catch (Zuzax::ZuzaxError) {
        Zuzax::showErrors();
        return -1;
    }
    /*
     *  Epetra (and trilinos) throw errors of type int!!!!
     *  This is kind of unusual. However, we can catch these errors
     *  and print out something useful, rather than letting them go uncaught
     *  to the operating system.
     */
    catch (int eNum) {
        cerr << "main:: Caught Trilinos integer error: " << eNum << "\n\n" << endl;
        return eNum;
    } catch (m1d::m1d_Error &mE) {
        cerr << "caught an error\n" << endl;
        m1d::showErrors(std::cerr);
        return -1;
    }

    /*
     * Release certain Zuzax resources
     */
    Zuzax::appdelete();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

}

//==============================================================================
void generateDomain1D(BatteryResidEval *ps)
{
    (ps->DL_ptr_)->generateDomain1D(ps);
}
//==================================================================================

//==================================================================================
//==================================================================================

