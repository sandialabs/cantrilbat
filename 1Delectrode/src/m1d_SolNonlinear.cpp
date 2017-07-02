/**
 *  @file m1d_SolNonlinear.cpp
 *    Damped Newton solver for 1D multi-domain problems
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "m1d_defs.h"

#include "m1d_Comm.h"
#include "m1d_SolNonlinear.h"
#include "cantera/base/clockWC.h"
#include "cantera/base/utilities.h"
#include "mdp_allo.h"
#include "cantera/base/stringUtils.h"

#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_DataAccess.h"

#include "m1d_LocalNodeIndices.h"
#include <stdio.h>
#include <math.h>


using namespace std;

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x,y)    (( (x) < (y) ) ? (x) : (y))
#endif
//----------------------------------------------------------------------------------------------------------------------------------
namespace m1d
{

//-----------------------------------------------------------
//                 Constants
//-----------------------------------------------------------

const double DampFactor = 4;
const int NDAMP = 10;

//-----------------------------------------------------------
//                 Static Functions
//-----------------------------------------------------------
//==================================================================================================================================
static void print_line(const char* const str, int n)
{
    for (int i = 0; i < n; i++) {
        printf("%s", str);
    }
    printf("\n");
}

//  Default value of a static quantities
bool  SolNonlinear::s_print_NumJac = false;

//==================================================================================================================================
SolNonlinear::SolNonlinear() :
    SolGlobalNonlinear(),
    solnType_(SteadyState_Solve), m_jacFormMethod(0), m_rowScaling(false), m_colScaling(false), colScaleUpdateFrequency_(1),
    m_matrixConditioning(false), m_reltol(1.0E-3), m_abstol(0), m_ewt(0),
    m_ewt_deltaDamping(0),
    m_absTol_deltaDamping(0),
    m_order(1), m_failure_counter(0),
    m_min_newt_its(0),
    maxNewtIts_(50),
    m_jacAge(0), m_curr_linearIts(0), m_curr_normLin(0.0), doResidSolnDamping_(true),
    doDeltaDamping_(true), doHighLowDamping_(false), m_NumLcEqns(0), m_NumLcOwnedEqns(0), m_NumGbEqns(0),
    m_y_curr(0), m_y_curr_owned(0), m_y_new(0), m_y_new_owned(0), m_y_nm1(0), m_y_pred_n(0), m_ydot_curr(0),
    m_ydot_curr_owned(0),
    m_ydot_new(0),
    m_ydot_new_owned(0),
    m_ydot_nm1(0), m_dumpJacobians(0),
    doTimeDependentResid_(false),
    time_n(0.0),
    delta_t_n(0.0),
    m_resid(0),
    m_resid_scaled(false), m_rhs(0), m_residWts(0), m_normResid0(0.0), m_normResidFRaw(0.0), m_normSolnFRaw(0.0),
    solnLowBound_(0), solnHighBound_(0), m_stp(0),
    m_step_2(0), m_func(0), m_rowScales(0), m_colScales(0),
    m_isAlgebraic(0),
    m_isArithmeticScaled(0),
    m_fbound(1.0), m_fdamp(1.0), tdjac_ptr(0), m_nfe(0), m_nJacEval(0),
    m_num_newt_its(0),
    m_numTotalNewtIts(0),
    m_numTotalLinearSolves(0), m_numTotalConvFails(0), m_numTotalTruncFails(0), num_failures(0),
    m_frequencyResidWtUpdatesTD(10),
    m_maxAge(5),
    m_elapsed(0.0), mypid_(0)
{
}
//==================================================================================================================================
SolNonlinear::~SolNonlinear()
{
    safeDelete(m_y_curr_owned);
    safeDelete(m_y_curr);
    safeDelete(m_y_new_owned);
    safeDelete(m_y_new);
    safeDelete(m_y_nm1);
    safeDelete(m_y_pred_n);
    safeDelete(m_ydot_curr);
    safeDelete(m_ydot_curr_owned);
    safeDelete(m_ydot_new);
    safeDelete(m_ydot_new_owned);
    safeDelete(m_ydot_nm1);
    safeDelete(m_resid);
    safeDelete(m_rhs);
    safeDelete(m_residWts);
    safeDelete(solnLowBound_);
    safeDelete(solnHighBound_);
    safeDelete(m_stp);
    safeDelete(m_step_2);
    safeDelete(m_rowScales);
    safeDelete(m_colScales);
    safeDelete(m_isAlgebraic);
    safeDelete(m_isArithmeticScaled);
    safeDelete(m_ewt);
    safeDelete(m_ewt_deltaDamping);
    safeDelete(m_absTol_deltaDamping);
    safeDelete(m_abstol);
}
//==================================================================================================================================
//    L2 Weighted Norm of a delta in the solution
/*
 *   The vector m_ewt[i]'s are always used to weight the solution errors in
 *   the calculation.
 *
 *   The second argument has a default of false. However,
 *   if true, then a table of the largest values is printed
 *   out to standard output.
 *
 *   @param delta_y  Norm of a delta of the solution vector
 *   @param printLargest if True a table is printed of the largest contributors.
 */
double
SolNonlinear::soln_error_norm(const Epetra_Vector_Owned& delta_y, const bool printLargest, const char* title,
                              const int typeYsoln,
                              const double dampFactor) const
{
    int i;
    int idLocalEqnMax = -1;
    double sum_norm = 0.0, error;
    stream0 ss;
    double gbSum = 0.0, gmax1;
    double soln_curr;
    double soln_new;
    for (i = 0; i < m_NumLcOwnedEqns; i++) {
        error = delta_y[i] / (*m_ewt)[i];
        sum_norm += (error * error);
    }
    Comm_ptr_->SumAll(&sum_norm, &gbSum, 1);
    sum_norm = sqrt(gbSum / m_NumGbEqns);
    if (printLargest) {
        if (typeYsoln != 1) {
            printf("not implemented\n");
            exit(-1);
        }
        if (m_print_flag >= 4 && m_print_flag <= 5) {
            if (!mypid_) {
                printf("\t  soln_error_norm(): ");
                if (title) {
                    printf("%s", title);
                } else {
                    printf(" Delta soln norm ");
                }
                printf(" = %-11.4E\n", sum_norm);
            }
        } else if (m_print_flag >= 6) {

            const int num_entries = 8;
            double dmax1, normContrib;
            int j;
            std::vector<int> imax(num_entries, -1);
            print0_sync_start(false, ss, *Comm_ptr_);
            if (!mypid_) {
                printf("\t  ");
                print_line("-", 90);
                printf("\t  soln_error_norm(): ");
                if (title) {
                    printf("%s", title);
                } else {
                    printf(" Delta soln norm ");
                }
                printf(" = %-11.4E\n", sum_norm);
                printf("\t\t      Printout of Largest Contributors:\n");
                printf("\t\t                                                      (damp = %g)\n", dampFactor);
                printf("\t\t      I             VarName  LcNode weightdeltaY/sqtN|     deltaY      ysolnOld     ysolnNew   Soln_Weights\n");
                printf("\t\t   ");
                print_line("-", 80);
            }
            print0_sync_end(false, ss, *Comm_ptr_);
            for (int jnum = 0; jnum < num_entries; jnum++) {
                dmax1 = -1.0;
                idLocalEqnMax = 0;
                // pick out the contribution for this processor
                for (i = 0; i < m_NumLcOwnedEqns; i++) {
                    bool used = false;
                    for (j = 0; j < jnum; j++) {
                        if (imax[j] == i) {
                            used = true;
                        }
                    }
                    if (!used) {
                        error = delta_y[i] / (*m_ewt)[i];
                        normContrib = sqrt(error * error);
                        if (normContrib > dmax1) {
                            idLocalEqnMax = i;
                            dmax1 = normContrib;
                        }
                    }
                }
                int procWinner = procChoice_Max(dmax1, *Comm_ptr_, gmax1);
                print0_sync_start(false, ss, *Comm_ptr_);
                if (procWinner == mypid_) {
                    imax[jnum] = idLocalEqnMax;
                    i = idLocalEqnMax;
                    int idGlobalEqnMax = m_func->LI_ptr_->IndexGbEqns_LcEqns[idLocalEqnMax];
                    if (i >= 0) {
                        int iLcNode;
                        int iGbNode;
                        int iNodeEqnNum;
                        VarType var;
                        VAR_TYPE_SUBNUM vtsub;
                        std::string vstring = m_func->variableID(i, iLcNode, iGbNode, iNodeEqnNum, var, vtsub);
                        string v16 = var.VariableName(16);
                        error = delta_y[i] / (*m_ewt)[i];
                        soln_curr = (*m_y_curr)[idLocalEqnMax];
                        soln_new = (*m_y_curr)[idLocalEqnMax] + delta_y[idLocalEqnMax] * dampFactor;
                        if (solnType_ == DAESystemInitial_Solve) {
                            if ((*m_isAlgebraic)[idLocalEqnMax] != 1) {
                                soln_curr = (*m_ydot_curr)[idLocalEqnMax];
                                soln_new = (*m_ydot_curr)[idLocalEqnMax] + delta_y[idLocalEqnMax] * dampFactor;
                            }
                        }
                        ss.print0("\t\t   %4d %24s-%-4d  %12.4e  | %12.4e %12.4e %12.4e %12.4e\n", idGlobalEqnMax, v16.c_str(),
                                  iLcNode, error / sqrt(m_NumGbEqns),
                                  delta_y[idLocalEqnMax], soln_curr, soln_new, (*m_ewt)[idLocalEqnMax]);
                    }
                }
                print0_sync_end(false, ss, *Comm_ptr_);
            }
            print0_sync_start(false, ss, *Comm_ptr_);
            if (!mypid_) {
                printf("\t\t   ");
                print_line("-", 80);
                printf("\t  ");
                print_line("-", 90);
            }
            print0_sync_end(false, ss, *Comm_ptr_);

        }
    }
    return sum_norm;
}
//==================================================================================================================================
//    L2 Weighted Norm of the residual
/*
 *   The vector m_residWts[i]'s are always used to weight the residual errors in
 *   the calculation.
 *
 *   The second argument has a default of false. However,
 *   if true, then a table of the largest values is printed
 *   out to standard output.
 *
 *   @param delta_y  Norm of a delta of the solution vector
 *   @param printLargest if True a table is printed of the largest contributors.
 */
double
SolNonlinear::res_error_norm(const Epetra_Vector_Owned& resid, const char* title, const int printLargest) const
{
    int i, idLocalEqnMax = -1;
    double sum_norm = 0.0, error;
    stream0 ss;
    double gbSum = 0.0, gmax1;
    if (m_resid_scaled) {
        for (i = 0; i < m_NumLcOwnedEqns; i++) {
            error = resid[i] / ((*m_rowScales)[i] * (*m_residWts)[i]);
            sum_norm += (error * error);
        }
    } else {
        for (i = 0; i < m_NumLcOwnedEqns; i++) {
            error = resid[i] / (*m_residWts)[i];
            sum_norm += (error * error);
        }
    }
    Comm_ptr_->SumAll(&sum_norm, &gbSum, 1);
    sum_norm = sqrt(gbSum / m_NumGbEqns);
    if (printLargest) {
        int num_entries = abs(printLargest);
        double dmax1, normContrib;
        int j;
        std::vector<int> imax(num_entries, -1);
        print0_sync_start(false, ss, *Comm_ptr_);
        if (!mypid_) {
            printf("\t  ");
            print_line("-", 90);
            printf("\t  res_error_norm():");
            if (title) {
                printf(" %s ", title);
            } else {
                printf("  residual L2 norm ");
            }
            printf("= %12.4E\n", sum_norm);

            if (m_print_flag >= 6) {
                printf("\t\tPrintout of Largest Contributors to norm:\n");
                printf("\t\t      I      EqnName       LcNode      | Resid/ResWt|     UnsclRes          ResWt\n");
                printf("\t\t   ");
                print_line("-", 80);
            }
        }
        print0_sync_end(false, ss, *Comm_ptr_);
        if (m_print_flag >= 6) {
            for (int jnum = 0; jnum < num_entries; jnum++) {
                dmax1 = -1.0;
                // pick out the contribution for this processor
                for (i = 0; i < m_NumLcOwnedEqns; i++) {
                    bool used = false;
                    for (j = 0; j < jnum; j++) {
                        if (imax[j] == i) {
                            used = true;
                        }
                    }
                    if (!used) {
                        if (m_resid_scaled) {
                            error = resid[i] / ((*m_rowScales)[i] * (*m_residWts)[i]);
                        } else {
                            error = resid[i] / (*m_residWts)[i];
                        }
                        normContrib = sqrt(error * error);
                        if (normContrib > dmax1) {
                            idLocalEqnMax = i;
                            //imax[jnum] = i;
                            dmax1 = normContrib;
                        }
                    }
                }
                int procWinner = procChoice_Max(dmax1, *Comm_ptr_, gmax1);
                print0_sync_start(false, ss, *Comm_ptr_);
                if (procWinner == mypid_) {
                    imax[jnum] = idLocalEqnMax;
                    i = idLocalEqnMax;
                    int idGlobalEqnMax = m_func->LI_ptr_->IndexGbEqns_LcEqns[idLocalEqnMax];
                    int iLcNode;
                    int iGbNode;
                    int iNodeEqnNum;
                    EqnType var;
                    EQ_TYPE_SUBNUM vtsub;
                    std::string vstring = m_func->equationID(i, iLcNode, iGbNode, iNodeEqnNum, var, vtsub);
                    string v24 = var.EquationName(24);
                    if (m_resid_scaled) {
                        error = resid[i] / ((*m_rowScales)[i] * (*m_residWts)[i]);
                    } else {
                        error = resid[i] / (*m_residWts)[i];
                    }
                    if (m_resid_scaled) {
                        error = resid[i]/ ((*m_rowScales)[i] * (*m_residWts)[i]);
                        if (i >= 0) {
                            ss.print0("\t\t   %4d  %24s-%-4d |%12.4e   %12.4e     %12.4e\n", idGlobalEqnMax,  v24.c_str(),
                                      iLcNode, fabs(error),
                                      (resid[idLocalEqnMax] / (*m_rowScales)[idLocalEqnMax]), (*m_residWts)[idLocalEqnMax]);
                        }
                    } else {
                        error = resid[i] / (*m_residWts)[i];
                        if (i >= 0) {
                            ss.print0("\t\t   %4d  %24s-%-4d |%12.4e   %12.4e     %12.4e\n", idGlobalEqnMax, v24.c_str(),
                                      iLcNode, fabs(error),
                                      resid[idLocalEqnMax], (*m_residWts)[idLocalEqnMax]);
                        }
                    }
                }
                print0_sync_end(false, ss, *Comm_ptr_);
            }
        }
        print0_sync_start(false, ss, *Comm_ptr_);
        if (!mypid_ && m_print_flag >= 6) {
            printf("\t\t   ");
            print_line("-", 80);
            printf("\t  ");
            print_line("-", 90);
        }
        print0_sync_end(false, ss, *Comm_ptr_);
    }
    return sum_norm;
}
//=====================================================================================================================
int
SolNonlinear::get_jac(EpetraJac& jac,
                      Epetra_Vector_Owned* res,
                      const bool doTimeDependentResid,
                      double time_curr,
                      double rdelta_t,
                      const Epetra_Vector_Ghosted* solnBase_ptr,
                      const Epetra_Vector_Ghosted* solnDotBase_ptr)
{
    m_nfe++;

    jac.matrixResEval(doTimeDependentResid, solnBase_ptr, solnDotBase_ptr, res, time_curr, rdelta_t, solnType_);
    // Only print the residual and matrix if print flag is high and static bool is set to true
    if (m_print_flag >= 7 && s_print_NumJac) {
        string ss = "Solution Values";
        m_func->showSolutionVector(ss, time_curr, 1.0/rdelta_t, *solnBase_ptr);
        ss = "Solution Time Derivative";
        m_func->showSolutionVector(ss, time_curr, 1.0/rdelta_t, *solnDotBase_ptr);
        ss = "Residual";
        m_func->showSolutionVector(ss, time_curr, 1.0/rdelta_t, *res);
        print0_epIntVector(*m_isAlgebraic, "IsAlgebraic");
        print0_epVbrMatrix(*(jac.A_), "Nonlinear matrix");
#ifdef DEBUG_INIT_CALCULATION
        if (m_nfe == 1) {
            FILE* ff = fopen("FirstJacDAE.txt", "w");
            print0_epVbrMatrix(*(jac.A_), "Nonlinear matrix", ff);
            fclose(ff);
        }
#endif
    }
    return 0;
}
//==================================================================================================================================
void
SolNonlinear::get_res(const double time_curr,
                      const double rdelta_t,
                      const Epetra_Vector_Ghosted* solnBase_ptr,
                      const Epetra_Vector_Ghosted* solnDotBase_ptr)
{
    m_func->residEval(m_resid, doTimeDependentResid_, solnBase_ptr, solnDotBase_ptr, time_curr, rdelta_t, Base_ResidEval,
                      solnType_);
    m_nfe++;
    m_resid_scaled = false;
}

//===================================================================================================================================
void
SolNonlinear::scaleMatrix(Epetra_Vector_Owned& delta_soln,
                          const Epetra_Vector_Ghosted& y_curr,
                          const Epetra_Vector_Ghosted& ydot_curr,
                          const double time_curr,
                          const double rdelta_t,
                          int loglevel)
{
    EpetraJac& jac = *tdjac_ptr;

    /*
     * Column scaling -> We scale the columns of the Jacobian
     * by the nominal important change in the solution vector
     */
    if (m_colScaling) {
        if (!jac.m_columnScaled) {
            /*
             * Go get new scales, put the scales in the vector,
             * m_colScales
             */
            if (colScaleUpdateFrequency_ >= 2) {
                setDefaultColumnScaleVector();
            }

            // There is a canned routine to do column scaling from the
            // right of matrices.
            jac.columnScale(m_colScales);
        }
    }

    /*
     * row sum scaling -> Note, this is an unequivocal success
     *      at keeping the small numbers well balanced and nonnegative.
     */
    if (m_rowScaling) {
        if (!jac.m_rowScaled) {
            /*
             * Go get new scales, put the scales in the vector,
             * m_rowScales
             */
            jac.getRowScales(m_rowScales);
            /*
             * Apply the scales to the matrix
             */
            jac.rowScale(m_rowScales);
        }
        /*
         * Apply the row scales to the right hand side
         */
        if (!m_resid_scaled) {
            for (int irow = 0; irow < m_NumLcOwnedEqns; irow++) {
                (*m_resid)[irow] *= (*m_rowScales)[irow];
            }
            m_resid_scaled = true;
        }
        int noww = (m_num_newt_its-1) % m_frequencyResidWtUpdatesTD;
        if (!noww) {
            computeResidWts();
        }
    }
}
//===================================================================================================================================
//  Compute the undamped Newton step.
/*
 * The residual function is
 * evaluated at the current time, time_curr, at the current values of the
 * solution vector, y_curr, and the solution time derivative, ydot_curr,
 * but the Jacobian is not recomputed.
 */
void
SolNonlinear::doNewtonSolve(Epetra_Vector_Owned& delta_soln,
                            const Epetra_Vector_Ghosted& y_curr,
                            const Epetra_Vector_Ghosted& ydot_curr,
                            const double time_curr,
                            const double rdelta_t,
                            int loglevel)
{
    EpetraJac& jac = *tdjac_ptr;
    int irow;

    for (int n = 0; n < m_NumLcOwnedEqns; n++) {
        (*m_rhs)[n] = -(*m_resid)[n];
    }

    // HKM I Missed row scaling the rhs !!!

    if (m_matrixConditioning) {

    }

#ifdef DEBUG_MODE
    // bool printJacContributions = true;
    // int numRows = 10;
    //  int focusGbEqn = 0;
    // int focusLcEqn = m_func->GbEqnToLcEqn(focusGbEqn);
#endif
    /*
     * Solve the system -> This also involves inverting the
     * matrix
     */
    bool doResCheck = false;
    if (m_print_flag >= 2) {
        doResCheck = true;
    }
    /*
     * Because we have previously row and column scaled the matrix and the rhs, the
     * norm computed by AX-b will be properly scaled.
     */
    int retn = jac.solve(m_rhs, &delta_soln, m_curr_linearIts, m_curr_normLin, doResCheck);
    if (retn) {
        throw m1d_Error("nonlinearSolve", "jac.solve returned an error");
    }

#ifdef DEBUG_MODE
    if (1) {
        for (int j = 0; j < m_NumLcOwnedEqns; j++) {
            ZZCantera::checkFinite(delta_soln[j]);
        }
    }
#endif
    /*
     * reverse the column scaling if there was any.
     */
    if (m_colScaling) {
        for (irow = 0; irow < m_NumLcOwnedEqns; irow++) {
            delta_soln[irow] *= (*m_colScales)[irow];
        }
    }

#ifdef DEBUG_INIT_CALCULATION
    /*
     *  HKM - This is an example of how to extract a formatted solution component and
     *        dump it to a file, when you want to dump it to a file
     */
    if (m_print_flag >= 7 && s_print_NumJac) {
        string ss = "Raw Solution Delta";
        FILE* ff = stdout;
        if (m_numTotalLinearSolves == 0) {
            ff = fopen("FirstSolDelta.txt", "w");
        }
        m_func->showSolutionVector(ss, time_curr, 1.0/rdelta_t, delta_soln, ff);
        if (ff != stdout) {
            fclose(ff);
        }
    }
#endif

#ifdef DEBUG_MODE_NOT
    if (printJacContributions) {
        for (int iNum = 0; iNum < numRows; iNum++) {
            if (iNum > 0) {
                focusRow++;
            }
            double dsum = 0.0;
            vector_fp& Jdata = jacBack.data();
            double dRow = Jdata[m_neq * focusRow + focusRow];
            printf("\n Details on delta_Y for row %d \n", focusRow);
            printf("  Value before = %15.5e, delta = %15.5e,"
                   "value after = %15.5e\n", y_curr[focusRow], delta_y[focusRow], y_curr[focusRow] + delta_y[focusRow]);
            if (!freshJac) {
                printf("    Old Jacobian\n");
            }
            printf("     col          delta_y            aij     "
                   "contrib   \n");
            printf("--------------------------------------------------"
                   "---------------------------------------------\n");
            printf(" Res(%d) %15.5e  %15.5e  %15.5e  (Res = %g)\n", focusRow, delta_y[focusRow], dRow, RRow[iNum] / dRow,
                   RRow[iNum]);
            dsum += RRow[iNum] / dRow;
            for (int ii = 0; ii < m_neq; ii++) {
                if (ii != focusRow) {
                    double aij = Jdata[m_neq * ii + focusRow];
                    double contrib = aij * delta_y[ii] * (-1.0) / dRow;
                    dsum += contrib;
                    if (fabs(contrib) > Pcutoff) {
                        printf("%6d  %15.5e  %15.5e  %15.5e\n", ii, delta_y[ii], aij, contrib);
                    }
                }
            }
            printf("--------------------------------------------------"
                   "---------------------------------------------\n");
            printf("        %15.5e                   %15.5e\n", delta_y[focusRow], dsum);
        }
    }

#endif

    m_numTotalLinearSolves++;
}
//===================================================================================================================================
int
SolNonlinear::doHardBounds(const Epetra_Vector_Ghosted& y_old, Epetra_Vector_Owned& step, double& fbound)
{
    double fbound_new = 1.0;
    bool lowered = false;
    if (m_print_flag > 3 && !mypid_) {
        printf("\t  ");
        print_line("-", 90);
        printf("\t  doHardBounds(): Checking for bounds on the step size\n");
    }
    // Compute the multiplier to keep all components in bounds
    // A value of one indicates that there is no limitation
    // on the current step size in the nonlinear method due to
    // bounds constraints (either negative values of delta
    // bounds constraints.
    if (doHighLowDamping_) {
        fbound_new = highLowBoundStep(y_old, step, m_print_flag);
    }
    if (fbound_new < fbound) {
        lowered = true;
        fbound = fbound_new;
    }
    //   Do the delta bounds constraint algorithm
    if (doDeltaDamping_) {
        fbound_new = deltaBoundStep(y_old, step);
        if (fbound_new < fbound) {
            lowered = true;
            fbound = fbound_new;
        }
    }
    if (m_print_flag > 3 && !mypid_) {
        if (!lowered && (fbound < 1.0)) {
            printf("\t              ... fbound increased maximally to %10.3g\n", fbound);
        }
    }
    //  If fbound is very small, then y0 is already close to the boundary and step0 points out of the allowed domain.
    //  In this case, the Newton algorithm fails, so return an error condition.
    /*
     *  This algorithm is just not sufficient, and I expect that there will be a lot more work that will be inserted into
     *  this routine, here.
     */
    if (fbound < 1.e-20) {
        if (m_print_flag > 1 && !mypid_) {
            printf("\t  doHardBounds(): ERROR! At limits.\n");
        }
        return -3;
    }
    // Scale the step size to the correct value to ensure that the bounds are satisfied.
    if (fbound < 1.0) {
        step.Scale(fbound);
    }
    // Print final line
    if (m_print_flag > 3 && !mypid_) {
        printf("\t  ");
        print_line("-", 90);

    }
    if (fbound < 1.0) {
        return 0;
    }
    return 1;
}
//==================================================================================================================================
/*
 *  Return the factor by which the undamped Newton step 'step0' must be multiplied in order to keep all solution components in
 *  all domains not-appreciably changing so much that the jacobian isn't representative.
 *  The idea behind these is that the Jacobian couldn't possibly be representative, if the  variable is changed by a lot. 
 *  (true for nonlinear systems, false for linear systems)
 *
 *    For variables which have a strict minimum of zero: 
 *      Maximum increase in variable in any one newton iteration: 
 *          a)   factor of 2 when above the value of the fabs(change) is above m_ewt_deltaDamping[i]
 *          b)   Equal to the change if the fabs(change) is below m_ewt_deltaDamping[i].
 *
 *      Maximum decrease in variable in any one newton iteration:
 *          a)   factor of 5 when above the value of the fabs(change) is above m_ewt_deltaDamping[i]
 *          b)   Equal to the change if the fabs(change) is below m_ewt_deltaDamping[i].
 *
 *    For arithmetically scaled variables, the maximum increase or decrease in an iteration is given by the value of 
 *    m_ewt_deltaDamping[i].
 */
double
SolNonlinear::deltaBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0)
{
    int i_fbounds = 0;
    int ifbd = 0;

    double fbound = 1.0;
    double f_delta_bounds = 1.0;
    double ff_alt;
    Epetra_IntVector& isS = *m_isArithmeticScaled;

    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        double y_curr = y[i];
        double y_new = y[i] + step0[i];
        if (solnType_ == DAESystemInitial_Solve) {
            if ((*m_isAlgebraic)[i] != 1) {
                y_curr = (*m_ydot_curr)[i];
                y_new = (*m_ydot_curr)[i] + step0[i];
            }
        }
        /*
         * Now do a delta bounds
         * Increase variables by a factor of 2 only
         * decrease variables by a factor of 5 only
         */
        double ff = 1.0;
        if (isS[i] == 1) {
            if ((y_new - y_curr) > (*m_ewt_deltaDamping)[i]) {
                ff = ((*m_ewt_deltaDamping)[i]) / (y_new - y_curr);
                ifbd = 1;
            }
            if ((y_new - y_curr) < - (*m_ewt_deltaDamping)[i]) {
                ff = (- (*m_ewt_deltaDamping)[i]) / (y_new - y_curr);
                ifbd = 0;
            }
        } else {
            const double DELTA_FAC_INCREASE = 2.0;

            if ((fabs(y_new) > DELTA_FAC_INCREASE * fabs(y_curr)) && (fabs(y_new - y_curr) > (*m_ewt_deltaDamping)[i])) {
                ff = fabs(y_curr / (y_new - y_curr));
                ff_alt = fabs((*m_ewt_deltaDamping)[i] / (y_new - y_curr));
                ff = MAX(ff, ff_alt);
                ifbd = 1;
            }
            const double DELTA_FAC_DECREASE = 5.0;
            if ((fabs(DELTA_FAC_DECREASE * y_new) < fabs(y_curr)) && (fabs(y_new - y_curr) > (*m_ewt_deltaDamping)[i])) {
                ff = y_curr / (y_new - y_curr) * (1.0 - DELTA_FAC_DECREASE) / DELTA_FAC_DECREASE;
                ff_alt = fabs((*m_ewt_deltaDamping)[i] / (y_new - y_curr));
                ff = MAX(ff, ff_alt);
                ifbd = 0;
            }
        }
        if (ff < f_delta_bounds) {
            f_delta_bounds = ff;
            i_fbounds = i;
        }
    }
    // Communicate between processors to find the lowest values
    double flocal = f_delta_bounds;
    int i_fbounds_gb = m_func->LI_ptr_->IndexGbEqns_LcEqns[i_fbounds];
    int iproc = procChoice_Min_Brcst1Int(flocal, *Comm_ptr_, f_delta_bounds, i_fbounds_gb);
    fbound = f_delta_bounds;

    int iLcNode;
    int iGbNode;
    int iNodeEqnNum;
    VarType var;
    VAR_TYPE_SUBNUM vtsub;
    std::string vstring = m_func->variableID(i_fbounds, iLcNode, iGbNode, iNodeEqnNum, var, vtsub);
    std::string v24 = var.VariableName(24);

    /*
     * Report on any corrections
     */
    if (fbound < 1.0 - 1.0E-13) {
        stream0 ss;
        if (m_print_flag == 3) {
            print0_sync_start(false, ss, *Comm_ptr_);
            if (iproc == mypid_) {
                if (ifbd) {
                    if ((solnType_ == DAESystemInitial_Solve) && ((*m_isAlgebraic)[i_fbounds] != 1)) {
                        printf("\t   ... delta damping: Increase of variableTimeDeriv  %24s-%-4d  (iG=%d) causing "
                               "delta damping of %12.5E: origVal = %12.5E, undampedNew = %12.5E, dampedNew = %12.5E\n",
                               v24.c_str(), iLcNode, i_fbounds_gb, f_delta_bounds,
                               (*m_ydot_curr)[i_fbounds], (*m_ydot_curr)[i_fbounds] + step0[i_fbounds] ,
                               (*m_ydot_curr)[i_fbounds] + f_delta_bounds*step0[i_fbounds]);
                    } else {
                        printf("\t   ... delta damping: Increase of variable  %24s-%-4d  (iG=%d) causing "
                               "delta damping of %12.5E: origVal = %12.5E, undampedNew = %12.5E, dampedNew = %12.5E\n",
                               v24.c_str(), iLcNode, i_fbounds_gb, f_delta_bounds,
                               y[i_fbounds], y[i_fbounds]+step0[i_fbounds] , y[i_fbounds]+f_delta_bounds*step0[i_fbounds]);
                    }
                } else {
                    if ((solnType_ == DAESystemInitial_Solve) && ((*m_isAlgebraic)[i_fbounds] != 1)) {
                        printf("\t   ... delta damping: Decrease of variableTimeDeriv  %24s-%-4d  (iG=%d) causing "
                               "delta damping of %12.5E: origVal = %12.5E, undampedNew = %12.5E, dampedNew = %12.5E\n",
                               v24.c_str(), iLcNode, i_fbounds_gb, f_delta_bounds,
                               (*m_ydot_curr)[i_fbounds], (*m_ydot_curr)[i_fbounds] + step0[i_fbounds] ,
                               (*m_ydot_curr)[i_fbounds] + f_delta_bounds*step0[i_fbounds]);
                    } else {
                        printf("\t   ... delta damping: Decrease of variable  %24s-%-4d  (iG=%d) causing "
                               "delta damping of %12.5E: origVal = %12.5E, undampedNew = %12.5E, dampedNew = %12.5E\n",
                               v24.c_str(), iLcNode, i_fbounds_gb, f_delta_bounds,
                               y[i_fbounds], y[i_fbounds]+step0[i_fbounds] , y[i_fbounds]+f_delta_bounds*step0[i_fbounds]);
                    }
                }
            }
            print0_sync_end(false, ss, *Comm_ptr_);
        } else if (m_print_flag > 3) {
            Comm_ptr_->Broadcast(&ifbd, 1, iproc);
            stream0 ss;
            print0_sync_start(false, ss, *Comm_ptr_);
            if (!mypid_) {
                if (ifbd) {
                    if ((solnType_ == DAESystemInitial_Solve) && ((*m_isAlgebraic)[i_fbounds] != 1)) {
                        printf("\t\tboundStep:  Increase of variableTimeDeriv  %24s-%-4d  (iG=%d) causing "
                               "delta damping of %12.5E: origVal = %12.5E, undampedNew = %12.5E, dampedNew = %12.5E\n",
                               v24.c_str(), iLcNode, i_fbounds_gb, f_delta_bounds,
                               (*m_ydot_curr)[i_fbounds], (*m_ydot_curr)[i_fbounds] + step0[i_fbounds] ,
                               (*m_ydot_curr)[i_fbounds] + f_delta_bounds*step0[i_fbounds]);
                    } else {
                        printf("\t\tboundStep:  Increase of variable  %24s-%-4d  (iG=%d) causing "
                               "delta damping of %12.5E: origVal = %12.5E, undampedNew = %12.5E, dampedNew = %12.5E\n",
                               v24.c_str(), iLcNode, i_fbounds_gb, f_delta_bounds,
                               y[i_fbounds], y[i_fbounds]+step0[i_fbounds] , y[i_fbounds]+f_delta_bounds*step0[i_fbounds]);
                    }
                } else {
                    if ((solnType_ == DAESystemInitial_Solve) && ((*m_isAlgebraic)[i_fbounds] != 1)) {
                        printf("\t\tboundStep:  Decrease of variableTimeDeriv  %24s-%-4d  (iG=%d) causing "
                               "delta damping of %12.5E: origVal = %12.5E, undampedNew = %12.5E, dampedNew = %12.5E\n",
                               v24.c_str(), iLcNode, i_fbounds_gb, f_delta_bounds,
                               (*m_ydot_curr)[i_fbounds], (*m_ydot_curr)[i_fbounds]+step0[i_fbounds] ,
                               (*m_ydot_curr)[i_fbounds]+f_delta_bounds*step0[i_fbounds]);
                    } else {
                        printf("\t\tboundStep:  Decrease of variable  %24s-%-4d  (iG=%d) causing "
                               "delta damping of %12.5E: origVal = %12.5E, undampedNew = %12.5E, dampedNew = %12.5E\n",
                               v24.c_str(), iLcNode, i_fbounds_gb, f_delta_bounds,
                               y[i_fbounds], y[i_fbounds]+step0[i_fbounds] , y[i_fbounds]+f_delta_bounds*step0[i_fbounds]);
                    }
                }
            }
            print0_sync_end(false, ss, *Comm_ptr_);
        }
    }

    return fbound;
}
//==================================================================================================================================
double
SolNonlinear::highLowBoundStep(const Epetra_Vector_Ghosted& y, const Epetra_Vector_Owned& step0, int loglevel)
{
    int i_lower = 0;
    int i_higher = 0;

    double f_bounds = 1.0;
    double f_low_bounds = 1.0;
    double f_high_bounds = 1.0;
    double ff, y_new;
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        y_new = y[i] + step0[i];

        double y_lbound = (*solnLowBound_)[i] - 0.01 * (*m_ewt)[i];
        if ((y_new < y_lbound) && y[i] >= (*solnLowBound_)[i]) {
            ff = 0.9 * (y[i] / (y[i] - y_new));
            if (ff < f_low_bounds) {
                f_low_bounds = ff;
                i_lower = i;
            }
        }

        double y_hbound = (*solnHighBound_)[i] + 0.01 * (*m_ewt)[i];
        if ((y_new > y_hbound) && y[i] <= (*solnHighBound_)[i]) {
            ff = 0.9 * (y[i] / (y[i] - y_new));
            if (ff < f_high_bounds) {
                f_high_bounds = ff;
                i_higher = i;
            }
        }

    }
    double flocal = f_high_bounds;
    int i_higher_gb = m_func->LI_ptr_->IndexGbNode_LcNode[i_higher];
    procChoice_Min_Brcst1Int(flocal, *Comm_ptr_, f_high_bounds, i_higher_gb);

    flocal = f_low_bounds;
    int i_lower_gb = m_func->LI_ptr_->IndexGbNode_LcNode[i_lower];
    procChoice_Min_Brcst1Int(flocal, *Comm_ptr_, f_low_bounds, i_lower_gb);

    f_bounds = MIN(f_low_bounds, f_high_bounds);
    /*
     * Report on any corrections
     */
    if (loglevel > 3) {
        stream0 ss;
        print0_sync_start(false, ss, *Comm_ptr_);
        if (!mypid_) {
            if (f_bounds < 1.0 - 1.0E-13) {
                if (f_low_bounds < f_bounds) {
                    printf("\t\tboundStep: Global Variable %d causing lower bounds "
                           "damping of %g\n", i_lower_gb, f_low_bounds);
                } else {
                    printf("\t\tboundStep: Global Variable %d causing high bounds "
                           "damping of %g\n", i_higher_gb, f_high_bounds);
                }
            }
        }
        print0_sync_end(false, ss, *Comm_ptr_);
    }
    return f_bounds;
}
//=====================================================================================================================
/*
 * On entry, step0 must contain an undamped Newton step for the
 * solution x0. This method attempts to find a damping coefficient
 * such that the next undamped step would have a norm smaller than
 * that of step0. If successful, the new solution after taking the
 * damped step is returned in y_new, and the undamped step at y_new is
 * returned in step1.
 *
 * @return   1 Successful step was taken: Next step was less than previous step.
 *                                        s1 is calculated
 *           2 Successful step: Next step's norm is less than 0.8
 *           3 Success:  The final residual is less than 1.0
 *                        A predicted deltaSoln is not produced however. s1 is estimated.
 *           4 Success:  The final residual is less than the residual
 *                       from the previous step.
 *                        A predicted deltaSoln is not produced however. s1 is estimated.
 *           0 Uncertain Success: s1 is about the same as s0
 *          -2 Unsuccessful step.
 */
int
SolNonlinear::dampStep(double time_curr,  const Epetra_Vector_Ghosted& y0,  const Epetra_Vector_Ghosted* ydot0_ptr,
                       double& s1, int& loglevel, int& num_backtracks)
{
    Epetra_Vector_Owned& step0 = *m_stp;
    Epetra_Vector_Owned& step1 = *m_step_2;
    int retnTrial = 0;

    // Compute the weighted norm of the undamped step size step0
    double s0 = m_normSolnFRaw;
    std::string ResS;

    double rdelta_t = 0.0;
    if (delta_t_n > 1.0E-300) {
        rdelta_t = 1.0 / delta_t_n;
    }
    //--------------------------------------------
    //           Attempt damped step
    //--------------------------------------------

    // damping coefficient starts at 1.0
    m_fdamp = 1.0;
    int m;
    double ff;
    num_backtracks = 0;
    double fbound = m_fbound;
    for (m = 0; m < NDAMP; m++) {
        ff = m_fdamp;
        // step the solution by the damped step size
        /*
         * Whenever we update the solution, we must also always update the time derivative.
         */
        updateSoln(y0, ydot0_ptr, ff, step0);
        /*
         *  Compute the next residual that would result if m_y_new[] were accepted.
         *  This is computed and storred in the vector m_resid[].
         *  The norm of this new residual is storred in m_normResidTrial and in m_normResidFRaw if
         *  this is the first iteration and therefore represents the raw Newton step.
         */
        get_res(time_curr, rdelta_t, m_y_new, m_ydot_new);

        if (m_print_flag >= 6) {
            ResS = "Residual For Damping Trial " + ZZCantera::int2str(m) + " with Damping Coeff " + ZZCantera::fp2str(ff);
            m_normResidTrial = res_error_norm(*m_resid, ResS.c_str(), 10);
        } else if (m_print_flag == 4 || m_print_flag == 5) {
            ResS = "Residual For Damp Trial " + ZZCantera::int2str(m) + " with Damping Coeff " + ZZCantera::fp2str(ff);
            m_normResidTrial = res_error_norm(*m_resid, ResS.c_str(), true);
        } else {
            m_normResidTrial = res_error_norm(*m_resid);
        }

        if (m == 0) {
            m_normResidFRaw = m_normResidTrial;
        }

        if (m_normResidTrial < 1.0 || m_normResidTrial < m_normResid0 || !doResidSolnDamping_) {
            if (loglevel >= 5 && !mypid_) {
                if (m_normResidTrial < 1.0) {
                    printf("\t  dampStep(): Current trial step and damping"
                           " coefficient accepted because residTrial test step < 1:\n");
                    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
                } else if (m_normResidTrial < m_normResid0) {
                    printf("\t  dampStep(): Current trial step and damping"
                           " coefficient accepted because resid0 > residTrial:\n");
                    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
                } else {
                    printf("\t  dampStep(): Current trial step and damping"
                           " coefficient accepted because residual solution damping is turned off:\n");
                    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
                }
            }
            /*
             *  We aren't going to solve the system if we don't need to. Therefore, return an estimate
             *  of the next solution update based on the ratio of the residual reduction.
             */
            if (m_normResid0 < 1.0E-10) {
                if (m_normResidTrial < 1.0E-5) {
                    s1 = 0.1 * s0;
                }
            } else {
                s1 = s0 * m_normResidTrial / m_normResid0;
            }
            if (m_normResidTrial < 1.0) {
                retnTrial = 3;
            } else {
                retnTrial = 4;
            }
            break;
        }
        /*
         *  Compute the next undamped step, step1[], that would result if y1[] were accepted.
         *  This is computed and storred in the vector m_resid
         */
        doNewtonSolve(step1, *m_y_new, *m_ydot_new, time_curr, rdelta_t, loglevel);

        doHardBounds(*m_y_new, step1, fbound);

        // compute the weighted norm of step1
        s1 = soln_error_norm(step1, false);

        // write log information
        if (loglevel >= 4) {
            print_solnDelta_norm_contrib(step0, "DeltaSoln", step1, "DeltaSolnTrialTest",
                                         "dampNewt: Important Entries for Weighted Soln Updates:", y0, *m_y_new, ff, 5);
        }
        if (loglevel >= 4 && !mypid_) {
            printf("\t\t\tdampStep(): s0 = %g, s1 = %g, "
                   "damp = %g\n", s0, s1, m_fdamp);
        }

        // if the norm of s1 is less than the norm of s0, then
        // accept this damping coefficient. Also accept it if this
        // step would result in a converged solution. Otherwise,
        // decrease the damping coefficient and try again.

        if (s1 < 0.8 || s1 < s0) {

            if (s1 > s0) {
                if (loglevel >= 4 && !mypid_) {
                    printf("\t\t\tdampStep(): current trial step and damping coefficient accepted because test step < 1\n");
                    printf("\t\t\t          s1 = %g, s0 = %g\n", s1, s0);
                }
                retnTrial = 2;
            } else {
                retnTrial = 1;
            }
            break;
        } else {
            if (loglevel >= 4 && !mypid_) {
                printf("\t\t\tdampStep(): current step rejected: (s1 = %g > s0 = %g)", s1, s0);
                if (m < (NDAMP - 1)) {
                    printf(" Decreasing damping factor and retrying");
                } else {
                    printf(" Giving up!!!");
                }
                printf("\n");
            }
        }
        num_backtracks++;
        m_fdamp /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the
    // solution after stepping by the damped step would represent
    // a converged solution, and return 0 otherwise. If no damping
    // coefficient could be found, return -2.
    if (m < NDAMP) {
        if (loglevel >= 4 && !mypid_) {
            printf("\t  dampStep(): current trial step accepted retnTrial = %d, its = %d, damp = %g\n", retnTrial, m+1, m_fdamp);
        }
        return retnTrial;
    } else {
        if (s1 < 0.5 && (s0 < 0.5)) {
            if (loglevel >= 4 && !mypid_) {
                printf("\t  dampStep(): current trial step accepted kindof retnTrial = %d, its = %d, damp = %g\n", 2, m+1, m_fdamp);
            }
            return 2;
        }
        if (s1 < 1.0) {
            if (loglevel >= 4 && !mypid_) {
                printf("\t  dampStep(): current trial step accepted and soln converged retnTrial = %d, its = %d, damp = %g\n",
                       0, m+1, m_fdamp);
            }
            return 0;
        }
        if (loglevel >= 4 && !mypid_) {
            printf("\t  dampStep(): current direction is rejected! retnTrial = %d, its = %d, damp = %g\n", -2, m+1, m_fdamp);
        }
        return -2;
    }
}
//===================================================================================================================================
/*
 * On entry, step0 must contain an undamped Newton step for the solution x0. This method attempts to find a damping coefficient
 * such that the residual is sufficiently reduced.
 *
 * @return   1 Successful step was taken: Next step was less than previous step.
 *                                        s1 is calculated
 *           2 Successful step: Next step's norm is less than 0.8
 *           3 Success:  The final residual is less than 1.0
 *                        A predicted deltaSoln is not produced however. s1 is estimated.
 *           4 Success:  The final residual is less than the residual
 *                       from the previous step.
 *                        A predicted deltaSoln is not produced however. s1 is estimated.
 *           0 Uncertain Success: s1 is about the same as s0
 *          -2 Unsuccessful step.
 */
int
SolNonlinear::dampStep_alt(double time_curr,  const Epetra_Vector_Ghosted& y0, const Epetra_Vector_Ghosted* ydot0_ptr,
                           double& s1, int& loglevel, int& num_backtracks)
{
    Epetra_Vector_Owned& step0 = *m_stp;
    //Epetra_Vector_Owned& step1 = *m_step_2;
    int retnTrial = 0;

    // Compute the weighted norm of the undamped step size step0
    double s0 = m_normSolnFRaw;
    string ResS;

    double rdelta_t = 0.0;
    if (delta_t_n > 1.0E-300) {
        rdelta_t = 1.0 / delta_t_n;
    }
    //--------------------------------------------
    //           Attempt damped step
    //--------------------------------------------

    // damping coefficient starts at 1.0
    m_fdamp = 1.0;
    int m;
    double ff;
    num_backtracks = 0;
    bool raccepted = false;
    for (m = 0; m < NDAMP; m++) {

        ff = m_fdamp * m_fbound;

        // step the solution by the damped step size, updating the time derivative
        /*
         *  New solution is put into m_y_new and  m_ydot_new
         */
        updateSoln(y0, ydot0_ptr, m_fdamp, step0);

        /*
         *  Compute the next residual that would result if m_y_new[] were accepted.
         *  This is computed and storred in the vector m_resid[].
         *  The norm of this new residual is storred in m_normResidTrial and in m_normResidFRaw if
         *  this is the first iteration and therefore represents the raw Newton step.
         */
        get_res(time_curr, rdelta_t, m_y_new, m_ydot_new);

        if (m_print_flag >= 6) {
            ResS = "Residual For Damping Trial " + ZZCantera::int2str(m) + " with Damping Coeff " + ZZCantera::fp2str(m_fdamp);
            m_normResidTrial = res_error_norm(*m_resid, ResS.c_str(), 10);
        } else if (m_print_flag == 4 || m_print_flag == 5) {
            ResS = "Residual For Damp Trial " + ZZCantera::int2str(m) + " with Damping Coeff " + ZZCantera::fp2str(m_fdamp);
            m_normResidTrial = res_error_norm(*m_resid, ResS.c_str(), true);
        } else {
            m_normResidTrial = res_error_norm(*m_resid);
        }

        if (m == 0) {
            m_normResidFRaw = m_normResidTrial;
        }
        //
        // Calculate an estimate of the next solution update
        //
        // s1 = s0 * m_normResidTrial / m_normResid0;
        s1 = s0;
        if (m_normResidTrial / m_normResid0 < 1.0) {
            s1 = s0 * m_normResidTrial / m_normResid0;
        }
        //
        //  We accept the step if the Residual is less than one, or if the residual is less than the initial residual
        //
        double rtest =   m_normResid0 * (0.2 * (1.0 - ff) * (1.0 - ff) * (1.0 - ff) * (1.0 - ff) + 0.8);
        bool steepEnough = (m_normResidTrial < rtest);

        bool atEnd = false;
        if (m_normResidTrial < 1.0) {
            if (s0 < 0.5 && s1 < 0.5) {
                atEnd = true;
            }
        }
        //
        // We allow the residual to go up in certain end cases. Basically the residual has to stay below the cutoff, and
        // the total step has to be a fraction of the step acceptance tolerance.
        // Basically, I've seen the damping step get stuck with residual < 1 and s0 < 1
        //
        if (m_normResidTrial < 0.3) {
            double s0damp = s0 *  m_fdamp;
            double fac = m_normResidTrial / m_normResid0;
            if (fac < 10.0) {
                double fac1 = fac * s0damp;
                if (fac1 < 1.0) {
                    atEnd = true;
                }
            }
        }

        if (atEnd || steepEnough || !doResidSolnDamping_) {
            raccepted = true;
            if (loglevel >= 5 && !mypid_) {
                if (m_normResidTrial < m_normResid0) {
                    printf("\t  dampStep(): Current trial step and damping"
                           " coefficient accepted because resid0 > residTrial < resid0:\n");
                    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
                } else if (m_normResidTrial < 1.0) {
                    printf("\t  dampStep(): Current trial step and damping"
                           " coefficient accepted because residTrial test step < 1:\n");
                    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
                } else {
                    printf("\t  dampStep(): Current trial step and damping"
                           " coefficient accepted because residual solution damping is turned off:\n");
                    printf("\t              resid0 = %g, residTrial = %g\n", m_normResid0, m_normResidTrial);
                }
            }
            /*
             *  We aren't going to solve the system if we don't need to. Therefore, return an estimate
             *  of the next solution update based on the ratio of the residual reduction.
             */
            if (m_normResid0 < 1.0E-10) {
                if (m_normResidTrial < 1.0E-5) {
                    s1 = 0.1 * s0;
                }
            }
            if (m_normResidTrial < 1.0) {
                retnTrial = 3;
            } else {
                retnTrial = 4;
            }
        }

        // xxxxx
        if (loglevel==3) {
            printf("              %11.4E      |                  |%11.4E | %6.2E | %10.4E  | %10.2E %2d |                  %11.4E\n",
                   m_normResid0, s0, m_fbound,  m_normResidFRaw, m_fdamp, m, m_normResidTrial);
        }


        // write log information
        if (loglevel >= 4) {
            print_solnDelta_norm_contrib(step0, "DeltaSoln", step0, "DeltaSolnTrialTest",
                                         "dampNewt: Important Entries for Weighted Soln Updates:", y0, *m_y_new, ff, 5);
        }
        if (loglevel >= 4 && !mypid_) {
            printf("\t\t\tdampStep(): s0 = %g, s1 = %g, damp = %g\n", s0, s1, m_fdamp);
        }



        if (raccepted) {
            if (m_normResidTrial < 1.0) {
                retnTrial = 3;
            } else {
                retnTrial = 4;
            }
            break;
        } else {
            if (s1 > s0) {
                retnTrial = 2;
            } else {
                retnTrial = 1;
            }
        }
        if (raccepted) {
            break;
        }

        if (loglevel >= 4 && !mypid_) {
            printf("\t\t\tdampStep(): current step rejected: (ResNorm = %g ResNorm0 = %g s0 = %g)", m_normResidTrial,  m_normResid0,
                   s0);
            if (m < (NDAMP - 1)) {
                printf(" Decreasing damping factor and retrying");
            } else {
                printf(" Giving up!!!");
            }
            printf("\n");
        }

        num_backtracks++;
        m_fdamp /= DampFactor;
    }

    // If a damping coefficient was found, return 1 if the
    // solution after stepping by the damped step would represent
    // a converged solution, and return 0 otherwise. If no damping
    // coefficient could be found, return -2.
    if (m < NDAMP) {
        if (loglevel >= 4 && !mypid_) {
            printf("\t  dampStep(): current trial step accepted retnTrial = %d, its = %d, damp = %g\n", retnTrial, m+1, m_fdamp);
        }
        return retnTrial;
    } else {
        if (s1 < 0.5 && (s0 < 0.5)) {
            if (loglevel >= 4 && !mypid_) {
                printf("\t  dampStep(): current trial step accepted kindof retnTrial = %d, its = %d, damp = %g\n", 2, m+1, m_fdamp);
            }
            return 2;
        }
        if (s1 < 1.0) {
            if (loglevel >= 4 && !mypid_) {
                printf("\t  dampStep(): current trial step accepted and soln converged retnTrial = %d, its = %d, damp = %g\n",
                       0, m+1, m_fdamp);
            }
            return 0;
        }
        if (loglevel >= 4 && !mypid_) {
            printf("\t  dampStep(): current direction is rejected! retnTrial = %d, its = %d, damp = %g\n", -2, m+1, m_fdamp);
        }
        return -2;
    }
}
//===================================================================================================================================
//  Update the solution vector using the step change that was just computed.
/*
 *    We update the solution vector and the solution dot vector (this is the time derivative vector),
 *    putting the answer into the fixed location,
 *       m_y_new[]   and m_ydot_new[]
 *    given the previous solution vector,
 *          y0[]     and ydot0_ptr[]
 *    and the update step vector with a damping factor
 *          ff           step_1[]
 *
 *   @param  y0         INPUT     Input solution vector     - Ghosted Epectra_Vector reference
 *   @param  ydot0_ptr  INPUT     Input solution dot vector - Ghosted Epectra_Vector ptr
 *   @param  ff         INPUT     Damping factor - double
 *   @param  step_1     INPUT     Input step vector         - Ghosted Epectra_Vector reference
 *
 *    OUTPUT
 * ----------------
 *    This routine changes
 *      m_y_new         OUTPUT     New solution vector     - Ghosted Epectra_Vector ptr
 *      m_ydot_new      OUTPUT     New solution dot vector - Ghosted Epectra_Vector ptr
 *
 *   DISCUSSION
 * ----------------
 *
 *    This routine will update the locally owned values. Then, it will call updateGhostEqns()
 *    for both the m_y_new and  m_ydot_new vectors to update the ghost unknowns on processors.
 *
 *    For the  DAESystemInitial_Solve problem this routine will scatter the step vector unknowns
 *    into the correct locations in the m_y_new and *m_ydot_new vectors according to whether the
 *    the value of (*m_isAlgebraic)[j] is equal to one or not. For DAE unknowns, theoretically
 *    it doesn't matter what the value of the time derivative is, since it doesn't enter into
 *    the equation set. However, here we set (*m_ydot_new)[j] = 0.0 for DAE unknowns to avoid
 *    doing nothing with the entry.
 *
 */
void SolNonlinear::updateSoln(const Epetra_Vector_Ghosted& y0, const Epetra_Vector_Ghosted* ydot0_ptr,
                              double ff, const Epetra_Vector_Ghosted& step_1)
{
    int j;
    if (solnType_ != DAESystemInitial_Solve) {

        // step the solution by the damped step size
        /*
         * Whenever we update the solution, we must also always
         * update the time derivative.
         */
        for (j = 0; j < m_NumLcOwnedEqns; j++) {
            (*m_y_new)[j] = y0[j] + ff * step_1[j];
        }

        m_func->updateGhostEqns(m_y_new, m_y_new_owned);
        if (solnType_ != SteadyState_Solve) {
            calc_ydot(m_order, *m_y_new, *m_ydot_new);
            m_func->updateGhostEqns(m_ydot_new, m_ydot_new_owned);
        }
    } else {
        for (j = 0; j < m_NumLcOwnedEqns; j++) {
            if ((*m_isAlgebraic)[j] == 1) {
                (*m_y_new)[j] = y0[j] + ff * step_1[j];
                // Theoretically, it shouldn't matter what the DAE time derivatives are
                (*m_ydot_new)[j] = 0.0;
            } else {
                (*m_ydot_new)[j] = (*ydot0_ptr)[j] + ff * step_1[j];
                (*m_y_new)[j] = y0[j];
            }
        }
        m_func->updateGhostEqns(m_y_new, m_y_new_owned);
        m_func->updateGhostEqns(m_ydot_new, m_ydot_new_owned);
    }
}
//=====================================================================================================================
/*
 * setColumnScales():
 *
 * Set the column scaling vector at the current time
 */
void
SolNonlinear::setDefaultColumnScaleVector()
{
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*m_colScales)[i] = (*m_ewt)[i];
    }
}
//=====================================================================================================================
// Return the column scales
/*
 *   Note, if there are no column scaling, then 1's are returned in the vector.
 *
 * @param colScales
 */
bool
SolNonlinear::getColumnScaleVector(Epetra_Vector_Owned& colScales) const
{
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        colScales[i] = (*m_colScales)[i];
    }
    return m_colScaling;
}
//=====================================================================================================================
// Set the column scales
/*
 *   Note, if there are no column scaling, then 1's are returned in the vector.
 *
 * @param colScales
 */
void
SolNonlinear::setColumnScaleVector(const Epetra_Vector_Owned& colScales)
{
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*m_colScales)[i] = colScales[i];
    }
    if (colScaleUpdateFrequency_ < 2) {
        if (m_print_flag > 1 && !mypid_) {
            printf("setColumnScaleVector: Warning: colScaleUpdateFrequency_ is set low!\n");
            exit(-1);
        }
    }
}

//=====================================================================================================================
// Setup the problem for solution.
void
SolNonlinear::setup_problem(Solve_Type_Enum solnType, const Epetra_Vector_Ghosted* const y_init,
                            const Epetra_Vector_Ghosted* const ydot_init, double time_curr,
                            ProblemResidEval& problem, EpetraJac& jac)
{
    tdjac_ptr = &jac;
    Epetra_VbrMatrix* A = jac.A_;
    m_func = &problem;
    Comm_ptr_ = &A->Comm();
    mypid_ = Comm_ptr_->MyPID();

    /*
     *   In the problems that we solve, the range map, the domain map, and the row map
     *   are all equal to one another.  This is the number of owned nodes or equations.
     *   The column map is different and larger than the others. The column map includes the
     *   ghosted nodes and equations.
     */
    const Epetra_BlockMap& rangeMap = A->RangeMap();
    const Epetra_BlockMap& domainMap = A->DomainMap();
    const Epetra_BlockMap& colMap = A->ColMap();
    const Epetra_BlockMap& rowMap = A->RowMap();

    int rangeMyNum = rangeMap.NumMyPoints();
    int rangeGbNum = rangeMap.NumGlobalPoints();

    int domainMyNum = domainMap.NumMyPoints();
    int domainGbNum = domainMap.NumGlobalPoints();

    int colMyNum = colMap.NumMyPoints();

    int rowMyNum = rowMap.NumMyPoints();
    int rowGbNum = rowMap.NumGlobalPoints();

    // Fill in some size arrays. These are obtained from the matrix size
    m_NumLcEqns = colMyNum;
    m_NumLcOwnedEqns = rangeMyNum;
    m_NumGbEqns = rangeGbNum;
    // Check some things that should be true
    AssertTrace(m_NumLcOwnedEqns == rangeMyNum);
    AssertTrace(m_NumLcOwnedEqns == domainMyNum);
    AssertTrace(m_NumLcOwnedEqns == rowMyNum);
    AssertTrace(m_NumGbEqns == rangeGbNum);
    AssertTrace(m_NumGbEqns == domainGbNum);
    AssertTrace(m_NumGbEqns == rowGbNum);

    // Check that the solution vector is ghosted.
    AssertTrace(y_init->MyLength() == m_NumLcEqns);
    if (ydot_init) {
        AssertTrace(ydot_init->MyLength() == m_NumLcEqns);
    }

    // Check block maps
    if (!colMap.SameAs(y_init->Map())) {
        printf("error");
        exit(-1);
    }
    if (!colMap.SameAs(ydot_init->Map())) {
        printf("error");
        exit(-1);
    }
    if (!rangeMap.SameAs(domainMap)) {
        printf("error");
        exit(-1);
    }
    if (!rangeMap.SameAs(rowMap)) {
        printf("error");
        exit(-1);
    }

    if (!rangeMap.SameAs(*(problem.LI_ptr_->GbBlockNodeEqnstoOwnedLcBlockNodeEqnsRowMap))) {
        printf("error");
        exit(-1);
    }

    // Malloc a bunch of internal memory

    safeDelete(m_y_curr);
    m_y_curr = new Epetra_Vector(*y_init);

    double* V;
    m_y_curr->ExtractView(&V);
    Epetra_DataAccess eee = View;
    safeDelete(m_y_curr_owned);
    m_y_curr_owned = new Epetra_Vector(eee, rangeMap, V);

    safeDelete(m_y_new);
    m_y_new = new Epetra_Vector(*y_init);

    m_y_new->ExtractView(&V);
    safeDelete(m_y_new_owned);
    m_y_new_owned = new Epetra_Vector_Owned(eee, rangeMap, V);

    safeDelete(m_y_nm1);
    m_y_nm1 = new Epetra_Vector_Ghosted(*y_init);

    /*
     *  We seed the predicted vector here with ones.
     *  The predicted vector is used to set error tolerances. Therefore, we are
     *  setting the expected values of all variables to one, until further notice.
     */
    safeDelete(m_y_pred_n);
    m_y_pred_n = new Epetra_Vector_Ghosted(*y_init);
    m_y_pred_n->PutScalar(1.0);

    safeDelete(m_ydot_curr);
    m_ydot_curr = new Epetra_Vector(*y_init);

    m_ydot_curr->ExtractView(&V);
    safeDelete(m_ydot_curr_owned);
    m_ydot_curr_owned = new Epetra_Vector_Owned(eee, rangeMap, V);

    safeDelete(m_ydot_new);
    m_ydot_new = new Epetra_Vector(*y_init);

    m_ydot_new->ExtractView(&V);
    safeDelete(m_ydot_new_owned);
    m_ydot_new_owned = new Epetra_Vector_Owned(eee, rangeMap, V);

    safeDelete(m_ydot_nm1);
    m_ydot_nm1 = new Epetra_Vector(*y_init);

    safeDelete(m_resid);
    m_resid = new Epetra_Vector(domainMap);

    safeDelete(m_rhs);
    m_rhs = new Epetra_Vector(domainMap);

    safeDelete(m_residWts);
    m_residWts = new Epetra_Vector(domainMap);
    m_residWts->PutScalar(1.0);

    safeDelete(m_stp);
    m_stp = new Epetra_Vector(domainMap);

    safeDelete(m_step_2);
    m_step_2 = new Epetra_Vector(domainMap);

    safeDelete(m_rowScales);
    m_rowScales = new Epetra_Vector(domainMap);
    m_rowScales->PutScalar(1.0);

    safeDelete(m_colScales);
    m_colScales = new Epetra_Vector(domainMap);
    m_colScales->PutScalar(1.0);

    safeDelete(m_ewt);
    m_ewt = new Epetra_Vector(domainMap);

    safeDelete(m_abstol);
    m_abstol = new Epetra_Vector(domainMap);
    m_abstol->PutScalar(1.0E-9);

    safeDelete(m_ewt_deltaDamping);
    m_ewt_deltaDamping = new Epetra_Vector(domainMap);

    /*
     *  Set up the delta damping to be the same size as
     */
    safeDelete(m_absTol_deltaDamping);
    m_absTol_deltaDamping = new Epetra_Vector(domainMap);
    m_absTol_deltaDamping->PutScalar(1.0E-9);

    safeDelete(m_isAlgebraic);
    m_isAlgebraic = new Epetra_IntVector(domainMap);


    safeDelete(m_isArithmeticScaled);
    m_isArithmeticScaled = new Epetra_IntVector(domainMap);
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*m_isArithmeticScaled)[i] = 0;
    }


    problem.fillIsAlgebraic(*m_isAlgebraic);

    problem.fillIsArithmeticScaled(*m_isArithmeticScaled);

    // Initialize a bunch of counters
    solnType_ = solnType;
    m_nfe = 0;
    m_nJacEval = 0;
    m_numTotalNewtIts = 0;
    m_numTotalLinearSolves = 0;
    m_numTotalConvFails = 0;
    m_numTotalTruncFails = 0;
    num_failures = 0;
    m_maxAge = 5;
    m_elapsed = 0.0;
}
//==================================================================================================================================
void
SolNonlinear::setPredicted_soln(const Epetra_Vector& y_pred)
{
    mdpUtil::mdp_copy_dbl_1(&(*m_y_pred_n)[0], &(y_pred[0]), m_NumLcEqns);
}
//==================================================================================================================================
/*
 * solve_nonlinear_problem():
 *
 *   Find the solution to F(X) = 0 by damped Newton iteration.  On
 *   entry, x0 contains an initial estimate of the solution.  On
 *   successful return, x1 contains the converged solution.
 *
 * SolnType = TRANSIENT -> we will assume we are relaxing a transient
 *        equation system for now. Will make it more general later,
 *        if an application comes up.
 *          = DAEINIT Solve the initial conditions problem
 */
int
SolNonlinear::solve_nonlinear_problem(Solve_Type_Enum solnType,
                                      Epetra_Vector_Ghosted* y_comm,
                                      Epetra_Vector_Ghosted* ydot_comm,
                                      double CJ,
                                      double time_curr,
                                      int& num_newt_its_comm,
                                      int& num_linear_solves,
                                      int& num_backtracks)
{
    ZZCantera::clockWC t0;
    if (solnType == SteadyState_Solve) {
        doTimeDependentResid_ = false;
    } else {
        doTimeDependentResid_ = true;
        //delta_t_n = timeStep;
        if (delta_t_n <= 0) {
            throw m1d_Error("SolNonlinear::solve_nonlinear_problem", "delta_t_n is not greater than zero");
        }
    }
    int m = 0;
    bool forceNewJac = true;
    double s1 = 1.e30;
    int retn = 0;
    solnType_ = solnType;

    mdpUtil::mdp_copy_dbl_1(&(*m_y_curr)[0], &(*y_comm)[0], m_NumLcEqns);
    if (solnType != SteadyState_Solve) {
        mdpUtil::mdp_copy_dbl_1(&(*m_ydot_curr)[0], &(*ydot_comm)[0], m_NumLcEqns);
    }
    setDefaultSolnWeights();
    /*
     *  Set the default column scales
     */
    if (m_colScaling && colScaleUpdateFrequency_ >= 1) {
        setDefaultColumnScaleVector();
    }

#ifdef DEBUG_MODE_NOT
    print_SVector("Input Soln Vector in nonlinear solver", *m_y_curr);
    print_SVector("atol in nonlinear solver", *m_abstol);
    print_SVector("atol_deltaDamping in nonlinear solver", *m_absTol_deltaDamping);
    print_SVector("Weighting Vector in nonlinear solver", *m_ewt);
    print_SVector("Weighting Vector for deltaDamping in nonlinear solver", *m_ewt_deltaDamping);
    print_IntSVector("Arithmetic Scaling id", *m_isArithmeticScaled);
    print_IntSVector("Algebraic Unknowns", *m_isAlgebraic);
#endif

    m_num_newt_its = 0;
    num_linear_solves = -m_numTotalLinearSolves;
    num_backtracks = 0;
    int i_backtracks;
    int convRes = 0;

    if (m_print_flag == 2 || m_print_flag == 3) {
        if (!mypid_) {
            if (solnType_ == DAESystemInitial_Solve) {
                printf("\tSolve_Nonlinear_Problem: DAESystemInitial\n\n");
            } else {
                printf("\tSolve_Nonlinear_Problem:\n\n");
            }
            printf("\t    Iter Resid0 NewJac |  LinearIts  Ax-b |DeltaSolnRaw|  FBound  |ResidAftBound| Fdamp DampIts |   DeltaSolnFinal  ResidFinal \n");
            printf("\t--------------------------------------------------------------------------------------------------------------------------------\n");
        }
    }

    while (1 > 0) {

        /*
         * Increment Newton Solve counter
         */
        m_numTotalNewtIts++;
        m_num_newt_its++;
        num_newt_its_comm = m_num_newt_its;
        if (m_print_flag >= 4 && !mypid_) {
            printf("\t");
            print_line("=", 90);
            printf("\tSolve_Nonlinear_Problem: iteration %d:\n", m_num_newt_its);
            printf("\t");
            print_line("=", 90);
        }

        /*
         * Make sure the ghost unknowns are up to date
         */
        m_func->updateGhostEqns(m_y_curr, m_y_curr_owned);
        // May need to do this for m_ydot_curr as well, do to new problem

        /*
         *  Either go get a new jacobian and residual
         *  or go get a new residual and rely on the old jacobian to solve
         *  the new right hand side
         */
        if (forceNewJac) {
            if (m_print_flag > 3 && !mypid_) {
                printf("\t  Getting a new Jacobian and solving system\n");
            }
            get_jac(*tdjac_ptr, m_resid, doTimeDependentResid_, time_curr, CJ, m_y_curr, m_ydot_curr);
            m_jacAge = 0;
        } else {
            get_res(time_curr, CJ, m_y_curr, m_ydot_curr);
            if (m_print_flag > 3 && !mypid_) {
                printf("\t  Solving system with old jacobian\n");
            }
            m_jacAge++;
        }

        /*
         * Scale the matrix and the rhs, if they aren't already scaled
         * Figure out and store the residual scaling factors.
         */
        scaleMatrix(*m_stp, *m_y_curr, *m_ydot_curr, time_curr, CJ, m_print_flag);

        /*
         *  Optional print out the initial residual
         */
        if (m_print_flag >= 6) {
            m_normResid0 = res_error_norm(*m_resid, "Initial norm of the residual", 10);
        } else if (m_print_flag == 4 || m_print_flag == 5) {
            m_normResid0 = res_error_norm(*m_resid, "Initial norm of the residual", true);
        } else {
            m_normResid0 = res_error_norm(*m_resid, "Initial norm of the residual");
        }

        // Compute the undamped Newton step
        /*
         * On return, the full step length is returned in the member vector m_stp
         */
        doNewtonSolve(*m_stp, *m_y_curr, *m_ydot_curr, time_curr, CJ, m_print_flag);

        if (m_print_flag > 3) {
            m_normSolnFRaw = soln_error_norm(*m_stp, 10, "Initial Undamped Step of the iteration");
        } else {
            m_normSolnFRaw = soln_error_norm(*m_stp, 0, "Initial Undamped Step of the iteration");
        }

        /*
         *  Check against the hard bounds in the problem
         *   retn =  1 : No bounds was violated fbound = 1.0;
         *   retn =  0 : Bounds was hit, fbound < 1.0; Don't consider the calculation as converged
         *   retn = -3 : We are at the bounds and the step size has become so short that we
         *               are heading in the forbidden zone. End the calculation.
         */
        /* Hewson suggests that when the hard bounds limit was previously
         * very small, we want to increase it more rapidly.  We get very
         * small hard bounds limits when we have variables with absolute
         * values close to teh absolute tolerances.
         * Changed from m_fbound *= 2.0;
         * to m_fbound *= 2.0 * ( 1 + fabs( log10 ( m_fbound ) ) );
         * This is just an attempt and will be adjusted in doHardBounds() below.
         */
        if (m_fbound < 1.0) {
            m_fbound *= 2.0 * (1 + fabs(log10(m_fbound)));
            if (m_fbound > 1.0) {
                m_fbound = 1.0;
            }
        }
        //
        //  On return m_stp is reduced by the factor, m_fbound
        //
        retn = doHardBounds(*m_y_curr, *m_stp, m_fbound);
        if (retn < 0) {
            m = retn;
            goto done;
        }

        // Damp the Newton step
        /*
         *  On return the recommended new solution and derivative is located in:
         *          m_y_new
         *          m_y_dot_new
         *  The estimate of the solution update norm for the next step is located in
         *          s1
         * @return   1 Successful step was taken: Next step was less than previous step.
         *                                        s1 is calculated
         *           2 Successful step: Next step's norm is less than 0.8
         *           3 Success:  The final residual is less than 1.0
         *                        A predicted deltaSoln is not produced however. s1 is estimated.
         *           4 Success:  The final residual is less than the residual
         *                       from the previous step.
         *                        A predicted deltaSoln is not produced however. s1 is estimated.
         *           0 Uncertain Success: s1 is about the same as s0
         *          -2 Unsuccessful step.
         */
//#define DOOLD
#ifdef DOOLD
        m = dampStep(time_curr, *m_y_curr, m_ydot_curr, s1, m_print_flag, i_backtracks);
#else
        m = dampStep_alt(time_curr, *m_y_curr, m_ydot_curr, s1, m_print_flag, i_backtracks);
#endif
        num_backtracks += i_backtracks;
        /*
         *  If we had to damp the step because of a hard bounds, then we can not be considered as converged
         */
        if (retn == 0) {
            if (m > 0) {
                //	if (m_print_flag > 2 && !mypid_) {
                //  printf("\t     ... Bounds damping was carried out before ... \n");
                // }
                m = 0;
            }
        }


        /*
         * Impose the minimum number of newton iterations critera
         */
        bool minNewItsCondition = false;
        if (m_num_newt_its < m_min_newt_its) {
            if (m == 3 || m == 2) {
                minNewItsCondition = true;
                if (m_print_flag > 3 && !mypid_) {
                    printf("\t   ... Minimum newton iterations not attained ...\n");
                }
            }
            if (m > 0) {
                m = 0;
            }
        }
        /*
         * Impose max newton iteration
         */
        if (m_num_newt_its > maxNewtIts_) {
            m = -1;
            if (m_print_flag > 1 && !mypid_) {
                printf("\t      ... Damped Newton unsuccessful (max newts exceeded) final soln_norm = %-12.4E\n", s1);
            }
        }

        if (m_print_flag > 3) {
            string ResS = "Accepted Step " + ZZCantera::int2str(m_num_newt_its) + " Residual Norm";
            res_error_norm(*m_resid, ResS.c_str(), 10);
        }

        convRes = 0;
        if (m > 0) {
            convRes = convergenceCheck(m, s1);
        }

        if (m_print_flag >= 4 && !mypid_) {
            if (convRes > 0) {
                printf("\t  Damped Newton iteration successful, nonlin "
                       "converged, final estimate of the next solution update norm = %-12.4E\n", s1);
            } else if (m == 0) {
                printf("\t  Damped Newton iteration successful, "
                       "final estimate of the next solution update norm = %-12.4E\n", s1);
            } else {
                printf("\t  Damped Newton unsuccessful, final estimate of the next solution update norm = %-12.4E\n", s1);
            }
        }

        // If we are converged, then let's use the best solution possible
        // for an end result. We did a resolve in dampStep(). Let's update
        // the solution to reflect that.
        // HKM 5/16 -> Took this out, since if the last step was a
        //             damped step, then adding stp1[j] is undamped, and
        //             may lead to oscillations. It kind of defeats the
        //             purpose of dampStep() anyway.
        // if (m == 1) {
        //  for (int j = 0; j < m_neq; j++) {
        //   y_new[j] += stp1[j];
        //                HKM setting intermediate y's to zero was a tossup.
        //                    slightly different, equivalent results
        // #ifdef DEBUG_HKM
        //	      y_new[j] = MAX(0.0, y_new[j]);
        // #endif
        //  }
        // }

        bool m_filterIntermediate = false;
        if (m_filterIntermediate) {
            // if (m == 0) {
            //    (void) filterNewStep(time_n, *y_new, *ydot_new);
            // }
        }
        // Write new solution into the curr solutions
        if (m >= 0) {
            if (solnType_ == DAESystemInitial_Solve) {
                mdpUtil::mdp_copy_dbl_1(&(*m_y_curr)[0], &(*m_y_new)[0], m_NumLcEqns);
                mdpUtil::mdp_copy_dbl_1(&(*m_ydot_curr)[0], &(*m_ydot_new)[0], m_NumLcEqns);
            } else {
                mdpUtil::mdp_copy_dbl_1(&(*m_y_curr)[0], &(*m_y_new)[0], m_NumLcEqns);
                if (solnType_ != SteadyState_Solve) {
                    calc_ydot(m_order, *m_y_curr, *m_ydot_curr);
                }
            }
        }
        //
        // Print out the Table values
        //
        if (m_print_flag == 2 || m_print_flag == 3) {
            if (!mypid_) {
                //printf("\t    Iter Resid0 NewJac |  LinearIts  Ax-b  | Fbound | ResidBound | Fdamp DampIts |   DeltaSolnF     ResidFinal\n");
                printf("\t%4d  %11.4E", m_num_newt_its, m_normResid0);
                if (!m_jacAge) {
                    printf("  Y   ");
                } else {
                    printf("  N   ");
                }
                // xxxxx
                printf("|%5d %11.4E | %10.4E |", m_curr_linearIts, m_curr_normLin,  m_normSolnFRaw);
                printf("%9.2E |", m_fbound);
                printf(" %10.4E  ", m_normResidFRaw);
                printf("| %10.2E %2d |    %11.4E   %11.4E", m_fdamp, i_backtracks, m_normSolnFRaw * m_fbound * m_fdamp,
                       m_normResidTrial);
                printf("\n");
                if (minNewItsCondition) {
                    if (m_print_flag > 2) {
                        printf("\t     ... Minimum newton iterations not attained ...\n");
                    }
                }
            }
        }

        // convergence
        if (convRes > 0) {
            goto done;
        }
        // If dampStep fails, first try a new Jacobian if an old
        // one was being used. If it was a new Jacobian, then
        // return -1 to signify failure.

        else if (m < 0) {
            goto done;
        }

    }
    //-------------------------------------------------------------------------------------
done:
    s1 = s1;

    if (m_print_flag == 2 || m_print_flag == 3) {
        if (convRes > 0) {
            if (convRes == 3) {
                printf("\t                       |                  |(%10.4E)|          |             | converged = 3 |    (%10.4E) \n",
                       s1, s1);
            } else {
                printf("\t                       |                  |(%10.4E)|          |             | converged = %1d |     %10.4E %11.4E \n",
                       s1, convRes,
                       s1, m_normResidTrial);
            }
        }
        if (!mypid_) {
            printf("\t  --------------------------------------------------------------------------------------------\n");
        }
    }
    /*
     *  we come here at the end of the calculation
     *  if m >0 then we have a successful calculation
     *
     *   m = -1  : We ran out of the allowed newton iterations. Nothing happened
     *             that was a show stopper. We just ran out of iterations
     *             Return the current solution as the best case.
     *   m = -2
     *   m = -3  : We ran into a hard bounds problem. This is a show stopper, because
     *             it is a hard bounds. This usually indicates that the either the nonlinear
     *             solution path veers into prohibited territory before settling down into
     *             valid territory, or, it may mean that there is sufficient round off error
     *             from ill-conditioning that the step goes into prohibited territory.
     */
    if (m >= -1) {
        mdpUtil::mdp_copy_dbl_1(&((*y_comm)[0]), &((*m_y_curr)[0]), m_NumLcEqns);
        if (solnType_ != SteadyState_Solve) {
            mdpUtil::mdp_copy_dbl_1(&(*ydot_comm)[0], &(*m_ydot_curr)[0], m_NumLcEqns);
        }
    }

    /*
     *  If we are successful, save the old solution as a prediction for
     *  use the next time the solver is called.
     */
    if (m == 1) {
        mdpUtil::mdp_copy_dbl_1(&((*m_y_pred_n)[0]), &((*m_y_curr)[0]), m_NumLcEqns);
    }

    num_linear_solves += m_numTotalLinearSolves;

    //  double time_elapsed = t0.secondsWC();
    if (m_print_flag > 1 && !mypid_) {
        if (m > 0) {
            printf("\t  Nonlinear problem solved successfully in %d iterations\n", m_num_newt_its);
            // printf("\t\tNonlinear problem solved successfully in "
            //  "%d its, time elapsed = %g sec\n", num_newt_its, time_elapsed);
        }
    }
    num_newt_its_comm = m_num_newt_its;

#ifdef DEBUG_HKM_NOT
    double rdelta_t = 0.0;
    if (delta_t_n > 1.0E-300) {
        rdelta_t = 1.0 / delta_t_n;
    }
    get_res(time_curr, rdelta_t, m_y_new, m_ydot_new);
    m_normResid0 = res_error_norm(*m_resid, "Final SolNonlinear norm of the residual", 10);
#endif

    return m;
}
//=====================================================================================================================
/* Function to calculate the acceleration vector ydot for the first or
 * second order predictor/corrector time integrator.  This routine can be
 * called by a first order - forward Euler / backward Euler predictor /
 * corrector or for a second order Adams - Bashforth / Trapezoidal Rule
 * predictor / corrector.  See Nachos documentation Sand86-1816 and Gresho,
 * Lee, Sani LLNL report UCRL - 83282 for more information.
 *
 *  variables:
 *
 *    on input:
 *
 *       N          - number of local unknowns on the processor
 *                    This is equal to internal plus border unknowns.
 *       order      - indicates order of method
 *                    = 1 -> first order forward Euler/backward Euler
 *                           predictor/corrector
 *                    = 2 -> second order Adams-Bashforth/Trapezoidal Rule
 *                           predictor/corrector
 *
 *      delta_t_n   - Magnitude of the current time step at time n
 *                    (i.e., = t_n - t_n-1)
 *      y_curr[]    - Current Solution vector at time n
 *      y_nm1[]     - Solution vector at time n-1
 *      ydot_nm1[] - Acceleration vector at time n-1
 *
 *   on output:
 *
 *      ydot_curr[]   - Current acceleration vector at time n
 *
 * Note we use the current attribute to denote the possibility that
 * y_curr[] may not be equal to m_y_n[] during the nonlinear solve
 * because we may be using a look-ahead scheme.
 */
void
SolNonlinear::calc_ydot(int order, const Epetra_Vector& y_curr, Epetra_Vector& ydot_curr)
{
    int i;
    double c1;
    const Epetra_Vector& y_nm1 = *m_y_nm1;
    switch (order) {
    case 0:
    case 1: /* First order forward Euler/backward Euler */
        c1 = 1.0 / delta_t_n;
        for (i = 0; i < m_NumLcEqns; i++) {
            ydot_curr[i] = c1 * (y_curr[i] - y_nm1[i]);
        }
        return;
    case 2: /* Second order Adams-Bashforth / Trapezoidal Rule */
        c1 = 2.0 / delta_t_n;
        for (i = 0; i < m_NumLcEqns; i++) {
            ydot_curr[i] = c1 * (y_curr[i] - y_nm1[i]) - (*m_ydot_nm1)[i];
        }
        return;
    }
}
//===================================================================================================================================

void
SolNonlinear::print_solnDelta_norm_contrib(const Epetra_Vector& solnDelta0, const char* const s0,
                                           const Epetra_Vector& solnDelta1, const char* const s1,
                                           const char* const title,
                                           const Epetra_Vector& y0, const Epetra_Vector& y1,
                                           double damp, int num_entries)
{
    int i, j, jnum;
    bool used;
    double dmax0, dmax1, error, rel_norm;
    stream0 ss;
    std::vector<int> imax(num_entries, -1);

    print0_sync_start(false, ss, *Comm_ptr_);
    if (!mypid_) {
        printf("\t\t%s currentDamp = %g\n", title, damp);
        printf("\t\t      I   ysoln     %10s  ysolnTrial %10s weight RelDeltaSoln RelDeltaSolnTrial\n", s0, s1);
        printf("\t\t   ");
        print_line("-", 100);
    }
    print0_sync_end(false, ss, *Comm_ptr_);

    for (jnum = 0; jnum < num_entries; jnum++) {
        int idLocalEqnMax = -1;
        dmax1 = -1.0;
        imax[jnum] = -1;
        for (i = 0; i < m_NumLcOwnedEqns; i++) {
            used = false;
            for (j = 0; j < jnum; j++) {
                if (imax[j] == i) {
                    used = true;
                }
            }
            if (!used) {
                error = solnDelta0[i] / (*m_ewt)[i];
                rel_norm = sqrt(error * error);
                error = solnDelta1[i] / (*m_ewt)[i];
                rel_norm += sqrt(error * error);
                if (rel_norm > dmax1) {
                    idLocalEqnMax = i;
                    dmax1 = rel_norm;
                }
            }
        }
        double flocal = dmax1;
        double gmax1 = -1.0;
        int idGlobalEqnMax = -1;
        if (idLocalEqnMax > 0) {
            idGlobalEqnMax = m_func->LI_ptr_->IndexGbEqns_LcEqns[idLocalEqnMax];
        }
        int iproc = procChoice_Max_Brcst1Int(flocal, *Comm_ptr_, gmax1, idGlobalEqnMax);
        if (iproc == mypid_) {
            imax[jnum] = idGlobalEqnMax;
        }
        print0_sync_start(false, ss, *Comm_ptr_);
        if (imax[jnum] >= 0) {
            i = idLocalEqnMax;
            error = solnDelta0[i] / (*m_ewt)[i];
            dmax0 = sqrt(error * error);
            error = solnDelta1[i] / (*m_ewt)[i];
            dmax1 = sqrt(error * error);
            ss.print0("\t\t   %4d %12.4e %12.4e %12.4e  %12.4e "
                      "%12.4e %12.4e %12.4e\n", idGlobalEqnMax, y0[i], solnDelta0[i], y1[i], solnDelta1[i], (*m_ewt)[i], dmax0, dmax1);
        }
        print0_sync_end(false, ss, *Comm_ptr_);
    }
    print0_sync_start(false, ss, *Comm_ptr_);
    if (!mypid_) {
        printf("\t\t   ");
        print_line("-", 100);
    }
    print0_sync_end(false, ss, *Comm_ptr_);
}
//=====================================================================================================================
// Check to see if the nonlinear problem has converged
/*
 *
 * @return integer is returned. If positive, then the problem has converged
 *           1 Successful step was taken: Next step's norm is less than 1.0.
 *                                        The final residual norm is less than 1.0.
 *           2 Successful step: Next step's norm is less than 0.8.
 *                              This step's norm is less than 1.0.
 *                              The residual norm can be anything.
 *           3 Success:  The final residual is less than 1.0
 *                        The predicted deltaSoln is below 1.0.
 *           4 Success:  The final residual is less than 1.0
 *                        The predicted deltaSoln is above 1.0.
 *           0 Not converged yet
 */
int
SolNonlinear::convergenceCheck(int dampCode, double s1)
{
    int retn = 0;
    if (m_fbound < 0.9999) {
        return retn;
    }
    if (m_fdamp < 0.9999) {
        return retn;
    }
    if (dampCode <= 0) {
        return retn;
    }
    if (dampCode == 3) {
        if (s1 < 1.0) {
            if (m_normResidTrial < 1.0) {
                return 3;
            }
        }
    }
    if (dampCode == 4) {
        if (s1 < 1.0) {
            if (m_normResidTrial < 1.0) {
                return 3;
            }
        }
    }
    if (s1 < 0.8) {
        if (m_normSolnFRaw < 1.0) {
            return 2;
        }
    }
    if (s1 < 1.0) {
        if (m_normSolnFRaw < 1.0) {
            return 1;
        }
    }
    return retn;
}
//==================================================================================================================================
void
SolNonlinear::setTolerances(double reltol, int n, const double* const abstol)
{
    if (!m_abstol) {
        m_abstol = new Epetra_Vector(m_y_curr->Map());
    }
    if (n != m_NumLcEqns) {
        printf("ERROR n is wrong\n");
        exit(-1);
    }
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*m_abstol)[i] = abstol[i];
    }
    m_reltol = reltol;
}
//==================================================================================================================================
void
SolNonlinear::setTolerances_deltaDamping(double reltol_dd, int n, const double* const abstol_dd)
{
    if (!m_absTol_deltaDamping) {
        m_absTol_deltaDamping = new Epetra_Vector(m_y_curr->Map());
    }
    if (n != m_NumLcOwnedEqns) {
        printf("ERROR n is wrong\n");
        exit(-1);
    }
    if (abstol_dd) {
        for (int i = 0; i < m_NumLcOwnedEqns; i++) {
            (*m_absTol_deltaDamping)[i] = abstol_dd[i];
        }
    } else {
        for (int i = 0; i < m_NumLcOwnedEqns; i++) {
            (*m_absTol_deltaDamping)[i] = (*m_abstol)[i];
        }
    }
}
//==================================================================================================================================
// Set the value of the maximum # of newton iterations
/*
 *  @param maxNewtIts   Maximum number of newton iterations
 */
void SolNonlinear:: setMaxNewtIts(const int maxNewtIts)
{
    maxNewtIts_ = maxNewtIts;
}
//==================================================================================================================================
void
SolNonlinear::setProblemType(int jacFormMethod)
{
    m_jacFormMethod = jacFormMethod;
}
//=====================================================================================================================
void
SolNonlinear::setDefaultSolnWeights()
{
    Epetra_Vector_Owned& ewt = *m_ewt;
    Epetra_Vector_Owned& ewt_deltaDamping = *m_ewt_deltaDamping;
    Epetra_IntVector& isS = *m_isArithmeticScaled;
    if (solnType_ == DAESystemInitial_Solve) {
        Epetra_IntVector& isA = *m_isAlgebraic;


        for (int i = 0; i < m_NumLcOwnedEqns; i++) {
            if (isA[i] == 1) {
                if (isS[i] == 1) {
                    ewt[i] =  1.0E4 * m_reltol * (*m_abstol)[i];
                    ewt_deltaDamping[i] = (*m_absTol_deltaDamping)[i];
                } else {
                    ewt[i] = (*m_abstol)[i] + m_reltol * fabs((*m_y_curr)[i]);
                    ewt_deltaDamping[i] = (*m_absTol_deltaDamping)[i] + m_reltol * fabs((*m_y_curr)[i]);
                }
            } else {
                if (isS[i] == 1) {
                    ewt[i] =  1.0E4 * m_reltol * (*m_abstol)[i];
                    ewt_deltaDamping[i] = (*m_absTol_deltaDamping)[i];
                } else {
                    ewt[i] = (*m_abstol)[i] + m_reltol * fabs((*m_ydot_curr)[i]);
                    ewt_deltaDamping[i] = (*m_absTol_deltaDamping)[i] + m_reltol * fabs((*m_ydot_curr)[i]);
                }
            }
        }
    } else {
        for (int i = 0; i < m_NumLcOwnedEqns; i++) {
            if (isS[i] == 1) {
                ewt[i] = 1.0E4 * m_reltol * (*m_abstol)[i];
                ewt_deltaDamping[i] = (*m_absTol_deltaDamping)[i];
            } else {
                ewt[i] = (*m_abstol)[i] + m_reltol * 0.5 * (fabs((*m_y_curr)[i]) + fabs((*m_y_pred_n)[i]));
                ewt_deltaDamping[i] =
                    (*m_absTol_deltaDamping)[i] + m_reltol * 0.5 * (fabs((*m_y_curr)[i]) + fabs((*m_y_pred_n)[i]));
            }
        }
    }
}
//=====================================================================================================================

//=====================================================================================================================
/*
 *
 * @param dumpJacobians Dump jacobians to disk.
 *
 *                   default = false
 */
void
SolNonlinear::setPrintSolnOptions(bool dumpJacobians)
{
    m_dumpJacobians = dumpJacobians;
}
//=====================================================================================================================
void
SolNonlinear::setNonLinOptions(int min_newt_its, bool matrixConditioning, bool colScaling, bool rowScaling,
                               int colScaleUpdateFrequency)
{
    m_min_newt_its = min_newt_its;
    m_matrixConditioning = matrixConditioning;
    setColScaling(colScaling, colScaleUpdateFrequency);
    setRowScaling(rowScaling);
}
//=====================================================================================================================
// set the previous time step
/*
 *
 * @param jac
 */
void
SolNonlinear::setPreviousTimeStep(const double timeStep_comm, const Epetra_Vector& y_nm1, const Epetra_Vector& ydot_nm1)
{
    delta_t_n = timeStep_comm;
    mdpUtil::mdp_copy_dbl_1(&(*m_y_nm1)[0], &(y_nm1[0]), m_NumLcEqns);
    mdpUtil::mdp_copy_dbl_1(&(*m_ydot_nm1)[0], &(ydot_nm1[0]), m_NumLcEqns);
}
//===================================================================================================================================
// Compute the Residual Weights
/*
 *  The residual weights are defined here to be equal to the inverse of the row scaling factors used to
 *  row scale the matrix, after column scaling is used. They are multiplied by 10-3 because the column
 *  weights are also multiplied by that same quantity.
 *
 *  The basic idea is that a change in the solution vector on the order of the convergence tolerance
 *  multiplied by  [RJC] which is of order one after row scaling should give you the relative weight
 *  of the row. Values of the residual for that row can then be normalized by the value of this weight.
 *  When the tolerance in delta x is achieved, the tolerance in the residual is also achieved.
 */
void
SolNonlinear::computeResidWts()
{
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*m_residWts)[i] = 1.0 / (*m_rowScales)[i];
    }
    if (m_print_flag >= 3 && m_num_newt_its > 1) {
        printf("           ... Computing New residual weights ...\n");
    }
}
//==================================================================================================================================
void
SolNonlinear::getResidWts(Epetra_Vector_Owned& residWts)
{
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        residWts[i] = (*m_residWts)[i];
    }
}
//====================================================================================================================
void
SolNonlinear::print_SVector(std::string header, const Epetra_MultiVector& v) const
{
    if (!mypid_) {
        printf("%s\n", header.c_str());
    }

    Epetra_VbrMatrix* A = (*tdjac_ptr).A_;
    const Epetra_BlockMap& domainMap = A->DomainMap();
    Epetra_Vector_Owned* vv = new Epetra_Vector(domainMap);

    for (int i = 0; i <  m_NumLcOwnedEqns; i++) {
        double* a = v[0];
        (*vv)[i] = a[i];
    }


    m_func->showSolutionVector(header, time_n, delta_t_n, *vv);

    delete vv;

// stream0 ss;
//Print0_epMultiVector(ss, v);
}
//====================================================================================================================
void
SolNonlinear::print_IntSVector(std::string header, const Epetra_IntVector& v) const
{
    if (!mypid_) {
        printf("%s\n", header.c_str());
    }
    m_func->showSolutionIntVector(header, time_n, delta_t_n, v);

// stream0 ss;
//Print0_epMultiVector(ss, v);
}
//=====================================================================================================================
void
SolNonlinear::setRowScaling(bool onoff)
{
    m_rowScaling = onoff;
}
//=====================================================================================================================
void
SolNonlinear::setColScaling(bool onoff, int colScaleUpdateFrequency)
{
    m_colScaling = onoff;
    colScaleUpdateFrequency_ = colScaleUpdateFrequency;
}
//=====================================================================================================================
// Set the toggles for solution damping
/*
 *
 * @param residSolnDamping
 * @param deltaDamping
 * @param highLowDamping
 */
void
SolNonlinear::setDampingToggles(const bool residSolnDamping, const bool deltaDamping, const bool highLowDamping)
{
    doResidSolnDamping_ = residSolnDamping;
    doDeltaDamping_ = deltaDamping;
    doHighLowDamping_ = highLowDamping;
}
//=====================================================================================================================
// Set the vectors for lower and upper boundaries.
/*
 * @param lowBounds
 * @param highBounds
 */
void
SolNonlinear::setSolutionBounds(const Epetra_Vector_Owned& lowBounds, const Epetra_Vector_Owned& highBounds)
{
    doHighLowDamping_ = true;
    if (!solnLowBound_) {
        solnLowBound_ = new Epetra_Vector(*m_stp);
        solnHighBound_ = new Epetra_Vector(*m_stp);
    }
    for (int i = 0; i < m_NumLcOwnedEqns; i++) {
        (*solnLowBound_)[i] = lowBounds[i];
        (*solnHighBound_)[i] = highBounds[i];
    }
}
//=====================================================================================================================
}

