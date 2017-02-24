/**
 * @file m1d_SDD_AnodeCollector.cpp
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this
 * work by or on behalf of the U.S. Government. Export of this program
 * may require a license from the United States Government.
 */

#include "m1d_SDD_AnodeCollector.h"
#include "m1d_SurDomain_AnodeCollector.h"

#include "m1d_ProblemStatementCell.h"
#include "m1d_BC_Battery.h"
#include "m1d_CanteraElectrodeGlobals.h"

//=====================================================================================================================
namespace m1d
{
//=====================================================================================================================
SDD_AnodeCollector::SDD_AnodeCollector(DomainLayout *dl_ptr, int pos, const char *domainName) :
    SDD_Mixed(dl_ptr,domainName),
    m_position(pos)
{
    /*
     * For this implementation, we are not adding a separate equation for this surface domain.
     * The current equation for the electrode phase from the bulk domain is used as the current in the anode
     * collector.
     *
     * This object will set an anode voltage of zero on that equation.
     */
    voltageVarBCType_ = PSCinput_ptr->anodeBCType_;
    anodeTempBCType_ = PSCinput_ptr->anodeTempBCType_;
    anodeTempCollector_ = PSCinput_ptr->anodeTempRef_;
    anodeHeatTranCoeff_ =  PSCinput_ptr->anodeHeatTranCoeff_;
    anodeCCThickness_ = PSCinput_ptr->anodeCCThickness_;

}
//=====================================================================================================================
SDD_AnodeCollector::SDD_AnodeCollector(const SDD_AnodeCollector &r) :
    SDD_Mixed(r.DL_ptr_),
    m_position(0)
{
    SDD_AnodeCollector::operator=(r);
}
//=====================================================================================================================
SDD_AnodeCollector::~SDD_AnodeCollector()
{
}
//=====================================================================================================================
SDD_AnodeCollector &
SDD_AnodeCollector::operator=(const SDD_AnodeCollector &r)
{
    if (this == &r) {
	return *this;
    }
    SDD_Mixed::operator=(r);
    m_position = r.m_position;
    voltageVarBCType_   = r.voltageVarBCType_;
    anodeTempBCType_    = r.anodeTempBCType_;
    anodeTempCollector_ = r.anodeTempCollector_;
    anodeHeatTranCoeff_ = r.anodeHeatTranCoeff_;
    anodeCCThickness_   = r.anodeCCThickness_;

    return *this;
}
//=====================================================================================================================
// Set the equation description
/*
 *  This routine is responsible for setting the variables:
 *    - NumEquationsPerNode
 *    - VariableNameList
 *    - EquationNameList
 *    - EquationIndexStart_EqName
 */
void
SDD_AnodeCollector::SetEquationDescription()
{
   /*
    * Set the policy for connecting bulk domains
    * This really isn't set yet.
    */
   setRLMapping(0);
   /*
    * Fill in the rest of the information
    */
   SurfDomainDescription::SetEquationDescription();

 
   EqnType e1(Current_Conservation, 1, "Anode Current Conservation");
   VarType v1(Voltage, 1, "AnodeVoltage");


   if (voltageVarBCType_ == 0) {
       /*
	*  Add the Dirichlet condition onto the bulk domain equation for the electrode current.
	*      Set the solid voltage at the anode to zero
	*/
       addDirichletCondition(e1, v1, 0.0);
   } else {
       if (anodeCCThickness_ > 0.0) {

	   BoundaryCondition* BC_timeDep = new BC_anodeCC(anodeCCThickness_ , 0.0);
	   addRobinCondition(e1, v1, BC_timeDep);
       } else {
	   
       }
   }


   /*
    *  All of the other boundary conditions default to zero flux at the interface
    *      This includes:
    *
    *           axial velocity
    *           flux of Li+
    *           flux of K+
    *           flux of Cl-
    *           flux of current
    */

   /*
    *  All of the other boundary conditions default to zero flux at the interface
    *      This includes:
    *
    *           flux of Li+
    *           flux of K+
    *           flux of Cl-
    *           flux of current
    *
    *  Because of the staggered grid, a dirichlet condition must be set on the axial velocity. Or else the
    *  last axial velocity unknown is not represented in the solution vector, leading to a singular matrix.
    *
    */
   EqnType e2(Continuity, 0, "Continuity: Bulk Velocity");
   VarType v2(Velocity_Axial, 0, "Axial velocity");
   addFluxCondition(e2, v2, 0.0);

   /*
    *  Temperature boundary condition
    */
   EqnType et(Enthalpy_Conservation, 0, "");
   VarType vt(Temperature, 0, "temperature");
   if (anodeTempBCType_ != -1) {
       if (anodeTempBCType_ == 0) {
	   double tref =  PSCinput_ptr->anodeTempRef_;
	   addDirichletCondition(et, vt, tref);
       }
       if (anodeTempBCType_ == 10) {
           double ht = PSCinput_ptr->anodeHeatTranCoeff_;
	   double tref =  PSCinput_ptr->anodeTempRef_;
           double area = PSCinput_ptr->crossSectionalArea_;
	   BoundaryCondition* BC_timeDep = new BC_heatTransfer(ht, tref, area);
	   addRobinCondition(et, vt, BC_timeDep);
       }
   }

}
//=====================================================================================================================
SurDomain1D *
SDD_AnodeCollector::mallocDomain1D()
{
  SurDomain_AnodeCollector * s1d = new SurDomain_AnodeCollector(*this, 1);
  SurDomain1DPtr_ = s1d;
  return s1d;
}
//=====================================================================================================================
} /* End of Namespace */
//=====================================================================================================================

