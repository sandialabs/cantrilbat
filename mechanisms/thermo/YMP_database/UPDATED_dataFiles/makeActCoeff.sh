#!/bin/sh
#
#
# Copywrite 2004 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
# retains certain rights in this software.
# See file License.txt for licensing information.
#
# $Id: makeActCoeff.sh,v 1.1 2009/01/12 21:12:09 hkmoffa Exp $
#
#

if test !  -e beta_cmx.xml
then

  ln -s  data0.ypf.R2_betas_ca_rev2.xml beta_cmx.xml
fi

if test ! -e  lambda_nj.xml 
then
  ln -s  data0.ypf.R2_extract_Pitzer_params_lambdas_rev3.xml lambda_nj.xml 
fi


if test ! -e    psi_cca_caa.xml
then
  ln -s data0.ypf.R2_extract_Pitzer_params_psi_rev3.xml  psi_cca_caa.xml
fi


if test ! -e  theta_cc_aa.xml
then
  ln -s data0.ypf.R2_Thetas_cc_aa_rev2.xml theta_cc_aa.xml
fi


if test ! -e zeta_nca.xml 
then
  ln -s data0.ypf.R2_extract_Pitzer_params_zetas_rev2.xml zeta_nca.xml 
fi




echo >actCoeffSection.xml


cat  >actCoeffSection.xml <<+
<activityCoefficients model="Pitzer" TempModel="complex1">
   <A_Debye model="water" />

+

cat beta_cmx.xml | cat >>actCoeffSection.xml 

cat theta_cc_aa.xml | cat >>actCoeffSection.xml

cat psi_cca_caa.xml | cat >>actCoeffSection.xml

cat lambda_nj.xml | cat >>actCoeffSection.xml

cat zeta_nca.xml | cat >> actCoeffSection.xml


cat  >>actCoeffSection.xml <<+

</activityCoefficients> 
+

