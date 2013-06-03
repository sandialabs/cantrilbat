open(DAT, "cantera_SPEQ06_HKFT_params_rev1.csv") || die("Could not open file!");
@raw_data=<DAT>;
close(DAT);

foreach $HKFT_Params (@raw_data)
{
	 chop($HKFT_Params);
	  ($species,$atomArray_species,$charge,$thermoModel,$HKFT_Pref,$HKFT_Tmax,$HKFT_Tmin,
		  $DG0_f_Pr_Tr_Units,$DG0_f_Pr_Tr_Value,$DH0_f_Pr_Tr_Units,$DH0_f_Pr_Tr_Value,$S0_Pr_Tr_Units,$S0_Pr_Tr_Value,
		  $a1_Param_Units,$a1_Param_Value,$a2_Param_Units,$a2_Param_Value,$a3_Param_Units,$a3_Param_Value,$a4_Param_Units,$a4_Param_Value,
		  $c1_Param_Units,$c1_Param_Value,$c2_Param_Units,$c2_Param_Value,$omega_Pr_Tr_Units,$omega_Pr_Tr_Value,$Source)=split(/\,/,$HKFT_Params);
	   print "               <species name=\"$species\">
	                           <atomArray> $atomArray_species </atomArray>
	                           <charge> $charge </charge>
	                           <thermo model=\"$thermoModel\"> 
                                     <HKFT Pref=\"$HKFT_Pref\" Tmax=\"$HKFT_Tmax\" Tmin=\"$HKFT_Tmin\"> 
				       <DG0_f_Pr_Tr units=\"$DG0_f_Pr_Tr_Units\"> $DG0_f_Pr_Tr_Value </DG0_f_Pr_Tr>
				       <DH0_f_Pr_Tr units=\"$DH0_f_Pr_Tr_Units\"> $DH0_f_Pr_Tr_Value </DH0_f_Pr_Tr>
				       <S0_Pr_Tr units=\"$S0_Pr_Tr_Units\"> $S0_Pr_Tr_Value </S0_Pr_Tr>
				     </HKFT>
				   </thermo>
				   <standardState model=\"$thermoModel\">
				      <a1 units=\"$a1_Param_Units\"> $a1_Param_Value </a1>
				      <a2 units=\"$a2_Param_Units\"> $a2_Param_Value </a2>
				      <a3 units=\"$a3_Param_Units\"> $a3_Param_Value </a3>
				      <a4 units=\"$a4_Param_Units\"> $a4_Param_Value </a4>
				      <c1 units=\"$c1_Param_Units\"> $c1_Param_Value </c1>
				      <c2 units=\"$c2_Param_Units\"> $c2_Param_Value </c2>
				      <omega_Pr_Tr units=\"$omega_Pr_Tr_Units\"> $omega_Pr_Tr_Value </omega_Pr_Tr>
				   </standardState>
				   <source>
				       $Source
				   </source>
				 </species> \n";
           print "\n";
      }
