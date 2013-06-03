
//
// beta - metal_Li_LiCl_electron
//
#include <stdio.h>
#include <math.h>

int main () {
 // Li(liquid) - up to 700K
  double aLi1 =    32.46972;
  double bLi1 =   -2.635975;
  double cLi1 =   -6.327128;
  double dLi1 =    4.230359;
  double eLi1 =    0.005686;
  double fLi1 =   -7.117319;
  double gLi1 =    74.29361;
  
  //  Li(l) 700 - 2700
  double aLi2 =    26.00896;
  double bLi2 =    5.632375;
  double cLi2 =   -4.013227;
  double dLi2 =    0.873686;
  double eLi2 =    0.344150;
  double fLi2 =   -4.199690;
  double gLi2 =    66.36284;

  // LiCl at all temperatures (strictly greater than 700)
  double aLiSi1 =   73.18025;
  double bLiSi1 =  -9.047232;
  double cLiSi1 =   -0.316390;
  double dLiSi1 =  0.079587;
  double eLiSi1 =    0.013594;
  double fLiSi1 =   -417.1314;
  double gLiSi1 =   157.6711;


  double aLiSi2 = aLiSi1   ;
  double bLiSi2 = bLiSi1   ;
  double cLiSi2 = cLiSi1 ;
  double dLiSi2 = dLiSi1   ;
  double eLiSi2 = eLiSi1    ;
  double fLiSi2 = fLiSi1  ;
  double gLiSi2 = gLiSi1 ;




  double a =  -aLiSi1  + aLi1;

  double b = - bLiSi1  + bLi1;

  double c = - cLiSi1  + cLi1;

  double d = - dLiSi1  + dLi1;

  double e = - eLiSi1  + eLi1;

  double f = - fLiSi1  + fLi1;

  double g = - gLiSi1  + gLi1;

  // f += sdiff * (-4.481143E4 / 1000.);

  g += 8.314472 * log(2.0);

  double a2 = - aLiSi2  + aLi2;

  double b2 = - bLiSi2  +  bLi2;

  double c2 = - cLiSi2  +  cLi2;

  double d2 = - dLiSi2  +  dLi2;

  double e2 = - eLiSi2  + eLi2;

  double f2 = - fLiSi2  + fLi2;

  double g2 = - gLiSi2  + gLi2;

 
  g2 +=  8.314472 * log(2.0);

  //======================================================================================================
  bool printMain = true;

  printf(" <?xml version=\"1.0\"?>\n");
  printf(" <ctml>\n");
  printf("   <validate reactions=\"yes\" species=\"yes\"/>\n");
  printf("    \n");

  if (printMain) {
    printf("   <!-- phase metal_Li_LiCl_electrons    -->\n");
    printf("   <phase dim=\"3\" id=\"metal_Li_LiCl_electrons\">\n");
    printf("       <!--      -->\n");
     
    printf("      <elementArray datasrc=\"elements.xml\">\n");
    printf("        E \n");
    printf("      </elementArray>\n");
    printf("      <speciesArray datasrc=\"#species_electrode\"> electron_Li_LiCl </speciesArray>\n");
    printf("      <thermo model=\"StoichSubstance\">\n");
    printf("      <density units=\"g/cm3\"> 100. </density>\n");
    printf("     </thermo>\n");
    printf("     <transport model=\"None\"/>\n");
    printf("     <kinetics model=\"none\"/>\n");
    printf("   </phase>\n");
  }
  printf("           \n");


  printf("          \n");

  printf("   <!-- species definitions     qq-->\n");
  printf("   <speciesData id=\"species_electrode\">\n");
  printf("        \n");

  printf("      <!-- species electron_Li_LiCl -->\n");
  printf("      <! --   temperature 400 to 2700 K -->\n");
  printf("      <! --   stoichiometric quantity in electrode rxn = %20.16g -->\n", log(2.0) * 8.314472);
  printf("      <species name=\"electron_Li_LiCl\">\n");
  printf("       <atomArray> E:1.0  </atomArray>\n");
  printf("       <charge> -1 </charge>\n"); 
  printf("       <thermo>\n");
  printf("        <Shomate Pref=\"1 bar\" Tmax=\"700.0\" Tmin=\"400.0\">\n");
  printf("           <floatArray size=\"7\">\n");
  printf("               %-11.8g,\n", a);
  printf("               %-11.8g,\n", b);
  printf("               %-11.8g,\n", c);
  printf("               %-11.8g,\n", d);
  printf("               %-11.8g,\n", e);
  printf("               %-11.8g,\n", f);
  printf("               %-11.8g \n", g);
  printf("           </floatArray>\n");
  printf("        </Shomate>\n");
  printf("        <Shomate Pref=\"1 bar\" Tmax=\"2700.0\" Tmin=\"700.0\">\n");
  printf("           <floatArray size=\"7\">\n");
  printf("               %-11.8g,\n", a2);
  printf("               %-11.8g,\n", b2);
  printf("               %-11.8g,\n", c2);
  printf("               %-11.8g,\n", d2);
  printf("               %-11.8g,\n", e2);
  printf("               %-11.8g,\n", f2);
  printf("               %-11.8g \n", g2);
  printf("           </floatArray>\n");
  printf("        </Shomate>\n");
  printf("       </thermo>\n");
  printf("        <density units=\"g/cm3\"> 100. </density>\n");
  printf("      </species>\n");
  printf("                \n");


  printf("                \n");
  printf("    </speciesData>\n");
  printf("        \n");
  printf("  </ctml>\n");
  printf("        \n");

  return 0;

}


