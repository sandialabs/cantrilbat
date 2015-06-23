#!/bin
#
#
#  Name of the program to run
#
PROGRAM=LiIon_PorousBat
#
#  Program input file
#
PROGRAM_INPUT=LiIon_PorousBat.inp
#
#  Any other program options
#
PROGRAM_OPTS=''
#
#
BLESSED_DATA_FILES=" good_out.txt timeDep_IV_blessed.dat  BulkDomain1D_0_blessed.dat BulkDomain1D_1_blessed.dat  BulkDomain1D_2_blessed.dat SurDomain1D_0_blessed.dat SurDomain1D_1_blessed.dat  SurDomain1D_3_blessed.dat" 
DATA_FILES="         out.txt      timeDep_IV.dat          BulkDomain1D_0.dat         BulkDomain1D_1.dat          BulkDomain1D_2.dat         SurDomain1D_0.dat         SurDomain1D_1.dat          SurDomain1D_3.dat "
DIFF_NAMES="         diff_out.txt diff_IV.txt             diff_b0.txt                diff_b1.txt                 diff_b2.txt                diff_s0.txt               diff_s1.txt                diff_s3.txt "
DIFF_REQ="           True         True                    True                       True                        True                       True                      True                       True "

#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
