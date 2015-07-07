#!/bin
#    runtest_Electrode input file:
#
#  Name of the program to run
#
PROGRAM=Fes2_4_MPNewman_cc
#
#  Program input file
#
PROGRAM_INPUT=
#
#  Any other program options
#
PROGRAM_OPTS=''
#
#
BLESSED_DATA_FILES=" good_out.txt FeS2_intResults_blessed.csv"
DATA_FILES="         out.txt      FeS2_intResults.csv        "
DIFF_NAMES="         diff_out.txt diff_csv.txt             "
DIFF_REQ="           True         True                    " 

#
#  Extra files to be removed before the test starts
#
EXTRA_WHACKED_FILES=''
