#!/bin/bash
#
# This shows how to run the planar reconstuction
# using the command line tool. 

PATH=$PATH:../bin

#run the planar example:
CDI_reconstruction.exe planar_example.config

#########################################################

#to run fresnel reconstruction:
#reconstruct the white-field then the sample
#uncomment the lines below.

#CDI_reconstruction.exe fresnel_example.config fresnel_wf
#CDI_reconstruction.exe fresnel_example.config fresnel

#########################################################

#to run multiple times with a different starting seed do
#something like the following:

#for a in `seq 10`
#do
#  CDI_reconstruction.exe planar_example.config planar $a &> log_$a
#  mv planar.cplx planar_result_${a}_.cplx
#done  

#a tool for merging the output can then be used (not written yet).

