#!/bin/bash
# Script to transform surfaces from native space EPI position
# Olivia Viessmann and Jon Polimeni, March 2018
set -e
#parse input arguments
Subject="$1"               #subject ID
orig="$2"                  #orig.mgz -> surfaces are reconstructed from here
ACPC_T1="$3"               #ac-pc aligned T1, this has a different tkr than the orig.mgz
EPI_dc_jac="$4"            #EPI, undistorted reference image
fMRI2str_fsl="$5"          #fsl .mat file that transforms EPI to ac-pc aligned T1
Surf="$6"                  #Path to surface
hemi="$7"   
NameOffMRI="$8"            #name of fMRI run, so that we can name the transformed surfaces accordingly
WriteTransformsToHere="$9" #location to write the .lta transforms to 
#The surfaces are reconstructed from the orig.mgz. 
#The orig.mgz is the conformed version of the input T1 (ac-pc aligned T1) 
#and has a different orientation and is flipped (the third dimension is A>>P 
#as if it is a coronal scan, the ac-pc is H>>F as if it is axial). 
#The fsl registration matrix, fMRI2str_fsl.mat, is between the ac-pc T1 
#and the undistorted EPI reference image. 
#So we cannot easily transform the surfaces using this transformation, 
#we first have to transform from orig to ac-pc and concatenate this with 
#the fMRI2str_fsl. Here we use lta_convert
#Use lta to do transformations
#Name of the output registration from EPI to ac-pc T1
lta_EPI2ACPC=${WriteTransformsToHere}/${NameOffMRI}2str.lta
lta_EPI2orig=${WriteTransformsToHere}/${NameOffMRI}2orig.lta
lta_ACPC2orig=${WriteTransformsToHere}/ACPC2orig.lta
#convert fsl to .lta
lta_convert --subject ${Subject} --infsl ${fMRI2str_fsl} --src ${EPI_dc_jac} --trg ${ACPC_T1} --outlta ${lta_EPI2ACPC}

#we create and identity lta file for the registration from AC-PC T1 to orig 
lta_convert --subject ${Subject} --inlta identity.nofile --src ${ACPC_T1} --trg ${orig} --outlta ${lta_ACPC2orig}

#Concatenate the two (out_type flag :1 - RAS2RAS, 0 - VOX2VOX (default))
mri_concatenate_lta -out_type 1 ${lta_EPI2ACPC} ${lta_ACPC2orig} ${lta_EPI2orig} 

#Then call mris_surf2surf to register the surface
/usr/local/freesurfer/stable6/bin/mri_surf2surf \
    --s ${Subject} \
    --reg ${lta_EPI2orig} \
    --sval-xyz ${Surf} \
    --tval-xyz ${EPI_dc_jac} \
    --hemi ${hemi} \
    --surfreg ${Surf} \
    --tval ${hemi}.${Surf}_2$NameOffMRI
echo "END: Transformed surface to EPI position"
