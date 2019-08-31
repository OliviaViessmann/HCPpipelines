#!/bin/bash 
# Script to create intermittent surfaces between pial and white
# Olivia Viessmann, March 2018

#Wed 12 Sep 2018 10:25:47 AM EDT
# OV: edited script to further construct surfaces below the white surface

set -e


# parse arguments
inputfolder="$1"
spring="$2"
outputfolder="$3"


#****** Reconstruct intermittent surfaces******
#
#go to folder with surface
cd ${inputfolder}

hemis=(lh rh)

for hemi in ${hemis[*]}; do
echo $hemi
 	mris_curvature -W ${hemi}.white
    	mris_curvature -W ${hemi}.pial
done

hemis=(lh rh)
for hemi in ${hemis[*]}; do
   	mris_expand -a 1 -s ${spring} -n 10 -thickness ${hemi}.white 1 ${outputfolder}/${hemi}.midgray
        #also expand to negative surface distances below the white surface
        mris_expand -a 1 -s ${spring} -n 10 -thickness ${hemi}.white -1 ${outputfolder}/${hemi}.subwhite 
done

rename midgray0 midgray. *midgray0*
rename subwhite0 subwhite. *subwhite0*

#Calculate the principal and gaussian curvature of each surface 
 for hemi in ${hemis[*]}; do
  	for depth in `seq -f %02.0f 0 10`; do
         	mris_curvature -w ${hemi}.midgray.${depth}
		mris_curvature -w ${hemi}.subwhite.${depth}
     	done
 done
 

echo "END: created surfaces"
