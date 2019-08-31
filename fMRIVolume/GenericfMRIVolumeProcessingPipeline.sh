#!/bin/bash 
# Derivative of the GenericfMRIVolumeProcessingPipeline.sh
# This script follows the minimal preprocessing of fmri data as in the HCP, but does not warp to MNI, but to the subject's native T1w space instead
# Edited by Olivia Viessmann, March 2018
# requires OneStepResampling2Native.sh
set -e

# Requirements for this script
#  installed versions of: FSL (version 5.0.6), FreeSurfer (version 5.3.0-HCP) , gradunwarp (HCP version 1.0.2) 
#  environment: use SetUpHCPPipeline.sh  (or individually set FSLDIR, FREESURFER_HOME, HCPPIPEDIR, PATH - for gradient_unwarp.py)

########################################## PIPELINE OVERVIEW ########################################## 

# TODO

########################################## OUTPUT DIRECTORIES ########################################## 

# TODO

# --------------------------------------------------------------------------------
#  Load Function Libraries
# --------------------------------------------------------------------------------

source $HCPPIPEDIR/global/scripts/log.shlib  # Logging related functions
source $HCPPIPEDIR/global/scripts/opts.shlib # Command line option functions

################################################ SUPPORT FUNCTIONS ##################################################

# --------------------------------------------------------------------------------
#  Usage Description Function
# --------------------------------------------------------------------------------

show_usage() {
    echo "Usage information To Be Written"
    exit 1
}

# --------------------------------------------------------------------------------
#   Establish tool name for logging
# --------------------------------------------------------------------------------
#log_SetToolName "GenericfMRIVolumeProcessingPipeline.sh"

################################################## OPTION PARSING #####################################################

opts_ShowVersionIfRequested $@

if opts_CheckForHelpRequest $@; then
    show_usage
fi

log_Msg "Parsing Command Line Options"

# parse arguments
Path=`opts_GetOpt1 "--path" $@`
log_Msg "Path: ${Path}"

Subject=`opts_GetOpt1 "--subject" $@`
log_Msg "Subject: ${Subject}"

NameOffMRI=`opts_GetOpt1 "--fmriname" $@`
log_Msg "NameOffMRI: ${NameOffMRI}"

fMRITimeSeries=`opts_GetOpt1 "--fmritcs" $@`
log_Msg "fMRITimeSeries: ${fMRITimeSeries}"

fMRIScout=`opts_GetOpt1 "--fmriscout" $@`
log_Msg "fMRIScout: ${fMRIScout}"

SpinEchoPhaseEncodeNegative=`opts_GetOpt1 "--SEPhaseNeg" $@`
log_Msg "SpinEchoPhaseEncodeNegative: ${SpinEchoPhaseEncodeNegative}"

SpinEchoPhaseEncodePositive=`opts_GetOpt1 "--SEPhasePos" $@`
log_Msg "SpinEchoPhaseEncodePositive: ${SpinEchoPhaseEncodePositive}"

MagnitudeInputName=`opts_GetOpt1 "--fmapmag" $@`  # Expects 4D volume with two 3D timepoints
log_Msg "MagnitudeInputName: ${MagnitudeInputName}"

PhaseInputName=`opts_GetOpt1 "--fmapphase" $@`  
log_Msg "PhaseInputName: ${PhaseInputName}"

GEB0InputName=`opts_GetOpt1 "--fmapgeneralelectric" $@`
log_Msg "GEB0InputName: ${GEB0InputName}"

DwellTime=`opts_GetOpt1 "--echospacing" $@`  
log_Msg "DwellTime: ${DwellTime}"

deltaTE=`opts_GetOpt1 "--echodiff" $@`  
log_Msg "deltaTE: ${deltaTE}"

UnwarpDir=`opts_GetOpt1 "--unwarpdir" $@`  
log_Msg "UnwarpDir: ${UnwarpDir}"

FinalfMRIResolution=`opts_GetOpt1 "--fmrires" $@`  
log_Msg "FinalfMRIResolution: ${FinalfMRIResolution}"

# FIELDMAP, SiemensFieldMap, GeneralElectricFieldMap, or TOPUP
# Note: FIELDMAP and SiemensFieldMap are equivalent
DistortionCorrection=`opts_GetOpt1 "--dcmethod" $@`
log_Msg "DistortionCorrection: ${DistortionCorrection}"

BiasCorrection=`opts_GetOpt1 "--biascorrection" $@`
# Convert BiasCorrection value to all UPPERCASE (to allow the user the flexibility to use NONE, None, none, legacy, Legacy, etc.)
BiasCorrection="$(echo ${BiasCorrection} | tr '[:lower:]' '[:upper:]')"
log_Msg "BiasCorrection: ${BiasCorrection}"

GradientDistortionCoeffs=`opts_GetOpt1 "--gdcoeffs" $@`  
log_Msg "GradientDistortionCoeffs: ${GradientDistortionCoeffs}"

TopupConfig=`opts_GetOpt1 "--topupconfig" $@`  # NONE if Topup is not being used
log_Msg "TopupConfig: ${TopupConfig}"

dof=`opts_GetOpt1 "--dof" $@`
dof=`opts_DefaultOpt $dof 6`
log_Msg "dof: ${dof}"

RUN=`opts_GetOpt1 "--printcom" $@`  # use ="echo" for just printing everything and not running the commands (default is to run)
log_Msg "RUN: ${RUN}"

#NOTE: the jacobian option only applies the jacobian of the distortion corrections to the fMRI data, and NOT from the nonlinear T1 to template registration
UseJacobian=`opts_GetOpt1 "--usejacobian" $@`
# Convert UseJacobian value to all lowercase (to allow the user the flexibility to use True, true, TRUE, False, False, false, etc.)
UseJacobian="$(echo ${UseJacobian} | tr '[:upper:]' '[:lower:]')"
log_Msg "UseJacobian: ${UseJacobian}"

MotionCorrectionType=`opts_GetOpt1 "--mctype" $@`  # use = "FLIRT" to run FLIRT-based mcflirt_acc.sh, or "MCFLIRT" to run MCFLIRT-based mcflirt.sh
MotionCorrectionType=`opts_DefaultOpt $MotionCorrectionType MCFLIRT` #use mcflirt by default

#error check
case "$MotionCorrectionType" in
    MCFLIRT|FLIRT)
        #nothing
    ;;
    
    *)
		log_Err_Abort "--mctype must be 'MCFLIRT' (default) or 'FLIRT'"
    ;;
esac

JacobianDefault="true"
if [[ $DistortionCorrection != "TOPUP" ]]
then
    #because the measured fieldmap can cause the warpfield to fold over, default to doing nothing about any jacobians
    JacobianDefault="false"
    #warn if the user specified it
    if [[ $UseJacobian == "true" ]]
    then
        log_Msg "WARNING: using --jacobian=true with --dcmethod other than TOPUP is not recommended, as the distortion warpfield is less stable than TOPUP"
    fi
fi
log_Msg "JacobianDefault: ${JacobianDefault}"

UseJacobian=`opts_DefaultOpt $UseJacobian $JacobianDefault`
log_Msg "After taking default value if necessary, UseJacobian: ${UseJacobian}"

if [[ -n $HCPPIPEDEBUG ]]
then
    set -x
fi

#sanity check the jacobian option
if [[ "$UseJacobian" != "true" && "$UseJacobian" != "false" ]]
then
	log_Err_Abort "the --usejacobian option must be 'true' or 'false'"
fi

# Setup PATHS
PipelineScripts=${HCPPIPEDIR_fMRIVol}
GlobalScripts=${HCPPIPEDIR_Global}

#Naming Conventions
T1wImage="T1w_acpc_dc"
T1wRestoreImage="T1w_acpc_dc_restore"
T1wRestoreImageBrain="T1w_acpc_dc_restore_brain"
T1wFolder="T1w" #Location of T1w images
AtlasSpaceFolder="MNINonLinear"
ResultsFolder="Results"
#OV this has not been changed to original version 
BiasField="BiasField_acpc_dc"
BiasFieldMNI="BiasField"
T1wAtlasName="T1w_restore"
MovementRegressor="Movement_Regressors" #No extension, .txt appended
MotionMatrixFolder="MotionMatrices"
MotionMatrixPrefix="MAT_"
FieldMapOutputName="FieldMap"
MagnitudeOutputName="Magnitude"
MagnitudeBrainOutputName="Magnitude_brain"
ScoutName="Scout"
OrigScoutName="${ScoutName}_orig"
OrigTCSName="${NameOffMRI}_orig"
FreeSurferBrainMask="brainmask_fs"
fMRI2strOutputTransform="${NameOffMRI}2str"
RegOutput="Scout2T1w"
AtlasTransform="acpc_dc2standard"
OutputfMRI2StandardTransform="${NameOffMRI}2standard"
Standard2OutputfMRITransform="standard2${NameOffMRI}"
QAImage="T1wMulEPI"
JacobianOut="Jacobian"
SubjectFolder="$Path"/"$Subject"
#note, this file doesn't exist yet, gets created by ComputeSpinEchoBiasField.sh during DistortionCorrectionAnd...
sebasedBiasFieldMNI="$SubjectFolder/$AtlasSpaceFolder/Results/$NameOffMRI/${NameOffMRI}_sebased_bias.nii.gz"

fMRIFolder="$Path"/"$Subject"/"$NameOffMRI"

#error check bias correction opt
case "$BiasCorrection" in
    NONE)
        UseBiasFieldMNI=""
		;;
    LEGACY)
        UseBiasFieldMNI="${fMRIFolder}/${BiasFieldMNI}.${FinalfMRIResolution}"
		;;    
    SEBASED)
        if [[ "$DistortionCorrection" != "TOPUP" ]]
        then
            log_Err_Abort "SEBASED bias correction is only available with --dcmethod=TOPUP"
        fi
        
	UseBiasFieldMNI="$sebasedBiasFieldMNI"

		;;
    "")
        log_Err_Abort "--biascorrection option not specified"
		;;
    *)
        log_Err_Abort "unrecognized value for bias correction: $BiasCorrection"
		;;
esac


########################################## DO WORK ########################################## 

T1wFolder="$Path"/"$Subject"/"$T1wFolder"
AtlasSpaceFolder="$Path"/"$Subject"/"$AtlasSpaceFolder"
#OV: change to T1 folder as our results are all in native space
ResultsFolder="$T1wFolder"/"$ResultsFolder"/"$NameOffMRI"

mkdir -p ${T1wFolder}/Results/${NameOffMRI}

if [ ! -e "$fMRIFolder" ] ; then
  log_Msg "mkdir ${fMRIFolder}"
  mkdir "$fMRIFolder"
fi
cp "$fMRITimeSeries" "$fMRIFolder"/"$OrigTCSName".nii.gz

#Create fake "Scout" if it doesn't exist
if [ $fMRIScout = "NONE" ] ; then
  ${RUN} ${FSLDIR}/bin/fslroi "$fMRIFolder"/"$OrigTCSName" "$fMRIFolder"/"$OrigScoutName" 0 1
else
  cp "$fMRIScout" "$fMRIFolder"/"$OrigScoutName".nii.gz
fi

#Gradient Distortion Correction of fMRI
log_Msg "Gradient Distortion Correction of fMRI"
if [ ! $GradientDistortionCoeffs = "NONE" ] ; then
    log_Msg "mkdir -p ${fMRIFolder}/GradientDistortionUnwarp"
    mkdir -p "$fMRIFolder"/GradientDistortionUnwarp
    ${RUN} "$GlobalScripts"/GradientDistortionUnwarp.sh \
		   --workingdir="$fMRIFolder"/GradientDistortionUnwarp \
		   --coeffs="$GradientDistortionCoeffs" \
		   --in="$fMRIFolder"/"$OrigTCSName" \
		   --out="$fMRIFolder"/"$NameOffMRI"_gdc \
		   --owarp="$fMRIFolder"/"$NameOffMRI"_gdc_warp
	
    log_Msg "mkdir -p ${fMRIFolder}/${ScoutName}_GradientDistortionUnwarp"	
    mkdir -p "$fMRIFolder"/"$ScoutName"_GradientDistortionUnwarp
    ${RUN} "$GlobalScripts"/GradientDistortionUnwarp.sh \
		   --workingdir="$fMRIFolder"/"$ScoutName"_GradientDistortionUnwarp \
		   --coeffs="$GradientDistortionCoeffs" \
		   --in="$fMRIFolder"/"$OrigScoutName" \
		   --out="$fMRIFolder"/"$ScoutName"_gdc \
		   --owarp="$fMRIFolder"/"$ScoutName"_gdc_warp
	
	if [[ $UseJacobian == "true" ]]
	then
	    ${RUN} ${FSLDIR}/bin/fslmaths "$fMRIFolder"/"$NameOffMRI"_gdc -mul "$fMRIFolder"/"$NameOffMRI"_gdc_warp_jacobian "$fMRIFolder"/"$NameOffMRI"_gdc
	    ${RUN} ${FSLDIR}/bin/fslmaths "$fMRIFolder"/"$ScoutName"_gdc -mul "$fMRIFolder"/"$ScoutName"_gdc_warp_jacobian "$fMRIFolder"/"$ScoutName"_gdc
	fi
else
    log_Msg "NOT PERFORMING GRADIENT DISTORTION CORRECTION"
    ${RUN} ${FSLDIR}/bin/imcp "$fMRIFolder"/"$OrigTCSName" "$fMRIFolder"/"$NameOffMRI"_gdc
    ${RUN} ${FSLDIR}/bin/fslroi "$fMRIFolder"/"$NameOffMRI"_gdc "$fMRIFolder"/"$NameOffMRI"_gdc_warp 0 3
    ${RUN} ${FSLDIR}/bin/fslmaths "$fMRIFolder"/"$NameOffMRI"_gdc_warp -mul 0 "$fMRIFolder"/"$NameOffMRI"_gdc_warp
   ${RUN} ${FSLDIR}/bin/imcp "$fMRIFolder"/"$OrigScoutName" "$fMRIFolder"/"$ScoutName"_gdc
    #make fake jacobians of all 1s, for completeness
    ${RUN} ${FSLDIR}/bin/fslmaths "$fMRIFolder"/"$OrigScoutName" -mul 0 -add 1 "$fMRIFolder"/"$ScoutName"_gdc_warp_jacobian
    ${RUN} ${FSLDIR}/bin/fslroi "$fMRIFolder"/"$NameOffMRI"_gdc_warp "$fMRIFolder"/"$NameOffMRI"_gdc_warp_jacobian 0 1
    ${RUN} ${FSLDIR}/bin/fslmaths "$fMRIFolder"/"$NameOffMRI"_gdc_warp_jacobian -mul 0 -add 1 "$fMRIFolder"/"$NameOffMRI"_gdc_warp_jacobian
fi

log_Msg "mkdir -p ${fMRIFolder}/MotionCorrection"
mkdir -p "$fMRIFolder"/MotionCorrection
${RUN} "$PipelineScripts"/MotionCorrection.sh \
       "$fMRIFolder"/MotionCorrection \
       "$fMRIFolder"/"$NameOffMRI"_gdc \
       "$fMRIFolder"/"$ScoutName"_gdc \
       "$fMRIFolder"/"$NameOffMRI"_mc \
       "$fMRIFolder"/"$MovementRegressor" \
       "$fMRIFolder"/"$MotionMatrixFolder" \
       "$MotionMatrixPrefix" \
       "$MotionCorrectionType"

# EPI Distortion Correction and EPI to T1w Registration
log_Msg "EPI Distortion Correction and EPI to T1w Registration"

DCFolderName=DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased
DCFolder=${fMRIFolder}/${DCFolderName}

if [ -e ${DCFolder} ] ; then
    ${RUN} rm -r ${DCFolder}
fi
log_Msg "mkdir -p ${DCFolder}"
mkdir -p ${DCFolder}

${RUN} ${PipelineScripts}/DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased.sh \
       --workingdir=${DCFolder} \
       --scoutin=${fMRIFolder}/${ScoutName}_gdc \
       --t1=${T1wFolder}/${T1wImage} \
       --t1restore=${T1wFolder}/${T1wRestoreImage} \
       --t1brain=${T1wFolder}/${T1wRestoreImageBrain} \
       --fmapmag=${MagnitudeInputName} \
       --fmapphase=${PhaseInputName} \
       --fmapgeneralelectric=${GEB0InputName} \
       --echodiff=${deltaTE} \
       --SEPhaseNeg=${SpinEchoPhaseEncodeNegative} \
       --SEPhasePos=${SpinEchoPhaseEncodePositive} \
       --echospacing=${DwellTime} \
       --unwarpdir=${UnwarpDir} \
       --owarp=${T1wFolder}/xfms/${fMRI2strOutputTransform} \
       --biasfield=${T1wFolder}/${BiasField} \
       --oregim=${fMRIFolder}/${RegOutput} \
       --freesurferfolder=${T1wFolder} \
       --freesurfersubjectid=${Subject} \
       --gdcoeffs=${GradientDistortionCoeffs} \
       --qaimage=${fMRIFolder}/${QAImage} \
       --method=${DistortionCorrection} \
       --topupconfig=${TopupConfig} \
       --ojacobian=${fMRIFolder}/${JacobianOut} \
       --dof=${dof} \
       --fmriname=${NameOffMRI} \
       --subjectfolder=${SubjectFolder} \
       --biascorrection=${BiasCorrection} \
       --usejacobian=${UseJacobian}

#OV
log_Msg "Distortion correction finished"
#OVstart: change to my script
#One Step Resampling (changed --t1, --freesurferbrainmask to T1,--biasfield
log_Msg "One Step Resampling to native space"
log_Msg "mkdir -p ${fMRIFolder}/OneStepResampling"
#OV change bias fiel to bias field in native space below in --biasfield option
mkdir -p ${fMRIFolder}/OneStepResampling
${RUN} ${PipelineScripts}/OneStepResampling2Native.sh \
       --workingdir=${fMRIFolder}/OneStepResampling \
      --infmri=${fMRIFolder}/${OrigTCSName}.nii.gz \
       --t1=${T1wFolder}/${T1wImage} \
       --fmriresout=${FinalfMRIResolution} \
       --fmrifolder=${fMRIFolder} \
       --fmri2structin=${T1wFolder}/xfms/${fMRI2strOutputTransform} \
       --struct2std=${AtlasSpaceFolder}/xfms/${AtlasTransform} \
       --owarp=${T1wFolder}/xfms/${NameOffMRI}2struct \
       --oiwarp=${T1wFolder}/xfms/struct2${NameOffMRI} \
       --motionmatdir=${fMRIFolder}/${MotionMatrixFolder} \
       --motionmatprefix=${MotionMatrixPrefix} \
       --ofmri=${fMRIFolder}/${NameOffMRI}_nonlin \
       --freesurferbrainmask=${T1wFolder}/${FreeSurferBrainMask} \
       --biasfield=${DCFolder}/ComputeSpinEchoBiasField/${NameOffMRI}_sebased_bias \
       --gdfield=${fMRIFolder}/${NameOffMRI}_gdc_warp \
       --scoutin=${fMRIFolder}/${OrigScoutName} \
       --scoutgdcin=${fMRIFolder}/${ScoutName}_gdc \
       --oscout=${fMRIFolder}/${NameOffMRI}_SBRef_nonlin \
       --ojacobian=${fMRIFolder}/${JacobianOut}_Native.${FinalfMRIResolution} #OV changed naming


log_Msg "mkdir -p ${ResultsFolder}"
mkdir -p ${ResultsFolder}

#now that we have the final Native fMRI space, resample the T1w-space sebased bias field related outputs
#the alternative is to add a bunch of optional arguments to OneStepResampling that just do the same thing
#we need to do this before intensity normalization, as it uses the bias field output
if [[ ${DistortionCorrection} == "TOPUP" ]]
then
#OVstart
echo "Resampling bias field to fMRI resolution"
    
    #create MNI space corrected fieldmap images
    #${FSLDIR}/bin/applywarp --rel --interp=spline --in=${DCFolder}/PhaseOne_gdc_dc_unbias -w ${AtlasSpaceFolder}/xfms/${AtlasTransform} -r ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin -o ${ResultsFolder}/${NameOffMRI}_PhaseOne_gdc_dc
    #${FSLDIR}/bin/applywarp --rel --interp=spline --in=${DCFolder}/PhaseTwo_gdc_dc_unbias -w ${AtlasSpaceFolder}/xfms/${AtlasTransform} -r ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin -o ${ResultsFolder}/${NameOffMRI}_PhaseTwo_gdc_dc
    #OV: the bias fields are already in T1 space so let's just copy them to not mess up naming later
#    cp ${DCFolder}/PhaseOne_gdc_dc_unbias.nii.gz ${ResultsFolder}/${NameOffMRI}_PhaseOne_gdc_dc.nii.gz
#    cp ${DCFolder}/PhaseTwo_gdc_dc_unbias.nii.gz ${ResultsFolder}/${NameOffMRI}_PhaseTwo_gdc_dc.nii.gz

    #create MNINonLinear final fMRI resolution bias field outputs
    if [[ ${BiasCorrection} == "SEBASED" ]]
    then
echo "Skipping the bias field warping to MNI"
#        ${FSLDIR}/bin/applywarp --interp=trilinear -i ${DCFolder}/ComputeSpinEchoBiasField/sebased_bias_dil.nii.gz -r ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin -w ${AtlasSpaceFolder}/xfms/${AtlasTransform} -o ${ResultsFolder}/${NameOffMRI}_sebased_bias.nii.gz
#       ${FSLDIR}/bin/fslmaths ${ResultsFolder}/${NameOffMRI}_sebased_bias.nii.gz -mas ${fMRIFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz ${ResultsFolder}/${NameOffMRI}_sebased_bias.nii.gz
        
#        ${FSLDIR}/bin/applywarp --interp=trilinear -i ${DCFolder}/ComputeSpinEchoBiasField/sebased_reference_dil.nii.gz -r ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin -w ${AtlasSpaceFolder}/xfms/${AtlasTransform} -o ${ResultsFolder}/${NameOffMRI}_sebased_reference.nii.gz
#       cp ${DCFolder}/ComputeSpinEchoBiasField/sebased_reference_dil.nii.gz ${ResultsFolder}/${NameOffMRI}_sebased_reference.nii.gz
#       ${FSLDIR}/bin/fslmaths ${ResultsFolder}/${NameOffMRI}_sebased_reference.nii.gz -mas ${fMRIFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz ${ResultsFolder}/${NameOffMRI}_sebased_reference.nii.gz
        
#        ${FSLDIR}/bin/applywarp --interp=trilinear -i ${DCFolder}/ComputeSpinEchoBiasField/${NameOffMRI}_dropouts.nii.gz -r ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin -w ${AtlasSpaceFolder}/xfms/${AtlasTransform} -o ${ResultsFolder}/${NameOffMRI}_dropouts.nii.gz
#       cp ${DCFolder}/ComputeSpinEchoBiasField/${NameOffMRI}_dropouts.nii.gz ${ResultsFolder}/${NameOffMRI}_dropouts.nii.gz
#OV Need to downsample the spin echo based field maps to the fMRI resolution
# Step 0: calculate the registratio between high res bias map and rfMRI SBRef, store .mat file
  ${FSLDIR}/bin/flirt \
	-in ${DCFolder}/ComputeSpinEchoBiasField/${NameOffMRI}_sebased_bias.nii.gz \
	-ref ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin\
	-omat ${fMRIFolder}/Highres_sebased_bias2${NameOffMRI}.mat

#Step 1: dilate the highres bias field to the full FOV to be on the safe side in case the subject moved a lot
  ${FSLDIR}/bin/fslmaths ${DCFolder}/ComputeSpinEchoBiasField/${NameOffMRI}_sebased_bias.nii.gz \
                        -dilall ${DCFolder}/ComputeSpinEchoBiasField/${NameOffMRI}_sebased_bias_dilated.nii.gz  

#Step 2: Apply registration to rfMRI reference (also does the downsampling) 
  ${FSLDIR}/bin/flirt \
	-in ${DCFolder}/ComputeSpinEchoBiasField/${NameOffMRI}_sebased_bias_dilated.nii.gz \
	-ref ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin\
        -applyxfm -init ${fMRIFolder}/Highres_sebased_bias2${NameOffMRI}.mat \
	-out ${fMRIFolder}/${NameOffMRI}_sebased_bias_dilated.${FinalfMRIResolution}.nii.gz

#STep 3: And apply brain mask in fMIR resolution
  ${FSLDIR}/bin/fslmaths ${fMRIFolder}/${NameOffMRI}_sebased_bias_dilated.${FinalfMRIResolution}.nii.gz \
	-mas ${fMRIFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz \
         ${fMRIFolder}/${NameOffMRI}_sebased_bias.${FinalfMRIResolution}.nii.gz
 

#OVend
   fi
fi

#Intensity Normalization and Bias Removal
#OV changed naming of jacobian to native one and swapped MNI bias field for native acpc bias field
#OV changed to bias field in native space, but at rfMRI resolution
log_Msg "Intensity Normalization and Bias Removal"
${RUN} ${PipelineScripts}/IntensityNormalization.sh \
       --infmri=${fMRIFolder}/${NameOffMRI}_nonlin \
       --biasfield=${fMRIFolder}/${NameOffMRI}_sebased_bias.${FinalfMRIResolution} \
       --jacobian=${fMRIFolder}/${JacobianOut}_Native.${FinalfMRIResolution}  \
       --brainmask=${fMRIFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution} \
       --ofmri=${fMRIFolder}/${NameOffMRI}_nonlin_norm \
       --inscout=${fMRIFolder}/${NameOffMRI}_SBRef_nonlin \
       --oscout=${fMRIFolder}/${NameOffMRI}_SBRef_nonlin_norm \
       --usejacobian=${UseJacobian}

# MJ QUERY: WHY THE -r OPTIONS BELOW?
# TBr Response: Since the copy operations are specifying individual files
# to be copied and not directories, the recursive copy options (-r) to the
# cp calls below definitely seem unnecessary. They should be removed in 
# a code clean up phase when tests are in place to verify that removing them
# has no unexpected bad side-effect.
${RUN} cp -r ${fMRIFolder}/${NameOffMRI}_nonlin_norm.nii.gz ${ResultsFolder}/${NameOffMRI}.nii.gz
${RUN} cp -r ${fMRIFolder}/${MovementRegressor}.txt ${ResultsFolder}/${MovementRegressor}.txt
${RUN} cp -r ${fMRIFolder}/${MovementRegressor}_dt.txt ${ResultsFolder}/${MovementRegressor}_dt.txt
${RUN} cp -r ${fMRIFolder}/${NameOffMRI}_SBRef_nonlin_norm.nii.gz ${ResultsFolder}/${NameOffMRI}_SBRef.nii.gz
${RUN} cp -r ${fMRIFolder}/${JacobianOut}_Native.${FinalfMRIResolution}.nii.gz ${ResultsFolder}/${NameOffMRI}_${JacobianOut}.nii.gz #OV changed naming to Native
${RUN} cp -r ${fMRIFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz ${ResultsFolder}/${FreeSurferBrainMask}.${FinalfMRIResolution}.nii.gz
#OV also copies the bias field at fMRI resolution
${RUN} cp -r ${fMRIFolder}/${NameOffMRI}_sebased_bias.${FinalfMRIResolution}.nii.gz ${ResultsFolder}/${NameOffMRI}_sebased_bias.${FinalfMRIResolution}.nii.gz
###Add stuff for RMS###
${RUN} cp -r ${fMRIFolder}/Movement_RelativeRMS.txt ${ResultsFolder}/Movement_RelativeRMS.txt
${RUN} cp -r ${fMRIFolder}/Movement_AbsoluteRMS.txt ${ResultsFolder}/Movement_AbsoluteRMS.txt
${RUN} cp -r ${fMRIFolder}/Movement_RelativeRMS_mean.txt ${ResultsFolder}/Movement_RelativeRMS_mean.txt
${RUN} cp -r ${fMRIFolder}/Movement_AbsoluteRMS_mean.txt ${ResultsFolder}/Movement_AbsoluteRMS_mean.txt
###Add stuff for RMS###

#OV: Run registration of surfaces from native space (ac_pc) to EPI space
log_Msg "Running registration of surfaces from native to EPI space"
#Set the correct subjects directory where the freesurfer econstruction lives
SUBJECTS_DIR=${T1wFolder}
export SUBJECTS_DIR
#undistorted EPI reference volume
EPI_dc_jac=${DCFolder}/FieldMap/SBRef_dc_jac.nii.gz
fMRI2str_fsl=${DCFolder}/fMRI2str.mat
ACPC_T1=${T1wFolder}/${T1wImage}.nii.gz
# mris_convert version of the ac-pc T1, this has a different tkr, but surfaces are constructed from here, 
#we will need this to convert from orig.mgz to ac-pc.nii to EPI because the fsl fMRI2str.mat file is a 
#registration between ac-pc and EPI and not from irig.mgz to EPI
orig=${T1wFolder}/${Subject}/mri/orig.mgz
NSurf=10
#Loop over all intermittent surfaces that we created in the PostFreeSurferPipelineBatchInterSurf.sh script
hemis=(lh rh)
WriteTransformsToHere=${T1wFolder}/${Subject}/surf #Location to store surfaces
for hemi in ${hemis[*]}; do
    for surf in `seq -f %02.0f 0 $NSurf`; do
        SurfName=midgray.${surf}
        ${PipelineScripts}/Surf2EPI.sh ${Subject} \
                                                            ${orig} \
                                                            ${ACPC_T1} \
                                                            ${EPI_dc_jac} \
                                                            ${fMRI2str_fsl} \
                                                            ${SurfName} \
                                                            ${hemi} \
                                                            ${NameOffMRI} \
                                                            ${WriteTransformsToHere}
        SurfName=subwhite.${surf}
        ${PipelineScripts}/Surf2EPI.sh ${Subject} \
                                                            ${orig} \
                                                            ${ACPC_T1} \
                                                            ${EPI_dc_jac} \
                                                            ${fMRI2str_fsl} \
                                                            ${SurfName} \
                                                            ${hemi} \
                                                            ${NameOffMRI} \
                                                            ${WriteTransformsToHere}
    done
done

#******************* WARP SURFACE OVERLAYS ************************************************************ 
#Also create surfae overlays of the warp field of the original distorted EPI. 
#We can then use the overlays later in Matlab to calculate the actual warp of the surfaces.
#Take the warp from EPI to T1-acpc and invert it, this is in relative conventions i.e. x'=x+w(x), but 
#for the matlab functio we need absolute x'=w(x)
#fsl.mat affine to get from EPI space to topup spin echo space
EPI_nodc_2TopUp=${DCFolder}/FieldMap/SBRef2WarpField.mat
EPI_nodc_2TopUp_lta=${DCFolder}/FieldMap/SBRef2WarpField.lta
TopUp2EPI_nodc_lta=${DCFolder}/FieldMap/WarpField2EPI.lta
#Get distortion direction (either LR or RL)
Distortion_dir=$(echo $NameOffMRI |cut -d "_" -f3)

#The WarpField files have opposite signs, the number 01 is for LR and 04 for RL
if [[ $Distortion_dir=="LR" ]]
 then 
   TopUpWarp=${DCFolder}/FieldMap/WarpField_01.nii.gz 
 else
   TopUpWarp=${DCFolder}/FieldMap/WarpField_04.nii.gz 
fi

TopUpWarp_surf_abs=${T1wFolder}/xfms/TopUpWarp_${Distortion_dir}_surf_abs.nii.gz
TopUpWarp_surf_abs_2EPI_nodc=${T1wFolder}/xfms/TopUpWarp_${Distortion_dir}_surf_abs2EPI_nodc.nii.gz
EPI_nodc_SBRef=${SubjectFolder}/unprocessed/3T/${NameOffMRI}/${Subject}_3T_${NameOffMRI}_SBRef.nii.gz

#Use convertwarp to convert the warp to absolute mm displacements
${FSLDIR}/bin/convertwarp -r ${TopUpWarp} --rel -w ${TopUpWarp} --absout -o ${TopUpWarp_surf_abs}

##convert the fsl .mat affine to an lta - we need this for mri_vol2surf
#OVout#lta_convert --infsl ${EPI_nodc_2TopUp} --src ${EPI_nodc_SBRef} --trg ${TopUpWarp} --outlta ${EPI_nodc_2TopUp_lta}
#OVout#lta_convert --inlta ${EPI_nodc_2TopUp_lta} --outlta ${TopUp2EPI_nodc_lta} --invert 

##use mri_vol2vol to transform the warp field into the epi space
#mri_vol2vol --mov ${TopUpWarp_surf_abs} --targ ${EPI_nodc_SBRef} --o ${TopUpWarp_surf_abs_2EPI_nodc} --lta ${TopUp2EPI_nodc_lta}

#Create surface overlays of the warp field for each surface
#OVout#for hemi in ${hemis[*]}; do
#OVout#    for surf in `seq -f %03.0f 0 $NSurf`; do
#OVout#        SurfName=midgray.${surf}_2${NameOffMRI}
#OVout#        TopUpWarp_surf=${T1wFolder}/${Subject}/surf/${hemi}.${SurfName}_WarpOverlay.mgz
#OVout#        mri_vol2surf \
#OVout#                      --srcsubject ${Subject} \
#OVout#                      --srcreg ${TopUp2EPI_nodc_lta} \
#OVout#                      --src ${TopUpWarp_surf_abs_2EPI_nodc} \
#OVout#                      --out ${TopUpWarp_surf} \
#OVout#                      --hemi ${hemi} \
#OVout#                      --surf ${SurfName} \
#OVout#                      --interp trilinear 
#OVout#
#OVout#        SurfName=subwhite.${surf}_2${NameOffMRI}
#OVout#        TopUpWarp_surf=${T1wFolder}/${Subject}/surf/${hemi}.${SurfName}_WarpOverlay.mgz
#OVout#        mri_vol2surf \
#OVout#                      --srcsubject ${Subject} \
#OVout#                      --srcreg ${TopUp2EPI_nodc_lta} \
#OVout#                     --src ${TopUpWarp_surf_abs_2EPI_nodc} \
#OVout#                    --out ${TopUpWarp_surf} \
#OVout#                      --hemi ${hemi} \
#OVout#                      --surf ${SurfName} \
#OVout#                      --interp trilineardone
#OVout#   done 
#OVout#done



#OV ***************************************Warp field to EPI space registration      *****
#	First we tranform the topup warpfield DistortionCorrectionAndEPIToT1wReg_FLIRTBBRAndFreeSurferBBRbased/FieldMap/WarpField01.nii.gz
#       into the position of the EPI to get each voxels displacement value. We do this with the inverse affine of SBRef2WarpField.mat.
#OVout#SBRef2WarpFieldTransform=${DCFolder}/FieldMap/SBRef2WarpField.mat
#OVout#WarpField2SBRefTransform=${DCFolder}/FieldMap/WarpField2SBRef.mat
#OVout#WarpFieldInEPI=${T1wFolder}/Results/${NameOffMRI}/WarpField2${NameOffMRI}.nii.gz

#OVout#${FSLDIR}/bin/convert_xfm -omat ${WarpField2SBRefTransform} -inverse ${SBRef2WarpFieldTransform}
#OVout#${FSLDIR}/bin/flirt -in ${TopUpWarp} -ref ${EPI_nodc_SBRef} -applyxfm -init ${WarpField2SBRefTransform} -out ${WarpFieldInEPI} 

#Also transform warp field into T1 acpc aligned space
#OVout#WarpFieldIn_acpc=${T1wFolder}/Results/${NameOffMRI}/WarpField2acpc_${NameOffMRI}.nii.gz
#OVout#${FSLDIR}/bin/flirt -in ${WarpFieldInEPI} -ref ${EPI_dc_jac} -applyxfm -init ${fMRI2str_fsl} -out ${WarpFieldIn_acpc}

log_Msg "Completed"

