#!/bin/bash

# Usage:
# ./align_UTE_to_Allen.sh UTE.nii.gz T2.nii.gz output_directory

source /etc/profile.d/modules.sh  
# Initialize conda for the shell
source /home/daisuke/anaconda3/etc/profile.d/conda.sh

conda activate antsenv

UTE=$1
T2=$2
OUTDIR=$3


echo "Running rigid registration with FLIRT..."

# minimization ... NOT VERY ACCURATE
flirt \
  -in $UTE \
  -ref $T2 \
  -out $OUTDIR/ute_in_t2\
  -dof 6 \
  -cost mutualinfo \
  -interp trilinear \
  -setbackground 0

echo "UTE aligned to T2"

module load freesurfer

# Compute global min/max using fslstats
UTE_FILE=${OUTDIR}/ute_in_t2.nii.gz
chmod 777 $UTE_FILE
MIN=$(fslstats $UTE_FILE -R | awk '{print $1}')
MAX=$(fslstats $UTE_FILE -R | awk '{print $2}')

fslmaths $OUTDIR/ute_in_t2.nii.gz \
-sub $MIN -div $(echo "$MAX - $MIN" | bc -l) \
-mul 65535 -sub 32768 \
  $OUTDIR/ute_in_t2_int16 \
  -odt short

echo "converted to int16"

# apply m2k
`dirname $0`/m2k.py $OUTDIR/ute_in_t2_int16.nii.gz $OUTDIR/ute_in_t2_m2k.nii.gz

echo "m2k applied"


# transform UTE in T2 space to Allen space, using the result of Atlas_T2_coreg_DS.sh

cd $OUTDIR

3drefit -xyzscale 10.0 ute_in_t2_m2k.nii.gz

 3dresample -dxyz 1 1 1 -rmode Cubic -inset ute_in_t2_m2k.nii -prefix r_ute_in_t2_m2k.nii.gz

antsApplyTransforms \
  -i r_ute_in_t2_m2k.nii.gz \
  -r Brain_template.nii \
  -o ute_in_allen.nii.gz \
  -t [Tem_to_T20GenericAffine.mat,1] \
  -t Tem_to_T21InverseWarp.nii.gz \
  -n Linear

echo "UTE aligned to Allen"