#!/bin/bash

source /etc/profile.d/modules.sh  
# Initialize conda for the shell
source /home/daisuke/anaconda3/etc/profile.d/conda.sh

# Activate your environment (replace antsenv with your env name)
conda activate antsenv

INPUTNII=$1
OUTDIR=$2

cd $OUTDIR

if [ ! -f "$OUTDIR/$INPUTNII" ]
then
    echo "Input nifti not found $1"
    exit
fi

#`dirname $0`/m2k.py $OUTDIR/TexVolSmooth.nii $OUTDIR/TexVolSmooth_m2k.nii

basename="${INPUTNII%.nii}"
OUTPUTNII="${basename}_T2.nii"



# forward transform from CCF to T2
antsApplyTransforms \
 -i $OUTDIR/$INPUTNII \
 -r T2w_brain_us.nii \
 -o $OUTDIR/$OUTPUTNII \
 -t Tem_to_T21Warp.nii.gz \
 -t Tem_to_T20GenericAffine.mat \
 --interpolation BSpline[1]
