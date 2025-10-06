#!/bin/bash

source /etc/profile.d/modules.sh  
# Initialize conda for the shell
source /home/daisuke/anaconda3/etc/profile.d/conda.sh

# Activate your environment (replace antsenv with your env name)
conda activate antsenv

INPUTNII=$1
OUTDIR=$2

if [ ! -f "$INPUTNII" ]
then
    echo "Input nifti not found $1"
    exit
fi

mkdir -p $OUTDIR

# split the inputs
`dirname $0`/m2k.py $INPUTNII $OUTDIR/input_swapped_m2k.nii.gz

module load freesurfer

fslmaths $OUTDIR/input_swapped_m2k.nii.gz -Tmean $OUTDIR/average_swapped_m2k.nii.gz -odt short

cd $OUTDIR

if [ -f T2w_refit.nii ]; then
    rm T2w_refit.nii
fi
if [ -f T2w_resample.nii ]; then
    rm T2w_resample.nii
fi
3dcopy average_swapped_m2k.nii T2w_refit.nii
3drefit -xyzscale 10.0 T2w_refit.nii
3dresample -dxyz 1 1 1 -rmode Cubic -inset T2w_refit.nii -prefix T2w_resample.nii
if [ -d brain ]; then
    rmdir -rf brain
fi
mkdir brain
cp T2w_resample.nii brain
cd brain
fslmaths T2w_resample.nii -Tmean T2w_a
mri_convert T2w_a.nii.gz T2w_a_conv.mnc
nu_correct T2w_a_conv.mnc T2w_a_conv.nu.mnc -iterations 20
imp2field -like T2w_a_conv.nu.mnc T2w_a_conv.nu.imp field.mnc
mri_convert field.mnc field.nii
mri_convert T2w_a_conv.nu.mnc T2w_a_nu.nii
rm -f T2w_a_conv.nu.mnc
fslmaths T2w_a_nu -thrP 20 a
bet a b -R -f 0.15 -g 0
3dAutomask -prefix T2w_brain_mask.nii.gz -apply_prefix c.nii.gz -clfrac 0.15 b.nii.gz

cp c.nii.gz ../
cp T2w_brain_mask.nii.gz ../
cd ..
if [ -f T2w_brain.nii ]; then
    rm T2w_brain.nii
fi
3dcopy c.nii.gz T2w_brain.nii
rm -f c.nii.gz
rm -rf brain

# Perform registration from Allen space to native T2w space
antsRegistrationSyN.sh -d 3 -f T2w_brain.nii -o Tem_to_T2 -m Brain_template.nii -t s

# Apply the forward transformation to the atlas annotation to bring it into T2w_brain.nii space
antsApplyTransforms \
 -i Allen_annotation_modified.nii \
 -r T2w_brain.nii \
 -o Atlas_anno_to_T2.nii.gz \
 -t Tem_to_T21Warp.nii.gz \
 -t Tem_to_T20GenericAffine.mat \
 -n NearestNeighbor

# Apply the inverse transformation to bring T2w_brain.nii into Brain_template.nii space
# Split 4D input into 3D volumes
fslsplit input_swapped_m2k.nii.gz vol_

# Apply the inverse transforms to each 3D volume and collect outputs
for v in vol_*.nii.gz; do
    3drefit -xyzscale 10.0 "$v"

    3dresample -dxyz 1 1 1 -rmode Cubic -inset "$v" -prefix "r_$v"

    antsApplyTransforms -i "r_$v" \
      -r Brain_template.nii -o "out_$v" \
      -t [Tem_to_T20GenericAffine.mat,1] \
      -t Tem_to_T21InverseWarp.nii.gz \
      -n Linear
done

# Merge transformed 3D volumes back into a 4D file
fslmerge -t input_swapped_m2k_Allen.nii.gz out_vol_*.nii.gz

# Clean up intermediate files
rm vol_*.nii.gz out_vol_*.nii.gz r_vol_*.nii.gz

# Inverse transform the brain mask to Allen space
antsApplyTransforms \
 -i T2w_brain_mask.nii.gz \
 -r Brain_template.nii \
 -o T2w_brain_mask_Allen.nii \
 -t [Tem_to_T20GenericAffine.mat,1] \
 -t Tem_to_T21InverseWarp.nii.gz \
 -n NearestNeighbor
 #-n Linear


if [ -f Atlas_anno_to_T2.nii ]; then
    rm Atlas_anno_to_T2.nii
fi
3dcopy Atlas_anno_to_T2.nii.gz Atlas_anno_to_T2.nii 

if [ -f $OUTDIR/pills_labels_Allen.nii ]; then
    rm $OUTDIR/pills_labels_Allen.nii
fi

# Run the pill finding script in Allen space 
python3 `dirname $0`/FindPillsExp_Allen.py $OUTDIR/input_swapped_m2k_Allen.nii.gz \
    $OUTDIR/T2w_brain_mask_Allen.nii $OUTDIR/Allen_pills_mask.nii $OUTDIR/pills_labels_Allen.nii

if [ -f $OUTDIR/pills_labels.nii ]; then
    rm $OUTDIR/pills_labels.nii
fi

# Forward transform the pill labels back to the original T2w space
antsApplyTransforms \
 -i pills_labels_Allen.nii \
 -r T2w_brain.nii \
 -o pills_labels.nii \
 -t Tem_to_T21Warp.nii.gz \
 -t Tem_to_T20GenericAffine.mat \
 -u uchar \
 -n GenericLabel \
# -n NearestNeighbor \ #No difference
 