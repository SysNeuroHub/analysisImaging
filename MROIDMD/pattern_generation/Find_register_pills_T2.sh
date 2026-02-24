#!/bin/bash
# detect pills in the first input.nii, then warp to T2 space

source /etc/profile.d/modules.sh  
# Initialize conda for the shell
source /home/daisuke/anaconda3/etc/profile.d/conda.sh

# Activate your environment (replace antsenv with your env name)
conda activate antsenv

# either T2* or UTE, warped to Allen space
OUTDIR=$1

# default value if second argument is missing
if [ $# -lt 2 ]; then
    PILL_IN_ALLEN="input_swapped_m2k_Allen.nii.gz"
else
    PILL_IN_ALLEN=$2
fi
 

# get directory of the current shell script
SCRIPT_DIR=$(dirname "$(realpath "$0")")

# Run the pill finding script in Allen space 
python3 "$SCRIPT_DIR/FindPillsExp_Allen.py" "$OUTDIR/$PILL_IN_ALLEN" \
    "$OUTDIR/T2w_brain_mask_Allen.nii" "$OUTDIR/Allen_pills_mask.nii" "$OUTDIR/pills_labels_Allen.nii"


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
 