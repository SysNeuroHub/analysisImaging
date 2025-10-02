#!/usr/bin/env python3
# created from change average_orientation.py

import nibabel
import numpy
import os
import sys


NII = nibabel.load(sys.argv[1])
A = NII.affine
IMG = numpy.int16(NII.get_fdata())
#if IMG.ndim == 4:
#    IMG = IMG.take([0], axis = 3).squeeze()
#B = numpy.vstack((
#    A[0, :],
#    A[2, :],
#    -A[1, :],
#    A[3, :]))
pixdim = NII.header.get_zooms() 
C = numpy.zeros((4,4))
C[0, 0] = pixdim[0]
C[1, 1] = pixdim[1]
C[2, 2] = pixdim[2] 
C[3, 3] = 1.0

# Equivalent to MATLAB operations: rot90(permute(rot90(v,1),[2 3 1]),2)
# rot90(v, 1): rotate 90 degrees counterclockwise in the first two axes
v1 = numpy.rot90(IMG, k=1, axes=(0, 1))
# permute axes [2, 3, 1] in MATLAB is [1, 2, 0] in Python (0-based)
v2 = numpy.transpose(v1, (1, 2, 0, 3))
# rot90(..., 2): rotate 180 degrees counterclockwise in the first two axes
IMG_M2K = numpy.rot90(v2, k=2, axes=(0, 1))
    
OutNII = nibabel.Nifti1Image(IMG_M2K, C)
nibabel.save(OutNII, sys.argv[2])
