#!/usr/bin/env python3

import os
import sys
import nibabel
import numpy

debug = False
if debug:
    print(os.getcwd())

import Otsu
import scipy.ndimage
import numpy.linalg
import Utils

nVoxels_th = 17 #200
EigFactorThresh = 4 #5

# segments the pills in Daisuke's mouse images

ALLENSpaceImage = sys.argv[1]
BrainMaskImage = sys.argv[2]
PillMaskImage = sys.argv[3]
OutSegImage = sys.argv[4]

head, tail = os.path.split(OutSegImage)

ALLENSpaceNII = nibabel.load(ALLENSpaceImage)
ALLENSpaceIMG = numpy.single(ALLENSpaceNII.get_fdata())
ALLENSpaceIMGFlat = numpy.reshape(ALLENSpaceIMG, [numpy.prod(ALLENSpaceIMG.shape[0:3]), ALLENSpaceIMG.shape[3]], order='F')

PillMaskNII = nibabel.load(PillMaskImage)
PillMaskIMG = numpy.uint8(PillMaskNII.dataobj)

OrigShape = PillMaskIMG.shape


BrainMaskNII = nibabel.load(BrainMaskImage)
#BrainMaskIMG = numpy.array(BrainMaskNII.dataobj) > 0
BrainMaskIMG = numpy.asanyarray(BrainMaskNII.dataobj) != 0

CentreMask = scipy.ndimage.binary_dilation(PillMaskIMG == 1, iterations=1)
#CentreMask = numpy.logical_and(BrainMaskNotIMG, scipy.ndimage.binary_dilation(CentreMask, iterations=7))
RightMask = scipy.ndimage.binary_dilation(PillMaskIMG == 2, iterations=1)
#RightMask = numpy.logical_and(BrainMaskNotIMG, scipy.ndimage.binary_dilation(RightMask, iterations=7))
LeftMask = scipy.ndimage.binary_dilation(PillMaskIMG == 3, iterations=1)
#LeftMask = numpy.logical_and(BrainMaskNotIMG, scipy.ndimage.binary_dilation(LeftMask, iterations=7))

#PillMaskIMG = PillMaskIMG.ravel(order='F')
#CentreMask = CentreMask.ravel(order='F')
#RightMask = RightMask.ravel(order='F')
#LeftMask = LeftMask.ravel(order='F')

def getCurves(IMG, IMGFlat, MaskFlat):
    Y = numpy.take(IMGFlat, numpy.where(MaskFlat)[0], axis = 0)
    E = numpy.zeros((Y.shape[0]))
    R = numpy.zeros((Y.shape[0]))
    
    X = numpy.arange(Y.shape[1])
    for z in range(Y.shape[0]):
        F = numpy.polyfit(X, numpy.log(Y[z, :] + 1), 1, w=numpy.sqrt(Y[z, :] + 1), full = False)
        E[z] = F[0]
        T = numpy.exp(F[1]) * numpy.exp(F[0] * X) - Y[z, :]
        R[z] = numpy.sum(T * T)
    T = numpy.zeros(IMGFlat.shape[0])
    T[MaskFlat] = E * R
    
    return numpy.reshape(T, IMG.shape[0:3], order='F')

#G = scipy.ndimage.median_filter(ALLENSpaceIMG, size = 9)
#OutNII = nibabel.Nifti1Image(M, PillMaskNII.affine)
#%nibabel.save(OutNII, os.path.join(head, 'M.nii.gz'))

ECentre = getCurves(ALLENSpaceIMG, ALLENSpaceIMGFlat, CentreMask.ravel(order='F'))
ELeft = getCurves(ALLENSpaceIMG, ALLENSpaceIMGFlat, LeftMask.ravel(order='F'))
ERight = getCurves(ALLENSpaceIMG, ALLENSpaceIMGFlat, RightMask.ravel(order='F'))

OrigE = numpy.minimum(ECentre + ELeft + ERight, 0)

if debug:
    ELeft_img = nibabel.Nifti1Image(ELeft, BrainMaskNII.affine)
    nibabel.save(ELeft_img, os.path.join(head, 'ELeft.nii.gz'))

G = scipy.ndimage.grey_opening(ALLENSpaceIMG[:, :, :, 0], size=(5, 5, 5))#size=(15, 15, 15)
GG = ALLENSpaceIMG[:, :, :, 0] - G
OrigE = OrigE * GG

if debug:
    OrigE_img = nibabel.Nifti1Image(OrigE, BrainMaskNII.affine)
    nibabel.save(OrigE_img, os.path.join(head, 'OrigE.nii.gz'))

def makeLabelsForMask(E, Mask):
    M = E[Mask]
    M = numpy.median(M[M < 0])
    L, numLabelsCentre = scipy.ndimage.label(numpy.logical_and(E < M * 10, Mask))
    H = numpy.bincount(L[L > 0])
    LToKeep = Utils.ismember(L, numpy.where(H > nVoxels_th)[0])
    L, numLabels = scipy.ndimage.label(LToKeep)
    return (L, numLabels)

LeftL, numLabelsLeft = makeLabelsForMask(OrigE, LeftMask)
CentreL, numLabelsCentre = makeLabelsForMask(OrigE, CentreMask)
RightL, numLabelsRight = makeLabelsForMask(OrigE, RightMask)

if debug:
    print("numLabelsLeft " + str(numLabelsLeft))
    print("numLabelsCentre " + str(numLabelsCentre))
    print("numLabelsRight " + str(numLabelsRight))
    CentreL_img = nibabel.Nifti1Image(CentreL, BrainMaskNII.affine)
    nibabel.save(CentreL_img, os.path.join(head, 'CentreL.nii.gz'))
    RightL_img = nibabel.Nifti1Image(RightL, BrainMaskNII.affine)
    nibabel.save(RightL_img, os.path.join(head, 'RightL.nii.gz'))

def EigensAllLabels(L):
    numLabels = numpy.max(L)
    eigenValuesVectors = list()
    for z in range(1, numLabels + 1):
        I = numpy.where(L == z)
        XYZ = numpy.vstack(I)
        if XYZ.shape[1] <= 3:
            eigenValuesVectors.append(None)
        else:
            C = numpy.cov(XYZ)
            eigenValuesVectors.append(numpy.linalg.eig(C))
    return eigenValuesVectors

outLabels = numpy.zeros(PillMaskIMG.shape, dtype = numpy.uint8)

def eigenFactorsAll(L, DominantAxis):
    E = EigensAllLabels(L)
    EigFactors = list()

    for z in range(len(E)):
        if E[z] is None:
            EigFactors.append(-1)
        else:
            # left-right, the first eigenvector's dominant element must be the dominant axis

            OtherIDX = numpy.setdiff1d([0, 1, 2], [DominantAxis])
            S = numpy.argsort(E[z][0])

            EigValues = E[z][0][S]
            EigVectors = numpy.take(E[z][1], S, axis = 1)
            #print(EigValues)
            #print(EigVectors)
            if numpy.all(numpy.abs(EigVectors[DominantAxis, 2]) > numpy.abs(EigVectors[OtherIDX, 2])):
                EigFactors.append(EigValues[2] / EigValues[1])
                # = numpy.max(numpy.abs(EigValues[OtherIDX, 2]))
                #if M == 0:
                #    EigFactors.append(-1)
                #else:
                    #EigFactors.append(numpy.abs(EigVectors[DominantAxis, 2]) / M)
                #    EigFactors.append(numpy.abs(EigVectors[DominantAxis, 2]) / M)
            else:
                EigFactors.append(-1)
    return numpy.array(EigFactors)

def heuristicChoice(IMG, Mask):
    P = numpy.min(IMG[Mask]) + 0.8 * numpy.ptp(IMG[Mask])
    return scipy.ndimage.binary_dilation(numpy.logical_and(IMG > P, Mask))



if numLabelsCentre > 50 or numLabelsCentre == 0:
    outLabels[heuristicChoice(ALLENSpaceIMG[:, :, :, 0], CentreMask)] = 1
elif numLabelsCentre > 1:
    # looking for left-right shape, dominant axis is 0
    EigFactorsCentre = eigenFactorsAll(CentreL, 0)
    H = numpy.bincount(CentreL[CentreL > 0])
    print(EigFactorsCentre)
    CentreToChoose = numpy.where(EigFactorsCentre > EigFactorThresh)[0] + 1
    print("Centre Choosing " + str(CentreToChoose + 1))
    outLabels[Utils.ismember(CentreL, CentreToChoose)] = 1
elif numLabelsCentre == 1:
    outLabels[CentreL == 1] = 1

if numLabelsLeft > 50 or numLabelsLeft == 0:
    outLabels[heuristicChoice(ALLENSpaceIMG[:, :, :, 0], LeftMask)] = 2
elif numLabelsLeft > 1:
    # looking for anterior-posterior shape, dominant axis is 2? instead of 1
    EigFactorsLeft = eigenFactorsAll(LeftL, 2)
    print(EigFactorsLeft)
    LeftToChoose = numpy.where(EigFactorsLeft > EigFactorThresh)[0] + 1
    print("Left Choosing " + str(LeftToChoose))
    T = Utils.ismember(LeftL, LeftToChoose)
    outLabels[Utils.ismember(LeftL, LeftToChoose)] = 2
else:
    outLabels[LeftL == 1] = 2


if numLabelsRight > 50 or numLabelsRight == 0:
    outLabels[heuristicChoice(ALLENSpaceIMG[:, :, :, 0], RightMask)] = 3
elif numLabelsRight > 1:
    # looking for anterior-posterior shape, dominant axis is 2? instead of 1
    EigFactorsRight = eigenFactorsAll(RightL, 2)
    print(EigFactorsRight)
    RightToChoose = numpy.where(EigFactorsRight > EigFactorThresh)[0] + 1
    print("Right Choosing " + str(RightToChoose))
    outLabels[Utils.ismember(RightL, RightToChoose)] = 3
else:
    outLabels[RightL == 1] = 3



if debug:
    outLabels_img = nibabel.Nifti1Image(outLabels, BrainMaskNII.affine)
    nibabel.save(outLabels_img, os.path.join(head, 'outLabels.nii.gz'))

# remove components that intersect with the brain mask
ALLENSpaceIMGOtsu = Otsu.robustOtsu(ALLENSpaceIMG[:, :, :, 0], [0.05, 0.95], NumberClasses=3)
L, numLabels = scipy.ndimage.label(ALLENSpaceIMGOtsu > 1)
H = numpy.bincount(L[L > 0])
ALLENSpaceIMGOtsuLargest = (L == numpy.argmax(H))

#OutNII = nibabel.Nifti1Image(numpy.int16(ALLENSpaceIMGOtsuLargest), PillMaskNII.affine)
#nibabel.save(OutNII, os.path.join(head, 'largest.nii.gz'))
#OutNII = nibabel.Nifti1Image(numpy.int16(outLabels), PillMaskNII.affine)
#nibabel.save(OutNII, os.path.join(head, 'outLabels.nii.gz'))

T = numpy.array(outLabels)
T[outLabels > 0] = 1
outLabelsL, ff = scipy.ndimage.label(T)

outLabelsOverLargest = outLabelsL[ALLENSpaceIMGOtsuLargest]
outLabelsOverLargest = numpy.unique(outLabelsOverLargest[outLabelsOverLargest > 0])

if outLabelsOverLargest.size > 0:
    print("removing " + str(outLabelsOverLargest))
    outLabels[Utils.ismember(outLabelsL, outLabelsOverLargest)] = 0

# remove labels less than 1mm from brain

D = scipy.ndimage.distance_transform_edt(numpy.logical_not(BrainMaskIMG), sampling=BrainMaskNII.header.get_zooms())

#print(BrainMaskNII.header.get_zooms())

T = numpy.array(outLabels)
T[outLabels > 0] = 1
outLabelsL, ff = scipy.ndimage.label(T)

outLabelsNearBrain = outLabelsL[D <= 0.5]
outLabelsNearBrain = numpy.unique(outLabelsNearBrain[outLabelsNearBrain > 0])

if outLabelsNearBrain.size > 0:
    print("removing " + str(outLabelsNearBrain))
    outLabels[Utils.ismember(outLabelsL, outLabelsNearBrain)] = 0

OutNII = nibabel.Nifti1Image(numpy.uint8(outLabels), PillMaskNII.affine)
nibabel.save(OutNII, OutSegImage)

#OutNII = nibabel.Nifti1Image(numpy.int16(LeftL), PillMaskNII.affine)
#nibabel.save(OutNII, os.path.join(head, 'LeftL.nii.gz'))
#OutNII = nibabel.Nifti1Image(numpy.int16(CentreL), PillMaskNII.affine)
#nibabel.save(OutNII, os.path.join(head, 'CentreL.nii.gz'))
#OutNII = nibabel.Nifti1Image(numpy.int16(RightL), PillMaskNII.affine)
#nibabel.save(OutNII, os.path.join(head, 'RightL.nii.gz'))

#OutNII = nibabel.Nifti1Image(Markers, PillMaskNII.affine)
#nibabel.save(OutNII, OutSegImage)
#OutNII = nibabel.Nifti1Image(HeightMap, PillMaskNII.affine)
#nibabel.save(OutNII, os.path.join(head, 'aaa.nii.gz'))
#OutNII = nibabel.Nifti1Image(OrigE, PillMaskNII.affine)
#nibabel.save(OutNII, os.path.join(head, 'OrigE.nii.gz'))

#OutNII = nibabel.Nifti1Image(numpy.uint8(outLabels), PillMaskNII.affine)
#nibabel.save(OutNII, OutSegImage)
