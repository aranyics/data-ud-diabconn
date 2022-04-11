#!/bin/bash

# The goal of this script is to extract preprocessed fMRI images
# from the nipype output directory ($PROCESSDIR)
# to a directory, where later analysis looks for imaging data ($GIFTINPUTDIR),
# and arrange/rename them in a simply automatized format:
# ($GIFTINPUTDIR/$s/rfmri_filtered.nii)

PROJDIR=/path/to/diabetes                # directory for the project
PROCESSDIR=$PROJDIR/process              # results from nipype preprocessing pipeline
PREFIX=d0                                # subject prefix (o0 for obese group)
GIFTINPUTDIR=$PROJDIR/analysis/ICA/input # subject directories for later analysis
POSTFIX=uni-1                            # miscellaneous filename components

for s in `ls $PROCESSDIR | grep $PREFIX`; do
    echo $s;
    gunzip -f ${PROCESSDIR}/${s}/1/filter/${s}_${POSTFIX}_filtered.nii.gz;
    mkdir -p ${GIFTINPUTDIR}/${s};
    cd ${GIFTINPUTDIR}/${s};
    #echo "ln -s ${PROCESSDIR}/${s}/1/filter/${s}_${POSTFIX}_filtered.nii ./rfmri_filtered.nii;"
    ln -s ${PROCESSDIR}/${s}/1/filter/${s}_${POSTFIX}_filtered.nii ./rfmri_filtered.nii;
done
