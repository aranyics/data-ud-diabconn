
def clusterize( data, minval = 3.0, minsize = 100, orthops = 2 ):
    '''based on: https://gist.github.com/ofgulban/27c4491592126dce37e97c578cbf307b'''

    import numpy as np
    from skimage.measure import label

    #intensity thresholding
    #data[data > minval] = 1.0
    data = np.where( data < minval, 0.0, 1.0 )
    data = data.astype("int")

    #connected clusters
    data = label(data, connectivity = orthops)
    labels, counts = np.unique(data, return_counts = True)

    #cluster size thresholding
    for i, (i_label, i_count) in enumerate(zip(labels[1:], counts[1:])):
        if i_count < minsize:
            data[data == i_label] = 0
    labels, counts = np.unique(data, return_counts = True)
    print('{0} clusters are found\n\t- sizes: {1}\n\t- labels: {2}'.format(labels.size, counts, labels))

    return data, labels, counts


def clustersToROI( SPM_nii, Tthres, Cthres, outdir, name ):

    import os
    import numpy as np
    import nibabel as nib

    outinfo = os.path.join( outdir, 'info' )

    try:
        os.makedirs( outdir )
    except OSError:
        pass
    try:
        os.makedirs( outinfo )
    except OSError:
        pass

    SPM_img = nib.load( SPM_nii )
    SPM = SPM_img.get_fdata()

    clusters, labels, counts = clusterize( abs(SPM), minval=Tthres, minsize=Cthres )
    outfn = ( os.path.join( outinfo, 'allROI.nii' ) )
    out = nib.Nifti1Image( clusters, header = SPM_img.header, affine = SPM_img.affine )
    nib.save( out, outfn )

    for l in labels:
        if (l == 0):
            continue
        tmp = np.copy(clusters)
        tmp[tmp != l] = 0
        tmp[tmp == l] = 1
        outfn = ( os.path.join( outdir, ''.join(['ROI', str(l), '_', name, '.nii']) ) )
        out = nib.Nifti1Image( tmp, header = SPM_img.header, affine = SPM_img.affine )
        nib.save( out, outfn )



if __name__=='__main__':
    import os, sys

    SPMmap = sys.argv[1]
    Tthres = float(sys.argv[2])
    Cthres = int(sys.argv[3])
    outdir = sys.argv[4]
    lab    = sys.argv[5]

    print("Image:\n\t{}".format(sys.argv[1]))
    print("T threshold:\n\t{}".format(sys.argv[2]))
    print("Minimum cluster size:\n\t{}".format(sys.argv[3]))

    clustersToROI( SPMmap, Tthres, Cthres, outdir, lab )
