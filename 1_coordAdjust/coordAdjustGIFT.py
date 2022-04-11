
def findNetworksGIFT( giftdir, giftprefix, NWfile ):

    import os
    import re
    import numpy as np
    import nibabel as nib
    from roitools import world2voxel as w2v

    nwpath = os.path.dirname(NWfile)
    nwmodule = os.path.splitext(os.path.basename(NWfile))[0]
    sys.path.append(nwpath)
    NW = __import__(nwmodule)

    Tthresh = 3.0

    groupComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_tmap_component_ica_s1_.nii', f) and f.startswith(giftprefix)]
    groupTs = [f for f in os.listdir(giftdir) if re.match(r'(.*)_tmap_timecourses_ica_s1_.nii', f) and f.startswith(giftprefix)]

    if len(groupComp) == 0: #individual ICA's only got 1 subject level component
        groupComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_component_ica_s1_.nii', f) and f.startswith(giftprefix)][0]
        groupTs = [f for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_timecourses_ica_s1_.nii', f) and f.startswith(giftprefix)][0]
        Tthresh = 1.4
    else:
        groupComp = groupComp[0]
        groupTs = groupTs[0]

    gift_img = nib.load( os.path.join(giftdir, groupComp) )
    gift = gift_img.get_fdata()

    alltable = list()
    besttable = list()
    for n in NW.loaded_networks:
        #print(n)
        best = [0, 0.0, 0, n]
        for t in range(gift.shape[3]):
            current = [0, 0.0, t, n]
            for r in getattr(NW, n).keys():
                vcoord = w2v( os.path.join(giftdir, groupComp), getattr(NW,n)[r][0:3] )
                vcoordr = tuple( [int(round(v)) for v in vcoord] )
                val = gift[vcoordr][t]
                if val > Tthresh:
                    current[0] = current[0] + 1
                    current[1] = current[1] + val / len(getattr(NW,n).keys()) # unweighted
                    #current[1] = current[1] + val * getattr(NW,n)[r][4]       # weighted
            #print(current)
            if current[1] > best[1]:
                best = current
            #elif current[0] == best[0] and current[1] > best[1]:
            #    best = current
            alltable.append(current)
        #print(best)
        besttable.append(best)

    return alltable, besttable


def createNetworkROI( giftdir, giftprefix, NWfile, outdir ):
    import os
    import re
    import subprocess as sp
    from roitools import world2voxel as w2v
    from roitools import draw_sphere as sph

    nwpath = os.path.dirname(NWfile)
    nwmodule = os.path.splitext(os.path.basename(NWfile))[0]
    sys.path.append(nwpath)
    NW = __import__(nwmodule)

    try:
        os.makedirs( os.path.join(outdir, 'ROI') )
    except OSError:
        pass

    groupComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_tmap_component_ica_s1_.nii', f) and f.startswith(giftprefix)]

    if len(groupComp) == 0:
        groupComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_component_ica_s1_.nii', f) and f.startswith(giftprefix)][0]
    else:
        groupComp = groupComp[0]

    for n in NW.loaded_networks:
        for r in getattr(NW, n).keys():
            vcoord = w2v( os.path.join(giftdir, groupComp), getattr(NW,n)[r][0:3] )
            vcoordr = tuple( [int(round(v)) for v in vcoord] )
            radius = getattr(NW,n)[r][3]
            roi_file = os.path.join(outdir, 'ROI', n + '_' + r + '.nii.gz')
            print( '  ' + roi_file )
            sphere = sph( os.path.join(giftdir, groupComp), getattr(NW,n)[r][0:3], getattr(NW,n)[r][3] )
            sphere.to_filename( roi_file )
            #sp.check_call(['fslmaths', os.path.join(giftdir, groupComp), '-roi', str(vcoordr[0]), '1', str(vcoordr[1]), '1', str(vcoordr[2]), '1', '0', '1', '-kernel', 'sphere', str(radius), '-fmean', '-add', '100', '-bin', roi_file ])


def coordAdjustGIFT( giftdir, giftprefix, NWfile, comptable, roidir, outdir ):

    import os
    import re
    import numpy as np
    import nibabel as nib
    import pandas as pd
    from shutil import copyfile
    from roitools import img_roi_summary
    from roitools import world2voxel as w2v
    from roitools import draw_sphere as sph

    nwpath = os.path.dirname(NWfile)
    nwmodule = os.path.splitext(os.path.basename(NWfile))[0]
    sys.path.append(nwpath)
    NW = __import__(nwmodule)
    copyfile( NWfile, os.path.join(outdir, 'network_regions.py') )

    subjects = [re.match(r'(.*)_(sub.*)_component_ica_s1_.nii', f).group(2) for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_component_ica_s1_.nii', f) and f.startswith(giftprefix)]
    subjectsComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_component_ica_s1_.nii', f) and f.startswith(giftprefix)]
    subjectsTs = [f for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_timecourses_ica_s1_.nii', f) and f.startswith(giftprefix)]

    networks = list( comptable['network'] )
    comps = list( comptable['component'] )

    for s in range(len(subjects)):
        print(subjects[s])
        roitable = list()
        gift_img = nib.load( os.path.join(giftdir, subjectsComp[s]) )
        gift = gift_img.get_fdata()

        try:
            os.makedirs( os.path.join(outdir, 'ROI_' + subjects[s]) )
        except OSError:
            pass

        for n in range(len(networks)):

            for r in getattr(NW, networks[n]).keys():

                roi_file = os.path.join(roidir, 'ROI', networks[n] + '_' + r + '.nii.gz')
                roi_file_subject = os.path.join(outdir, 'ROI_' + subjects[s], networks[n] + '_' + r + '.nii.gz')

                roi_summary = img_roi_summary( gift, roi_file, comps[n] )

                sphere = sph( os.path.join(giftdir, subjectsComp[s]), [roi_summary['max_Xmm'], roi_summary['max_Ymm'], roi_summary['max_Zmm']], getattr(NW, networks[n])[r][5] ) #sphere with fix roi radius (8mm) around adjusted roi centre
                sphere.to_filename( roi_file_subject )

                roitable.append( [subjects[s], networks[n], r,
                                  roi_summary['max_Xmm'], roi_summary['max_Ymm'], roi_summary['max_Zmm'],
                                  roi_summary['max'], roi_summary['min'], roi_summary['median'], roi_summary['mean'], roi_summary['sd']] )

        roitable_df = pd.DataFrame( roitable, columns=['giftID', 'network', 'region', 'max_Xmm', 'max_Ymm', 'max_Zmm', 'max', 'min', 'median', 'mean', 'sd'] )
        roitable_df.to_csv( os.path.join(outdir, 'ROI_' + subjects[s], 'networks.csv') )

    return 0


def coordAdjustGroupGIFT( giftdir, giftprefix, NWfile, comptable, roidir, outdir ):

    import os
    import re
    import numpy as np
    import nibabel as nib
    import pandas as pd
    from shutil import copyfile
    from roitools import img_roi_summary
    from roitools import world2voxel as w2v
    from roitools import draw_sphere as sph

    nwpath = os.path.dirname(NWfile)
    nwmodule = os.path.splitext(os.path.basename(NWfile))[0]
    sys.path.append(nwpath)
    NW = __import__(nwmodule)
    copyfile( NWfile, os.path.join(outdir, 'network_regions.py') )

    groupComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_mean*_component_ica_s1_.nii', f) and f.startswith(giftprefix)]
    groupTs = [f for f in os.listdir(giftdir) if re.match(r'(.*)_mean*_timecourses_ica_s1_.nii', f) and f.startswith(giftprefix)]

    if len(groupComp) == 0:
        groupComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_component_ica_s1_.nii', f) and f.startswith(giftprefix)][0]
        groupTs = [f for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_timecourses_ica_s1_.nii', f) and f.startswith(giftprefix)][0]
    else:
        groupComp = groupComp[0]
        groupTs = groupTs[0]

    networks = list( comptable['network'] )
    comps = list( comptable['component'] )

    print(giftprefix)
    roitable = list()
    gift_img = nib.load( os.path.join(giftdir, groupComp) )
    gift = gift_img.get_fdata()

    try:
        os.makedirs( os.path.join(outdir, 'ROI', 'adjusted') )
    except OSError:
        pass

    for n in range(len(networks)):

        for r in getattr(NW, networks[n]).keys():

            roi_file = os.path.join(roidir, 'ROI', networks[n] + '_' + r + '.nii.gz')
            roi_file_adjusted = os.path.join(outdir, 'ROI', 'adjusted', networks[n] + '_' + r + '.nii.gz')

            roi_summary = img_roi_summary( gift, roi_file, comps[n] )

            sphere = sph( os.path.join(giftdir, groupComp), [roi_summary['max_Xmm'], roi_summary['max_Ymm'], roi_summary['max_Zmm']], getattr(NW, networks[n])[r][5] ) #sphere with fix roi radius (8mm) around adjusted roi centre
            sphere.to_filename( roi_file_adjusted )

            roitable.append( [giftprefix, networks[n], r,
                              roi_summary['max_Xmm'], roi_summary['max_Ymm'], roi_summary['max_Zmm'],
                              roi_summary['max'], roi_summary['min'], roi_summary['median'], roi_summary['mean'], roi_summary['sd']] )

    roitable_df = pd.DataFrame( roitable, columns=['giftID', 'network', 'region', 'max_Xmm', 'max_Ymm', 'max_Zmm', 'max', 'min', 'median', 'mean', 'sd'] )
    roitable_df.to_csv( os.path.join(outdir, 'ROI', 'adjusted', 'networks.csv') )

    return 0


def topICARegions( giftdir, giftprefix, outdir, df_alltable, top=8, clusthr=500 ):
    import os
    import re
    import subprocess as sp
    import pandas as pd
    import nibabel as nib
    import numpy as np
    from roitools import world2voxel as w2v
    from roitools import draw_sphere as sph
    from clusterize import clusterize

    comps = df_alltable.sort_values(by=['Tmean'], ascending=False).component[0:8]
    print(comps)

    outdir = os.path.join(outdir, 'topROI')
    try:
        os.makedirs( outdir )
    except OSError:
        pass

    groupComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_tmap_component_ica_s1_.nii', f) and f.startswith(giftprefix)]

    if len(groupComp) == 0:
        groupComp = [f for f in os.listdir(giftdir) if re.match(r'(.*)_sub[0-9]*_component_ica_s1_.nii', f) and f.startswith(giftprefix)][0]
    else:
        groupComp = groupComp[0]

    roitable = list()
    gift_img = nib.load( os.path.join(giftdir, groupComp) )
    gift = gift_img.get_fdata()
    print(gift.shape)

    for c in comps:
        data = gift[:,:,:,c]
        clusters, labels, sizes = clusterize( data, minsize = clusthr )
        print(labels)
        print(sizes)
        print(data[clusters == 1])
        exit()




if __name__ == '__main__':
    import os, sys
    import pandas as pd

    giftdir = sys.argv[1]
    giftprefix = sys.argv[2]
    NWfile = sys.argv[3]
    outdir = sys.argv[4]
    onlyFindNW = False

    outdir = os.path.join(outdir, giftprefix + '_network')

    try:
        os.makedirs(outdir)
    except OSError:
        pass
#    try:
#	os.makedirs(os.path.join(targetDir, 'network/icats'))
#    except OSError:
#	pass

    # 1. Find ICA components most corresponding to RSNs defined in NWfile (rsn_dict.py)
    alltable, besttable = findNetworksGIFT( giftdir, giftprefix, NWfile )
    df_besttable = pd.DataFrame(besttable, columns=['match', 'Tmean', 'component', 'network'])
    print('\nRSN components')
    print(df_besttable)
    df_besttable.to_csv( os.path.join(outdir, 'network_components.csv') )
    df_alltable = pd.DataFrame(alltable, columns=['match', 'Tmean', 'component', 'network'])
    df_alltable.to_csv( os.path.join(outdir, 'network_components_all.csv') )

    #topICARegions( giftdir, giftprefix, outdir, df_alltable )

    if not onlyFindNW:
        print('Create group ROIs')
        # 2. Create spherical ROIs of original regions
        createNetworkROI( giftdir, giftprefix, NWfile, outdir )

        print('Adjust coordinates')
        # 3. Create sperical ROIs of spatially adjusted regions
        # Adjust group-wise region coordinates
        coordAdjustGroupGIFT( giftdir, giftprefix, NWfile, df_besttable, outdir, outdir )
        # Adjust region coordinates for each subject
        #coordAdjustGIFT( giftdir, giftprefix, NWfile, df_besttable, outdir, outdir )

