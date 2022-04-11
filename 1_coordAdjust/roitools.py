def img_roi_summary(img_file, mask_file, frame=0):
    import os, sys
    import numpy                        as np
    import nibabel                      as nib

    scriptsrc = os.path.realpath('__file__')
    scriptpath = os.path.dirname(scriptsrc)
    sys.path.append(scriptpath)

    #from convert_coordinates import voxel2world as v2w

    mask = nib.load(mask_file)
    mask_data = mask.get_data()
    if type(img_file) == type('str'):
        img = nib.load(img_file)
        img_data = img.get_fdata()
        if len(img_data.shape) == 4:
            img_data = img_data[:, :, :, frame]
    if type(img_file).__module__ == np.__name__ or type(img_file).__module__ == np.core.memmap.__module__:
        img  = img_file
        img_data = img[:, :, :, frame]
    if type(img_file).__module__ == 'nibabel.nifti1':
        img  = img_file
        img_data = img.get_fdata()[:, :, :, frame]

    if img_data.shape != mask.shape:
        print(img_data.shape)
        print(mask.shape)
        print( "volume mismatch" )
        exit()

    maxcoord = [0, 0, 0]
    maxval = -10000.0
    roivalues = []

    for x in range(mask.shape[0]):
        for y in range(mask.shape[1]):
            for z in range(mask.shape[2]):
                if mask_data[x][y][z] == 1:
                    roivalues.append(img_data[x][y][z])
                    if img_data[x][y][z] > maxval:
                        maxval = img_data[x][y][z]
                        maxcoord = voxel2world(mask_file, [x, y, z])

    nproivalues = np.array( roivalues )

    #print( "min:{0}".format( str(np.amin(nproivalues)) ) )
    #print( "max:{0}".format( str(np.amax(nproivalues)) ) )
    #print( "sd:{0}".format( str(np.nanstd(nproivalues)) ) )
    #print( "mean:{0}".format( str(np.nanmean(nproivalues)) ) )
    #print( "median:{0}".format( str(np.nanmedian(nproivalues)) ) )
    #print( "max_Xmm:{0}".format( str(maxcoord[0] )) )
    #print( "max_Ymm:{0}".format( str(maxcoord[1] )) )
    #print( "max_Zmm:{0}".format( str(maxcoord[2] )) )

    img_summary={'min' : np.amin(nproivalues),
                 'max' : np.amax(nproivalues),
                 'sd' : np.nanstd(nproivalues),
                 'mean' : np.nanmean(nproivalues),
                 'median' : np.nanmedian(nproivalues),
                 'max_Xmm' : maxcoord[0],
                 'max_Ymm' : maxcoord[1],
                 'max_Zmm' : maxcoord[2]}
    return( img_summary )


def voxel2world(img_file, coord):
    import nibabel as nib

    img = nib.load(img_file)

    out = list( nib.affines.apply_affine(img.get_affine(), coord) )

    return out


def world2voxel(img_file, coord):
    import nibabel  as nib
    import numpy    as np
    from   math import floor

    img = nib.load(img_file)

    out = list( nib.affines.apply_affine(np.linalg.inv(img.get_affine()), coord) )

    return out


def draw_sphere(img_file, wcoord, r):
  import nibabel  as nib
  import numpy    as np

  img = nib.load(img_file)
  data = img.get_fdata()
  hdr = img.header
  aff = img.get_affine()

  sphere_data = data[:, :, :, 0]
  sphere_data[:] = 0
  sphere_data = sphere_data.astype(int)

  voxel = world2voxel(img_file, wcoord)

  if ( r < min( hdr['pixdim'][1:4] ) ):
    sphere_data[int(voxel[0]), int(voxel[1]), int(voxel[2])] = 1
    sphere = nib.Nifti1Image( sphere_data, aff, hdr )
    return sphere

  X_voxelsize = hdr['pixdim'][1]
  Y_voxelsize = hdr['pixdim'][2]
  Z_voxelsize = hdr['pixdim'][3]

  sx = round(float(r) / X_voxelsize) * 2 + 1
  sy = round(float(r) / Y_voxelsize) * 2 + 1
  sz = round(float(r) / Z_voxelsize) * 2 + 1

  dx2 = (X_voxelsize) ** 2
  dy2 = (Y_voxelsize) ** 2
  dz2 = (Z_voxelsize) ** 2

  for z in np.arange(-(sz/2), (sz/2+1), 1):
    for y in np.arange(-(sy/2), (sy/2+1), 1):
      for x in np.arange(-(sx/2), (sx/2+1), 1):
        if (x*x*dx2 + y*y*dy2 + z*z*dz2) <= (r) ** 2:
          x_idx = int(x + sx/2 + voxel[0] - round(r/2))
          y_idx = int(y + sy/2 + voxel[1] - round(r/2))
          z_idx = int(z + sz/2 + voxel[2] - round(r/2))
          if x_idx > 0 and y_idx > 0 and z_idx > 0:
            sphere_data[x_idx, y_idx, z_idx] = 1

  sphere = nib.Nifti1Image( sphere_data, aff, hdr )
  return sphere





if __name__=='__main__':
    import os, sys
    from roi_mask import roi_sphere as sph

    group_img_file = sys.argv[1]
    subj_img_file = sys.argv[2]
    coord = [float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])]
    radius = 8
    name = 'tmpmask'
    mask_file = sph(group_img_file, name, coord, radius)
    img_roi_summary(subj_img_file, mask_file)
    os.remove(mask_file)
