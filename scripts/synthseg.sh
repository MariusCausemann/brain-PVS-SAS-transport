docker run --rm -ti \
  --user $(id -u):$(id -g) \
  -v "$(pwd)"/results/freesurfer:/usr/local/freesurfer/subjects \
  -v "$(pwd)"/freesurfer_license:/usr/local/freesurfer/licence \
  -v "$(pwd)"/data:/usr/local/freesurfer/data \
  -e FS_LICENSE=/usr/local/freesurfer/licence/license.txt \
  freesurfer/freesurfer:7.4.1 \
mri_synthseg --i /usr/local/freesurfer/data/T1.nii.gz --o /usr/local/freesurfer/subjects --threads 4 --robust