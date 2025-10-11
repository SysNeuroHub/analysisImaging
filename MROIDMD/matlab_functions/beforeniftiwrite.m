function properNifti = beforeniftiwrite(loadedNifti)
%properNifti = preprocessAfterniftiread(loadedNifti)

properNifti=flip(flip(permute(loadedNifti,[2 3 1]),3),1);

%now the niftidata should be
