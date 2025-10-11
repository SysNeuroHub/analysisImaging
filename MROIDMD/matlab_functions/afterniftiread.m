function properNifti = afterniftiread(loadedNifti)
%properNifti = preprocessAfterniftiread(loadedNifti)

% c = permute(loadedNifti,[1 3 2]);  
% properNifti = fliplr(rot90(c));
% properNifti=flip(rot90(permute(loadedNifti,[3 1 2]),2),3);
properNifti=rot90(permute(loadedNifti,[3 1 2]),2);

%now the niftidata should be
