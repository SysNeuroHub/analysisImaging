function mrimg_brain = getSurfaceData(mr_brain)
%mrimg_brain = getSurfaceData(mr_brain)
aa_brain=permute(mr_brain,[1 3 2]);  
[rows, cols, slices] = size(aa_brain);
surfZ = zeros(rows, cols);

%detect depth of the surface in z
for r = 1:rows
    for c = 1:cols
        idx = find(aa_brain(r,c,:) > 0, 1, 'last'); % first non-zero along z
        if ~isempty(idx)
            surfZ(r,c) = idx;   % store slice index of surface
        end
    end
end
mrimg_brain= double(fliplr(rot90(surfZ)));
