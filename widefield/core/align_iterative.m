function [AlignNanThresh, ErrorInitialAlign, dsprealign, mimg] = align_iterative(data, ops)
% [AlignNanThresh, ErrorInitialAlign, dsprealign, mimg] = align_iterative(data, ops)
% iteratively finds mean image across frames (mimg) and image shift (dsprealign) 
%
% Inputs:
%   data: Y x X x nframes
%   ops must include:
%       NiterPrealign: number of iterations
%       
% Outputs:
%   AlignNanThrsh: putative threshold of correlation between target and frames after registration
%   ErrorInitialAlign: change in shift size between iteration step (i-1) and (i)
%   dsprealign: shift size after all iterations
%   mimg: average image across good frames after registration (after all iterations)
%
% called by determineTargetFrame.m

mimg = pick_reg_init(data);

dsold = zeros(size(data,3), 2);
err = zeros(ops.NiterPrealign, 1);

ops.SubPixel = Inf;
for i = 1:ops.NiterPrealign    
    fprintf(1, 'iteration %d/%d\n', i, ops.NiterPrealign);
    [dsnew, Corr]  = registration_offsets(data, ops, mimg, 1);
    dreg  = register_movie(data, ops, dsnew);
    [~, igood] = sort(Corr, 'descend');
    if i<floor(ops.NiterPrealign/2)        
        igood = igood(1:100);  
    else
        igood = igood(1:round(size(data,3)/2));  
    end
    mimg = mean(dreg(:,:,igood),3);
    
    err(i) = mean(sum((dsold - dsnew).^2,2)).^.5;
        
    dsold = dsnew;
end

AlignNanThresh = median(Corr) - 4*std(Corr);
ErrorInitialAlign = err;
dsprealign = dsnew;

end 
