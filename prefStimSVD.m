function [xMap, yMap] = prefStimSVD(U, trigV_grid, screen_resize_scale, n_boot, use_method)
%[xMap, yMap] = prefStimSVD(U, trigV_grid)
% returns cortical map of stimulus preference from U and stim-triggered V
% trigV_grid{nstim}(npresentation, nSV)
%[...] = prefStimSVD(U, trigV_grid, screen_resize_scale)
% lets you to upsample in stimulus space(visual field) to estimate pref stim in higher
% resolution
%[...] = prefStimSVD(U, trigV_grid, screen_resize_scale, n_boot)
% estimates preferred stimulus n_boot times using boot strap (default=1)
%[...] = prefStimSVD(U, trigV_grid, screen_resize_scale, n_boot, use_method)
% lets you choose algorithm to estimate pref stim : 'max' (default), 'com'
%
% 7/7/20 created from AP_sparsenoise
% TODO: artefact near mask when U is resized

if nargin < 3
    screen_resize_scale = 1;
end
if nargin < 4
    n_boot = 1;
end
if nargin < 5
    use_method = 'max';
end

nx = size(trigV_grid,2); %stimulus dimension in x
ny = size(trigV_grid,1); %stimulus dimension in y

filter_sigma = (screen_resize_scale*2);

if n_boot > 1
    response_mean_bootstrap = cellfun(@(x) bootstrp(n_boot,@mean,x')',trigV_grid,'uni',false);%{visY,visX}(nSV, n_boot)
else
    response_mean_bootstrap = cellfun(@(x) mean(x,2), trigV_grid, 'uni', false);
end

% (to split trials instead of bootstrap)
%split_trials = cellfun(@(x) shake(discretize(1:size(x,2),round(linspace(1,size(x,2),n_boot+1)))),response_grid,'uni',false);
%response_mean_boostrap = cellfun(@(x,y) grpstats(x',y','mean')',response_grid,split_trials,'uni',false);
gauss_filt = fspecial('gaussian',[ny,nx],filter_sigma);
xMap = nan(size(U,1),size(U,2),n_boot);
yMap = nan(size(U,1),size(U,2),n_boot);
for curr_boot = 1:n_boot
%     keyboard
if n_boot > 1
    response_mean = cell2mat(cellfun(@(x) x(:,curr_boot),response_mean_bootstrap(:),'uni',false)');%[nSV x visPos]
else
    response_mean = cell2mat(cellfun(@(x) x,response_mean_bootstrap(:),'uni',false)');
end

    %brain response to each stimulus %[brainPosY,brainPosX,visPosition]
    stim_im = svdFrameReconstruct(U,response_mean);
    
    stim_im_px = reshape(permute(stim_im,[3,1,2]),ny,nx,[]); %[visPosY, visPosX, brainPos]    
    %filtering in visual space not in brain space. 
    stim_im_smoothed = imfilter(imresize(stim_im_px,screen_resize_scale,'bilinear'),gauss_filt);

    switch use_method
        case 'max'
            % Upsample each pixel's response map and find maximum            
            [~,mi] = max(reshape(stim_im_smoothed,[],size(stim_im_px,3)),[],1);
            [m_y,m_x] = ind2sub(size(stim_im_smoothed),mi);
            yMap(:,:,curr_boot) = reshape(m_y/screen_resize_scale,size(U,1),size(U,2));
            xMap(:,:,curr_boot) = reshape(m_x/screen_resize_scale,size(U,1),size(U,2));
            
        case 'com'
            % Conversely, do COM on original^2
            [xx,yy] = meshgrid(1:size(stim_im_smoothed,2),1:size(stim_im_smoothed,1));
            xMap = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,xx),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(U,1),size(U,2))/screen_resize_scale;
            yMap = reshape(sum(sum(bsxfun(@times,stim_im_smoothed.^2,yy),1),2)./sum(sum(stim_im_smoothed.^2,1),2),size(U,1),size(U,2))/screen_resize_scale;
    end
    
    disp(curr_boot);    
end
