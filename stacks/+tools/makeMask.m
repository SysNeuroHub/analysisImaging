function mask = makeMask( image, shape )
% mask = makeMask( image);
% mask = makeMask( image, shape )

%   walther May 2010
% 2013-10-09 DS return mask = 0 when roi not specified

switch nargin
    case 0
        error('Input error in subroutine RATIO OFFSET');
    case 1
        shape = 'box';
end

figure('visible','on')
% image = imshow(mat2gray(image), 'InitialMagnification', 'fit');
image = imagesc(image, prctile(image(:),[10 90]));%,'InitialMagnification', 'fit');
axis equal tight;
colormap(gray);
msgbox('Go to image, define a region of interest, then close the image window', 'modal'); 
if strcmp(shape,'box')
    try
        mask = createMask(imrect(gca), image);
    catch err
        mask = 0;
    end
elseif strcmp(shape,'elipse')
    try
        mask = createMask(imellipse(gca), image);
    catch err
        mask = 0;
    end
end

% check for keyboard input
uiwait(gcf)
if ishandle(gcf), close(gcf), end

end

