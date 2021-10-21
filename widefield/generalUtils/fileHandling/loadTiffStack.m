function stack = loadTiffStack(tiffFilename, mode, verbose)
% stack = loadTiffStack(tiffFilename, mode, verbose)
% loads tiff stacks with a specified mode
%
% Inputs:
%   tiffFilename: fullpath the tiff file
%   mode: specifies how load tiff {'imread', 'tiffobj'}
%       tiffobj is faster for USB-attached hard drive. 
%       imread is faster for local hard drive (I think!)
%   verbose: if true, show progress on the command line
%
% Outputs:
%   stack: Y x X x nFrames

%TODO: try
%https://au.mathworks.com/matlabcentral/fileexchange/32025-dylanmuir-tiffstack
%for faster loading

InfoImage=imfinfo(tiffFilename);
nImagesThisFile=length(InfoImage);

switch mode
    case 'imread'
    case 'tiffobj'
        t = Tiff(tiffFilename,'r');
end

stack = zeros(InfoImage(1).Height, InfoImage(1).Width, nImagesThisFile, 'uint16');

for i=1:nImagesThisFile
    if verbose && mod(i,100)==0
        fprintf(1, '  image %d\n', i);
    end
    
    switch mode
        case 'imread'
            stack(:,:,i) = imread(tiffFilename,i);
        case 'tiffobj'
            t.setDirectory(i);
            stack(:,:,i)=t.read();
    end
    w = warning ('off','all');
end

switch mode
    case 'imread'
    case 'tiffobj'
        t.close();
end


warning(w);