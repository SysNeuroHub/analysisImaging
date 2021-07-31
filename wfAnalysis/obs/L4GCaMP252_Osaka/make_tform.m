

 blue = loadTiffStack('\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\L4GCaMP6s_252\2021-01-09\blue.tif',...
     'imread',0);
baseImage = imresize(double(blue),0.5);

blue = loadTiffStack('\\vault-v2.erc.monash.edu.au\MNHS-dshi0006\Subjects\L4GCaMP6s_252\2020-11-05\blue focus.tif',...
    'imread',0);
inputImage = imresize(double(blue), 0.5);

[t_concord, input_points, base_points, fliplrCam1, fliplrCam2, ...
    flipudCam1, flipudCam2,cam1_registered, resultFig, scale_recovered, theta_recovered] = ...
    rotate_flip_cams(inputImage, baseImage);
save('obs_to_new','t_concord');
