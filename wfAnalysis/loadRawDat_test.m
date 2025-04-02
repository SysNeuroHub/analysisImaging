datPath = 'E:\svdinput\data\izanami\2025-01-30\vid1raw.dat';
imageSize = [900 1000];%results.imageSize;
nFr = 2e4;%results.nFrames;
batchSize = round(0.1*nFr);%ops.nRegisterBatchLimit; % images at once


fid = fopen(datPath);
imstack = fread(fid,  imageSize(1)*imageSize(2)*batchSize, '*uint16');
imstack = single(imstack);
imstack = reshape(imstack, imageSize(1), imageSize(2), []);
fclose(fid);
imstack = imresize(imstack, .1);

dFF = imstack./mean(imstack,3);

for ii = 1:size(imstack,3)
    imagesc(100*(dFF(:,:,ii)-1)); caxis([-1 1])
    drawnow;
    pause(.01);
end