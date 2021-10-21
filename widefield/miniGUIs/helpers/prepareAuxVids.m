function auxVid = prepareAuxVids(mouseName, thisDate, expNum)
%auxVid = prepareAuxVids(mouseName, thisDate, expNum)
% called in quckMovieWithVids

%movieDir = fileparts(dat.expFilePath(mouseName, thisDate, expNum, 'eyetracking', 'master'));
movieDir = sprintf('//zserver/Data/EyeCamera/%s/%s/%d',...
        mouseName, thisDate, expNum); %4/6/18

vidNum = 1;
auxVid = [];

facePath = fullfile(movieDir, 'face.mj2');
if exist(facePath, 'file')
    try
        faceT = fullfile(movieDir, 'face_timeStamps.mat');
        if ~exist(faceT, 'file')
            alignVideo(mouseName, thisDate, expNum, 'face');
        end
        vr = VideoReader(facePath);
        load(faceT);
        auxVid(vidNum).data = {vr, tVid};
        auxVid(vidNum).f = @plotMJ2frame;
        auxVid(vidNum).name = 'face';
        vidNum = vidNum+1;
    catch err
        disp(err);
    end
end

eyePath = fullfile(movieDir, 'eye.mj2');
if exist(eyePath, 'file')
    try
        eyeT = fullfile(movieDir, 'eye_timeStamps.mat');
        if ~exist(eyeT, 'file')
            alignVideo(mouseName, thisDate, expNum, 'eye');
        end
        load(eyeT);
        vr = VideoReader(eyePath);
        auxVid(vidNum).data = {vr, tVid};
        auxVid(vidNum).f = @plotMJ2frame;
        auxVid(vidNum).name = 'eye';
    catch err
        disp(err)
    end
end