function repeatList = availableRepeats(protocol, stimList, FileString, sizeTh)
% repeatList = availableRepeats(protocol, stimList, FileString)
% returns repeat indexes by looking in the directory of the
% RAW data. If stimList is a vector, repeatList is a cell of the lists for
% each stimulus. 
%
% repeatList = availableRepeats(protocol, [], FileString)
% returns repeatList of all stimuli
%
% repeatList = availableRepeats(protocol, stimList, FileString, sizeTh)
% lets you specify the threshold for fraction of number of actual frames and expected frames.
% (default: 0.95)
%
% Currently only applicable to mpep experiment with pco.edge camera

% 2014-10-29 DS created
% 2015-08-14 DS moved under Stackset/+tools/

if nargin < 4
    sizeTh = 0.95;
end

if isempty(stimList)
    stimList = 1:protocol.nstim;
end

if ~strcmp(FileString, '_ratio')
    repeatList = validRepeats(protocol, stimList, FileString, sizeTh);
else
    repeatList_cam1 = validRepeats(protocol, stimList, '_cam1', sizeTh);
    repeatList_cam2 = validRepeats(protocol, stimList, '_cam2', sizeTh);
    
    for jjj = 1:length(stimList)
        repeatList{jjj} = intersect(repeatList_cam1{jjj}, repeatList_cam2{jjj});
    end
end



    function repeatList = validRepeats(protocol, stimList, FileString, sizeTh)
        
        
        ServerDir = tools.getServerDirectory(protocol.animal);
        AnimalDir   = fullfile(ServerDir, [protocol.animal FileString]); % 13.10.18 DS
        if isstr(protocol.iseries)
            SeriesDir   = fullfile(AnimalDir, sprintf('%s', protocol.iseries));
        else
            SeriesDir   = fullfile(AnimalDir, sprintf('%d', protocol.iseries));
        end
        ExpDir      = fullfile(SeriesDir, sprintf('%d', protocol.iexp));
        
        disp(['Detecting available repeats in ' ExpDir]);
        
        prefix = [protocol.animal FileString '_'];
        dd = dir([ExpDir '\' prefix '*.mat']);
        nfl = length(dd);
        iii = 1;
        for ifl = 1:nfl
            [name, exten] = strread(dd(ifl).name,'%s %s','delimiter','.');
            
            %if (length(name{1}) > length(prefix)) && (length(name{1}) < length(prefix)+4)
            serialNumber(iii) = str2num(name{1}(length(prefix)+1:end));
            bytes(iii) = dd(ifl).bytes;
            iii = iii + 1;
            %end
        end
        
        if nfl<1
            error('availableRepeats: no raw data is detected.');
        else
            %% load timestamps from raw data
            for iii = 1:length(stimList)
                iStim = stimList(iii);
                seqnums = protocol.seqnums(iStim,:);
                [~, repeatList_iStim, serialNumIdx] = intersect(seqnums, serialNumber);
                
                duration = protocol.pars(1,iii)/10;%[s]
                
                nframes = zeros(length(repeatList_iStim), 1);
                srate = zeros(length(repeatList_iStim), 1);
                for ifl = 1:length(repeatList_iStim)
                    try
                    [~, ~, TimeVecMM] = tools.LoadPCO( fullfile(ExpDir, dd(serialNumIdx(ifl)).name));%07/11/15
                    catch err %added on 24/12/15
                        nframes(ifl) = 0;
                        srate(ifl) = 0;
                        continue
                    end
                    nframes(ifl) = length(TimeVecMM);
                    srate(ifl) = 1/nanmedian(diff(TimeVecMM));%[Hz]
                end
                
                repeatList{iii} = repeatList_iStim(nframes > sizeTh * srate*duration);%';
                
                if size(repeatList{iii},1) > size(repeatList{iii},2)
                    repeatList{iii} = repeatList{iii}';
                end
            end
        end
    end


end

