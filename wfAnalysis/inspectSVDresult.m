function inspectSVDresult(results, nSVD)
% inspectSVDresult(results, nSVD)
% produces three figures per vid (U, Ulast, V) for sanity check of SV
% decomposition
%
% 15/5/20 DS created from inspectSVDresult(Exps, suffix, nSVD)

if nargin< 3
    nSVD = [];
end
scrsz = get(0,'screensize');
% % %
% % % %where is ops info used to build SVD?
% %
% % thisDate = ['20' num2str(Exps.iseries(3:4)) '-' num2str(Exps.iseries(5:6)) ...
% %     '-' num2str(Exps.iseries(7:8))];
% % thisPath = fullfile('\\zserver\Data\Subjects', Exps.animal, thisDate);
% % %thisPath = '\\zserver\Data\Subjects\M010101_87\2016-06-27';
% % mkdir([thisPath '\' num2str(Exps.iexp)]);


for ss = 1:length(results)
    
    %% load cam2or1 SVD
    %load( [thisPath '\dataSummary_cam' num2str(cc) suffix{ss} '.mat'] );
    %Fs = 1/median(diff(dataSummary.timeStampsFromStamp)); %FIX THIS
    Fs = 15;
    %U_sm = readNPY([thisPath '\svdSpatialComponents_cam' num2str(cc) suffix{ss} '.npy']);
    %         try
    %             V = readNPY([thisPath '\' num2str(Exps.iexp) '\svdTemporalComponents_cam' num2str(cc) suffix{ss} '.npy']);
    %         catch err
    %             V = readNPY([thisPath '\svdTemporalComponents_cam' num2str(cc) suffix{ss} '.npy']);
    %         end
    %         V_sm = V';
    U_sm = results(ss).U;
    V_sm = results(ss).V;
    
    if ~isempty(nSVD)
        U_sm = U_sm(:,:,1:min(nSVD, size(U_sm,3)));
        V_sm = V_sm(1:min(size(V_sm,1), nSVD),:);
    end
    
    Sv_sm = results(ss).Sv;
    timeStamps_sm = results(ss).timeStampsFromStamp;
    meanImage_sm = imresize(results(ss).meanImage, [size(U_sm,1) size(U_sm,2)]);%hack this is due to image registration
    
    %  svdViewer(U_sm, Sv_sm{1}, V_sm{1}, 1/median(diff(dataSummary.timeStampsFromStamp)));
    %  movieWithTracesSVD(U_sm,V_sm{1},1:size(V_sm{1},1),[],[])
    
    figure('position',scrsz);
    for ii=1:24
        subplot(4,6,ii);
        U_c = U_sm(:,:,ii);
        imagesc(U_c);axis equal tight off
        caxis(prctile(U_c(:),[1 99]));
        if ii==1
            title(['U' results(ss).name]);
        else
            title(['comp: ' num2str(ii)]);
        end
    end
    %screen2png([thisPath '\' num2str(Exps.iexp) '\U cam' num2str(cc) suffix{ss}]);
    screen2png(['U' results(ss).name]);
    close
    
    figure('position',scrsz);
    for ii=1:24
        subplot(4,6,ii);
        
        jj = ii + size(U_sm,3) -24;
        U_c = U_sm(:,:,jj);
        imagesc(U_c);axis equal tight off
        caxis(prctile(U_c(:),[1 99]));
        if ii==1
            title(['U' results(ss).name]);
        else
            title(['comp: ' num2str(jj)]);
        end
    end
    %screen2png([thisPath '\' num2str(Exps.iexp) '\U cam' num2str(cc) suffix{ss} '_last']);
    screen2png(['U' results(ss).name '_last']);
    close
    
    figure('position',scrsz);
    subplot(221);
    imagesc(timeStamps_sm, 1:24, log(abs(V_sm(:,1:24)')));
    xlabel('time');
    ylabel('time course of log(abs(V)) initial');
    
    [Pxx1,F] = pwelch(V_sm(1:24,:)',[],[],[],Fs);
    [Pxx2,F] = pwelch(V_sm(end-24:end,:)',[],[],[],Fs);
    %[Pxx,F] = pwelch(V_sm',[],[],2^nextpow2(size(V_sm,2)),Fs);
    subplot(223);
    Pc = log(Pxx1');
    Pc(isinf(Pc))=0;
    Pc(isnan(Pc))=0;
    imagesc(F,1:24,Pc);
    caxis(prctile(Pc(:), [1 99]));
    title('pspec of V initial');
    colorbar;
    
    subplot(222);
    imagesc(timeStamps_sm, size(U_sm,3)+1-24:size(U_sm,3), ...
        log(abs(V_sm(:,size(U_sm,3)+1-24:size(U_sm,3))')));
    xlabel('time');
    ylabel('time course of log(abs(V)) last');
    
    subplot(224);
    Pc = log(Pxx2');
    Pc(isinf(Pc))=0;
    Pc(isnan(Pc))=0;
    imagesc(F,size(U_sm,3)+1-24:size(U_sm,3),Pc);
    caxis(prctile(Pc(:), [1 99]));
    title('pspec of V last');
    colorbar
    
    screen2png(['Vspec' results(ss).name]);
    close
    
    
end