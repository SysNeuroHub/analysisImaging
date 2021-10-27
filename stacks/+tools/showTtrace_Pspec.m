function [h, axtime, axspec] = showTtrace_Pspec(ttrace, taxis, visible, superimpose)
% h = showTtrace_Pspec(ttrace, frameRate)
% return figure handle for 1D time trace and corresponding powerspectrum,
% where ttrace is frames x conditions
%
% showTtrace_Pspec(ttrace, frameRate, 'on')
% shows the figure in the current window

%FrameRate in Hz

% 2014-03-31 DS
% 2014-06-22 ttrace is converted from a cell to matrix
% 2014-10-13 ttrace has to be a matrix rather than a cell
% 2015-08-13 2nd input is changed from frameRate to taxis
% 2016-09-16 powerspectrum using welch's method

if nargin < 4
    superimpose = false;
end
if nargin < 3
    visible = 'off';
end

color = [0 0 1;1 0 0; .1 1 .1];

numTraces = size(ttrace, 2);

frameRate = 1/median(diff(taxis));

%% timeVector
for pp = 1:numTraces
    nsteps = length(ttrace(:,pp));
    taxis_cell{pp} = taxis;%(1:nsteps)/frameRate; %[s]
end

%% powerspectrum
for pp = 1:numTraces
    %     [ faxis{pp}, pspec{pp}] = tools.signalPowerSpectrum( ttrace(:,pp), frameRate);
    L = size(ttrace,1);
    NFFT = 2^nextpow2(L);
    [pspec{pp},faxis{pp}] = pwelch(ttrace(:,pp),[],[],NFFT,frameRate);
end

%% range in y-axis
for pp = 1:numTraces
    y_cache_t(pp,:) = prctile(ttrace(:,pp), [1 99]);
    y_cache_f(pp,:) = prctile(pspec{pp}, [1 99]);
end
ymin_t = min(y_cache_t(:,1));
ymax_t = max(y_cache_t(:,2));
ymin_f = min(y_cache_f(:,1));
ymax_f = max(y_cache_f(:,2));

%% make figure
scrsz = get(0,'ScreenSize');
h = figure('Position',scrsz,'visible',visible);

if ~superimpose
    for pp = 1:numTraces
        
        subplot(numTraces, 2, 2*pp - 1);
        plot(taxis_cell{pp}, ttrace(:,pp)); %, 'color',color(:,pp));
        
        grid on;
        set(gca,'xticklabel','')
        xlim([taxis_cell{1}(1) taxis_cell{1}(end)]);
        ylim([ymin_t ymax_t]);
        
        if pp == numTraces
            xlabel('Time[s]');
            set(gca,'xticklabelmode','auto');
            axis on;
        end
        
        
        subplot(numTraces, 2, 2*pp);
        semilogy(faxis{pp}, pspec{pp}); %, 'color',color(:,pp));
        
        grid on;
        set(gca,'xticklabel','')
        xlim([faxis{1}(1) faxis{1}(end)]);
        ylim([ymin_f ymax_f]);
        
        if pp == numTraces
            xlabel('Frequency[Fz]');
            set(gca,'xticklabelmode','auto');
            axis on;
        end
    end
else
    for pp = 1:numTraces
        
        axtime=subplot(1, 2, 1);
        plot(taxis_cell{pp}, ttrace(:,pp)); %, 'color',color(:,pp));
        hold on
        grid on;
        set(gca,'xticklabel','')
        xlim([taxis_cell{1}(1) taxis_cell{1}(end)]);
        ylim([ymin_t ymax_t]);
        
        if pp == numTraces
            xlabel('Time[s]');
            set(gca,'xticklabelmode','auto');
            axis on;
        end
        
        
        axspec=subplot(1, 2, 2);
        semilogy(faxis{pp}, pspec{pp}); %, 'color',color(:,pp));
        hold on
        grid on;
        set(gca,'xticklabel','')
        xlim([faxis{1}(1) faxis{1}(end)]);
        ylim([ymin_f ymax_f]);
        
        if pp == numTraces
            xlabel('Frequency[Fz]');
            set(gca,'xticklabelmode','auto');
            axis on;
        end
    end
end
    
    
    
    %% example
    % [MyStack, ratioInfo, registerInfo] = StackSet.LoadStacks( DIRS.Stacks, p(iExp), ...
    % ResizeFac,iStim, iTr, Exps(iExp).Cam, registerInfo, 'ratioInfo',ratioInfo, 'registerInfo',registerInfo);
    % MyStack = MyStack.Crop([50 100],[50 100]);
    % MyStack = MyStack.Trim([1 12]);
    % ttrace = MyStack.SpaceAverages(1);
    % tools.showTtrace_Pspec(ttrace,MyStack.FrameRate,'on')
    
