function h = drawFTresult(ComplexMaps, Bright, visible, freqList)  
% h = drawFTresult(ComplexMaps, Bright, visible, freqList)  
% inputs:
% ComplexMaps: 3D, 3rd dimension can be idx of frequency or stimulus
%use matteobox

% 2016-6-29 show amplitude with jet, show phase in hsv

% TODO
% 3rd input "visible" is not working.

if nargin < 2
    Bright = 1.5;
end
if nargin < 3
    visible = 'on';
end


%cmaps = cat(3,ComplexMaps{:});
cmaps = ComplexMaps;
nConds = size(cmaps,3);

cmaps(isinf(cmaps)) = NaN;
lims = prctile(abs(cmaps(:)),[10 90]);
lims = round(lims*10000)/10000;

AbsMaps = abs(cmaps);
AngleMaps = angle(cmaps);


h = figure('visible',visible); 
%clf; 
%ax = gridplot( 3,nConds);
for icond = 1:nConds
    
    %% ABSMAP
    %axes(ax(1+3*(icond-1)));
    subplot(3,nConds,icond);
    imagesc(AbsMaps(:,:,icond));
    ax1 = gca;
    colormap(ax1,'jet');
    mcolorbar;

    AbsMaps_cache = AbsMaps(:,:,icond);
    caxis(prctile(AbsMaps_cache(:), [10 95]));
    axis equal tight
    set(gca,'xticklabel','','yticklabel','');
    if nargin > 3
        title([num2str(freqList(icond)) ' [Hz]']);
    end
    if icond == 1
        ylabel('amplitude');
        %title([num2str(firstTimeBase) ' - ' num2str(lastTimeBase) '[s], ' num2str(freqList(icond)) ' [Hz]']);
    end
    
    %% ANGLEMAP
    %axes(ax(2+3*(icond-1)));
    subplot(3,nConds,nConds+icond);
    imagesc(AngleMaps(:,:,icond));
    ax2 = gca;
    colormap(ax2,'hsv');
    caxis([-pi pi]);
    axis equal tight
    set(gca,'xticklabel','','yticklabel','');
    if icond == 1
        ylabel('phase [-pi ~ +pi]');
    end
    
    
    %% ABS+ANGLEMAP
    %axes(ax(3+3*(icond-1)));
    %ShowComplexImage( cmaps(:,:,icond), lims, ax(3+3*(icond-1)), Bright );
    subplot(3,nConds,2*nConds+icond);
    imagesc(AngleMaps(:,:,icond),'alphadatamapping','scaled','alphadata',AbsMaps(:,:,icond));
    ax3 = gca;
    colormap(ax3,'hsv');
    caxis([-pi pi]);
    axis equal tight
    set(gca,'xticklabel','','yticklabel','');
    if icond == 1
        ylabel('phase + amp');
    end
    
    colormap(jet);
end
            