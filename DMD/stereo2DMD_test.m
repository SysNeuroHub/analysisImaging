%% prepare DMD image
imageData = rgb2gray(imread('/home/daisuke/Documents/git/analysisImaging/MROIDMD/matlab_functions/star_800x500.png'));

% scale = 5/9;%
% width = scale*1440;
% height = scale*900;
% brainImage = zeros(height,width);
% MmPerPixel_t = 0.0104 / scale; %measured w scale 27/1/25 from getMmPerPix.m
% 
% f=figure;
% f.InnerPosition = [1 1 width height];
% image(zeros(size(brainImage)));colormap(gray); 
% hold on;
% line([width/2-scale*50 width/2-scale*50+1/MmPerPixel_t], scale*[800 800],'linewidth',2,'color','w');
% 
% %% add reference star
% ax = gca;
% mark=plot(width/2+1, height/2+1,'wp','MarkerSize',0.6*height,'MarkerFaceColor','none','MarkerEdgeColor','w');
% vline(width/2+1, ax, '-','w');
% hline(height/2+1, ax, '-','w');
% text(width/4, height/4,'TL','FontSize',30,'color','w')
% ax.Position = [0 0 1 1];
% 
% 
% %% superimpose grid on CCF
% xfrombregma = 3.5;%-2.7; %[mm]
% yfrombregma = -3.6; %A>0, P<0
% radiusmm = MmPerPixel_t*[1.5 5 10]; %[mm];
% 
% xfrombregmapix = 1/MmPerPixel_t * xfrombregma + bregma(2);
% yfrombregmapix = -1/MmPerPixel_t * yfrombregma + bregma(1);
% radiuspix = 1/MmPerPixel_t * radiusmm;
% 
% % hline(yfrombregmapix, ax,'-','w');
% % vline(xfrombregmapix, ax,'-','w');
% % 
% % screen2png([saveName '_grid'],f);%


% %% show only patches along the grid
% patchNumber = 0;
% for ii = 1%:2
%     for yy = 1:numel(yfrombregma)
%         for xx = 1:numel(xfrombregma)
%             for rr = 1:numel(radiusmm)
%                 patchNumber = patchNumber + 1;
% 
%                 fpatch=figure;
%                 fpatch.InnerPosition = [1 1 width height];
% 
%                 image(zeros(size(brainImage)));colormap(gray);
%                 if ii == 2
%                     suffix = '_wCCF';
%                     addAllenCtxOutlines(bregma, lambda, 'w', MmPerPixel_t);%this looks at lambda and shrinks the CCF
%                 else
%                     suffix = '';
%                 end
% 
%                 hold on;
%                 %         viscircles([xfrombregmapix(xx) yfrombregmapix(yy)], radiuspix(rr), 'color','w');
%                 rectangle('position',[xfrombregmapix(xx)-radiuspix(rr) yfrombregmapix(yy)-radiuspix(rr) ...
%                     2*radiuspix(rr) 2*radiuspix(rr)],'curvature',[1 1], ...
%                     'facecolor','w');
%                 %          plot(xfrombregmapix(xx), yfrombregmapix(yy),'rp','MarkerSize',40,'MarkerEdgeColor','w')
%                 axis ij image off
%                 xlim([1 size(brainImage,2)]);
%                 ylim([1 size(brainImage,1)]);
% 
%                 axpatch = gca;
%                 axpatch.Position = [0 0 1 1];
%                 screen2png([num2str(patchNumber) suffix],fpatch);
%                 close(fpatch);
%             end
%         end
%     end
% end

%% prepare allen brain volume
MRIdir = '/home/daisuke/Documents/git/analysisImaging/MROIDMD';
v=niftiread(fullfile(MRIdir, 'pattern_generation/Allen_annotation_modified.nii'));
allenVolume = (v~=0).*(v~=2000);

fv = isosurface(allenVolume, 0.5);

% Normalize x,y of surface vertices to image coordinates
xv = fv.vertices(:,2);
yv = fv.vertices(:,3);

% Scale to image dimensions
u = round( rescale(xv, 1, size(imageData,2)) );
v = round( rescale(yv, 1, size(imageData,1)) );

% Get color values from image
textureColors = zeros(size(fv.vertices,1), 3);  % for RGB

for k = 1:length(u)
    textureColors(k,:) = squeeze(imageData(v(k), u(k), :));  % sample color
end

%% visualize
figure
p = patch(fv, 'FaceVertexCData', textureColors, ...
              'FaceColor', 'interp', 'EdgeColor', 'none');
axis equal; camlight; lighting gouraud

% TODO
% register CCF2D contour from UCL to Kim
% DMD image 
