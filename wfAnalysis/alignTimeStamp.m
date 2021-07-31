function [Values_align, idx_align] = alignTimeStamp(Values, TimeStamp, TimeStamp_align, method)
% alignTimeStamp(Values, TimeStamp, TimeStamp_align, method)
% Aligns N-dimentional tensor in time to the specified time stamps.
% Supported methods are 'nearest' or 'interp'
%
% Values_align = alignTimeStamp(Values, TimeStamp, TimeStamp_align) returns
% a 1~4-dimentional matrix (Nth dim corresponds time), in which time is aligned to TimeStamp_align.
% Length of the resulting matrix in time is the same as length(TimeStamp_align).
%
% [Stack_align, idx] = alignTimeStamp(..) also returns time index from
% original TimeStamp to the aligned TimeStamp


%2014-16-7 DS allow up to 4 dimension as Values
%2016-2-26 DS added option when dataDim=2
%2020-06-10 DS copied from stacks/+tools

if nargin < 4 || isempty(method)
    method = 'nearest';
end

dataDim = length(size(Values));

switch method
    case 'nearest'
        for tt = 1:length(TimeStamp_align)
            [~,idx_align(tt)] = min(abs(TimeStamp_align(tt) - TimeStamp));
        end
        
        switch dataDim
            case 1
                Values_align = Values(idx_align);
            case 2
                Values_align = Values(:,idx_align);
            case 3
                Values_align = Values(:,:,idx_align);
            case 4
                Values_align = Values(:,:,idx_align,:);
        end
        
    case 'interp'
        %         [xi, yi, zi] = meshgrid(1:size(Stack,2),1:size(Stack,1),TimeStamp);
        %         [xi_align, yi_align, zi_align] = meshgrid(1:size(Stack,2),1:size(Stack,1),TimeStamp_align);
        %zi_align = zi_align(~isnan(zi_align));
        
        %Stack_align = interp3(xi,yi,zi,Stack,xi_align,yi_align,zi_align);
        
        %          TimeVec_original = S.TimeVec;
        %     TimeVec_upsample = S.TimeVec(1):1/srate_upsample:S.TimeVec(end);
        
        
        switch dataDim
            case 1
            case 2
                 Values_align = interp1(TimeStamp, Values', TimeStamp_align, 'spline', 'extrap')';
                 %extrapolation does bad ...
            case {3, 4}
                [xo yo zo] = meshgrid(1:size(Values,2),1:size(Values,1),TimeStamp);
                [xq yq zq] = meshgrid(1:size(Values,2),1:size(Values,1),TimeStamp_align);
                
                if dataDim == 3
                    Values_align = single(interp3(xo,yo,zo,Values,xq,yq,zq,'spline'));
                elseif dataDim == 4
                    for ddd = 1:size(Values,4)
                        Values_align(:,:,:,ddd) = single(interp3(xo,yo,zo,Values(:,:,:,ddd),xq,yq,zq,'spline'));
                    end
                end
        end
        
        %this is just a place holder. needs to be revised
        idx_align = 1:length(TimeStamp_align);
        
end
